# subspecies-phylogeny

A Nextflow DSL2 pipeline for whole-genome SNP phylogeny of closely related bacterial genomes. It uses [SKA2](https://github.com/bacpop/ska.rust) (split k-mer analysis) to build a SNP alignment, [Gubbins](https://github.com/nickjcroucher/gubbins) to mask recombinant regions, and [IQ-TREE](http://www.iqtree.org/) to infer a maximum-likelihood phylogeny. [FastANI](https://github.com/ParBLiSS/FastANI) provides an all-vs-all pairwise identity check to flag outliers before tree construction.

Input genomes are expected to be closely related — same species or subspecies. This is the appropriate regime for SKA2, which breaks down if genomes are too diverged for k-mers to align.

---

## Contents

- [Quick start](#quick-start)
- [Intended workflow](#intended-workflow)
- [Pipeline parameters](#pipeline-parameters)
- [Output structure](#output-structure)
- [Configuring tool arguments](#configuring-tool-arguments)
  - [SKA2](#ska2-arguments)
  - [Gubbins](#gubbins-arguments)
  - [IQ-TREE](#iq-tree-arguments)
- [Samplesheet format](#samplesheet-format)
- [Requirements](#requirements)

---

## Quick start

```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --outdir results
```

---

## Intended workflow

The pipeline is designed to be run in stages so that intermediate outputs can be inspected before committing to the computationally expensive steps.

Each stage has a ready-to-use example in the `examples/` directory. Copy the relevant `params.yml` and `modules.config`, edit them for your dataset, and pass them to Nextflow:

```bash
nextflow run main.nf -profile docker \
    -params-file examples/stage1_build/params.yml
```

Stages 2, 3, and 4 also take a `-c modules.config` to set tool-specific arguments:

```bash
nextflow run main.nf -profile docker \
    -params-file examples/stage3_trial/params.yml \
    -c examples/stage3_trial/modules.config
```

---

### Step 1 — Build the merged SKF

**Example:** [`examples/stage1_build/params.yml`](examples/stage1_build/params.yml)

Run with `--skip_alignment` to build per-sample `.skf` files and merge them. FastANI all-vs-all runs at the same time, producing a pairwise ANI table and a neighbour-joining tree that can be used to spot outliers.

```bash
nextflow run main.nf -profile docker \
    -params-file examples/stage1_build/params.yml
```

**Outputs to inspect:**
- `results/01_build/ska2/merged.skf` — merged split-kmer file for reuse in later steps
- `results/01_build/fastani/fastani.txt` — pairwise ANI table
- `results/01_build/nj_tree/fastani_nj.nwk` — quick NJ tree from ANI distances; check for any genomes that cluster far from the rest

If any genomes look like outliers in the FastANI tree, add their names to a plain-text file (one name per line) for use with `--ska_delete_samples` in the next step. Sample names correspond to the FASTA filename basename (without extension), which may differ from the samplesheet `sample` column — check `ska2/skf/` filenames to confirm exact names, or run `ska nk results/01_build/ska2/merged.skf`.

### Step 2 — Explore min-freq values and remove outliers

**Example:** [`examples/stage2_explore/params.yml`](examples/stage2_explore/params.yml), [`examples/stage2_explore/modules.config`](examples/stage2_explore/modules.config)

Use the pre-computed merged SKF to trial a range of `--min-freq` values for `ska align`. This step is fast because it skips the genome-building phase. Optionally remove outlier samples identified in step 1.

```bash
nextflow run main.nf -profile docker \
    -params-file examples/stage2_explore/params.yml \
    -c examples/stage2_explore/modules.config
```

`--min-freq` controls what fraction of samples must have a called base at a position for it to be included in the alignment. Lower values retain more positions at the cost of more missing data; higher values give a cleaner but smaller alignment. Examine the alignment lengths in `ska2/min_freq_*/alignment.fasta` alongside the Gubbins recombination maps in `gubbins/min_freq_*/` to choose a value.

**Outputs to inspect:**
- `ska2/min_freq_*/alignment.fasta` — SNP alignment at each frequency threshold; compare lengths
- `snpsites/min_freq_*/` — variable-sites counts
- `gubbins/min_freq_*/` — Gubbins recombination predictions at each threshold

### Step 3 — Trial Gubbins parameters

**Example:** [`examples/stage3_trial/params.yml`](examples/stage3_trial/params.yml), [`examples/stage3_trial/modules.config`](examples/stage3_trial/modules.config)

With a chosen min-freq value, run the full alignment block (with `--skip_iqtree`) to inspect the Gubbins outputs before committing to tree inference. Adjust Gubbins parameters in `modules.config` (see [Configuring tool arguments](#configuring-tool-arguments)).

```bash
nextflow run main.nf -profile docker \
    -params-file examples/stage3_trial/params.yml \
    -c examples/stage3_trial/modules.config
```

Re-run this stage (updating `modules.config`) until satisfied with the recombination predictions and filtered alignments.

**Outputs to inspect:**
- `gubbins/min_freq_0.9/*.recombination_predictions.gff` — predicted recombinant regions
- `gubbins/min_freq_0.9/*.filtered_polymorphic_sites.fasta` — recombination-free alignment that will feed IQ-TREE; check this is non-empty and has a reasonable number of variable sites

### Step 4 — Run IQ-TREE to produce final phylogenies

**Example:** [`examples/stage4_final/params.yml`](examples/stage4_final/params.yml), [`examples/stage4_final/modules.config`](examples/stage4_final/modules.config)

Once satisfied with the alignment and Gubbins parameters, copy the finalised settings from stage 3 into `examples/stage4_final/modules.config`, add IQ-TREE model and bootstrap arguments, then run the full pipeline. Both the unmasked (`no_gubbins`) and Gubbins-masked (`gubbins`) trees are inferred in parallel.

```bash
nextflow run main.nf -profile docker \
    -params-file examples/stage4_final/params.yml \
    -c examples/stage4_final/modules.config
```

**Outputs:**
- `results/04_final/iqtree/no_gubbins/min_freq_0.9/tree.treefile`
- `results/04_final/iqtree/gubbins/min_freq_0.9/tree.treefile`
- `results/04_final/multiqc/multiqc_report.html`

---

## Pipeline parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | `null` | Path to samplesheet CSV (required unless `--ska_merged_skf` is set) |
| `--outdir` | `./results` | Output directory |
| `--ska_k` | `31` | K-mer size for `ska build`. Larger values increase specificity at the cost of sensitivity in more diverged genomes |
| `--ska_align_min_freq` | `0.9` | Comma-separated list of `--min-freq` values for `ska align`, e.g. `"0.5,0.9,1.0"`. Each value is run as a separate analysis branch |
| `--skip_alignment` | `false` | Skip `ska align` and all downstream steps. Produces only the merged SKF and FastANI outputs |
| `--skip_gubbins` | `false` | Skip Gubbins recombination filtering. IQ-TREE runs on snp-sites output only (no gubbins track) |
| `--skip_iqtree` | `false` | Skip IQ-TREE. Alignment, snp-sites, and Gubbins still run on all branches. Useful for inspecting alignments before tree inference |
| `--ska_merged_skf` | `null` | Path to a pre-computed `merged.skf`. When set, `ska build`, `ska merge`, and FastANI are skipped |
| `--ska_delete_samples` | `null` | Path to a plain-text file with one sample name per line. Those samples are removed from the merged SKF before alignment |
| `--ska_distance` | `false` | Run `ska distance` to produce a pairwise SNP distance table and NJ tree |
| `--ska_lo` | `false` | Run `ska lo` to identify SNPs and INDELs left out of the split-kmer graph (proxy for ambiguous regions) |
| `--ska_lo_reference` | `null` | Optional reference FASTA to anchor `ska lo` coordinates |
| `--multiqc_title` | `null` | Custom title for the MultiQC report |
| `--max_cpus` | `16` | Maximum CPUs available per process |
| `--max_memory` | `128.GB` | Maximum memory available per process |
| `--max_time` | `240.h` | Maximum wall time per process |

---

## Output structure

```
results/
  ska2/
    skf/                          per-sample .skf files
    merged.skf                    merged split-kmer file (all samples)
    merged_deleted.skf            merged SKF after ska delete (if --ska_delete_samples)
    distances.tsv                 pairwise SNP distances (if --ska_distance)
    lo_output_snps.fas            left-out SNPs (if --ska_lo)
    lo_output_indels.vcf          left-out INDELs (if --ska_lo)
    min_freq_<f>/
      alignment.fasta             SNP alignment
  snpsites/
    min_freq_<f>/
      <id>.fas                    variable-sites FASTA
      <id>.sites.txt              constant-sites counts for IQ-TREE correction
  gubbins/
    min_freq_<f>/
      *.filtered_polymorphic_sites.fasta    recombination-free alignment
      *.recombination_predictions.gff       predicted recombinant regions
      *.summary_of_snp_distribution.vcf     SNP distribution summary
  iqtree/
    {no_gubbins,gubbins}/min_freq_<f>/
      tree.treefile               ML phylogeny (Newick)
      tree.log                    IQ-TREE log
  fastani/
    fastani.txt                   all-vs-all ANI table
  nj_tree/
    fastani_nj.nwk                NJ tree from FastANI distances
    ska2_distance_nj.nwk          NJ tree from ska distance (if --ska_distance)
  multiqc/
    multiqc_report.html
  pipeline_info/
    execution_trace_<datetime>.txt      per-task CPU/memory/time stats
    execution_report_<datetime>.html    interactive run summary
    execution_timeline_<datetime>.html  Gantt-style task timeline
    pipeline_dag_<datetime>.html        workflow DAG
```

---

## Configuring tool arguments

Tool-specific arguments are injected via `ext.args` in a Nextflow config file. Create a `custom.config` alongside your run command:

```bash
nextflow run main.nf -profile docker -c custom.config --input samplesheet.csv --outdir results
```

The sections below describe the most important arguments for each tool, and show how to set them.

### SKA2 arguments

The k-mer size (`--ska_k`) is the most important global SKA2 parameter and is set as a pipeline parameter rather than via `ext.args`.

#### `ska build`

```groovy
// custom.config
process {
    withName: 'SKA2_BUILD' {
        ext.args = { "-k ${params.ska_k} --single-strand" }
    }
}
```

| Argument | Default | Notes |
|---|---|---|
| `-k` | `31` | K-mer size. Set via `--ska_k`. Increase (e.g. 41, 51) to reduce spurious k-mer matches in more diverged genomes; decrease for very short reads or fragmented assemblies. Must be odd. |
| `--single-strand` | off | Only use the forward strand. Use for single-stranded data or if palindromic k-mers cause issues. |
| `--min-count` | `1` | Minimum count to include a k-mer (relevant for read input, not assemblies). |

#### `ska align`

The `--min-freq` value is set per analysis branch via `--ska_align_min_freq`. Additional arguments can be appended via `ext.args`.

```groovy
process {
    withName: 'SKA2_ALIGN' {
        // --min-freq is always set by the pipeline; append extra args here
        ext.args = { "--min-freq ${meta.min_freq} --filter no-filter" }
    }
}
```

| Argument | Default | Notes |
|---|---|---|
| `--min-freq` | — | **Set by the pipeline** via `--ska_align_min_freq`. Fraction of samples that must have a called base at a position for it to be included. Higher values give cleaner alignments with fewer positions; lower values retain more positions with more missing data. |
| `--filter` | `core` | Position filter: `no-filter`, `core` (present in all samples), or `bi-allelic`. `core` produces a full-core alignment; `no-filter` maximises the number of sites. |
| `--constant-sites` | off | Include invariant sites in the output. Usually not needed; snp-sites is used downstream to strip them. |

---

### Gubbins arguments

Gubbins detects recombinant regions and removes them from the alignment before tree inference. The pipeline always runs Gubbins on every alignment branch; the `no_gubbins` IQ-TREE track uses the snp-sites alignment directly without Gubbins masking.

```groovy
process {
    withName: 'GUBBINS' {
        ext.args = '--tree-builder fasttree --iterations 5 --min-snps 3'
    }
}
```

| Argument | Default (pipeline) | Notes |
|---|---|---|
| `--tree-builder` | `fasttree` | Tree inference engine used internally. Options: `raxml`, `raxmlng`, `fasttree`, `iqtree`, `iqtree-fast`, `rapidnj`. `fasttree` is the pipeline default because it handles small sample counts (<4) that cause RAxML to fail. Switch to `raxmlng` or `iqtree` for larger datasets. |
| `--iterations` | `5` | Maximum number of Gubbins iterations. More iterations refine recombination boundaries but increase runtime. |
| `--min-snps` | `3` | Minimum number of SNPs required to call a recombination event. Increase to reduce false positives in very similar genomes. |
| `--min-window-size` | `100` | Minimum length of a predicted recombination block (bp). |
| `--max-window-size` | `10000` | Maximum length of a recombination block to consider. |
| `--filter-percentage` | `25.0` | Sequences with more than this percentage of masked bases are excluded from recombination detection. |
| `--model` | `GTRGAMMA` | Substitution model for tree inference. `GTR` or `GTRGAMMA` are standard; use `JC` for very small alignments. |
| `--first-tree-builder` | same as `--tree-builder` | Tree builder for the first Gubbins iteration only. Use `rapidnj` for a fast initial tree on large datasets. |
| `--tree-args` | — | Extra arguments passed directly to the tree builder executable. |

**Note on alignment input:** SKA2 alignments may contain IUPAC ambiguity codes (e.g. `R`, `Y`). The pipeline pre-processes the alignment with `awk` to replace these with `N` before passing to Gubbins, since Gubbins only accepts `ACGTNacgtn-`.

---

### IQ-TREE arguments

IQ-TREE infers a maximum-likelihood phylogeny. The pipeline passes `-fconst` automatically on the `no_gubbins` track to correct for constant sites stripped by snp-sites. On the `gubbins` track this correction is omitted because Gubbins already provides a variable-sites-only alignment.

```groovy
process {
    withName: 'IQTREE' {
        ext.args = { (meta.constant_sites ? "-fconst ${meta.constant_sites}" : '') + ' -m GTR+G -B 1000' }
    }
}
```

**Important:** always include `meta.constant_sites ? "-fconst ${meta.constant_sites}" : ''` in your `ext.args` closure to preserve the ascertainment-bias correction on the `no_gubbins` track.

| Argument | Default (pipeline) | Notes |
|---|---|---|
| `-fconst A,C,G,T` | set automatically | Constant-site counts for ascertainment-bias correction. The pipeline provides these from snp-sites on the `no_gubbins` track; do not remove this from `ext.args`. |
| `-m` | `MFP` (auto) | Substitution model. `MFP` triggers ModelFinder to select the best model. Specify explicitly (e.g. `-m GTR+G`, `-m GTR+G+ASC`) to skip model selection and reduce runtime. Note: do not use `+ASC` on the `no_gubbins` track — the pipeline handles constant-sites correction via `-fconst` instead. |
| `-B` | — | Number of ultrafast bootstrap replicates (e.g. `-B 1000`). Adds support values to the tree. |
| `--ufboot` / `-b` | — | Alternative bootstrap methods: `-b 100` for standard (slower) bootstrap. |
| `-T` | `AUTO` | Number of threads. Leave as `AUTO` to let IQ-TREE choose, or set explicitly. Overridden by the Nextflow process `cpus` directive. |
| `--seed` | — | Random seed for reproducibility. |
| `-nt AUTO` | — | Synonym for `-T AUTO` in older IQ-TREE versions. |
| `--polytomy` | off | Collapse near-zero branches into polytomies. |

---

## Samplesheet format

A CSV file with at minimum two columns: `sample` and `fasta`.

```csv
sample,fasta
isolate_A,/path/to/isolate_A.fasta
isolate_B,/path/to/isolate_B.fasta
isolate_C,/path/to/isolate_C.fasta
```

- `sample` — used as the output label in per-sample files (e.g. `.skf` filenames).
- `fasta` — absolute path or path relative to the project directory. HTTP/FTP URLs are also accepted.

**Note:** the sample name stored inside the merged SKF (and used by `ska delete`) is derived from the **FASTA filename basename** (without extension), which may differ from the `sample` column. For example, `sample=BU_61` with `fasta=BU_61_NT5381.1.fa` will store the name `BU_61_NT5381.1` in the SKF. Run `ska nk merged.skf` to list the stored sample names.

---

## Requirements

- [Nextflow](https://nextflow.io/) ≥ 25.0
- Docker, Singularity/Apptainer, or Conda (select with `-profile docker|singularity|conda`)

No other software needs to be installed; all tools are pulled from container images automatically.

| Tool | Version | Container |
|---|---|---|
| SKA2 | 0.5.1 | `quay.io/biocontainers/ska2:0.5.1--h4349ce8_0` |
| FastANI | — | `quay.io/biocontainers/fastani:1.34--h0ffd775_2` |
| snp-sites | — | nf-core/snpsites |
| Gubbins | 3.4.3 | `quay.io/biocontainers/gubbins:3.4.3--py310h5140242_0` |
| IQ-TREE | — | nf-core/iqtree |
| R (ape) | ≥ 5.8 | `community.wave.seqera.io/library/r-ape:5.8--48d6804841ebe369` |
| MultiQC | — | nf-core/multiqc |
