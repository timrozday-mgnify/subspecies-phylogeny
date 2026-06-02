# subspecies-phylogeny

A Nextflow DSL2 pipeline for whole-genome SNP phylogeny of closely related bacterial genomes. It uses [SKA2](https://github.com/bacpop/ska.rust) (split k-mer analysis) to build a SNP alignment, [Gubbins](https://github.com/nickjcroucher/gubbins) to mask recombinant regions, and [IQ-TREE](http://www.iqtree.org/) to infer a maximum-likelihood phylogeny. [FastANI](https://github.com/ParBLiSS/FastANI) provides an all-vs-all pairwise identity check to flag outliers before tree construction.

Input genomes are expected to be closely related — same species or subspecies. This is the appropriate regime for SKA2, which breaks down if genomes are too diverged for k-mers to align.

---

## Contents

- [Quick start](#quick-start)
- [Intended workflow](#intended-workflow)
  - [Stage 1 — Explore](#stage-1--explore)
  - [Stage 2 — Build](#stage-2--build)
- [Genome QC](#genome-qc)
- [Pipeline parameters](#pipeline-parameters)
- [Output structure](#output-structure)
- [Configuring tool arguments](#configuring-tool-arguments)
  - [SKA2](#ska2-arguments)
  - [Gubbins](#gubbins-arguments)
  - [IQ-TREE](#iq-tree-arguments)
- [Samplesheet format](#samplesheet-format)
- [Reference databases](#reference-databases)
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

The pipeline runs in two stages. Stage 1 explores the dataset — running QC,
building the alignment at multiple thresholds, and assessing recombination —
without committing to tree inference. Stage 2 runs the final IQ-TREE phylogeny
with the parameters chosen from stage 1.

Ready-to-use examples live in `examples/`. Copy the `params.yml` and
`modules.config`, edit for your dataset, and pass them to Nextflow:

```bash
nextflow run main.nf -profile docker \
    -params-file examples/stage1_explore/params.yml \
    -c examples/stage1_explore/modules.config
```

---

### Stage 1 — Explore

**Example:** [`examples/stage1_explore/params.yml`](examples/stage1_explore/params.yml), [`examples/stage1_explore/modules.config`](examples/stage1_explore/modules.config)

Runs genome QC, builds the merged SKF, runs all-vs-all FastANI, trials a range
of `--min-freq` alignment thresholds in parallel, and runs Gubbins recombination
detection across all thresholds. IQ-TREE is skipped so the stage completes
quickly. Review the outputs, choose your final min-freq value, and identify any
genomes to remove before committing to stage 2.

```bash
nextflow run main.nf -profile docker \
    -params-file examples/stage1_explore/params.yml \
    -c examples/stage1_explore/modules.config
```

**Outputs to review** (`results/01_explore/`):

*Genome QC:*
- `qc/quast/<sample>/` — N50, contig count, GC%
- `qc/magpurify/<sample>/` — per-contig chimera flags
- `qc/checkm2/<sample>/` — completeness and contamination (ML-based)
- `qc/checkm/<sample>/` — completeness and contamination (HMM, legacy)
- `qc/gunc/<sample>/` — chimeric genome detection score
- `qc/busco/<sample>/` — lineage-specific ortholog completeness
- `fastani/fastani.txt` — pairwise ANI (flag samples < 95 % ANI to cohort)
- `nj_tree/fastani_nj.nwk` — quick NJ tree (check for long outlier branches)

*Alignment exploration:*
- `ska2/min_freq_*/alignment.fasta` — SNP alignment per threshold; compare lengths
- `snpsites/min_freq_*/` — variable-site counts per threshold
- `gubbins/min_freq_*/` — Gubbins recombination predictions per threshold

*Reuse in stage 2:*
- `ska2/merged.skf` — pass to `--ska_merged_skf` to skip rebuilding
- `ska2/skf/` — per-sample SKF basenames (for `--ska_delete_samples`)

**Decision guidance:**

| Signal | Action |
|---|---|
| CheckM2 completeness < 90 % or contamination > 5 % | Remove genome |
| GUNC maxCSS > 0.45 | Chimeric — inspect and consider removing |
| FastANI ANI < 95 % to cohort median | Wrong species — remove |
| Unusually high contig count (QUAST) | Fragmented assembly — review |
| Alignment length drops steeply above a threshold | Choose that threshold |

To remove problem genomes: add their FASTA basenames (without `.fa`/`.fasta`)
to a plain-text file (one per line) and uncomment `ska_delete_samples` in
`params.yml`. Confirm exact names with `ska nk results/01_explore/ska2/merged.skf`.

**Activate database-dependent QC tools** by uncommenting their paths in `params.yml`
(see [Reference databases](#reference-databases)):

| Tool | Database parameter |
|---|---|
| QUAST, MAGpurify, BUSCO | Always run (no database required by default) |
| CheckM2 | `checkm2_db` |
| CheckM v1 | `checkm_db` |
| GUNC | `gunc_db` |
| BUSCO (pre-downloaded, faster) | `busco_lineage_db` + `busco_lineage` |

---

### Stage 2 — Build

**Example:** [`examples/stage2_build/params.yml`](examples/stage2_build/params.yml), [`examples/stage2_build/modules.config`](examples/stage2_build/modules.config)

Runs the final IQ-TREE phylogeny using the parameters chosen from stage 1.
Reuses the merged SKF so the SKA2 build and FastANI steps are skipped.
Both the unmasked (`no_gubbins`) and Gubbins-masked (`gubbins`) trees are
inferred in parallel.

```bash
nextflow run main.nf -profile docker \
    -params-file examples/stage2_build/params.yml \
    -c examples/stage2_build/modules.config
```

**Outputs** (`results/02_build/`):
- `iqtree/no_gubbins/min_freq_0.95/tree.treefile` — ML tree without recombination masking
- `iqtree/gubbins/min_freq_0.95/tree.treefile` — ML tree with Gubbins masking
- `gubbins/min_freq_0.95/` — Gubbins recombination statistics
- `multiqc/multiqc_report.html`

---

## Genome QC

Before any phylogenetic analysis the pipeline runs a genome quality-control section that
assesses completeness, contamination, and chimeric content. Results are **informational
only** — no genomes are filtered from the downstream analysis. Skip the entire section with
`--skip_qc`.

| Tool | Purpose | Needs database? |
|---|---|---|
| [QUAST](https://quast.sourceforge.net/) | Assembly statistics (N50, contig count, GC%) | No |
| [MAGpurify](https://github.com/snayfach/MAGpurify) | Per-contig chimera detection (GC-content + tetranucleotide frequency) | No (optional for `phylo-markers` module) |
| [CheckM2](https://github.com/chklovski/CheckM2) | ML-based completeness and contamination | Yes — `--checkm2_db` |
| [CheckM](https://github.com/Ecogenomics/CheckM) v1 | HMM-based completeness and contamination (legacy) | Yes — `--checkm_db` |
| [GUNC](https://grp-bork.embl-community.io/gunc/) | Chimeric genome detection via gene taxonomy | Yes — `--gunc_db` |
| [BUSCO](https://busco.ezlab.org/) | Lineage-specific single-copy ortholog completeness | Yes — `--busco_lineage_db` |

QUAST and MAGpurify (gc-content + tetra-freq modules) run on every genome whenever
`--skip_qc` is not set, regardless of database availability. The remaining tools only run
when their database path is supplied. See [Reference databases](#reference-databases) for
download instructions and [Stage 1](#stage-1--explore) for the recommended workflow.

Outputs are published under `results/qc/<toolname>/<sample>/`.

---

## Pipeline parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | `null` | Path to samplesheet CSV (required unless `--ska_merged_skf` is set) |
| `--outdir` | `./results` | Output directory |
| `--ska_k` | `31` | K-mer size for `ska build`. Larger values increase specificity at the cost of sensitivity in more diverged genomes |
| `--ska_align_min_freq` | `0.9` | Comma-separated list of `--min-freq` values for `ska align`, e.g. `"0.5,0.9,1.0"`. Each value is run as a separate analysis branch |
| `--skip_phylo` | `false` | Skip all phylogenetic steps (FastANI, SKA2 build, alignment, Gubbins, IQ-TREE). Genome QC still runs. Use to run QC only |
| `--skip_alignment` | `false` | Skip `ska align` and all downstream steps. Produces only the merged SKF and FastANI outputs |
| `--skip_gubbins` | `false` | Skip Gubbins recombination filtering. IQ-TREE runs on snp-sites output only (no gubbins track) |
| `--skip_iqtree` | `false` | Skip IQ-TREE. Alignment, snp-sites, and Gubbins still run on all branches. Useful for inspecting alignments before tree inference |
| `--ska_merged_skf` | `null` | Path to a pre-computed `merged.skf`. When set, `ska build`, `ska merge`, and FastANI are skipped |
| `--ska_delete_samples` | `null` | Path to a plain-text file with one sample name per line. Those samples are removed from the merged SKF before alignment |
| `--ska_distance` | `false` | Run `ska distance` to produce a pairwise SNP distance table and NJ tree |
| `--ska_lo` | `false` | Run `ska lo` to identify SNPs and INDELs left out of the split-kmer graph (proxy for ambiguous regions) |
| `--ska_lo_reference` | `null` | Optional reference FASTA to anchor `ska lo` coordinates |
| `--skip_qc` | `false` | Skip the entire genome QC section (QUAST, MAGpurify, CheckM2, CheckM, GUNC, BUSCO) |
| `--checkm2_db` | `null` | Path to CheckM2 DIAMOND database file (`.dmnd`). CheckM2 only runs when this is set |
| `--checkm_db` | `null` | Path to CheckM v1 database root directory. CheckM only runs when this is set |
| `--gunc_db` | `null` | Path to GUNC DIAMOND database file (`.dmnd`). GUNC only runs when this is set |
| `--busco_lineage_db` | `null` | Path to pre-downloaded BUSCO lineage datasets directory. When not set, BUSCO auto-detects lineage and downloads at runtime (slower, needs internet in container) |
| `--busco_lineage` | `null` | BUSCO lineage name (e.g. `bacteria_odb10`). Defaults to `auto_prok` |
| `--magpurify_db` | `null` | Path to MAGpurify database directory for the `phylo-markers` module (optional) |
| `--multiqc_title` | `null` | Custom title for the MultiQC report |
| `--max_cpus` | `16` | Maximum CPUs available per process |
| `--max_memory` | `128.GB` | Maximum memory available per process |
| `--max_time` | `240.h` | Maximum wall time per process |

---

## Output structure

```
results/                              (example layout using stage1_explore → stage2_build)
  01_explore/                        stage 1 outputs
    qc/  fastani/  nj_tree/  ska2/   QC and build outputs (detail below)
    snpsites/  gubbins/              alignment exploration
  02_build/                          stage 2 outputs
    gubbins/  iqtree/  multiqc/

results/                              (QC detail, within 01_explore/)
  qc/
    quast/<sample>/
      report.html                     QUAST assembly report
      report.tsv.gz                   TSV version of the report
    magpurify/<sample>/
      *_gc_contaminants.txt           contigs flagged by GC-content check
      *_tetra_contaminants.txt        contigs flagged by tetranucleotide-frequency check
      *_phylo_contaminants.txt        contigs flagged by phylo-markers (if --magpurify_db)
    checkm2/<sample>/
      *_checkm2_report.tsv            completeness and contamination estimates
    checkm/<sample>/
      *_checkm_qa.tsv                 CheckM v1 lineage-workflow QA table
    gunc/<sample>/
      *.maxCSS_level.tsv              GUNC contamination assessment (max CSS level)
    busco/<sample>/
      *.batch_summary.txt             BUSCO batch summary
      short_summary.*.txt             per-lineage short summary
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

Tool-specific arguments are injected via `ext.args` in a Nextflow config file. The stage 1 example config at [`examples/stage1_explore/modules.config`](examples/stage1_explore/modules.config) contains commented-out overrides for each QC tool, SKA2, and Gubbins. The stage 2 config at [`examples/stage2_build/modules.config`](examples/stage2_build/modules.config) shows the IQ-TREE arguments. Create a `custom.config` for any run:

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

## Generating reports

Two Quarto reports are provided in `notebooks/` to explore and interpret pipeline results interactively. They are rendered with [Quarto](https://quarto.org/) using the `basicpython` conda kernel and require the `quarto` conda environment for rendering.

**Prerequisites:**

```bash
# quarto env needs papermill to inject -P parameters
pip install papermill          # or: conda install -c conda-forge papermill

# basicpython kernel needs: numpy pandas yaml toytree toyplot plotly
# Optional for static SVG/PNG export of Plotly figures:
mamba install -n basicpython -c conda-forge python-kaleido
```

### Highlights config

Strains to mark on trees and charts are defined in a YAML file rather than as command-line strings. Example files are in `notebooks/highlights/`:

```yaml
# notebooks/highlights/MyRun.yml
highlights:
  - name: AAYH02            # exact tip label as it appears in the tree
    label: "20HM reference" # shown in legend / tooltip
    color: "#e63946"        # optional; defaults to a built-in palette
  - name: BU_61_NT5381.1
    label: "Close relative"
```

Pass the **absolute path** to the file via `-P highlights_file:...`.

### Stage 1 — Explore report

```bash
conda run -n quarto quarto render \
    /path/to/subspecies-phylogeny/notebooks/stage1_explore.qmd \
    -P build_dir:/absolute/path/to/run/results/01_explore \
    -P explore_dir:/absolute/path/to/run/results/01_explore \
    -P qc_dir:/absolute/path/to/run/results/01_explore \
    -P run_name:MyRun \
    -P highlights_file:/absolute/path/to/subspecies-phylogeny/notebooks/highlights/MyRun.yml \
    --output-dir /absolute/path/to/run/notebooks/
```

Quarto resolves `-P` paths relative to the **notebook file**, not the working directory, so absolute paths are required.

### Stage 2 — Build report

```bash
conda run -n quarto quarto render \
    /path/to/subspecies-phylogeny/notebooks/stage2_build.qmd \
    -P build_dir:/absolute/path/to/run/results/01_explore \
    -P final_dir:/absolute/path/to/run/results/02_build \
    -P run_name:MyRun \
    -P highlights_file:/absolute/path/to/subspecies-phylogeny/notebooks/highlights/MyRun.yml \
    -P min_freq:0.95 \
    --output-dir /absolute/path/to/run/notebooks/
```

| Parameter | Description |
|---|---|
| `build_dir` | Stage 1 results directory (`01_explore`) — FastANI and SKA2 build outputs |
| `explore_dir` | Stage 1 results directory (`01_explore`) — alignment and Gubbins outputs |
| `qc_dir` | Stage 1 results directory (`01_explore`) — QC outputs |
| `final_dir` | Stage 2 results directory (`02_build`) — stage 2 report only |
| `run_name` | Label shown in figure titles |
| `highlights_file` | Path to highlights YAML (see above); omit to render without highlighting |
| `min_freq` | The `--min-freq` value used for the final run — stage 2 report only |

Rendered HTML files are self-contained (`embed-resources: true`) and can be shared without additional assets. Figures and tables are also saved to `results/<stage>/figs/stage<N>/` and `results/<stage>/tables/stage<N>/` respectively.

---

## Reference databases

The genome QC tools that require databases are not downloaded automatically. Download them
once to a stable location on your system and point the pipeline to them via the parameters
below.

### CheckM2

```bash
# Using the CheckM2 CLI (recommended — verifies the download)
checkm2 database --download --path /path/to/checkm2_db/

# Or download directly from Zenodo (~3.5 GB):
# https://zenodo.org/records/5571251/files/checkm2_database.tar.gz
```

Pass to the pipeline: `--checkm2_db /path/to/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd`

### CheckM v1

```bash
# Download the database archive (~1.4 GB):
# https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz

mkdir /path/to/checkm_db && tar xzf checkm_data_2015_01_16.tar.gz -C /path/to/checkm_db
```

Pass to the pipeline: `--checkm_db /path/to/checkm_db`

### GUNC

```bash
# Using the GUNC CLI (progenomes2 database, ~3.5 GB)
gunc download_db -db progenomes2 /path/to/gunc_db/

# GTDB database (~11 GB, higher resolution)
gunc download_db -db gtdb /path/to/gunc_db/
```

Pass to the pipeline: `--gunc_db /path/to/gunc_db/gunc_db_progenomes2.0.dmnd`

Reference: https://grp-bork.embl-community.io/gunc/install.html

### BUSCO

```bash
# Download a specific lineage (recommended for prokaryotes)
busco --download bacteria_odb10 --download_path /path/to/busco_lineages/

# Browse available lineages:
# https://busco-data.ezlab.org/v5/data/lineages/
```

Pass to the pipeline:
- `--busco_lineage_db /path/to/busco_lineages/`
- `--busco_lineage bacteria_odb10` (or another lineage name, or leave unset to use `auto_prok`)

### MAGpurify (phylo-markers module, optional)

```bash
# Download from Zenodo (~2 GB) and extract
wget https://zenodo.org/record/3688811/files/MAGpurify-db-v1.0.tar.bz2
tar xjf MAGpurify-db-v1.0.tar.bz2 -C /path/to/magpurify_db/
```

Pass to the pipeline: `--magpurify_db /path/to/magpurify_db/MAGpurify-db-v1.0`

The MAGpurify GC-content and tetranucleotide-frequency modules run without any database.

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
| QUAST | 5.2.0 | `quay.io/biocontainers/quast:5.2.0--py39pl5321heaaa4ec_4` |
| MAGpurify | 2.1.2 | `quay.io/biocontainers/magpurify:2.1.2--pyhdfd78af_2` |
| CheckM2 | 1.0.1 | `quay.io/biocontainers/checkm2:1.0.1--pyh7cba7a3_0` |
| CheckM v1 | 1.2.3 | `quay.io/biocontainers/checkm-genome:1.2.3--pyhdfd78af_1` |
| GUNC | 1.0.6 | `quay.io/biocontainers/gunc:1.0.6--pyhdfd78af_0` |
| BUSCO | 5.4.7 | `quay.io/biocontainers/busco:5.4.7--pyhdfd78af_0` |
