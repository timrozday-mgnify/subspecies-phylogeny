"""Pytest tests for Quarto report rendering.

Each notebook is rendered once per test session (session-scoped fixtures) against
the fixture data in tests/notebooks/results/.  Individual tests then make assertions
on the exit code, output file, and HTML content.

Running locally (requires the quarto conda env and basicpython kernel):
    pytest tests/notebooks/ -v

In CI the basicpython kernel is registered by the workflow before pytest runs.
"""
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent.parent
NOTEBOOKS_DIR = REPO_ROOT / "notebooks"
FIXTURE_DIR = Path(__file__).parent


def _render(notebook: str, params: dict, output_dir: Path) -> subprocess.CompletedProcess:
    cmd = ["conda", "run", "-n", "quarto", "quarto", "render", str(NOTEBOOKS_DIR / notebook)]
    for key, value in params.items():
        cmd += ["-P", f"{key}:{value}"]
    cmd += ["--output-dir", str(output_dir)]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=300)


# ---------------------------------------------------------------------------
# Session-scoped render fixtures — notebook is rendered once for all tests
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def stage1_render(tmp_path_factory):
    out = tmp_path_factory.mktemp("stage1")
    result = _render(
        "stage1_explore.qmd",
        {
            "results_dir": str(FIXTURE_DIR / "results" / "01_explore"),
            "run_name": "test",
        },
        out,
    )
    return result, out


@pytest.fixture(scope="session")
def stage1_skip_large_pairwise_render(tmp_path_factory):
    out = tmp_path_factory.mktemp("stage1_skip_large_pairwise")
    result = _render(
        "stage1_explore.qmd",
        {
            "results_dir": str(FIXTURE_DIR / "results" / "01_explore"),
            "run_name": "test",
            "skip_large_pairwise_plots": "true",
        },
        out,
    )
    return result, out


@pytest.fixture(scope="session")
def stage1_metadata_render(tmp_path_factory):
    out = tmp_path_factory.mktemp("stage1_metadata")
    metadata_file = tmp_path_factory.mktemp("stage1_metadata_file") / "genome_metadata.csv"
    metadata_file.write_text(
        "sample,type\n"
        "BU_61_NT5381,isolate\n"
        "BU_909_NT5401,mag\n",
        encoding="utf-8",
    )
    result = _render(
        "stage1_explore.qmd",
        {
            "results_dir": str(FIXTURE_DIR / "results" / "01_explore"),
            "run_name": "test",
            "metadata_file": str(metadata_file),
            "skip_large_pairwise_plots": "true",
        },
        out,
    )
    return result, out


@pytest.fixture(scope="session")
def stage1_malformed_gff_render(tmp_path_factory):
    run = tmp_path_factory.mktemp("stage1_malformed_gff_results")
    out = tmp_path_factory.mktemp("stage1_malformed_gff")
    gubbins_dir = run / "gubbins" / "min_freq_0.9"
    gubbins_dir.mkdir(parents=True)
    (gubbins_dir / "0.9.per_branch_statistics.csv").write_text(
        "\t".join([
            "Node",
            "Number of Recombination Blocks",
            "Total SNPs",
            "r/m",
            "Bases in Clonal Frame",
            "Genome Length",
            "Number of SNPs Inside Recombinations",
        ])
        + "\n"
        + "\t".join(["sampleA", "0", "5", "0", "1000", "1000", "0"])
        + "\n",
        encoding="utf-8",
    )
    (gubbins_dir / "0.9.recombination_predictions.gff").write_text(
        "\n".join([
            "##gff-version 3",
            "##sequence-region SEQUENCE 1 1000",
            "SEQUENCE\tGUBBINS\tCDS\tSEQUENCE\tGUBBINS\tCDS\t384\t512\t0.000\t.\t0\t"
            'node="Node_1->sampleA";taxa="sampleA";snp_count="27";',
            "",
        ]),
        encoding="utf-8",
    )
    result = _render(
        "stage1_explore.qmd",
        {
            "results_dir": str(run),
            "run_name": "malformed_gff",
            "skip_large_pairwise_plots": "true",
        },
        out,
    )
    return result, out


@pytest.fixture(scope="session")
def stage2_render(tmp_path_factory):
    out = tmp_path_factory.mktemp("stage2")
    result = _render(
        "stage2_build.qmd",
        {
            "explore_dir": str(FIXTURE_DIR / "results" / "01_explore"),
            "build_dir": str(FIXTURE_DIR / "results" / "02_build"),
            "run_name": "test",
            "min_freq": "0.9",
        },
        out,
    )
    return result, out


@pytest.fixture(scope="session")
def stage2_skip_large_pairwise_render(tmp_path_factory):
    out = tmp_path_factory.mktemp("stage2_skip_large_pairwise")
    result = _render(
        "stage2_build.qmd",
        {
            "explore_dir": str(FIXTURE_DIR / "results" / "01_explore"),
            "build_dir": str(FIXTURE_DIR / "results" / "02_build"),
            "run_name": "test",
            "min_freq": "0.9",
            "skip_large_pairwise_plots": "true",
        },
        out,
    )
    return result, out


@pytest.fixture(scope="session")
def stage2_metadata_render(tmp_path_factory):
    out = tmp_path_factory.mktemp("stage2_metadata")
    metadata_file = tmp_path_factory.mktemp("stage2_metadata_file") / "genome_metadata.csv"
    metadata_file.write_text(
        "sample,type\n"
        "BU_61_NT5381,isolate\n"
        "BU_909_NT5401,mag\n",
        encoding="utf-8",
    )
    result = _render(
        "stage2_build.qmd",
        {
            "explore_dir": str(FIXTURE_DIR / "results" / "01_explore"),
            "build_dir": str(FIXTURE_DIR / "results" / "02_build"),
            "run_name": "test",
            "min_freq": "0.9",
            "metadata_file": str(metadata_file),
            "skip_large_pairwise_plots": "true",
        },
        out,
    )
    return result, out


# ---------------------------------------------------------------------------
# Stage 1 tests
# ---------------------------------------------------------------------------

class TestStage1Explore:
    def test_renders_without_error(self, stage1_render):
        result, _ = stage1_render
        assert result.returncode == 0, f"quarto render failed:\n{result.stderr}"

    def test_html_output_exists(self, stage1_render):
        _, out = stage1_render
        assert (out / "stage1_explore.html").exists()

    def test_html_output_non_empty(self, stage1_render):
        _, out = stage1_render
        assert (out / "stage1_explore.html").stat().st_size > 10_000

    @pytest.mark.parametrize("text", [
        "test",                     # run_name injected via -P
        "FastANI",
        "Neighbour-joining tree",
        "Alignment statistics",
        "Gubbins",
    ])
    def test_html_contains(self, stage1_render, text):
        _, out = stage1_render
        html = (out / "stage1_explore.html").read_text(encoding="utf-8")
        assert text in html, f"Expected '{text}' in stage1_explore.html"

    def test_skip_large_pairwise_plots_renders(self, stage1_skip_large_pairwise_render):
        result, out = stage1_skip_large_pairwise_render
        assert result.returncode == 0, f"quarto render failed:\n{result.stderr}"
        assert (out / "stage1_explore.html").exists()

    @pytest.mark.parametrize("text", [
        "The Pairwise ANI heatmap was skipped",
        "The Pairwise SKA2 SNP distance matrix was skipped",
        "The ANI vs SNP distance scatterplot was skipped",
    ])
    def test_skip_large_pairwise_plots_callouts(self, stage1_skip_large_pairwise_render, text):
        _, out = stage1_skip_large_pairwise_render
        html = (out / "stage1_explore.html").read_text(encoding="utf-8")
        assert text in html, f"Expected '{text}' in skipped stage1_explore.html"

    def test_malformed_gubbins_gff_rows_are_skipped(self, stage1_malformed_gff_render):
        result, out = stage1_malformed_gff_render
        assert result.returncode == 0, f"quarto render failed:\n{result.stderr}"
        html = (out / "stage1_explore.html").read_text(encoding="utf-8")
        assert "Malformed Gubbins GFF rows skipped" in html
        assert "min_freq=0.9: 1 row(s)" in html

    def test_compact_unlabeled_nj_tree_is_rendered(self, stage1_render):
        _, out = stage1_render
        html = (out / "stage1_explore.html").read_text(encoding="utf-8")
        assert "Compact unlabeled NJ tree" in html

    def test_optional_genome_metadata_highlights_render(self, stage1_metadata_render):
        result, out = stage1_metadata_render
        assert result.returncode == 0, f"quarto render failed:\n{result.stderr}"
        html = (out / "stage1_explore.html").read_text(encoding="utf-8")
        assert "Genome metadata legend" in html
        assert "Isolate" in html
        assert "MAG" in html


# ---------------------------------------------------------------------------
# Stage 2 tests
# ---------------------------------------------------------------------------

class TestStage2Build:
    def test_renders_without_error(self, stage2_render):
        result, _ = stage2_render
        assert result.returncode == 0, f"quarto render failed:\n{result.stderr}"

    def test_html_output_exists(self, stage2_render):
        _, out = stage2_render
        assert (out / "stage2_build.html").exists()

    def test_html_output_non_empty(self, stage2_render):
        _, out = stage2_render
        assert (out / "stage2_build.html").stat().st_size > 10_000

    @pytest.mark.parametrize("text", [
        "test",                     # run_name injected via -P
        "IQ-TREE",
        "ML phylogeny",
        "Branch-length",
        "Gubbins",
    ])
    def test_html_contains(self, stage2_render, text):
        _, out = stage2_render
        html = (out / "stage2_build.html").read_text(encoding="utf-8")
        assert text in html, f"Expected '{text}' in stage2_build.html"

    def test_skip_large_pairwise_plots_renders(self, stage2_skip_large_pairwise_render):
        result, out = stage2_skip_large_pairwise_render
        assert result.returncode == 0, f"quarto render failed:\n{result.stderr}"
        assert (out / "stage2_build.html").exists()

    def test_skip_large_pairwise_plots_callout(self, stage2_skip_large_pairwise_render):
        _, out = stage2_skip_large_pairwise_render
        html = (out / "stage2_build.html").read_text(encoding="utf-8")
        assert "The Pairwise IQ-TREE ML distance heatmap was skipped" in html

    @pytest.mark.parametrize("text", [
        "Compact unlabeled IQ-TREE ML phylogeny — no recombination masking",
        "Compact unlabeled IQ-TREE ML phylogeny — Gubbins-masked",
    ])
    def test_compact_phylogenies_are_rendered(self, stage2_render, text):
        _, out = stage2_render
        html = (out / "stage2_build.html").read_text(encoding="utf-8")
        assert text in html

    def test_optional_genome_metadata_highlights_render(self, stage2_metadata_render):
        result, out = stage2_metadata_render
        assert result.returncode == 0, f"quarto render failed:\n{result.stderr}"
        html = (out / "stage2_build.html").read_text(encoding="utf-8")
        assert "Genome metadata legend" in html
        assert "Right-side ticks mark genome type (see legend above) and the FastANI medoid in red" in html
