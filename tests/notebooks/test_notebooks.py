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
    cmd = ["quarto", "render", str(NOTEBOOKS_DIR / notebook)]
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
