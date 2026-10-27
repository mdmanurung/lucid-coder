"""CLI tests for scripts/sc_markers.py: demo goldens, output contract, fail strings.

Run from skills/sc-markers with pytest;
see references/maintainer.md "Running and testing".
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "sc_markers.py"
DEMO = SKILL_DIR / "examples" / "demo.h5ad"
TEMPLATES = SKILL_DIR.parents[1] / "templates"

# Golden demo values: the planted markers, and nothing else, pass.
PLANTED = {
    "c1": {"G000", "G001", "G002", "G003", "G004"},
    "c2": {"G005", "G006", "G007", "G008", "G009"},
    "c3": {"G010", "G011", "G012", "G013", "G014"},
}

ARGS = {
    "--adata": str(DEMO),
    "--groupby": "cluster",
    "--method": "wilcoxon",
    "--padj-max": "0.05",
    "--lfc-min": "1",
    "--min-pct": "0.25",
    "--n-top": "5",
    "--limit": "2",
    "--reference": "synthetic demo genes, no genome build",
    "--seed": "0",
}


def run(*argv):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *argv],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )


def run_with(tmp_path, drop=(), **changes):
    """Demo arguments with `changes` ({"min_pct": "2"}) and `drop` applied."""
    args = {**ARGS, **{f"--{k.replace('_', '-')}": v for k, v in changes.items()}}
    argv = [x for k, v in args.items() if k not in drop for x in (k, v)]
    return run(*argv, "--outdir", str(tmp_path / "out"))


def contract_paths():
    """File paths listed in the SKILL.md output tree, relative to outdir."""
    text = (SKILL_DIR / "SKILL.md").read_text()
    tree = text.split("## Output structure")[1].split("```")[1]
    paths, stack = [], []
    for line in tree.splitlines():
        match = re.search(r"[├└]── (\S+)", line)
        if not match:
            continue
        depth = match.start() // 4
        stack[depth:] = [match.group(1)]
        if not match.group(1).endswith("/"):
            paths.append("".join(stack))
    return paths


def is_error(completed, message):
    return (
        completed.returncode == 2
        and f"ERROR: {message}" in completed.stderr
        and "Traceback" not in completed.stderr
    )


def edited_demo(tmp_path, edit):
    adata = ad.read_h5ad(DEMO)
    edit(adata)
    path = tmp_path / "edited.h5ad"
    adata.write_h5ad(path)
    return str(path)


@pytest.fixture(scope="module")
def demo_out(tmp_path_factory):
    out = tmp_path_factory.mktemp("demo") / "out"
    completed = run("--demo", "--outdir", str(out))
    assert completed.returncode == 0, completed.stderr
    return out


def test_demo_golden_values(demo_out):
    result = json.loads((demo_out / "result.json").read_text())
    assert (result["n_cells"], result["n_groups"], result["n_genes_tested"]) == (
        150,
        3,
        200,
    )
    assert result["n_pass"] == {"c1": 5, "c2": 5, "c3": 5}
    assert {g: set(genes) for g, genes in result["top_markers"].items()} == PLANTED
    assert result["seed"] == 0 and result["method"] == "wilcoxon"
    assert result["reference"] == "synthetic demo genes, no genome build"
    markers = pd.read_csv(demo_out / "tables" / "markers.tsv", sep="\t")
    assert len(markers) == 600
    passing = markers[markers["passes"]].groupby("group")["gene"].agg(set)
    assert passing.to_dict() == PLANTED
    top = pd.read_csv(demo_out / "tables" / "top_markers.tsv", sep="\t")
    assert top["rank"].tolist() == [1, 2, 3, 4, 5] * 3


def test_markers_sorted_by_group_then_score(demo_out):
    markers = pd.read_csv(demo_out / "tables" / "markers.tsv", sep="\t")
    resorted = markers.sort_values(
        ["group", "score", "gene"], ascending=[True, False, True]
    )
    assert markers.index.equals(resorted.index)


def test_t_test_recovers_planted(tmp_path):
    completed = run_with(tmp_path, method="t-test")
    assert completed.returncode == 0, completed.stderr
    result = json.loads((tmp_path / "out" / "result.json").read_text())
    assert {g: set(genes) for g, genes in result["top_markers"].items()} == PLANTED
    assert result["n_pass"] == {"c1": 5, "c2": 5, "c3": 5}


def test_output_contract(demo_out):
    paths = contract_paths()
    assert len(paths) == 13
    for rel in paths:
        assert (demo_out / rel).stat().st_size > 0, rel
    check = subprocess.run(
        ["sha256sum", "-c", "reproducibility/checksums.sha256"],
        cwd=demo_out,
        capture_output=True,
        text=True,
        check=False,
    )
    assert check.returncode == 0, check.stdout + check.stderr
    assert "--seed 0" in (demo_out / "reproducibility" / "commands.sh").read_text()


def test_missing_argument_fails_without_demo(tmp_path):
    completed = run_with(tmp_path, drop=["--limit"])
    assert is_error(completed, "missing required arguments: --limit")


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--seed", "1", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --seed")


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"min_pct": "2"}, "--min-pct must be in [0, 1], got 2.0"),
        ({"limit": "0"}, "--limit must be > 0, got 0.0"),
        ({"adata": "/no/such.h5ad"}, "input file not found: /no/such.h5ad"),
        ({"layer": "logcounts"}, "layer not found in the h5ad: logcounts"),
        ({"groupby": "leiden"}, "obs column not found: leiden"),
        (
            {"layer": "counts"},
            "layer counts does not look log-normalised: all integers (raw counts?)",
        ),
        ({"lfc_min": "100"}, "zero genes pass padj < 0.05, logfc > 100.0"),
    ],
)
def test_input_errors(tmp_path, changes, message):
    assert is_error(run_with(tmp_path, **changes), message)


def test_negative_values(tmp_path):
    def scale(adata):
        adata.X = adata.X.toarray() - 1

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, scale))
    assert is_error(completed, "X does not look log-normalised: negative values")


def test_one_group(tmp_path):
    def one_group(adata):
        adata.obs["cluster"] = "c1"

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, one_group))
    assert is_error(completed, "--groupby cluster needs at least 2 groups, found 1")


def test_small_group(tmp_path):
    def singleton(adata):
        adata.obs["cluster"] = adata.obs["cluster"].astype(str)
        adata.obs.iloc[0, 0] = "c4"

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, singleton))
    assert is_error(completed, "groups with fewer than 2 cells in cluster: c4")


@pytest.mark.parametrize(("axis", "what"), [("var", "gene names"), ("obs", "cell ids")])
def test_duplicate_ids(tmp_path, axis, what):
    def duplicate(adata):
        names = list(getattr(adata, f"{axis}_names"))
        names[1] = names[0]
        setattr(adata, f"{axis}_names", names)

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, duplicate))
    assert is_error(completed, f"duplicate {what}: ")


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
