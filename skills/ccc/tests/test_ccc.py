"""CLI tests for scripts/ccc.py: demo goldens, output contract, fail strings.

Run from skills/ccc with pytest; see references/maintainer.md "Running and testing".
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
SCRIPT = SKILL_DIR / "scripts" / "ccc.py"
DEMO = SKILL_DIR / "examples" / "demo.h5ad"
RESOURCE = SKILL_DIR / "examples" / "demo_resource.tsv"
TEMPLATES = SKILL_DIR.parents[1] / "templates"

# Golden demo values: the two planted sender -> receiver pairs, and nothing else.
PLANTED = [
    {"source": "sender", "target": "receiver", "ligand": "LIG1", "receptor": "REC1"},
    {"source": "sender", "target": "receiver", "ligand": "LIG2", "receptor": "REC2"},
]

ARGS = {
    "--adata": str(DEMO),
    "--resource": str(RESOURCE),
    "--groupby": "cell_type",
    "--methods": "rank_aggregate",
    "--expr-prop": "0.1",
    "--min-cells": "5",
    "--n-perms": "100",
    "--magnitude-max": "0.002",
    "--specificity-max": "0.01",
    "--sources": "all",
    "--targets": "all",
    "--n-top": "10",
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


def run_with(tmp_path, drop=(), out="out", **changes):
    """Demo arguments with `changes` ({"n_top": "0"}) and `drop` applied."""
    args = {**ARGS, **{f"--{k.replace('_', '-')}": v for k, v in changes.items()}}
    argv = [x for k, v in args.items() if k not in drop for x in (k, v)]
    return run(*argv, "--outdir", str(tmp_path / out))


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
    assert (result["n_cells"], result["n_groups"], result["n_resource_pairs"]) == (
        120,
        3,
        20,
    )
    assert (result["n_interactions"], result["n_pass"]) == (180, 2)
    assert (result["n_matched"], result["n_in_resource"]) == (40, 40)
    assert result["reference"] == "synthetic demo genes, no genome build"
    assert result["top_interactions"] == PLANTED
    assert result["seed"] == 0 and len(result["methods"]) == 5
    interactions = pd.read_csv(demo_out / "tables" / "interactions.tsv", sep="\t")
    head = interactions.head(2)[["source", "target", "ligand", "receptor"]]
    assert head.to_dict("records") == PLANTED
    assert interactions["passes"].tolist() == [True, True] + [False] * 178
    counts = pd.read_csv(demo_out / "tables" / "pair_counts.tsv", sep="\t")
    assert len(counts) == 9 and counts["n_pass"].sum() == 2
    planted_cell = (counts["source"] == "sender") & (counts["target"] == "receiver")
    assert counts.loc[planted_cell, "n_pass"].item() == 2


def test_interactions_sorted_by_ranks(demo_out):
    interactions = pd.read_csv(demo_out / "tables" / "interactions.tsv", sep="\t")
    keys = ["magnitude_rank", "specificity_rank", "source", "target", "ligand"]
    assert interactions.index.equals(interactions.sort_values(keys).index)


def test_rerun_is_identical(demo_out, tmp_path):
    completed = run("--demo", "--outdir", str(tmp_path / "again"))
    assert completed.returncode == 0, completed.stderr
    for name in ("interactions.tsv", "pair_counts.tsv"):
        again = (tmp_path / "again" / "tables" / name).read_bytes()
        assert again == (demo_out / "tables" / name).read_bytes(), name


def test_per_sample_subset(tmp_path):
    completed = run_with(
        tmp_path,
        sample="sample",
        methods="cellphonedb,natmi",
        sources="sender",
        magnitude_max="0.005",
        specificity_max="0.05",
    )
    assert completed.returncode == 0, completed.stderr
    out = tmp_path / "out"
    interactions = pd.read_csv(out / "tables" / "interactions.tsv", sep="\t")
    assert set(interactions["source"]) == {"sender"}
    passing = interactions[interactions["passes"]]
    for sample in ("s1", "s2"):
        rows = passing[passing["sample"] == sample]
        assert rows[["source", "target", "ligand", "receptor"]].to_dict("records") == (
            PLANTED
        ), sample
    counts = pd.read_csv(out / "tables" / "pair_counts.tsv", sep="\t")
    assert list(counts.columns) == ["sample", "source", "target", "n_pass"]
    assert len(counts) == 18 and counts["n_pass"].sum() == 4


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
    completed = run_with(tmp_path, drop=["--n-perms"])
    assert is_error(completed, "missing required arguments: --n-perms")


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--seed", "1", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --seed")


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"expr_prop": "2"}, "--expr-prop must be in [0, 1], got 2.0"),
        ({"n_perms": "0"}, "--n-perms must be an integer >= 1, got 0"),
        ({"methods": "cellchat"}, "--methods: unknown cellchat"),
        ({"methods": "logfc"}, "--methods logfc has no magnitude score"),
        ({"adata": "/no/such.h5ad"}, "input file not found: /no/such.h5ad"),
        ({"groupby": "leiden"}, "obs column not found: leiden"),
        ({"sample": "donor"}, "obs column not found: donor"),
        ({"targets": "receiver,tcell"}, "--targets: labels not in cell_type: tcell"),
        ({"magnitude_max": "1e-9"}, "zero interactions pass magnitude_rank <= 1e-09"),
    ],
)
def test_input_errors(tmp_path, changes, message):
    assert is_error(run_with(tmp_path, **changes), message)


def test_resource_without_receptor(tmp_path):
    path = tmp_path / "resource.tsv"
    path.write_text("ligand\tpartner\nLIG1\tREC1\n")
    completed = run_with(tmp_path, resource=str(path))
    assert is_error(completed, "resource column not found: receptor")


def test_resource_without_overlap(tmp_path):
    path = tmp_path / "resource.tsv"
    path.write_text("ligand\treceptor\nTGFB1\tTGFBR1_TGFBR2\n")
    completed = run_with(tmp_path, resource=str(path))
    assert is_error(completed, "no resource pair has its ligand and receptor")


def test_raw_counts(tmp_path):
    def counts(adata):
        adata.X = adata.layers["counts"]

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, counts))
    assert is_error(completed, "X does not look log-normalised: all integers")


def test_negative_values(tmp_path):
    def scale(adata):
        adata.X = adata.X.toarray() - 1

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, scale))
    assert is_error(completed, "X does not look log-normalised: negative values")


def test_one_group(tmp_path):
    def one_group(adata):
        adata.obs["cell_type"] = "sender"

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, one_group))
    assert is_error(completed, "--groupby cell_type needs at least 2 groups, found 1")


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
