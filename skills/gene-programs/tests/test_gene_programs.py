"""CLI tests for scripts/gene_programs.py: demo goldens, output contract, fail strings.

Run from skills/gene-programs with pytest;
see references/maintainer.md "Running and testing".
"""

import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "gene_programs.py"
DEMO = SKILL_DIR / "examples" / "demo.h5ad"
TEMPLATES = SKILL_DIR.parents[1] / "templates"

# Golden demo values: each program's top 10 genes are one planted set.
PLANTED = {
    "GP1": {f"G{i:03d}" for i in range(10, 20)},  # P2, used by g2
    "GP2": {f"G{i:03d}" for i in range(0, 10)},  # P1, used by g1
    "GP3": {f"G{i:03d}" for i in range(20, 30)},  # P3, used by g3
}
OWN_GROUP = {"GP1": "g2", "GP2": "g1", "GP3": "g3"}

ARGS = {
    "--adata": str(DEMO),
    "--groupby": "group",
    "--reference": "synthetic demo, no genome build",
    "--k": "3",
    "--n-hvg": "200",
    "--n-replicates": "20",
    "--density-threshold": "0.5",
    "--active-threshold": "0.1",
    "--n-top-genes": "10",
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
    """Demo arguments with `changes` ({"k": "1"}) and `drop` applied."""
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
    assert (result["n_cells"], result["n_genes"], result["n_hvg"], result["k"]) == (
        180,
        300,
        200,
        3,
    )
    assert (result["n_components_kept"], result["n_components_total"]) == (60, 60)
    assert result["n_genes_scored"] == 300
    assert result["reference"] == "synthetic demo, no genome build"
    assert {gp: set(genes) for gp, genes in result["top_genes"].items()} == PLANTED
    by_group = pd.read_csv(demo_out / "tables" / "program_by_group.tsv", sep="\t")
    best = by_group.loc[by_group.groupby("program")["auc"].idxmax()]
    assert dict(zip(best["program"], best["group"], strict=True)) == OWN_GROUP
    assert (best["auc"] == 1).all()
    stability = pd.read_csv(demo_out / "tables" / "stability.tsv", sep="\t")
    assert stability["n_components"].tolist() == [20, 20, 20]
    assert (stability["silhouette"] > 0.9).all()


def test_usage_sums_to_one_and_h5ad(demo_out):
    usage = pd.read_csv(demo_out / "tables" / "usage.tsv", sep="\t")
    assert len(usage) == 180 * 3
    assert np.allclose(usage.groupby("cell")["usage"].sum(), 1)
    adata = ad.read_h5ad(demo_out / "data" / "programs.h5ad")
    assert list(adata.obsm["usage"].columns) == ["GP1", "GP2", "GP3"]


def test_tables_sorted(demo_out):
    loadings = pd.read_csv(demo_out / "tables" / "gene_loadings.tsv", sep="\t")
    resorted = loadings.sort_values(
        ["program", "score", "gene"], ascending=[True, False, True]
    )
    assert loadings.index.equals(resorted.index)
    top = pd.read_csv(demo_out / "tables" / "top_genes.tsv", sep="\t")
    assert top["rank"].tolist() == list(range(1, 11)) * 3


def test_output_contract(demo_out):
    paths = contract_paths()
    assert len(paths) == 23
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


def test_rerun_identical_tables(demo_out, tmp_path):
    completed = run("--demo", "--outdir", str(tmp_path / "again"))
    assert completed.returncode == 0, completed.stderr
    for name in ("usage.tsv", "top_genes.tsv", "stability.tsv"):
        first = (demo_out / "tables" / name).read_bytes()
        assert (tmp_path / "again" / "tables" / name).read_bytes() == first, name


def test_commands_sh_reproduces_run(demo_out, tmp_path):
    line = (demo_out / "reproducibility" / "commands.sh").read_text().splitlines()[-1]
    argv = shlex.split(line)
    argv[argv.index("--outdir") + 1] = str(tmp_path / "replay")
    completed = subprocess.run(argv, capture_output=True, text=True, check=False)
    assert completed.returncode == 0, completed.stderr
    for name in ("usage.tsv", "gene_loadings.tsv", "program_by_group.tsv"):
        first = (demo_out / "tables" / name).read_bytes()
        assert (tmp_path / "replay" / "tables" / name).read_bytes() == first, name


def test_missing_reference_fails(tmp_path):
    completed = run_with(tmp_path, drop=["--reference"])
    assert is_error(completed, "missing required arguments: --reference")


@pytest.mark.parametrize(
    ("axis", "message"),
    [("var", "duplicate gene names: G001"), ("obs", "duplicate cell ids: cell001")],
)
def test_duplicate_ids(tmp_path, axis, message):
    def duplicate(adata):
        names = getattr(adata, f"{axis}_names").to_list()
        names[0] = names[1]
        setattr(adata, f"{axis}_names", names)

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, duplicate))
    assert is_error(completed, message)


def test_missing_argument_fails_without_demo(tmp_path):
    completed = run_with(tmp_path, drop=["--density-threshold"])
    assert is_error(completed, "missing required arguments: --density-threshold")


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--k", "4", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --k")


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"k": "1"}, "--k must be an integer >= 2, got 1"),
        ({"n_replicates": "3"}, "--n-replicates must be >= 4, got 3"),
        ({"active_threshold": "1"}, "--active-threshold must be in (0, 1), got 1.0"),
        ({"adata": "/no/such.h5ad"}, "input file not found: /no/such.h5ad"),
        ({"layer": "counts"}, "layer not found in the h5ad: counts"),
        ({"groupby": "leiden"}, "obs column not found: leiden"),
        ({"n_hvg": "400"}, "--n-hvg 400 exceeds the 300 expressed genes"),
        (
            {"density_threshold": "1e-6"},
            "0 of 60 replicate components pass --density-threshold 1e-06; "
            "need at least --k 3",
        ),
    ],
)
def test_input_errors(tmp_path, changes, message):
    assert is_error(run_with(tmp_path, **changes), message)


def test_not_counts(tmp_path):
    def lognorm(adata):
        adata.X = np.log1p(adata.X.toarray())

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, lognorm))
    assert is_error(completed, "X does not look like raw counts: non-integer values")


def test_one_group(tmp_path):
    def one_group(adata):
        adata.obs["group"] = "g1"

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, one_group))
    assert is_error(completed, "--groupby group needs at least 2 groups, found 1")


def test_zero_cells(tmp_path):
    def empty_cell(adata):
        x = adata.X.toarray()
        x[0] = 0
        adata.X = x

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, empty_cell))
    assert is_error(completed, "1 cells have zero counts over the 200 HVGs")


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
