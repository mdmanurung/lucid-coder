"""CLI tests for scripts/sc_normalize.py: demo goldens, output contract, fail strings.

Run from skills/sc-normalize with pytest;
see references/maintainer.md for the exact command.
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "sc_normalize.py"
EXAMPLES = SKILL_DIR / "examples"
TEMPLATES = SKILL_DIR.parents[1] / "templates"

# Golden demo values, computed once with scanpy 1.12.4, muon 0.2.0,
# scikit-learn 1.9.1. Pearson residuals must take the depth signal off PC1.
PEARSON_RHO_MAX = 0.3
LOG1P_RHO_MIN = 0.9
PEARSON_RHO = 0.2022
LOG1P_RHO = 0.9724
# dsb: each planted marker is above this value exactly in its positive group.
DSB_THRESHOLD = 3.5
POSITIVE = {"CD3": "g1", "CD19": "g2", "CD14": "g3"}

ARGS = {
    "--adata": str(EXAMPLES / "demo.h5ad"),
    "--method": "pearson",
    "--layer-out": "normalized",
    "--seed": "0",
    "--reference": "synthetic demo features, no genome build",
    "--context": "synthetic demo, no real tissue",
    "--theta": "100",
    "--clip": "25",
    "--gene-id-col": "gene_ids",
}
DSB_ARGS = {
    "--adata": str(EXAMPLES / "demo_adt.h5ad"),
    "--method": "dsb",
    "--layer-out": "normalized",
    "--seed": "0",
    "--reference": "synthetic demo features, no genome build",
    "--context": "synthetic demo, no real tissue",
    "--empty-droplets": str(EXAMPLES / "demo_empty.h5ad"),
    "--isotype-controls": "IgG1,IgG2a,IgG2b",
}


def run(*argv):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *argv],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )


def run_with(tmp_path, base=ARGS, drop=(), **changes):
    """`base` arguments with `changes` ({"theta": "0"}) and `drop` applied."""
    args = {**base, **{f"--{k.replace('_', '-')}": v for k, v in changes.items()}}
    argv = [x for k, v in args.items() if k not in drop for x in (k, v)]
    return run(*argv, "--outdir", str(tmp_path / "out"))


def contract_paths(method):
    """Output-tree file paths from SKILL.md; rows marked `dsb only` need dsb."""
    text = (SKILL_DIR / "SKILL.md").read_text()
    tree = text.split("## Output structure")[1].split("```")[1]
    paths, stack = [], []
    for line in tree.splitlines():
        match = re.search(r"[├└]── (\S+)", line)
        if not match:
            continue
        depth = match.start() // 4
        stack[depth:] = [match.group(1)]
        if match.group(1).endswith("/") or ("dsb only" in line and method != "dsb"):
            continue
        paths.append("".join(stack))
    return paths


def is_error(completed, message):
    return (
        completed.returncode == 2
        and f"ERROR: {message}" in completed.stderr
        and "Traceback" not in completed.stderr
    )


def demo(tmp_path_factory, method):
    out = tmp_path_factory.mktemp(method) / "out"
    completed = run("--demo", "--method", method, "--outdir", str(out))
    assert completed.returncode == 0, completed.stderr
    return out


@pytest.fixture(scope="module")
def pearson_out(tmp_path_factory):
    return demo(tmp_path_factory, "pearson")


@pytest.fixture(scope="module")
def log1p_out(tmp_path_factory):
    return demo(tmp_path_factory, "log1p")


@pytest.fixture(scope="module")
def dsb_out(tmp_path_factory):
    return demo(tmp_path_factory, "dsb")


def result(outdir):
    return json.loads((outdir / "result.json").read_text())


def test_pearson_removes_depth_from_pc1(pearson_out, log1p_out):
    pearson, log1p = result(pearson_out), result(log1p_out)
    assert pearson["pc1_depth_abs_spearman"] == PEARSON_RHO < PEARSON_RHO_MAX
    assert log1p["pc1_depth_abs_spearman"] == LOG1P_RHO > LOG1P_RHO_MIN
    assert pearson["n_cells"] == 600 and pearson["n_features"] == 300
    assert pearson["n_features_dropped"] == 0
    assert pearson["reference"] == "synthetic demo features, no genome build"
    assert 2 <= len(pearson["chat_summary_lines"]) <= 3
    assert pearson["context"].startswith("synthetic demo")
    assert pearson["chat_summary_lines"][0].startswith("Context: ")
    assert pearson["status_counts"] == {"tested": 300}
    assert pearson["gene_id_col"] == "gene_ids"
    table = pd.read_csv(pearson_out / "tables" / "gene_mean_variance.tsv", sep="\t")
    assert list(table.columns[:2]) == ["feature", "gene_id"]
    assert table["gene_id"].iloc[0] == "ENSG90000000000.1"


def test_layers(pearson_out, log1p_out):
    raw = ad.read_h5ad(EXAMPLES / "demo.h5ad")
    for out in (pearson_out, log1p_out):
        norm = ad.read_h5ad(out / "data" / "normalized.h5ad")
        assert (norm.layers["counts"] != raw.X).nnz == 0
        assert np.array_equal(norm.X, norm.layers["normalized"])
    # log1p of counts scaled to 1e4: expm1 of each row sums to 1e4.
    log1p = ad.read_h5ad(log1p_out / "data" / "normalized.h5ad")
    assert np.allclose(np.expm1(log1p.X).sum(axis=1), 1e4)


def test_dsb_separates_positive_populations(dsb_out):
    norm = ad.read_h5ad(dsb_out / "data" / "normalized.h5ad")
    group = norm.obs["true_group"].to_numpy()
    for protein, positive in POSITIVE.items():
        above = norm[:, protein].X.ravel() > DSB_THRESHOLD
        assert np.array_equal(above, group == positive), protein
    assert (norm[:, "CD16"].X < DSB_THRESHOLD).all()
    assert result(dsb_out)["n_features"] == 7


def test_rerun_is_identical(pearson_out, tmp_path):
    out = tmp_path / "again"
    assert run("--demo", "--outdir", str(out)).returncode == 0
    for name in ["per_cell_summary.tsv", "gene_mean_variance.tsv"]:
        assert (out / "tables" / name).read_bytes() == (
            pearson_out / "tables" / name
        ).read_bytes()


@pytest.mark.parametrize("method", ["pearson", "dsb"])
def test_output_contract(method, pearson_out, dsb_out):
    outdir = {"pearson": pearson_out, "dsb": dsb_out}[method]
    paths = contract_paths(method)
    assert len(paths) == {"pearson": 17, "dsb": 20}[method]
    for rel in paths:
        assert (outdir / rel).stat().st_size > 0, rel
    if method != "dsb":
        assert not (outdir / "figures" / "isotype_vs_background.pdf").exists()
    check = subprocess.run(
        ["sha256sum", "-c", "reproducibility/checksums.sha256"],
        cwd=outdir,
        capture_output=True,
        text=True,
        check=False,
    )
    assert check.returncode == 0, check.stdout + check.stderr
    assert "--seed 0" in (outdir / "reproducibility" / "commands.sh").read_text()
    cells = pd.read_csv(outdir / "tables" / "per_cell_summary.tsv", sep="\t")
    assert cells["cell"].is_monotonic_increasing


@pytest.mark.parametrize("flag", ["--seed", "--context"])
def test_missing_argument_fails_without_demo(tmp_path, flag):
    completed = run_with(tmp_path, drop=[flag])
    assert is_error(completed, f"missing required arguments: {flag}")


def test_missing_reference_fails(tmp_path):
    completed = run_with(tmp_path, drop=["--reference"])
    assert is_error(completed, "missing required arguments: --reference")


@pytest.mark.parametrize(
    ("source", "base", "flag", "axis", "what"),
    [
        ("demo.h5ad", ARGS, "adata", "var", "feature names in"),
        ("demo.h5ad", ARGS, "adata", "obs", "cell ids in"),
        ("demo_empty.h5ad", DSB_ARGS, "empty_droplets", "var", "feature names in"),
    ],
)
def test_duplicate_ids(tmp_path, source, base, flag, axis, what):
    def duplicate(adata):
        names = list(getattr(adata, f"{axis}_names"))
        names[1] = names[0]
        setattr(adata, f"{axis}_names", names)

    path = edited(tmp_path, source, duplicate)
    completed = run_with(tmp_path, base=base, **{flag: path})
    assert is_error(completed, f"duplicate {what} {path}: ")


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--theta", "1", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --theta")


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


@pytest.mark.parametrize(
    ("drop", "changes", "message"),
    [
        ((), {"theta": "0"}, "--theta must be > 0, got 0.0"),
        ((), {"layer_out": "counts"}, "--layer-out must be a layer name other than"),
        (("--clip",), {}, "--method pearson needs: --clip"),
        (("--gene-id-col",), {}, "--method pearson needs: --gene-id-col"),
        ((), {"gene_id_col": "ensembl"}, "var column not found in --adata: ensembl"),
        ((), {"target_sum": "1e4"}, "--method pearson does not use: --target-sum"),
        ((), {"adata": "/no/such.h5ad"}, "input file not found: /no/such.h5ad"),
    ],
)
def test_argument_errors(tmp_path, drop, changes, message):
    assert is_error(run_with(tmp_path, drop=drop, **changes), message)


def edited(tmp_path, source, edit, name="bad.h5ad"):
    adata = ad.read_h5ad(EXAMPLES / source)
    edit(adata)
    adata.write_h5ad(tmp_path / name)
    return str(tmp_path / name)


def test_non_integer_counts(tmp_path):
    def halve(adata):
        adata.X = adata.X * 0.5

    path = edited(tmp_path, "demo.h5ad", halve)
    assert is_error(run_with(tmp_path, adata=path), "counts in X of")


def test_zero_count_cell(tmp_path):
    def zero(adata):
        adata.X = adata.X.toarray()
        adata.X[0] = 0

    path = edited(tmp_path, "demo.h5ad", zero)
    assert is_error(run_with(tmp_path, adata=path), "1 cells have zero total counts")


def test_all_zero_gene_dropped(tmp_path):
    def zero(adata):
        adata.X = adata.X.toarray()
        adata.X[:, 0] = 0

    completed = run_with(tmp_path, adata=edited(tmp_path, "demo.h5ad", zero))
    assert completed.returncode == 0, completed.stderr
    res = result(tmp_path / "out")
    assert res["n_features_dropped"] == 1
    assert res["status_counts"] == {"tested": 299, "not_detected": 1}
    table = pd.read_csv(tmp_path / "out/tables/gene_mean_variance.tsv", sep="\t")
    row = table.set_index("feature").loc["G000"]
    assert row["status"] == "not_detected" and row["gene_id"] == "ENSG90000000000.1"
    assert row[["mean_counts", "var_normalised"]].isna().all()


def test_ambiguous_gene_ids(tmp_path):
    def duplicate(adata):
        adata.var.loc["G001", "gene_ids"] = adata.var.loc["G000", "gene_ids"]

    completed = run_with(tmp_path, adata=edited(tmp_path, "demo.h5ad", duplicate))
    assert is_error(completed, "duplicate gene ids: ENSG90000000000.1")


def test_no_gene_ids_writes_na(tmp_path):
    completed = run_with(tmp_path, gene_id_col="none")
    assert completed.returncode == 0, completed.stderr
    assert result(tmp_path / "out")["gene_id_col"] == "none"
    table = pd.read_csv(tmp_path / "out/tables/gene_mean_variance.tsv", sep="\t")
    assert table["gene_id"].isna().all()


def test_isotype_missing(tmp_path):
    completed = run_with(tmp_path, base=DSB_ARGS, isotype_controls="IgG1,IgG3")
    assert is_error(completed, "isotype controls not in the ADT features: IgG3")


def test_panel_mismatch(tmp_path):
    path = edited(tmp_path, "demo_empty.h5ad", lambda a: None)
    reordered = ad.read_h5ad(path)[:, ::-1].copy()
    reordered.write_h5ad(path)
    completed = run_with(tmp_path, base=DSB_ARGS, empty_droplets=path)
    assert is_error(completed, "--adata and --empty-droplets must have the same")


def test_barcode_overlap(tmp_path):
    def rename(adata):
        adata.obs_names = ["cell0000", *adata.obs_names[1:]]

    path = edited(tmp_path, "demo_empty.h5ad", rename)
    completed = run_with(tmp_path, base=DSB_ARGS, empty_droplets=path)
    assert is_error(completed, "1 barcodes are in both --adata and --empty-droplets")


def test_flat_empty(tmp_path):
    def flatten(adata):
        adata.X = adata.X.toarray()
        adata.X[:, 0] = 5

    path = edited(tmp_path, "demo_empty.h5ad", flatten)
    completed = run_with(tmp_path, base=DSB_ARGS, empty_droplets=path)
    assert is_error(completed, "features with zero variance in the empty droplets: CD3")


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
