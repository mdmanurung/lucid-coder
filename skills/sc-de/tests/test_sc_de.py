"""CLI tests for scripts/sc_de.py: demo goldens, output contract, fail strings.

Run from skills/sc-de with pytest;
see references/maintainer.md for the exact command and the R
environment. MAST goldens were set with R 4.5.2, MAST 1.36.0, lme4 2.1.0.
"""

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "sc_de.py"
H5AD = SKILL_DIR / "examples" / "demo.h5ad"
NPX = SKILL_DIR / "examples" / "demo_npx.tsv"
TEMPLATES = SKILL_DIR.parents[1] / "templates"
HAS_R = bool(os.environ.get("RSCRIPT") or shutil.which("Rscript"))
needs_r = pytest.mark.skipif(not HAS_R, reason="no Rscript; set RSCRIPT")

# Planted truth from examples/make_demo.py.
MAST_TRUTH = {f"G{i:03d}" for i in range(15)}
STATUS_TRUTH = {f"P{i:02d}" for i in range(1, 9)}
AGE_TRUTH = {f"P{i:02d}" for i in range(9, 13)}

LM = {
    "--mode": "lm",
    "--input": str(NPX),
    "--formula": "~ status + age + sex",
    "--term": "status,disease,healthy",
    "--padj-max": "0.05",
    "--lfc-min": "0.5",
    "--seed": "0",
    "--reference": "GRCh38, GENCODE 44",
    "--context": "synthetic demo: 40 samples, disease vs healthy NPX, no real cohort",
}
MAST = {
    "--mode": "mast",
    "--adata": str(H5AD),
    "--gene-id-col": "gene_ids",
    "--cell-type-col": "cell_type",
    "--cell-type": "T",
    "--donor-col": "donor",
    "--formula": "~ condition",
    "--contrast-term": "condition,stim,ctrl",
    "--n-cells-per-donor": "50",
    "--min-frac": "0.1",
    "--cdr": "yes",
    "--log-base": "e",
    "--padj-max": "0.05",
    "--lfc-min": "0.5",
    "--seed": "0",
    "--reference": "GRCh38, GENCODE 44",
    "--context": "synthetic demo: 6 donors, ctrl vs stim T cells, no real tissue",
}


def run(*argv, env=None):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *argv],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1", **(env or {})},
    )


def run_with(tmp_path, base, drop=(), env=None, **changes):
    """`base` arguments with `changes` ({"padj_max": "2"}) and `drop` applied."""
    args = {**base, **{f"--{k.replace('_', '-')}": v for k, v in changes.items()}}
    argv = [x for k, v in args.items() if k not in drop for x in (k, v)]
    return run(*argv, "--outdir", str(tmp_path / "out"), env=env)


def significant(outdir):
    table = pd.read_csv(outdir / "tables" / "de_results.tsv", sep="\t")
    key = "gene_symbol" if "gene_symbol" in table.columns else "feature"
    return set(table.loc[table["significant"], key])


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


def demo(tmp_path_factory, *mode):
    out = tmp_path_factory.mktemp("demo") / "out"
    completed = run("--demo", *mode, "--outdir", str(out))
    assert completed.returncode == 0, completed.stderr
    return out


@pytest.fixture(scope="module")
def mast_out(tmp_path_factory):
    return demo(tmp_path_factory)


@pytest.fixture(scope="module")
def lm_out(tmp_path_factory):
    return demo(tmp_path_factory, "--mode", "lm")


@needs_r
def test_mast_demo_golden(mast_out):
    result = json.loads((mast_out / "result.json").read_text())
    assert (result["n_cells"], result["n_donors"]) == (600, 6)
    assert (result["n_tested"], result["n_significant"]) == (200, 16)
    # All 15 planted genes plus G141: a cell-level false positive from donor noise.
    assert significant(mast_out) == MAST_TRUTH | {"G141"}
    assert result["formula_fitted"] == "~ condition + cngeneson"
    assert "overstate donor-level evidence" in result["chat_summary_lines"][1]
    assert result["chat_summary_lines"][0].startswith("Context: ")
    assert result["context"] == MAST["--context"]
    assert result["status_counts"] == {"tested": 200}
    assert result["gene_id_col"] == "gene_ids"
    table = pd.read_csv(mast_out / "tables" / "de_results.tsv", sep="\t")
    assert list(table.columns[:2]) == ["gene_symbol", "gene_id"]
    assert set(table.loc[table["gene_symbol"] == "G000", "gene_id"]) == {
        "ENSG90000000000.1"
    }
    assert result["reference"] == "synthetic demo, no genome"
    assert result["dropped_by_level"] == {}


@needs_r
def test_lm_demo_golden(lm_out):
    result = json.loads((lm_out / "result.json").read_text())
    assert (result["n_samples"], result["n_tested"]) == (40, 50)
    assert significant(lm_out) == STATUS_TRUTH
    assert result["dropped_by_level"] == {}
    assert result["status_counts"] == {"tested": 50}
    assert "gene_id_col" not in result
    table = pd.read_csv(lm_out / "tables" / "de_results.tsv", sep="\t")
    assert set(table["feature"].head(8)) == STATUS_TRUTH
    assert "gene_id" not in table.columns


@needs_r
def test_lm_numeric_term(tmp_path):
    completed = run_with(tmp_path, LM, term="age", lfc_min="0")
    assert completed.returncode == 0, completed.stderr
    assert significant(tmp_path / "out") == AGE_TRUTH


@needs_r
@pytest.mark.parametrize("fixture", ["mast_out", "lm_out"])
def test_output_contract(fixture, request):
    outdir = request.getfixturevalue(fixture)
    paths = contract_paths()
    assert len(paths) == 12
    for rel in paths:
        assert (outdir / rel).stat().st_size > 0, rel
    check = subprocess.run(
        ["sha256sum", "-c", "reproducibility/checksums.sha256"],
        cwd=outdir,
        capture_output=True,
        text=True,
        check=False,
    )
    assert check.returncode == 0, check.stdout + check.stderr
    assert "--seed 0" in (outdir / "reproducibility" / "commands.sh").read_text()


@needs_r
def test_mast_rerun_is_identical(mast_out, tmp_path):
    completed = run("--demo", "--outdir", str(tmp_path / "again"))
    assert completed.returncode == 0, completed.stderr
    check = subprocess.run(
        ["sha256sum", "-c", str(mast_out / "reproducibility" / "checksums.sha256")],
        cwd=tmp_path / "again",
        capture_output=True,
        text=True,
        check=False,
    )
    assert check.returncode == 0, check.stdout


def test_missing_argument_fails_without_demo(tmp_path):
    assert is_error(
        run_with(tmp_path, LM, drop=["--seed"]), "missing required arguments: --seed"
    )


@pytest.mark.parametrize("arg", ["--reference", "--context"])
def test_required_arguments(tmp_path, arg):
    completed = run_with(tmp_path, LM, drop=[arg])
    assert is_error(completed, f"missing required arguments: {arg}")


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--padj-max", "0.1", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --padj-max")


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path, LM), "outdir already exists")


@pytest.fixture(scope="module")
def bad_inputs(tmp_path_factory):
    """A raw-count h5ad and an NPX table with a repeated row."""
    import anndata as ad
    import numpy as np

    tmp = tmp_path_factory.mktemp("bad")
    adata = ad.read_h5ad(H5AD)
    adata.X = np.round(np.expm1(adata.X.toarray()))
    adata.write_h5ad(tmp / "counts.h5ad")
    npx = pd.read_csv(NPX, sep="\t", dtype=str)
    pd.concat([npx, npx.head(1)]).to_csv(tmp / "dup.tsv", sep="\t", index=False)
    # One disease sample loses its age: its rows drop under the "disease" level.
    lost = npx["sample"] == npx.loc[npx["status"] == "disease", "sample"].iloc[0]
    npx.assign(age=npx["age"].where(~lost, "")).to_csv(
        tmp / "missing.tsv", sep="\t", index=False
    )
    adata = ad.read_h5ad(H5AD)
    adata.var_names = ["G000", *adata.var_names[1:-1], "G000"]
    adata.write_h5ad(tmp / "dup_genes.h5ad")
    # Low detection: the last 20 genes are zeroed in every cell (never detected).
    adata = ad.read_h5ad(H5AD)
    x = adata.X.toarray()
    x[:, -20:] = 0
    adata.X = x
    adata.write_h5ad(tmp / "silent.h5ad")
    adata = ad.read_h5ad(H5AD)
    adata.obs["batch"] = adata.obs["condition"]
    adata.write_h5ad(tmp / "confounded.h5ad")
    return tmp


@pytest.mark.parametrize(
    ("base", "changes", "drop", "message"),
    [
        (LM, {"padj_max": "1.5"}, (), "--padj-max must be in (0, 1), got 1.5"),
        (LM, {}, ("--term",), "--mode lm needs --term"),
        (LM, {"cdr": "yes"}, (), "--mode lm does not use --cdr"),
        (LM, {"formula": "status"}, (), "--formula must start with ~"),
        (LM, {"term": "status,disease"}, (), "--term must be column or"),
        (LM, {"term": "status,disease,disease"}, (), "--term must be column or"),
        (
            LM,
            {"formula": "~ status + age", "term": "sex,M,F"},
            (),
            "--term column 'sex' is not in --formula",
        ),
        (LM, {"input": "/no/such.tsv"}, (), "input file not found: /no/such.tsv"),
        (LM, {"formula": "~ status + bmi"}, (), "input column not found: bmi"),
        (
            LM,
            {"term": "status,disease,control"},
            (),
            "level 'control' not found in column 'status'",
        ),
        (LM, {"term": "sex"}, (), "column 'sex' must be numeric"),
        (MAST, {"layer": "counts"}, (), "layer not found in h5ad: counts"),
        (MAST, {}, ("--gene-id-col",), "--mode mast needs --gene-id-col"),
        (LM, {"gene_id_col": "gene_ids"}, (), "--mode lm does not use --gene-id-col"),
        (
            MAST,
            {"gene_id_col": "gene_id"},
            (),
            "gene id column not found in var: gene_id",
        ),
        (MAST, {"cell_type": "B"}, (), "cell type 'B' not found in column 'cell_type'"),
        (
            MAST,
            {"n_cells_per_donor": "0"},
            (),
            "--n-cells-per-donor must be 'all' or an integer >= 1",
        ),
        (
            MAST,
            {"min_frac": "0.95"},
            (),
            "no gene is detected in at least 0.95 of the kept cells",
        ),
    ],
)
def test_input_errors(tmp_path, base, changes, drop, message):
    assert is_error(run_with(tmp_path, base, drop=drop, **changes), message)


def test_raw_counts_refused(tmp_path, bad_inputs):
    completed = run_with(tmp_path, MAST, adata=str(bad_inputs / "counts.h5ad"))
    assert is_error(completed, "X must hold log-normalised expression")


def test_duplicate_rows_refused(tmp_path, bad_inputs):
    completed = run_with(tmp_path, LM, input=str(bad_inputs / "dup.tsv"))
    assert is_error(completed, "2 rows repeat a (feature, sample) combination")


def test_ambiguous_gene_ids_refused(tmp_path, bad_inputs):
    completed = run_with(tmp_path, MAST, adata=str(bad_inputs / "dup_genes.h5ad"))
    assert is_error(
        completed,
        "1 gene symbol(s) map to several gene ids: "
        "G000 -> ENSG90000000000.1, ENSG90000000199.1",
    )


def test_duplicate_symbols_without_ids(tmp_path, bad_inputs):
    path = str(bad_inputs / "dup_genes.h5ad")
    completed = run_with(tmp_path, MAST, adata=path, gene_id_col="none")
    assert is_error(completed, "duplicate gene symbols: G000")


@needs_r
def test_gene_id_none_and_untested_status(tmp_path, bad_inputs):
    completed = run_with(
        tmp_path,
        MAST,
        adata=str(bad_inputs / "silent.h5ad"),
        gene_id_col="none",
        n_cells_per_donor="10",
    )
    assert completed.returncode == 0, completed.stderr
    result = json.loads((tmp_path / "out" / "result.json").read_text())
    assert result["gene_id_col"] == "none"
    assert result["status_counts"]["not_detected"] == 20
    table = pd.read_csv(tmp_path / "out" / "tables" / "de_results.tsv", sep="\t")
    assert table["gene_id"].isna().all() and len(table) == 200
    silent = table[table["status"] == "not_detected"]
    assert set(silent["gene_symbol"]) == {f"G{i:03d}" for i in range(180, 200)}
    assert silent["pvalue"].isna().all() and not silent["significant"].any()


def test_rank_deficient_refused(tmp_path, bad_inputs):
    completed = run_with(
        tmp_path,
        MAST,
        adata=str(bad_inputs / "confounded.h5ad"),
        formula="~ condition + batch",
    )
    assert is_error(completed, "design is rank-deficient")


@needs_r
def test_lm_drops_per_level(tmp_path, bad_inputs):
    completed = run_with(tmp_path, LM, input=str(bad_inputs / "missing.tsv"))
    assert completed.returncode == 0, completed.stderr
    result = json.loads((tmp_path / "out" / "result.json").read_text())
    assert result["dropped_by_level"] == {"disease": 50}


def test_no_rscript(tmp_path):
    completed = run_with(
        tmp_path, LM, env={"RSCRIPT": "/no/Rscript", "PATH": "/usr/bin:/bin"}
    )
    assert is_error(completed, "Rscript not found; set RSCRIPT or put Rscript on PATH")


@needs_r
def test_r_failure(tmp_path):
    completed = run_with(tmp_path, LM, formula="~ status + log(age")
    assert is_error(completed, "lm.R failed with exit 1")


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
