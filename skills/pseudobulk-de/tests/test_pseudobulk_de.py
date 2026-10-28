"""CLI tests for scripts/pseudobulk_de.py: demo goldens, output contract, fail strings.

Run from skills/pseudobulk-de with pytest;
see references/maintainer.md for the R engine environment.
"""

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "pseudobulk_de.py"
DEMO = SKILL_DIR / "examples" / "demo.h5ad"
TEMPLATES = SKILL_DIR.parents[1] / "templates"

# Golden demo values, computed once with pydeseq2 0.5.4 (python engine).
N_SAMPLES, N_GENES_TESTED, N_SIGNIFICANT = 12, 250, 28
# Golden dreamlet demo values (dreamlet 1.8.0, ~ condition + (1|patient)):
# 29 significant, 29 of them among the 30 true DE genes.
DREAMLET_SIGNIFICANT, DREAMLET_TRUE_HITS = 29, 29
RANDOM_DESIGN = "~ condition + (1|patient)"

ARGS = {
    "--adata": str(DEMO),
    "--gene-id-col": "gene_ids",
    "--cell-type-col": "cell_type",
    "--cell-type": "T",
    "--sample-col": "sample",
    "--donor-col": "patient",
    "--design": "~ patient + condition",
    "--contrast": "condition,stim,ctrl",
    "--min-cells": "10",
    "--min-frac": "0.1",
    "--padj-max": "0.05",
    "--lfc-min": "1",
    "--seed": "0",
    "--engine": "python",
    "--reference": "synthetic demo, no genome",
    "--context": "synthetic demo: 7 patients, ctrl vs stim T cells, no real tissue",
}


def run(*argv, env=None):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *argv],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1", **(env or {})},
    )


def run_with(tmp_path, env=None, drop=(), **changes):
    """Demo arguments with `changes` ({"--min-frac": "2"}) and `drop` applied."""
    args = {**ARGS, **{f"--{k.replace('_', '-')}": v for k, v in changes.items()}}
    argv = [x for k, v in args.items() if k not in drop for x in (k, v)]
    return run(*argv, "--outdir", str(tmp_path / "out"), env=env)


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


@pytest.fixture(scope="module")
def demo_out(tmp_path_factory):
    out = tmp_path_factory.mktemp("demo") / "out"
    completed = run("--demo", "--outdir", str(out))
    assert completed.returncode == 0, completed.stderr
    return out


def test_demo_golden_values(demo_out):
    result = json.loads((demo_out / "result.json").read_text())
    assert result["n_samples"] == N_SAMPLES
    assert result["n_genes_tested"] == N_GENES_TESTED
    assert result["n_significant"] == N_SIGNIFICANT
    assert result["dropped_samples"] == ["P7_ctrl"]
    assert result["dropped_by_level"] == {"stim": [], "ctrl": ["P7_ctrl"]}
    assert result["n_donors"] == 6
    assert result["reference"] == "synthetic demo, no genome"
    assert result["seed"] == 0 and result["engine"] == "python"
    assert 2 <= len(result["chat_summary_lines"]) <= 3
    truth = set(ad.read_h5ad(DEMO).var.query("true_de").index)
    de = pd.read_csv(demo_out / "tables" / "de_results.tsv", sep="\t")
    assert set(de.loc[de["significant"], "gene_symbol"]) <= truth


def test_biology_core_fields(demo_out):
    result = json.loads((demo_out / "result.json").read_text())
    assert result["context"] == ARGS["--context"]
    assert result["chat_summary_lines"][0].startswith("Context: ")
    assert result["gene_id_col"] == "gene_ids"
    counts = result["status_counts"]
    assert counts["tested"] == N_GENES_TESTED and sum(counts.values()) == 300
    assert set(counts) <= {"tested", "not_tested", "not_detected"}
    de = pd.read_csv(demo_out / "tables" / "de_results.tsv", sep="\t")
    assert list(de.columns[:2]) == ["gene_symbol", "gene_id"] and len(de) == 300
    untested = de[de["status"] != "tested"]
    assert untested["pvalue"].isna().all() and not untested["significant"].any()
    assert {f"G{i:03d}" for i in range(260, 300)} <= set(untested["gene_symbol"])
    g000 = de.loc[de["gene_symbol"] == "G000", "gene_id"]
    assert g000.tolist() == ["ENSG90000000000.1"]
    counts = pd.read_csv(demo_out / "tables" / "pseudobulk_counts.tsv", sep="\t")
    assert list(counts.columns[:2]) == ["gene_id", "gene_symbol"]
    assert len(counts) == N_GENES_TESTED


def test_gene_id_col_none_writes_na(tmp_path):
    completed = run_with(tmp_path, gene_id_col="none")
    assert completed.returncode == 0, completed.stderr
    result = json.loads((tmp_path / "out" / "result.json").read_text())
    assert result["gene_id_col"] == "none"
    de = pd.read_csv(tmp_path / "out" / "tables" / "de_results.tsv", sep="\t")
    assert de["gene_id"].isna().all()


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


def test_cell_type_col_none_tests_all_cells(tmp_path, demo_out):
    """One-population input: `--cell-type-col none` equals subsetting first."""
    t_cells = ad.read_h5ad(DEMO)
    t_cells = t_cells[t_cells.obs["cell_type"] == "T"].copy()
    del t_cells.obs["cell_type"]
    t_cells.write_h5ad(tmp_path / "t_only.h5ad")
    done = run_with(tmp_path, adata=str(tmp_path / "t_only.h5ad"), cell_type_col="none")
    assert done.returncode == 0, done.stderr
    table = "tables/de_results.tsv"
    assert (tmp_path / "out" / table).read_bytes() == (demo_out / table).read_bytes()


@pytest.mark.parametrize("arg", ["--seed", "--context", "--gene-id-col"])
def test_missing_argument_fails_without_demo(tmp_path, arg):
    assert is_error(
        run_with(tmp_path, drop=[arg]), f"missing required arguments: {arg}"
    )


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--seed", "1", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --seed")


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"min_frac": "2"}, "--min-frac must be in [0, 1], got 2.0"),
        ({"design": "condition"}, "--design must be '~ a + b'"),
        (
            {"contrast": "condition,stim"},
            "--contrast must be factor,numerator,denominator",
        ),
        ({"adata": "/no/such.h5ad"}, "input file not found: /no/such.h5ad"),
        ({"layer": "counts"}, "layer not found in the h5ad: counts"),
        ({"gene_id_col": "gene_id"}, "gene id column not found in var: gene_id"),
        ({"design": "~ batch + condition"}, "obs column not found: batch"),
        ({"donor_col": "donor"}, "obs column not found: donor"),
        (
            {"design": "~ condition"},
            "6 donor(s) in 'patient' have more than one sample, but the design has "
            "no 'patient' term",
        ),
        ({"cell_type": "NK"}, "cell type 'NK' not found in obs column 'cell_type'"),
        (
            {"design": "~ patient"},
            "contrast factor 'condition' is not a term of --design",
        ),
        ({"contrast": "condition,stim,mock"}, "contrast level 'mock' not found"),
        (
            {"sample_col": "patient", "design": "~ condition"},
            "obs column 'condition' varies within a sample",
        ),
        ({"min_cells": "41"}, "fewer than 2 samples for condition='stim'"),
        ({"min_frac": "1"}, "zero genes pass --min-frac 1.0"),
    ],
)
def test_input_errors(tmp_path, changes, message):
    assert is_error(run_with(tmp_path, **changes), message)


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


def test_non_integer_counts(tmp_path):
    adata = ad.read_h5ad(DEMO)
    adata.X = adata.X * 0.5
    adata.write_h5ad(tmp_path / "bad.h5ad")
    completed = run_with(tmp_path, adata=str(tmp_path / "bad.h5ad"))
    assert is_error(completed, "counts in X must be non-negative integers")


def test_missing_reference(tmp_path):
    completed = run_with(tmp_path, drop=["--reference"])
    assert is_error(completed, "missing required arguments: --reference")


@pytest.mark.parametrize(
    ("axis", "message"),
    [
        (
            "var",
            "1 gene symbol(s) map to several gene ids: "
            "G001 -> ENSG90000000001.1, ENSG90000000002.1",
        ),
        ("obs", "duplicate cell ids: cell0001"),
    ],
)
def test_duplicate_ids(tmp_path, axis, message):
    adata = ad.read_h5ad(DEMO)
    names = getattr(adata, f"{axis}_names").tolist()
    names[2] = names[1]
    setattr(adata, f"{axis}_names", names)
    adata.write_h5ad(tmp_path / "dup.h5ad")
    assert is_error(run_with(tmp_path, adata=str(tmp_path / "dup.h5ad")), message)


def test_rank_deficient_design(tmp_path):
    # batch is a copy of condition, so the two terms are aliased.
    adata = ad.read_h5ad(DEMO)
    adata.obs["batch"] = adata.obs["condition"].map({"ctrl": "b1", "stim": "b2"})
    adata.write_h5ad(tmp_path / "conf.h5ad")
    completed = run_with(
        tmp_path,
        adata=str(tmp_path / "conf.h5ad"),
        design="~ patient + batch + condition",
    )
    assert is_error(completed, "design is rank-deficient (7 < 8 columns)")


def test_no_rscript(tmp_path):
    env = {"RSCRIPT": "", "PATH": str(tmp_path)}
    completed = run_with(tmp_path, env=env, engine="r")
    assert is_error(completed, "Rscript not found; set RSCRIPT or put Rscript on PATH")


def test_r_failure(tmp_path):
    completed = run_with(tmp_path, env={"RSCRIPT": shutil.which("false")}, engine="r")
    assert is_error(completed, "deseq2.R failed with exit 1")


def test_random_effect_refuses_python_and_r(tmp_path):
    for engine in ("python", "r"):
        completed = run_with(tmp_path, design=RANDOM_DESIGN, engine=engine)
        assert is_error(
            completed,
            "--design has the random-effect term (1|patient): repeated measures "
            f"need --engine dreamlet, not --engine {engine}",
        )


def test_fixed_design_still_needs_engine(tmp_path):
    completed = run_with(tmp_path, drop=["--engine"])
    assert is_error(completed, "missing required arguments: --engine")


def test_random_effect_selects_dreamlet(tmp_path):
    # RSCRIPT=false stops the run at the R call, after the engine is chosen.
    completed = run_with(
        tmp_path,
        env={"RSCRIPT": shutil.which("false")},
        drop=["--engine"],
        design=RANDOM_DESIGN,
    )
    assert is_error(completed, "dreamlet.R failed with exit 1")


def test_random_contrast_factor(tmp_path):
    completed = run_with(
        tmp_path, design="~ patient + (1|condition)", engine="dreamlet"
    )
    assert is_error(completed, "contrast factor 'condition' is not a term of --design")


RSCRIPT = os.environ.get("RSCRIPT") or shutil.which("Rscript")


@pytest.mark.skipif(not RSCRIPT, reason="Rscript not found")
def test_r_engine_matches_python(tmp_path, demo_out):
    completed = run("--demo", "--engine", "r", "--outdir", str(tmp_path / "r"))
    assert completed.returncode == 0, completed.stderr
    result = json.loads((tmp_path / "r" / "result.json").read_text())
    assert result["n_significant"] == N_SIGNIFICANT
    truth = set(ad.read_h5ad(DEMO).var.query("true_de").index)
    de = pd.read_csv(tmp_path / "r" / "tables" / "de_results.tsv", sep="\t")
    assert set(de.loc[de["significant"], "gene_symbol"]) <= truth
    env = (tmp_path / "r" / "reproducibility" / "environment.txt").read_text()
    assert "DESeq2 " in env


@pytest.mark.skipif(not RSCRIPT, reason="Rscript not found")
def test_dreamlet_demo(tmp_path):
    out = tmp_path / "dl"
    completed = run("--demo", "--engine", "dreamlet", "--outdir", str(out))
    assert completed.returncode == 0, completed.stderr
    result = json.loads((out / "result.json").read_text())
    assert result["engine"] == "dreamlet"
    assert result["params"]["design"] == RANDOM_DESIGN
    assert result["n_samples"] == N_SAMPLES
    assert result["n_genes_tested"] == N_GENES_TESTED
    assert result["n_significant"] == DREAMLET_SIGNIFICANT
    truth = set(ad.read_h5ad(DEMO).var.query("true_de").index)
    de = pd.read_csv(out / "tables" / "de_results.tsv", sep="\t")
    assert list(de.columns) == [
        "gene_symbol",
        "gene_id",
        "base_mean",
        "log2_fc",
        "lfc_se",
        "stat",
        "pvalue",
        "padj",
        "status",
        "significant",
    ]
    hits = set(de.loc[de["significant"], "gene_symbol"])
    assert len(hits & truth) == DREAMLET_TRUE_HITS
    assert (de["lfc_se"] - de["log2_fc"] / de["stat"]).abs().max() < 1e-9
    for rel in contract_paths():
        assert (out / rel).stat().st_size > 0, rel
    check = subprocess.run(
        ["sha256sum", "-c", "reproducibility/checksums.sha256"],
        cwd=out,
        capture_output=True,
        text=True,
        check=False,
    )
    assert check.returncode == 0, check.stdout + check.stderr
    env = (out / "reproducibility" / "environment.txt").read_text()
    assert "dreamlet " in env and "variancePartition " in env
    assert "--engine dreamlet" in (out / "reproducibility" / "commands.sh").read_text()


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["plots.py", "theme.py", "palettes.py", "skillkit.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
