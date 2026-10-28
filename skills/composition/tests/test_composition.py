"""CLI tests for scripts/composition.py: demo goldens, output contract, fail strings.

Run from skills/composition with pytest;
see references/maintainer.md for the exact command.
Three tests call R (about 3 minutes each here): the demo, the counts-TSV rerun and
bad_coef. Every other failure exits in Python before R starts.
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
SCRIPT = SKILL_DIR / "scripts" / "composition.py"
DEMO = SKILL_DIR / "examples" / "demo.h5ad"
TEMPLATES = SKILL_DIR.parents[1] / "templates"

ARGS = {
    "--adata": str(DEMO),
    "--sample-col": "sample",
    "--label-col": "cell_type",
    "--donor-col": "donor",
    "--context": "synthetic demo",
    "--formula": "~ condition + (1|donor)",
    "--coef": "conditionstim",
    "--group-col": "condition",
    "--tau": "1",
    "--min-count": "50",
    "--padj-max": "0.05",
    "--seed": "0",
}


def run(*argv, env=None):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *argv],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1", **(env or {})},
    )


def run_with(tmp_path, drop=(), env=None, **changes):
    """Demo arguments with `changes` ({"tau": "0"}) and `drop` applied."""
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


@pytest.fixture(scope="module")
def tsv_inputs(tmp_path_factory):
    """The demo as --counts / --meta files: all 7 labels, one meta row per sample."""
    tmp = tmp_path_factory.mktemp("tsv")
    obs = ad.read_h5ad(DEMO).obs.astype(str)
    pd.crosstab(obs["sample"], obs["cell_type"]).to_csv(tmp / "counts.tsv", sep="\t")
    meta = obs.drop_duplicates("sample")[["sample", "donor", "condition"]]
    meta.to_csv(tmp / "meta.tsv", sep="\t", index=False)
    return tmp


def tsv_args(tsv_inputs, counts="counts.tsv", meta="meta.tsv"):
    return {"counts": str(tsv_inputs / counts), "meta": str(tsv_inputs / meta)}


H5AD_FLAGS = ("--adata", "--sample-col", "--label-col")


def test_demo_golden_values(demo_out):
    result = json.loads((demo_out / "result.json").read_text())
    assert result["n_samples"] == 16
    assert result["n_donors"] == 8
    pdc = ad.read_h5ad(DEMO).obs.query("cell_type == 'pDC'")
    expected = pdc["condition"].astype(str).value_counts().to_dict()
    assert result["dropped_cells_by_level"] == expected
    assert result["n_labels_tested"] == 6
    assert result["dropped_labels"] == ["pDC"]
    assert result["significant"] == {"Mono": "up", "B": "down"}
    assert result["paired_by"] == "donor"
    assert result["seed"] == 0 and result["engine"] == "crumblr"
    da = pd.read_csv(demo_out / "tables" / "da_results.tsv", sep="\t")
    assert da["label"].tolist()[:2] == ["Mono", "B"]
    assert da.columns.tolist() == [
        *["label", "logFC", "se", "t", "ci_low", "ci_high", "pvalue", "padj"],
        *["significant", "status"],
    ]
    # The label dropped by --min-count stays, with NA statistics.
    pdc_row = da.set_index("label").loc["pDC"]
    assert pdc_row["status"] == "not_tested" and pd.isna(pdc_row["padj"])
    assert result["status_counts"] == {"tested": 6, "not_tested": 1}
    assert result["context"].startswith("synthetic demo")
    assert result["chat_summary_lines"][0].startswith("Context: ")
    props = pd.read_csv(demo_out / "tables" / "proportions.tsv", sep="\t", index_col=0)
    assert props.sum(axis=1).round(5).eq(1).all()


def test_counts_tsv_matches_demo(demo_out, tsv_inputs, tmp_path):
    completed = run_with(tmp_path, drop=H5AD_FLAGS, **tsv_args(tsv_inputs))
    assert completed.returncode == 0, completed.stderr
    for table in ("counts.tsv", "da_results.tsv", "clr.tsv"):
        tsv = (tmp_path / "out" / "tables" / table).read_text()
        assert tsv == (demo_out / "tables" / table).read_text(), table


def test_output_contract(demo_out):
    paths = contract_paths()
    assert len(paths) == 33
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
    assert "crumblr 1.2.0" in (demo_out / "reproducibility/environment.txt").read_text()


def test_missing_argument_fails_without_demo(tmp_path):
    assert is_error(
        run_with(tmp_path, drop=["--tau"]), "missing required arguments: --tau"
    )


def test_missing_context_fails(tmp_path):
    assert is_error(
        run_with(tmp_path, drop=["--context"]), "missing required arguments: --context"
    )


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--tau", "2", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --tau")


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


@pytest.mark.parametrize(
    ("drop", "changes", "message"),
    [
        (
            (),
            {"counts": "x.tsv", "meta": "y.tsv"},
            "give --adata with --sample-col and --label-col, or --counts with --meta",
        ),
        (("--label-col",), {}, "give --adata with --sample-col and --label-col"),
        ((), {"tau": "0"}, "--tau must be > 0, got 0.0"),
        ((), {"min_count": "-1"}, "--min-count must be an integer >= 0, got -1"),
        ((), {"padj_max": "0"}, "--padj-max must be in (0, 1], got 0.0"),
        ((), {"seed": "-1"}, "--seed must be an integer >= 0, got -1"),
        (
            (),
            {"formula": "~ condition * donor"},
            "--formula must be '~ a + b' with sample-level column names",
        ),
        ((), {"adata": "/no/such.h5ad"}, "input file not found: /no/such.h5ad"),
        ((), {"group_col": "batch"}, "column not found: batch"),
        ((), {"label_col": "celltype"}, "column not found: celltype"),
        ((), {"donor_col": "subject"}, "column not found: subject"),
        (
            (),
            {"formula": "~ condition"},
            "8 donor(s) in 'donor' have more than one sample, but the design has "
            "no 'donor' term",
        ),
        (
            (),
            {"sample_col": "donor"},
            "column 'condition' varies within a sample",
        ),
        (
            (),
            {"min_count": "100000"},
            "0 label(s) keep >= --min-count 100000 cells; crumblr needs at least 2",
        ),
    ],
)
def test_input_errors(tmp_path, drop, changes, message):
    assert is_error(run_with(tmp_path, drop=drop, **changes), message)


def test_bad_counts(tmp_path, tsv_inputs):
    counts = pd.read_csv(tsv_inputs / "counts.tsv", sep="\t", index_col=0)
    counts.iloc[0, 0] = -1
    counts.to_csv(tmp_path / "bad.tsv", sep="\t")
    args = tsv_args(tsv_inputs) | {"counts": str(tmp_path / "bad.tsv")}
    completed = run_with(tmp_path, drop=H5AD_FLAGS, **args)
    assert is_error(completed, "--counts must hold non-negative integers")


def test_sample_mismatch(tmp_path, tsv_inputs):
    meta = pd.read_csv(tsv_inputs / "meta.tsv", sep="\t")
    meta.iloc[1:].to_csv(tmp_path / "meta.tsv", sep="\t", index=False)
    args = tsv_args(tsv_inputs) | {"meta": str(tmp_path / "meta.tsv")}
    completed = run_with(tmp_path, drop=H5AD_FLAGS, **args)
    assert is_error(completed, "--counts and --meta disagree on samples: D1_ctrl")


def test_duplicate_cell_ids(tmp_path):
    adata = ad.read_h5ad(DEMO)
    adata.obs_names = ["cell0"] * 2 + list(adata.obs_names[2:])
    adata.write_h5ad(tmp_path / "dup.h5ad")
    completed = run_with(tmp_path, adata=str(tmp_path / "dup.h5ad"))
    assert is_error(completed, "duplicate cell ids in --adata: cell0")


def test_duplicate_sample_ids(tmp_path, tsv_inputs):
    meta = pd.read_csv(tsv_inputs / "meta.tsv", sep="\t")
    pd.concat([meta, meta.iloc[:1]]).to_csv(
        tmp_path / "meta.tsv", sep="\t", index=False
    )
    args = tsv_args(tsv_inputs) | {"meta": str(tmp_path / "meta.tsv")}
    completed = run_with(tmp_path, drop=H5AD_FLAGS, **args)
    assert is_error(completed, "duplicate sample ids in --meta: D1_ctrl")
    counts = pd.read_csv(tsv_inputs / "counts.tsv", sep="\t", index_col=0)
    pd.concat([counts, counts.iloc[:1]]).to_csv(tmp_path / "counts.tsv", sep="\t")
    args = tsv_args(tsv_inputs) | {"counts": str(tmp_path / "counts.tsv")}
    completed = run_with(tmp_path, drop=H5AD_FLAGS, **args)
    assert is_error(completed, "duplicate sample ids in --counts: D1_ctrl")


def test_rank_deficient_design(tmp_path):
    adata = ad.read_h5ad(DEMO)
    adata.obs["batch"] = adata.obs["condition"].astype(str)
    adata.write_h5ad(tmp_path / "batch.h5ad")
    completed = run_with(
        tmp_path,
        adata=str(tmp_path / "batch.h5ad"),
        formula="~ condition + batch + (1|donor)",
    )
    assert is_error(completed, "design is rank-deficient (2 < 3 columns)")


def test_no_rscript(tmp_path):
    completed = run_with(tmp_path, env={"RSCRIPT": "/no/such/Rscript"})
    assert is_error(completed, "Rscript not found; set RSCRIPT")


def test_r_failed(tmp_path):
    completed = run_with(tmp_path, env={"RSCRIPT": "/bin/false"})
    assert is_error(completed, "crumblr.R failed with exit 1")


def test_bad_coef(tmp_path):
    assert is_error(
        run_with(tmp_path, coef="conditionctrl"),
        "--coef 'conditionctrl' is not a coefficient of the fit; "
        "choose one of: (Intercept), conditionstim",
    )


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
