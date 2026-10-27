"""CLI tests for scripts/sc_qc.py: demo goldens, output contract, fail strings.

Run from skills/sc-qc with pytest; see references/maintainer.md for the exact uv line.
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
from scipy import sparse

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "sc_qc.py"
DEMO = SKILL_DIR / "examples" / "demo.h5ad"
TEMPLATES = SKILL_DIR.parents[1] / "templates"

# Golden demo values (mad mode 5/5/3, auto threshold), scanpy 1.12.4.
N_BEFORE, N_AFTER, N_FAILED_QC, N_DOUBLETS = 1035, 877, 86, 83
AFTER_PER_SAMPLE = {"S1": 289, "S2": 293, "S3": 295}
# Fixed mode on the demo: genes 200-2500, counts >= 500, % mito <= 10.
FIXED_AFTER_PER_SAMPLE = {"S1": 300, "S2": 300, "S3": 300}

ARGS = {
    "--adata": str(DEMO),
    "--sample-col": "sample",
    "--mito-prefix": "MT-",
    "--mode": "mad",
    "--mads-counts": "5",
    "--mads-genes": "5",
    "--mads-mito": "3",
    "--expected-doublet-rate": "0.08",
    "--doublet-threshold": "auto",
    "--seed": "0",
    "--reference": "synthetic demo genes, no genome build",
}
FIXED = {
    "mode": "fixed",
    "min_genes": "200",
    "max_genes": "2500",
    "min_counts": "500",
    "max_pct_mito": "10",
}
MAD_FLAGS = ["--mads-counts", "--mads-genes", "--mads-mito"]


def run(*argv):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *argv],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )


def run_with(tmp_path, drop=(), **changes):
    """Demo arguments with `changes` ({"seed": "1"}) and `drop` applied."""
    args = {**ARGS, **{f"--{k.replace('_', '-')}": v for k, v in changes.items()}}
    argv = [x for k, v in args.items() if k not in drop for x in (k, v)]
    return run(*argv, "--outdir", str(tmp_path / "out"))


def write_h5ad(tmp_path, adata):
    path = tmp_path / "input.h5ad"
    adata.write_h5ad(path)
    return str(path)


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
    assert result["n_cells_before"] == N_BEFORE
    assert result["n_cells_after"] == N_AFTER
    assert result["n_failed_qc"] == N_FAILED_QC
    assert result["n_doublets"] == N_DOUBLETS
    assert result["cells_after_per_sample"] == AFTER_PER_SAMPLE
    assert result["reference"] == "synthetic demo genes, no genome build"
    for s, n in AFTER_PER_SAMPLE.items():
        assert n + result["failed_qc_per_sample"][s] <= N_BEFORE // 3
    assert sum(result["failed_qc_per_sample"].values()) == N_FAILED_QC
    assert sum(result["doublets_per_sample"].values()) == N_DOUBLETS
    assert result["seed"] == 0 and result["mode"] == "mad"
    assert 2 <= len(result["chat_summary_lines"]) <= 3
    # Every planted doublet and low-quality cell is removed; only good cells stay.
    planted = ad.read_h5ad(DEMO).obs["planted"]
    cells = pd.read_csv(demo_out / "tables" / "qc_per_cell.tsv", sep="\t")
    kept = planted[cells.loc[cells["kept"], "cell"]]
    assert set(kept) == {"good"}
    filtered = ad.read_h5ad(demo_out / "data" / "qc_filtered.h5ad")
    assert filtered.n_obs == N_AFTER and filtered.n_vars == 533


def test_fixed_mode(tmp_path):
    completed = run_with(tmp_path, drop=MAD_FLAGS, **FIXED)
    assert completed.returncode == 0, completed.stderr
    result = json.loads((tmp_path / "out" / "result.json").read_text())
    assert result["cells_after_per_sample"] == FIXED_AFTER_PER_SAMPLE


def test_output_contract(demo_out):
    paths = contract_paths()
    assert len(paths) == 20
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
    completed = run_with(tmp_path, drop=["--seed"])
    assert is_error(completed, "missing required arguments: --seed")


def test_missing_reference_fails(tmp_path):
    completed = run_with(tmp_path, drop=["--reference"])
    assert is_error(completed, "missing required arguments: --reference")


@pytest.mark.parametrize(
    ("axis", "what"), [("var", "gene names: MT-G00"), ("obs", "cell ids: cell0000")]
)
def test_duplicate_ids(tmp_path, axis, what):
    adata = ad.read_h5ad(DEMO)
    names = list(getattr(adata, f"{axis}_names"))
    names[1] = names[0]
    setattr(adata, f"{axis}_names", names)
    completed = run_with(tmp_path, adata=write_h5ad(tmp_path, adata))
    assert is_error(completed, f"duplicate {what}")


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--seed", "1", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --seed")


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    completed = run_with(tmp_path, drop=MAD_FLAGS, **FIXED)
    assert is_error(completed, "outdir already exists")


@pytest.mark.parametrize(
    ("drop", "changes", "message"),
    [
        (["--mads-mito"], {}, "--mode mad needs: --mads-mito"),
        ([], {"min_genes": "200"}, "--mode mad does not use: --min-genes"),
        ([], {"expected_doublet_rate": "1.5"}, "--expected-doublet-rate must be"),
        ([], {"doublet_threshold": "high"}, "--doublet-threshold must be 'auto'"),
        ([], {"mads_genes": "0"}, "--mads-genes must be > 0, got 0.0"),
        (
            MAD_FLAGS,
            {**FIXED, "max_genes": "100"},
            "--max-genes must be > --min-genes, got 100",
        ),
        ([], {"adata": "/no/such.h5ad"}, "input file not found: /no/such.h5ad"),
        ([], {"layer": "counts"}, "layer not found in the h5ad: counts"),
        ([], {"sample_col": "donor"}, "obs column not found: donor"),
        ([], {"mito_prefix": "mt-"}, "no gene name starts with --mito-prefix 'mt-'"),
        (
            MAD_FLAGS,
            {**FIXED, "min_counts": "1000000"},
            "zero cells pass QC and doublet removal",
        ),
    ],
)
def test_input_errors(tmp_path, drop, changes, message):
    assert is_error(run_with(tmp_path, drop=drop, **changes), message)


def test_non_integer_counts(tmp_path):
    adata = ad.read_h5ad(DEMO)
    adata.X = adata.X * 0.5
    completed = run_with(tmp_path, adata=write_h5ad(tmp_path, adata))
    assert is_error(completed, "counts in X must be non-negative integers")


def test_too_few_cells(tmp_path):
    adata = ad.read_h5ad(DEMO)
    s1 = np.flatnonzero(adata.obs["sample"] == "S1")
    small = adata[~np.isin(np.arange(adata.n_obs), s1[20:])].copy()
    completed = run_with(tmp_path, adata=write_h5ad(tmp_path, small))
    message = "sample 'S1' has 20 cells; Scrublet needs more than 30"
    assert is_error(completed, message)


def test_scrublet_failed(tmp_path):
    # Constant counts leave Scrublet no variable gene, so its PCA raises.
    names = [f"MT-{i}" for i in range(5)] + [f"G{i}" for i in range(50)]
    adata = ad.AnnData(
        X=sparse.csr_matrix(np.full((40, len(names)), 3, dtype=np.int32)),
        obs=pd.DataFrame({"sample": "A"}, index=[f"c{i}" for i in range(40)]),
        var=pd.DataFrame(index=names),
    )
    completed = run_with(tmp_path, adata=write_h5ad(tmp_path, adata))
    assert is_error(completed, "Scrublet failed:")


def test_no_threshold(monkeypatch):
    """Auto threshold failed for a sample: scanpy leaves no threshold in uns."""
    sys.path.insert(0, str(SCRIPT.parent))
    import sc_qc
    import scanpy as sc
    import skillkit

    def no_threshold(adata, **kwargs):
        adata.obs["doublet_score"] = 0.1
        adata.uns["scrublet"] = {"batches": {"A": {"parameters": {}}}}

    monkeypatch.setattr(sc.pp, "scrublet", no_threshold)
    cells = ad.AnnData(
        X=np.ones((40, 2), dtype=np.float32),
        obs=pd.DataFrame({"sample": "A"}, index=[f"c{i}" for i in range(40)]),
    )
    p = {"doublet_threshold": "auto", "expected_doublet_rate": 0.1, "seed": 0}
    with pytest.raises(skillkit.SkillError, match="found no doublet-score threshold"):
        sc_qc.scrublet(cells, p)


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
