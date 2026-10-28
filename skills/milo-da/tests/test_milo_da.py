"""CLI tests for scripts/milo_da.py: demo goldens, output contract, fail strings.

Run from skills/milo-da with pytest; see references/maintainer.md for the exact command.
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
SCRIPT = SKILL_DIR / "scripts" / "milo_da.py"
DEMO = SKILL_DIR / "examples" / "demo.h5ad"
TEMPLATES = SKILL_DIR.parents[1] / "templates"

# Golden demo values, computed once with pertpy 1.4.0, pydeseq2 0.5.4, scanpy 1.12.4.
# Mono is planted at 3x in condition B; 17 of 20 Mono nhoods come out up.
GOLDEN = {
    "n_cells": 1200,
    "n_samples": 6,
    "n_donors": 6,
    "n_nhoods": 75,
    "n_significant": 18,
    "n_up": 17,
    "n_down": 1,
    "dominant_annotation": "Mono",
    "significant_by_annotation": {"B": 1, "Mono": 17},
    "graph_rebuilt": False,
}

ARGS = {
    "--adata": str(DEMO),
    "--context": "synthetic demo",
    "--sample-col": "sample",
    "--donor-col": "donor",
    "--design": "~condition",
    "--model-contrasts": "conditionB-conditionA",
    "--prop": "0.1",
    "--use-rep": "X_pca",
    "--n-neighbors": "30",
    "--spatial-fdr-max": "0.1",
    "--annotation-col": "cell_type",
    "--min-annotation-frac": "0.6",
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
    """Demo arguments with `changes` ({"prop": "0"}) and `drop` applied."""
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
    adata.write_h5ad(tmp_path / "bad.h5ad")
    return str(tmp_path / "bad.h5ad")


@pytest.fixture(scope="module")
def demo_out(tmp_path_factory):
    out = tmp_path_factory.mktemp("demo") / "out"
    completed = run("--demo", "--outdir", str(out))
    assert completed.returncode == 0, completed.stderr
    return out


def test_demo_golden_values(demo_out):
    result = json.loads((demo_out / "result.json").read_text())
    assert {k: result[k] for k in GOLDEN} == GOLDEN
    assert result["seed"] == 0 and result["solver"] == "pydeseq2"
    assert len(result["chat_summary_lines"]) == 4
    assert result["chat_summary_lines"][0] == "Context: " + result["context"]
    nhoods = pd.read_csv(demo_out / "tables" / "nhoods.tsv", sep="\t")
    keys = ["SpatialFDR", "PValue", "nhood"]
    assert nhoods[keys].equals(nhoods.sort_values(keys)[keys].reset_index(drop=True))
    by_anno = pd.read_csv(demo_out / "tables" / "da_by_annotation.tsv", sep="\t")
    assert by_anno.set_index("annotation").loc["Mono", "n_up"] == 17


def test_rerun_is_identical(demo_out, tmp_path):
    out = tmp_path / "again"
    assert run("--demo", "--outdir", str(out)).returncode == 0
    for name in ["nhoods.tsv", "da_by_annotation.tsv"]:
        assert (out / "tables" / name).read_bytes() == (
            demo_out / "tables" / name
        ).read_bytes()


def test_output_contract(demo_out):
    paths = contract_paths()
    assert len(paths) == 17
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


def test_other_k_rebuilds_graph(tmp_path):
    out = tmp_path / "out"
    completed = run("--demo", "--n-neighbors", "20", "--outdir", str(out))
    assert completed.returncode == 0, completed.stderr
    result = json.loads((out / "result.json").read_text())
    assert result["graph_rebuilt"] is True
    assert "rebuilt with k = 20" in result["chat_summary_lines"][2]


def test_mixed_label(tmp_path):
    # Mono nhoods are 100% Mono, so a cutoff of 1 marks only impure nhoods Mixed.
    def blend(adata):
        adata.obs["cell_type"] = adata.obs["cell_type"].astype(str)
        adata.obs.iloc[::7, adata.obs.columns.get_loc("cell_type")] = "T"

    completed = run_with(
        tmp_path, adata=edited_demo(tmp_path, blend), min_annotation_frac="0.95"
    )
    assert completed.returncode == 0, completed.stderr
    nhoods = pd.read_csv(tmp_path / "out" / "tables" / "nhoods.tsv", sep="\t")
    mixed = nhoods["annotation"] == "Mixed"
    assert mixed.any()
    assert (nhoods.loc[mixed, "annotation_frac"] < 0.95).all()


def test_missing_argument_fails_without_demo(tmp_path):
    assert is_error(
        run_with(tmp_path, drop=["--seed"]), "missing required arguments: --seed"
    )


def test_missing_context_fails(tmp_path):
    assert is_error(
        run_with(tmp_path, drop=["--context"]), "missing required arguments: --context"
    )


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--seed", "1", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --seed")


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"prop": "0"}, "--prop must be in (0, 1], got 0.0"),
        ({"design": "condition"}, "--design must be a formula like ~x, got condition"),
        ({"spatial_fdr_max": "1"}, "--spatial-fdr-max must be in (0, 1), got 1.0"),
        ({"n_neighbors": "5000"}, "--n-neighbors must be < the 1200 cells, got 5000"),
        ({"adata": "/no/such.h5ad"}, "input file not found: /no/such.h5ad"),
        ({"sample_col": "patient"}, "obs column not found: patient"),
        ({"donor_col": "patient"}, "obs column not found: patient"),
        ({"design": "~status"}, "obs column not found: status"),
        ({"use_rep": "X_scVI"}, "obsm key not found: X_scVI"),
        ({"design": "~cell_type"}, "design column 'cell_type' takes more than one"),
        (
            {"model_contrasts": "conditionC-conditionA"},
            "pertpy da_nhoods refused the model: Contrast term 'conditionC'",
        ),
    ],
)
def test_input_errors(tmp_path, changes, message):
    completed = run_with(tmp_path, **changes)
    assert is_error(completed, message), completed.stderr
    assert not (tmp_path / "out").exists()


def test_column_na(tmp_path):
    def blank(adata):
        adata.obs["sample"] = adata.obs["sample"].astype(object)
        adata.obs.iloc[0, adata.obs.columns.get_loc("sample")] = None

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, blank))
    assert is_error(completed, "obs column 'sample' has missing values")


def test_duplicate_cell_ids(tmp_path):
    def dupe(adata):
        adata.obs_names = [adata.obs_names[0], *adata.obs_names[1:-1], "cell0000"]

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, dupe))
    assert is_error(completed, "duplicate cell ids: cell0000")


def test_rank_deficient_design(tmp_path):
    # batch equals condition: the two terms are aliased.
    def alias(adata):
        adata.obs["batch"] = adata.obs["condition"]

    completed = run_with(
        tmp_path, adata=edited_demo(tmp_path, alias), design="~batch + condition"
    )
    assert is_error(completed, "design is rank-deficient (2 < 3 columns)")


def test_donor_without_design_term(tmp_path):
    # D1 gives samples A1 and B1, but the design has no donor term.
    def share(adata):
        adata.obs["donor"] = adata.obs["donor"].astype(str).replace("D4", "D1")

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, share))
    assert is_error(completed, "1 donor(s) in 'donor' have more than one sample")


def test_zero_significant(tmp_path):
    # Sample labels shuffled over cells: every sample is a random draw, so no DA.
    def shuffle(adata):
        rng = np.random.default_rng(0)
        sample = rng.permutation(adata.obs["sample"].astype(str).to_numpy())
        adata.obs["sample"] = sample
        adata.obs["donor"] = ["D" + x for x in sample]
        adata.obs["condition"] = [x[0] for x in sample]

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, shuffle))
    assert completed.returncode == 0, completed.stderr
    out = tmp_path / "out"
    result = json.loads((out / "result.json").read_text())
    assert result["n_significant"] == 0 and result["dominant_annotation"] is None
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


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
