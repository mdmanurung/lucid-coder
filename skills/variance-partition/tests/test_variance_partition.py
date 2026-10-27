"""CLI tests for scripts/variance_partition.py: demo goldens, output contract, fails.

Run from skills/variance-partition with pytest. The R engine needs RSCRIPT,
R_LIBS_USER and LD_LIBRARY_PATH; see references/maintainer.md. Each R run loads
variancePartition, which takes about 2 minutes on the shared filesystem.
"""

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "variance_partition.py"
COUNTS = SKILL_DIR / "examples" / "demo_counts.tsv"
META = SKILL_DIR / "examples" / "demo_meta.tsv"
TEMPLATES = SKILL_DIR.parents[1] / "templates"
RSCRIPT = os.environ.get("RSCRIPT") or shutil.which("Rscript")
needs_r = pytest.mark.skipif(not RSCRIPT, reason="Rscript not found")

# Golden demo values (variancePartition 1.40.2, lme4 2.1.0). Planted blocks of 20:
# G000-G019 donor, G020-G039 condition, G040-G059 noise.
PLANTED = {
    "donor": {f"G{i:03d}" for i in range(0, 20)},
    "condition": {f"G{i:03d}" for i in range(20, 40)},
    "Residuals": {f"G{i:03d}" for i in range(40, 60)},
}

ARGS = {
    "--counts": str(COUNTS),
    "--meta": str(META),
    "--formula": "~ (1|donor) + (1|condition) + rin",
    "--min-frac": "0.1",
    "--n-top": "5",
    "--max-cancor": "0.9",
    "--reference": "synthetic demo genes, no genome build",
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


def run_with(tmp_path, env=None, drop=(), args=ARGS, **changes):
    """`args` with `changes` ({"min_frac": "2"}) and `drop` applied."""
    merged = {**args, **{f"--{k.replace('_', '-')}": v for k, v in changes.items()}}
    argv = [x for k, v in merged.items() if k not in drop for x in (k, v)]
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


def write_cells(path, cells_per_sample, rin_varies=False):
    """Split each demo sample's counts over its cells; the sums stay exact."""
    import anndata as ad

    counts = pd.read_csv(COUNTS, sep="\t", index_col=0)
    meta = pd.read_csv(META, sep="\t", index_col=0)
    rng = np.random.default_rng(0)
    blocks, obs = [], []
    for sample in counts.columns:
        n = cells_per_sample.get(sample, 4)
        share = np.full(n, 1 / n)
        blocks.append(
            np.stack([rng.multinomial(c, share) for c in counts[sample]], axis=1)
        )
        obs.append(meta.loc[[sample] * n].assign(sample=sample))
    obs = pd.concat(obs, ignore_index=True)
    if rin_varies:
        obs["rin"] = np.arange(len(obs), dtype=float)
    obs.index = [f"cell{i:04d}" for i in range(len(obs))]
    adata = ad.AnnData(
        X=np.vstack(blocks).astype(np.float32),
        obs=obs,
        var=pd.DataFrame(index=counts.index),
    )
    adata.write_h5ad(path)
    return path


ADATA_ARGS = {
    "--sample-col": "sample",
    "--donor-col": "donor",
    "--min-cells": "3",
    **{k: v for k, v in ARGS.items() if k not in ("--counts", "--meta")},
}


@pytest.fixture(scope="module")
def demo_out(tmp_path_factory):
    if not RSCRIPT:
        pytest.skip("Rscript not found")
    out = tmp_path_factory.mktemp("demo") / "out"
    completed = run("--demo", "--outdir", str(out))
    assert completed.returncode == 0, completed.stderr
    return out


def test_demo_golden_values(demo_out):
    result = json.loads((demo_out / "result.json").read_text())
    assert result["n_samples"] == 24 and result["n_genes"] == 120
    assert result["top_covariate"] == "donor"
    assert list(result["median_fractions"]) == [
        "donor",
        "Residuals",
        "rin",
        "condition",
    ]
    assert result["max_cancor"] == pytest.approx(0.761163, abs=1e-6)
    for term, planted in PLANTED.items():
        assert set(result["top_genes"][term]) <= planted, term
    assert result["seed"] == 0 and result["n_donors"] is None
    assert result["reference"] == "synthetic demo genes, no genome build"
    assert 2 <= len(result["chat_summary_lines"]) <= 3
    fractions = pd.read_csv(demo_out / "tables" / "variance_fractions.tsv", sep="\t")
    assert list(fractions.columns) == ["gene", "donor", "condition", "rin", "Residuals"]
    total = fractions.drop(columns="gene").sum(axis=1)
    assert np.allclose(total, 1)


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
    assert (
        "variancePartition "
        in (demo_out / "reproducibility" / "environment.txt").read_text()
    )


@needs_r
def test_adata_input_matches_counts(tmp_path, demo_out):
    """Cells that sum back to the demo counts give the same fractions."""
    adata = write_cells(tmp_path / "cells.h5ad", {"D01_ctrl": 2})
    completed = run_with(tmp_path, args=ADATA_ARGS, adata=str(adata), min_cells="1")
    assert completed.returncode == 0, completed.stderr
    result = json.loads((tmp_path / "out" / "result.json").read_text())
    assert (result["n_donors"], result["n_samples"]) == (12, 24)
    got = pd.read_csv(tmp_path / "out" / "tables" / "variance_fractions.tsv", sep="\t")
    want = pd.read_csv(demo_out / "tables" / "variance_fractions.tsv", sep="\t")
    pd.testing.assert_frame_equal(got, want, rtol=1e-6)


@needs_r
def test_collinear(tmp_path):
    completed = run_with(tmp_path, max_cancor="0.5")
    assert is_error(
        completed,
        "canonical correlation 0.761 between donor and rin exceeds --max-cancor 0.5; "
        "all pairs: donor/condition ",
    )
    assert not (tmp_path / "out").exists()


def test_missing_argument_fails_without_demo(tmp_path):
    assert is_error(
        run_with(tmp_path, drop=["--seed"]), "missing required arguments: --seed"
    )


def test_missing_adata_arguments(tmp_path):
    completed = run_with(tmp_path, drop=["--counts", "--meta"], adata="x.h5ad")
    assert is_error(
        completed, "missing required arguments: --sample-col --donor-col --min-cells"
    )


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--seed", "1", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --seed")


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"min_frac": "2"}, "--min-frac must be in [0, 1], got 2.0"),
        ({"formula": "donor + rin"}, "--formula must be '~ (1|a) + b'"),
        ({"formula": "~ " + " + ".join("abcdefghij")}, "--formula has 10 terms"),
        ({"adata": "x.h5ad"}, "give --counts with --meta, or --adata"),
        ({"min_cells": "3"}, "--min-cells applies only to --adata input"),
        ({"counts": "/no/such.tsv"}, "input file not found: /no/such.tsv"),
        ({"formula": "~ (1|batch)"}, "column not found in the sample metadata: batch"),
        ({"formula": "~ donor + rin"}, "covariate 'donor' is categorical"),
    ],
)
def test_input_errors(tmp_path, changes, message):
    assert is_error(run_with(tmp_path, **changes), message)


def test_bad_counts(tmp_path):
    counts = pd.read_csv(COUNTS, sep="\t", index_col=0) * 0.5
    counts.to_csv(tmp_path / "half.tsv", sep="\t")
    completed = run_with(tmp_path, counts=str(tmp_path / "half.tsv"))
    assert is_error(completed, f"counts in {tmp_path / 'half.tsv'} must be")


def test_no_genes(tmp_path):
    """One silent sample: no gene is detected in every sample."""
    counts = pd.read_csv(COUNTS, sep="\t", index_col=0)
    counts.iloc[:, 0] = 0
    counts.to_csv(tmp_path / "zero.tsv", sep="\t")
    completed = run_with(tmp_path, counts=str(tmp_path / "zero.tsv"), min_frac="1")
    assert is_error(completed, "zero genes pass --min-frac 1.0")


def test_too_few_count_samples(tmp_path):
    counts = pd.read_csv(COUNTS, sep="\t", index_col=0).iloc[:, :2]
    meta = pd.read_csv(META, sep="\t").iloc[:2]
    counts.to_csv(tmp_path / "c.tsv", sep="\t")
    meta.to_csv(tmp_path / "m.tsv", sep="\t", index=False)
    completed = run_with(
        tmp_path, counts=str(tmp_path / "c.tsv"), meta=str(tmp_path / "m.tsv")
    )
    assert is_error(completed, "2 samples in --counts; need at least 3")


def test_duplicate_ids(tmp_path):
    counts = pd.read_csv(COUNTS, sep="\t", index_col=0)
    dup_gene = pd.concat([counts, counts.head(1)])
    dup_gene.to_csv(tmp_path / "g.tsv", sep="\t")
    completed = run_with(tmp_path, counts=str(tmp_path / "g.tsv"))
    assert is_error(completed, "duplicate gene names: G000")
    lines = (COUNTS.read_text()).splitlines()
    header = lines[0].split("\t")
    header[2] = header[1]
    (tmp_path / "s.tsv").write_text("\n".join(["\t".join(header), *lines[1:]]) + "\n")
    completed = run_with(tmp_path, counts=str(tmp_path / "s.tsv"))
    assert is_error(completed, f"duplicate sample ids: {header[1]}")


def test_rank_deficient(tmp_path):
    meta = pd.read_csv(META, sep="\t").assign(rin2=lambda m: 2 * m["rin"])
    meta.to_csv(tmp_path / "m.tsv", sep="\t", index=False)
    completed = run_with(
        tmp_path, meta=str(tmp_path / "m.tsv"), formula="~ (1|donor) + rin + rin2"
    )
    assert is_error(completed, "design is rank-deficient")


def test_replication_needs_donor_term(tmp_path):
    """Two samples per donor and no donor term: pseudoreplication stops the run."""
    adata = str(write_cells(tmp_path / "cells.h5ad", {}))
    completed = run_with(
        tmp_path, args=ADATA_ARGS, adata=adata, formula="~ (1|condition) + rin"
    )
    assert is_error(completed, "12 donor(s) in 'donor' have more than one sample")


def test_samples_mismatch(tmp_path):
    meta = pd.read_csv(META, sep="\t").iloc[1:]
    meta.to_csv(tmp_path / "meta.tsv", sep="\t", index=False)
    completed = run_with(tmp_path, meta=str(tmp_path / "meta.tsv"))
    assert is_error(
        completed, "--counts columns and --meta rows name different samples: D01_ctrl"
    )


def test_adata_errors(tmp_path):
    adata = str(write_cells(tmp_path / "cells.h5ad", {}))
    cases = [
        ({"layer": "counts"}, "layer not found in the h5ad: counts"),
        ({"min_cells": "5"}, "0 samples left after --min-cells 5; need at least 3"),
    ]
    for i, (changes, message) in enumerate(cases):
        sub = tmp_path / str(i)
        completed = run_with(sub, args=ADATA_ARGS, adata=adata, **changes)
        assert is_error(completed, message), completed.stderr
    varies = str(write_cells(tmp_path / "varies.h5ad", {}, rin_varies=True))
    completed = run_with(tmp_path / "v", args=ADATA_ARGS, adata=varies)
    assert is_error(completed, "obs column 'rin' varies within a sample")


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


def test_no_rscript(tmp_path):
    completed = run_with(tmp_path, env={"RSCRIPT": "", "PATH": str(tmp_path)})
    assert is_error(completed, "Rscript not found; set RSCRIPT or put Rscript on PATH")


def test_r_failure(tmp_path):
    completed = run_with(tmp_path, env={"RSCRIPT": shutil.which("false")})
    assert is_error(completed, "varpart.R failed with exit 1")


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
