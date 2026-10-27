"""CLI tests for scripts/mashr_sharing.py: demo goldens, output contract, fail strings.

Run from skills/mashr-sharing with pytest;
see references/maintainer.md for the R environment.
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
SCRIPT = SKILL_DIR / "scripts" / "mashr_sharing.py"
EXAMPLES = SKILL_DIR / "examples"
TEMPLATES = SKILL_DIR.parents[1] / "templates"
CELL_TYPES = ["B", "CD4T", "CD8T", "Mono"]
RSCRIPT = os.environ.get("RSCRIPT") or shutil.which("Rscript")
needs_r = pytest.mark.skipif(not RSCRIPT, reason="Rscript not found")

# Planted truth from examples/make_demo.py.
SHARED = {f"G{i:04d}" for i in range(80)}
SPECIFIC = {
    ct: {f"G{80 + 20 * k + i:04d}" for i in range(20)}
    for k, ct in enumerate(CELL_TYPES)
}

DE = [x for ct in CELL_TYPES for x in ("--de", f"{ct}={EXAMPLES / f'demo_{ct}.tsv'}")]
ARGS = {
    "--gene-col": "gene",
    "--effect-col": "log2_fc",
    "--se-col": "lfc_se",
    "--align": "intersect",
    "--covariances": "canonical",
    "--strong-lfsr": "0.05",
    "--lfsr-max": "0.05",
    "--sharing-factor": "0.5",
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


def run_with(tmp_path, de=DE, env=None, drop=(), **changes):
    """Demo tables and arguments (canonical) with `changes` and `drop` applied."""
    args = {**ARGS, **{f"--{k.replace('_', '-')}": v for k, v in changes.items()}}
    argv = [x for k, v in args.items() if k not in drop for x in (k, v)]
    return run(*de, *argv, "--outdir", str(tmp_path / "out"), env=env)


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


def write_table(tmp_path, name, frame):
    path = tmp_path / f"{name}.tsv"
    frame.to_csv(path, sep="\t", index=False)
    return f"{name}={path}"


def demo_table(ct):
    return pd.read_csv(EXAMPLES / f"demo_{ct}.tsv", sep="\t")


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
    assert result["n_genes"] == 800 and result["cell_types"] == CELL_TYPES
    assert result["n_strong"] == 154
    assert result["n_significant_any"] == 153
    assert result["n_shared_all"] == 80
    assert result["n_specific"] == 73
    assert result["significant_per_cell_type"] == {
        "B": 96,
        "CD4T": 99,
        "CD8T": 100,
        "Mono": 98,
    }
    assert result["scatter_pair"] == ["B", "Mono"]
    assert result["genes_dropped_by_align"] == dict.fromkeys(CELL_TYPES, 0)
    assert result["reference"] == "synthetic demo genes, no genome build"
    assert result["seed"] == 0 and 2 <= len(result["chat_summary_lines"]) <= 3
    post = pd.read_csv(demo_out / "tables" / "posterior.tsv", sep="\t")
    n_sig = post[post["significant"]].groupby("gene")["cell_type"].agg(list)
    # Every shared gene is recovered in all four types; every one-type call is
    # a planted specific gene in its own cell type; no null gene passes.
    assert set(n_sig[n_sig.str.len() == 4].index) == SHARED
    for gene, types in n_sig[n_sig.str.len() == 1].items():
        assert gene in SPECIFIC[types[0]], gene
    assert set(n_sig.index) <= SHARED.union(*SPECIFIC.values())
    counts = pd.read_csv(demo_out / "tables" / "significant_counts.tsv", sep="\t")
    assert counts["n_genes"].tolist() == [73, 0, 0, 80]
    sharing = pd.read_csv(demo_out / "tables" / "pairwise_sharing.tsv", sep="\t")
    assert sharing["sharing"].between(0.67, 0.71).all() and len(sharing) == 6


def test_output_contract(demo_out):
    paths = contract_paths()
    assert len(paths) == 18
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
    commands = (demo_out / "reproducibility" / "commands.sh").read_text()
    assert "--seed 0" in commands and commands.count("--de ") == 4


def test_rerun_is_identical(demo_out, tmp_path):
    completed = run("--demo", "--outdir", str(tmp_path / "again"))
    assert completed.returncode == 0, completed.stderr
    check = subprocess.run(
        ["sha256sum", "-c", str(demo_out / "reproducibility" / "checksums.sha256")],
        cwd=tmp_path / "again",
        capture_output=True,
        text=True,
        check=False,
    )
    assert check.returncode == 0, check.stdout


@needs_r
def test_intersect_drops_per_cell_type(tmp_path):
    de = ["--de", write_table(tmp_path, "B", demo_table("B").iloc[50:]), *DE[2:]]
    completed = run_with(tmp_path, de=de)
    assert completed.returncode == 0, completed.stderr
    result = json.loads((tmp_path / "out" / "result.json").read_text())
    assert result["genes_dropped_by_align"] == {
        "B": 0,
        "CD4T": 50,
        "CD8T": 50,
        "Mono": 50,
    }


@needs_r
def test_refusal_leaves_no_outdir(tmp_path):
    completed = run_with(tmp_path, covariances="canonical+ed", strong_lfsr="1e-300")
    assert completed.returncode == 2 and not (tmp_path / "out").exists()


@needs_r
def test_union_keeps_missing_genes(tmp_path):
    """Genes absent from one table stay NA there, and sharing has no NA."""
    de = [
        x
        for ct in CELL_TYPES[1:]
        for x in ("--de", f"{ct}={EXAMPLES / f'demo_{ct}.tsv'}")
    ]
    b = demo_table("B")
    de = ["--de", write_table(tmp_path, "B", b.iloc[50:]), *de]
    completed = run_with(tmp_path, de=de, align="union")
    assert completed.returncode == 0, completed.stderr
    post = pd.read_csv(tmp_path / "out" / "tables" / "posterior.tsv", sep="\t")
    assert len(post) == 4 * 800 - 50
    assert post["pm"].notna().all() and post["lfsr"].notna().all()
    sharing = pd.read_csv(
        tmp_path / "out" / "tables" / "pairwise_sharing.tsv", sep="\t"
    )
    assert sharing["sharing"].notna().all()
    spearman = pd.read_csv(tmp_path / "out" / "tables" / "raw_spearman.tsv", sep="\t")
    assert spearman.loc[spearman["cell_type_a"] == "B", "n_genes"].eq(750).all()


def test_missing_argument_fails_without_demo(tmp_path):
    assert is_error(
        run_with(tmp_path, drop=["--seed"]), "missing required arguments: --seed"
    )


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--seed", "1", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --seed")


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


@pytest.mark.parametrize(
    ("de", "changes", "message"),
    [
        (DE, {"lfsr_max": "1"}, "--lfsr-max must be in (0, 1), got 1.0"),
        (["--de", "B"], {}, "--de must be name=path, got 'B'"),
        (DE[:2] * 2, {}, "--de name 'B' is given twice"),
        (DE[:2], {}, "need --de for at least 2 cell types, got 1"),
        (["--de", "B=/no/such.tsv", *DE[2:]], {}, "input file not found: /no/such.tsv"),
        (DE, {"se_col": "se"}, "column 'se' not found in the --de table for 'B'"),
    ],
)
def test_argument_errors(tmp_path, de, changes, message):
    assert is_error(run_with(tmp_path, de=de, **changes), message)


@pytest.mark.parametrize(
    ("edit", "message"),
    [
        (
            lambda t: pd.concat([t, t.head(2)]),
            "2 genes repeat in the --de table for 'B', e.g. 'G0000'",
        ),
        (
            lambda t: t.assign(
                log2_fc=t["log2_fc"].astype(object).where(t.index > 0, "x")
            ),
            "column 'log2_fc' in the --de table for 'B' must be numeric",
        ),
        (
            lambda t: t.assign(lfc_se=t["lfc_se"].where(t.index > 2, 0)),
            "3 rows of 'lfc_se' in the --de table for 'B' are <= 0",
        ),
        (
            lambda t: t.assign(gene="X" + t["gene"]),
            "zero genes are shared by all --de tables (--align intersect)",
        ),
    ],
)
def test_table_errors(tmp_path, edit, message):
    de = ["--de", write_table(tmp_path, "B", edit(demo_table("B"))), *DE[2:]]
    assert is_error(run_with(tmp_path, de=de), message)


def test_no_rscript(tmp_path):
    completed = run_with(tmp_path, env={"RSCRIPT": "", "PATH": str(tmp_path)})
    assert is_error(completed, "Rscript not found; set RSCRIPT or put Rscript on PATH")


def test_r_failed(tmp_path):
    completed = run_with(tmp_path, env={"RSCRIPT": shutil.which("false")})
    assert is_error(completed, "mash.R failed with exit 1")


@needs_r
def test_missing_r_package(tmp_path):
    # An empty R_ENVIRON_USER stops ~/.Renviron from restoring the library path.
    (tmp_path / "Renviron").write_text("")
    env = {"R_LIBS_USER": str(tmp_path), "R_ENVIRON_USER": str(tmp_path / "Renviron")}
    completed = run_with(tmp_path, env=env)
    assert is_error(completed, "R package not installed: mashr; ask for it")


@needs_r
def test_too_few_strong(tmp_path):
    completed = run_with(tmp_path, covariances="canonical+ed", strong_lfsr="1e-300")
    assert is_error(
        completed,
        "0 genes pass --strong-lfsr 1e-300 in any cell type, fewer than the 4",
    )


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
