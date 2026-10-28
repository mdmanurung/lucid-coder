"""CLI tests for scripts/group_tests.py: demo goldens, output contract, fail strings.

Run from skills/group-tests with pytest;
see references/maintainer.md for the exact command.
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "group_tests.py"
DEMO = SKILL_DIR / "examples" / "demo.tsv"
TEMPLATES = SKILL_DIR.parents[1] / "templates"

# Planted truth from examples/make_demo.py.
VISIT_TRUTH = {f"F{i:02d}" for i in range(1, 7)}
COHORT_TRUTH = {f"F{i:02d}" for i in range(7, 13)}
AGE_TRUTH = {f"F{i:02d}" for i in range(13, 19)}

ARGS = {
    "--context": "synthetic demo, pre vs post",
    "--input": str(DEMO),
    "--test": "paired-wilcoxon",
    "--group-col": "visit",
    "--levels": "pre,post",
    "--pair-col": "donor",
    "--value-unit": "a.u.",
    "--adjust": "bh",
    "--alpha": "0.05",
    "--n-plot": "4",
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
    """Demo arguments with `changes` ({"alpha": "2"}) and `drop` applied."""
    args = {**ARGS, **{f"--{k.replace('_', '-')}": v for k, v in changes.items()}}
    argv = [x for k, v in args.items() if k not in drop for x in (k, v)]
    return run(*argv, "--outdir", str(tmp_path / "out"))


def significant(outdir):
    tests = pd.read_csv(outdir / "tables" / "tests.tsv", sep="\t")
    return set(tests.loc[tests["significant"], "feature"])


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
def pre_only(tmp_path_factory):
    """Pre-visit rows only: one value per donor, for unpaired and correlation tests."""
    path = tmp_path_factory.mktemp("pre") / "pre.tsv"
    data = pd.read_csv(DEMO, sep="\t", dtype=str)
    data[data["visit"] == "pre"].to_csv(path, sep="\t", index=False)
    return str(path)


def test_demo_golden_values(demo_out):
    result = json.loads((demo_out / "result.json").read_text())
    assert result["n_tested"] == 30
    assert result["n_significant"] == 7
    # All six planted shifts plus F15 (an age feature, null for pre vs post):
    # BH bounds the FDR, not zero.
    assert significant(demo_out) == VISIT_TRUTH | {"F15"}
    assert result["seed"] == 0 and result["test"] == "paired-wilcoxon"
    assert 3 <= len(result["chat_summary_lines"]) <= 4
    assert result["chat_summary_lines"][0].startswith("Context: synthetic demo")
    assert result["context"].startswith("synthetic demo")
    tests = pd.read_csv(demo_out / "tables" / "tests.tsv", sep="\t")
    assert tests["feature"].head(6).tolist() == sorted(VISIT_TRUTH)
    assert (tests.loc[tests["feature"].isin(VISIT_TRUTH), "hodges_lehmann"] > 0).all()


@pytest.mark.parametrize(
    ("changes", "expected"),
    [
        (
            {"test": "wilcoxon", "group_col": "cohort", "levels": "healthy,disease"},
            COHORT_TRUTH | {"F17"},
        ),
        (
            {
                "test": "wilcoxon",
                "group_col": "cohort",
                "levels": "healthy,disease",
                "adjust": "bonferroni",
            },
            {"F07", "F08", "F09", "F10", "F12"},
        ),
        ({"test": "spearman", "covariate": "age"}, AGE_TRUTH),
        ({"test": "pearson", "covariate": "age"}, AGE_TRUTH),
    ],
)
def test_other_tests_golden(tmp_path, pre_only, changes, expected):
    drop = ["--group-col", "--levels", "--pair-col"]
    if "group_col" in changes:
        drop = ["--pair-col"]
    completed = run_with(tmp_path, drop=drop, input=pre_only, **changes)
    assert completed.returncode == 0, completed.stderr
    assert significant(tmp_path / "out") == expected


def test_drops_per_group(tmp_path):
    """A missing value is counted under its own group level, not only in a total."""
    path = tmp_path / "gap.tsv"
    data = pd.read_csv(DEMO, sep="\t", dtype=str)
    data.loc[(data["feature"] == "F20") & (data["subject"] == "D01_post"), "value"] = ""
    data.to_csv(path, sep="\t", index=False)
    completed = run_with(tmp_path, input=str(path))
    assert completed.returncode == 0, completed.stderr
    result = json.loads((tmp_path / "out" / "result.json").read_text())
    assert result["n_missing_dropped_by_group"] == {"pre": 0, "post": 1}


def test_duplicate_feature_subject(tmp_path):
    """Unique ids: one row per (feature, subject)."""
    path = tmp_path / "dup.tsv"
    data = pd.read_csv(DEMO, sep="\t", dtype=str)
    pd.concat([data, data.head(1)]).to_csv(path, sep="\t", index=False)
    completed = run_with(tmp_path, input=str(path))
    assert is_error(completed, "2 rows repeat a (feature, subject) combination")


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
    assert "--seed 0" in (demo_out / "reproducibility" / "commands.sh").read_text()


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


def test_missing_argument_fails_without_demo(tmp_path):
    assert is_error(
        run_with(tmp_path, drop=["--seed"]), "missing required arguments: --seed"
    )


def test_missing_context_fails(tmp_path):
    assert is_error(
        run_with(tmp_path, drop=["--context"]), "missing required arguments: --context"
    )


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--alpha", "0.1", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --alpha")


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


@pytest.mark.parametrize(
    ("drop", "changes", "message"),
    [
        ((), {"alpha": "1.5"}, "--alpha must be in (0, 1), got 1.5"),
        (
            ("--pair-col",),
            {},
            "--test paired-wilcoxon needs --pair-col",
        ),
        (
            (),
            {"test": "spearman", "covariate": "age"},
            "--test spearman does not use --group-col --levels --pair-col",
        ),
        ((), {"levels": "pre"}, "--levels must be two different values a,b"),
        ((), {"input": "/no/such.tsv"}, "input file not found: /no/such.tsv"),
        ((), {"group_col": "batch"}, "input column not found: batch"),
        (
            ("--group-col", "--levels", "--pair-col"),
            {"test": "pearson", "covariate": "cohort"},
            "column 'cohort' must be numeric",
        ),
        ((), {"levels": "pre,mid"}, "level 'mid' not found in column 'visit'"),
        (
            (),
            {"pair_col": "cohort"},
            "1200 rows repeat a (feature, visit, cohort) combination",
        ),
        ((), {"pair_col": "subject"}, "no feature could be tested"),
    ],
)
def test_input_errors(tmp_path, drop, changes, message):
    assert is_error(run_with(tmp_path, drop=drop, **changes), message)


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
