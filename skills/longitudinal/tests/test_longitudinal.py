"""CLI tests for scripts/longitudinal.py: demo goldens, output contract, fail strings.

Run from skills/longitudinal with pytest and the R environment;
see references/maintainer.md.
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
SCRIPT = SKILL_DIR / "scripts" / "longitudinal.py"
DEMO = SKILL_DIR / "examples" / "demo.tsv"
TEMPLATES = SKILL_DIR.parents[1] / "templates"
RSCRIPT = os.environ.get("RSCRIPT") or shutil.which("Rscript")
needs_r = pytest.mark.skipif(not RSCRIPT, reason="Rscript not found")

# Planted truth from examples/make_demo.py.
TIME_TRUTH = {f"F{i:02d}" for i in range(1, 7)}
GAP_SUBJECT = "S24"  # misses visit 2 in every feature

ARGS = {
    "--context": "synthetic demo, 3 visits",
    "--input": str(DEMO),
    "--formula": "value ~ time + age + (1|subject)",
    "--term": "time3",
    "--time-type": "factor",
    "--reference-time": "1",
    "--value-unit": "a.u.",
    "--time-unit": "visit",
    "--adjust": "bh",
    "--padj-max": "0.05",
    "--n-plot": "4",
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


def run_with(tmp_path, env=None, drop=(), **changes):
    """Demo arguments with `changes` ({"padj_max": "2"}) and `drop` applied."""
    args = {**ARGS, **{f"--{k.replace('_', '-')}": v for k, v in changes.items()}}
    argv = [x for k, v in args.items() if k not in drop for x in (k, v)]
    return run(*argv, "--outdir", str(tmp_path / "out"), env=env)


def results(outdir):
    return pd.read_csv(outdir / "tables" / "model_results.tsv", sep="\t")


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


def demo_variant(tmp_path, edit):
    """The demo TSV after `edit(frame)`, written to tmp_path; returns its path."""
    path = tmp_path / "variant.tsv"
    edit(pd.read_csv(DEMO, sep="\t", dtype=str)).to_csv(path, sep="\t", index=False)
    return str(path)


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
    assert result["n_features"] == 30 and result["n_fitted"] == 30
    assert set(result["significant_features"]) == TIME_TRUTH
    assert result["singular_features"] == ["F30"]
    assert result["n_not_converged"] == 0 and result["failed_features"] == []
    assert result["seed"] == 0 and result["term"] == "time3"
    assert 3 <= len(result["chat_summary_lines"]) <= 4
    assert result["chat_summary_lines"][0].startswith("Context: synthetic demo")
    table = results(demo_out)
    assert set(table["feature"].head(6)) == TIME_TRUTH
    assert (table["estimate"].head(6).between(1.3, 1.7)).all()
    assert (
        table.loc[table["feature"] == "F30", "fit_message"]
        .str.contains("singular")
        .all()
    )


@needs_r
def test_numeric_time_golden(tmp_path):
    completed = run_with(tmp_path, time_type="numeric", term="time")
    assert completed.returncode == 0, completed.stderr
    table = results(tmp_path / "out")
    assert set(table.loc[table["significant"], "feature"]) == TIME_TRUTH
    # Planted slope 0.75 per visit.
    assert table["estimate"].head(6).between(0.65, 0.85).all()


def test_gap_is_not_bridged(demo_out):
    """The drawn grid has an empty visit 2 for S24, and its Line2D holds a NaN there."""
    drawn = pd.read_csv(demo_out / "figures" / "top01.csv")
    gap = drawn[drawn["subject"] == GAP_SUBJECT].set_index("x")["value"]
    assert gap.isna().tolist() == [False, True, False]
    report = (demo_out / "report.md").read_text()
    gap_note = "break at a missed visit (subjects with a gap: 1; dropped, no values: 0)"
    assert gap_note in report

    sys.path.insert(0, str(SCRIPT.parent))
    import longitudinal

    rows = pd.DataFrame(
        {
            "subject": ["A", "A", "B", "B", "B"],
            "x": [0.0, 2.0, 0.0, 1.0, 2.0],
            "value": [1.0, 3.0, 1.0, 2.0, 3.0],
        }
    )
    p = {"time_unit": "visit"}
    fig, _ = longitudinal.draw_trajectories(rows, ["1", "2", "3"], p, "v (a.u.)")
    lines = {line.get_gid(): line for line in fig.axes[0].lines if line.get_gid()}
    y = np.asarray(lines["A"].get_ydata(), dtype=float)
    assert np.isnan(y[1]) and not np.isnan(y[[0, 2]]).any()
    assert not np.isnan(np.asarray(lines["B"].get_ydata(), dtype=float)).any()


def test_drops_per_time(tmp_path):
    def blank(d):
        d.loc[
            (d["feature"] == "F20") & (d["subject"] == "S01") & (d["time"] == "3"),
            "value",
        ] = ""
        return d

    completed = run_with(tmp_path, input=demo_variant(tmp_path, blank))
    assert completed.returncode == 0, completed.stderr
    result = json.loads((tmp_path / "out" / "result.json").read_text())
    assert result["n_rows_dropped_by_time"] == {"1": 0, "2": 0, "3": 1}


def test_rank_deficient(tmp_path):
    """A covariate that copies time is aliased with it: the run stops before R."""
    path = demo_variant(tmp_path, lambda d: d.assign(batch="b" + d["time"]))
    completed = run_with(
        tmp_path, input=path, formula="value ~ time + batch + (1|subject)"
    )
    assert is_error(completed, "design is rank-deficient")


def test_output_contract(demo_out):
    paths = contract_paths()
    assert len(paths) == 24
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
    assert "lmerTest" in (demo_out / "reproducibility" / "environment.txt").read_text()


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
    completed = run("--demo", "--term", "time2", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --term")


@needs_r
def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"padj_max": "1.5"}, "--padj-max must be in (0, 1], got 1.5"),
        (
            {"time_type": "numeric", "reference_time": "V1"},
            "--reference-time must be a number for --time-type numeric, got V1",
        ),
        ({"formula": "value ~ time + (time|subject)"}, "--formula needs the random"),
        ({"formula": "value ~ time + log(age) + (1|subject)"}, "--formula must be"),
        ({"formula": "value ~ age + (1|subject)"}, "--formula must have time"),
        ({"input": "/no/such.tsv"}, "input file not found: /no/such.tsv"),
        (
            {"formula": "value ~ time + bmi + (1|subject)"},
            "input column not found: bmi",
        ),
        ({"reference_time": "9"}, "reference time '9' not found in column 'time'"),
    ],
)
def test_argument_errors(tmp_path, changes, message):
    assert is_error(run_with(tmp_path, **changes), message)


def test_not_numeric(tmp_path):
    path = demo_variant(tmp_path, lambda d: d.assign(time="v" + d["time"]))
    completed = run_with(tmp_path, input=path, time_type="numeric", term="time")
    assert is_error(completed, "column 'time' must be numeric")


def test_duplicate_rows(tmp_path):
    path = demo_variant(tmp_path, lambda d: pd.concat([d, d.head(2)]))
    completed = run_with(tmp_path, input=path)
    assert is_error(completed, "4 rows repeat a (feature, subject, time) combination")


def test_no_rscript(tmp_path):
    completed = run_with(tmp_path, env={"RSCRIPT": "", "PATH": str(tmp_path)})
    assert is_error(completed, "Rscript not found; set RSCRIPT or put Rscript on PATH")


def test_r_failed(tmp_path):
    completed = run_with(tmp_path, env={"RSCRIPT": shutil.which("false")})
    assert is_error(completed, "lmer.R failed with exit 1")


@needs_r
def test_missing_r_package(tmp_path):
    # An empty user library hides lme4, which is not in R's own library; ~/.Renviron
    # would set R_LIBS_USER back, so it is skipped too.
    env = {"R_LIBS_USER": str(tmp_path), "R_LIBS": "", "R_ENVIRON_USER": "/dev/null"}
    completed = run_with(tmp_path, env=env)
    assert is_error(completed, "missing R package: lme4")


@needs_r
def test_term_absent(tmp_path):
    completed = run_with(tmp_path, term="timeV3")
    assert is_error(
        completed,
        "--term 'timeV3' is not a fixed-effect coefficient; "
        "choose one of: (Intercept), time2, time3, age",
    )


@needs_r
def test_no_tests(tmp_path):
    # Each subject seen once: lmer refuses a grouping factor with one row per level.
    path = demo_variant(tmp_path, lambda d: d.assign(subject=d["subject"] + d["time"]))
    completed = run_with(tmp_path, input=path)
    assert is_error(completed, "no feature could be fitted; first lmer message:")


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
