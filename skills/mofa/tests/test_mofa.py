"""CLI tests for scripts/mofa.py: demo goldens, output contract, fail strings.

Run from skills/mofa with pytest; see references/maintainer.md for the exact command.
"""

import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "mofa.py"
EXAMPLES = SKILL_DIR / "examples"
TEMPLATES = SKILL_DIR.parents[1] / "templates"

# Planted truth from examples/make_demo.py: A in all views, B in rna+adt, C in tf.
ACTIVE = {"Factor1": ["adt", "rna", "tf"], "Factor2": ["adt", "rna"], "Factor3": ["tf"]}
# R2 (%) of the demo run. Reruns on one machine were bit-identical; the tolerance
# allows for another BLAS build.
R2 = {
    ("Factor1", "adt"): 27.12,
    ("Factor1", "rna"): 33.95,
    ("Factor1", "tf"): 33.52,
    ("Factor2", "adt"): 25.88,
    ("Factor2", "rna"): 13.95,
    ("Factor3", "tf"): 25.35,
}
R2_TOL = 0.5

ARGS = {
    "--context": "synthetic demo, healthy vs disease",
    "--metadata": str(EXAMPLES / "demo_meta.tsv"),
    "--group-col": "group",
    "--likelihoods": "rna=gaussian,adt=gaussian,tf=gaussian",
    "--factors": "6",
    "--scale-views": "yes",
    "--convergence": "fast",
    "--min-var-explained": "5",
    "--n-top-features": "10",
    "--seed": "0",
    "--reference": "synthetic demo, no genome build",
}
VIEWS = {v: str(EXAMPLES / f"demo_{v}.tsv") for v in ("rna", "adt", "tf")}


def run(*argv):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *argv],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )


def run_with(tmp_path, drop=(), views=None, **changes):
    """Demo arguments with `changes` ({"seed": "-1"}), `views` and `drop` applied."""
    args = {**ARGS, **{f"--{k.replace('_', '-')}": v for k, v in changes.items()}}
    argv = [x for k, v in args.items() if k not in drop for x in (k, v)]
    views = VIEWS if views is None else views
    argv += [x for name, path in views.items() for x in ("--view", f"{name}={path}")]
    return run(*argv, "--outdir", str(tmp_path / "out"))


def write_view(tmp_path, view, edit):
    """A copy of a demo view after `edit(table)`; returns its path."""
    table = pd.read_csv(VIEWS[view], sep="\t", index_col=0)
    path = tmp_path / f"{view}.tsv"
    edit(table).to_csv(path, sep="\t")
    return str(path)


def write_meta(tmp_path, edit):
    meta = pd.read_csv(EXAMPLES / "demo_meta.tsv", sep="\t")
    path = tmp_path / "meta.tsv"
    edit(meta).to_csv(path, sep="\t", index=False)
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
    assert result["n_samples"] == 60 and result["n_factors"] == 6
    assert result["views"] == {"adt": 40, "rna": 200, "tf": 80}
    assert result["n_active"] == 3
    assert result["active_factors"] == ACTIVE
    assert result["seed"] == 0 and result["device"] == "cpu"
    assert "OMP_NUM_THREADS=1" in result["threads"]
    assert len(result["chat_summary_lines"]) == 4
    assert result["chat_summary_lines"][0].startswith("Context: synthetic demo")
    assert result["n_shared_ids"] == 60
    assert result["reference"] == "synthetic demo, no genome build"
    # Every metadata sample is in a view, so none is left out.
    assert result["dropped_per_group"] == {}
    r2 = pd.read_csv(demo_out / "tables" / "variance_explained.tsv", sep="\t")
    got = r2.set_index(["factor", "view"])["r2"]
    for key, value in R2.items():
        assert abs(got[key] - value) < R2_TOL, key
    # Flagged, not dropped: all six factors stay in every table.
    assert r2["factor"].nunique() == 6
    assert r2.groupby("factor")["active"].first().sum() == 3


def test_demo_group_association(demo_out):
    assoc = pd.read_csv(demo_out / "tables" / "factor_vs_metadata.tsv", sep="\t")
    assert len(assoc) == 12
    group = assoc[assoc["covariate"] == "group"].set_index("factor")
    assert (group["test"] == "kruskal-wallis").all()
    # Only the planted factor A (Factor1) separates the groups; age is null.
    assert group["padj"].idxmin() == "Factor1" and group.loc["Factor1", "padj"] < 0.01
    assert (assoc.loc[assoc["covariate"] == "age", "padj"] > 0.05).all()
    assert (assoc.loc[assoc["covariate"] == "age", "test"] == "spearman").all()


def test_output_tables(demo_out):
    factors = pd.read_csv(demo_out / "tables" / "factors.tsv", sep="\t")
    assert factors.shape == (60, 7) and factors["sample"].is_monotonic_increasing
    weights = pd.read_csv(demo_out / "tables" / "weights.tsv", sep="\t")
    assert list(weights.columns) == ["view", "feature", "factor", "weight"]
    assert len(weights) == (200 + 40 + 80) * 6


def test_output_contract(demo_out):
    paths = contract_paths()
    assert len(paths) == 34
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
    assert "--seed 0" in commands and commands.count("--view ") == 3
    assert "--reference 'synthetic demo, no genome build'" in commands


def test_commands_sh_reproduces_the_run(demo_out, tmp_path):
    """The resolved command in commands.sh, with a new --outdir, gives the same fit."""
    line = (demo_out / "reproducibility" / "commands.sh").read_text().splitlines()[-1]
    argv = shlex.split(line)
    argv[argv.index("--outdir") + 1] = str(tmp_path / "again")
    completed = subprocess.run(argv, capture_output=True, text=True, check=False)
    assert completed.returncode == 0, completed.stderr
    for name in ("factors.tsv", "variance_explained.tsv", "factor_vs_metadata.tsv"):
        assert (tmp_path / "again" / "tables" / name).read_bytes() == (
            demo_out / "tables" / name
        ).read_bytes(), name


def test_captions_say_how(demo_out):
    report = (demo_out / "report.md").read_text()
    assert "seed 0; no downsampling" in report
    assert "rows fixed in factor order" in report


def test_left_out_samples_counted_per_group(tmp_path):
    extra = pd.DataFrame({"sample": ["X1", "X2"], "group": ["disease", "disease"]})
    meta = write_meta(tmp_path, lambda m: pd.concat([m, extra]))
    completed = run_with(tmp_path, metadata=meta)
    assert completed.returncode == 0, completed.stderr
    result = json.loads((tmp_path / "out" / "result.json").read_text())
    assert result["dropped_per_group"] == {"disease": 2}


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


def test_input_order_does_not_matter(demo_out, tmp_path):
    shuffled = {
        **VIEWS,
        "rna": write_view(
            tmp_path, "rna", lambda t: t.sample(frac=1, random_state=1).iloc[:, ::-1]
        ),
    }
    completed = run_with(tmp_path, views=dict(reversed(shuffled.items())))
    assert completed.returncode == 0, completed.stderr
    for name in ("factors.tsv", "variance_explained.tsv"):
        assert (tmp_path / "out" / "tables" / name).read_bytes() == (
            demo_out / "tables" / name
        ).read_bytes(), name


def test_missing_argument_fails_without_demo(tmp_path):
    assert is_error(
        run_with(tmp_path, drop=["--seed"]), "missing required arguments: --seed"
    )
    assert is_error(
        run_with(tmp_path, drop=["--reference"]),
        "missing required arguments: --reference",
    )
    assert is_error(
        run_with(tmp_path, drop=["--context"]),
        "missing required arguments: --context",
    )


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--factors", "3", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --factors")


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"factors": "1"}, "--factors must be an integer >= 2, got 1"),
        ({"min_var_explained": "100"}, "--min-var-explained must be a percentage"),
        ({"n_top_features": "0"}, "--n-top-features must be >= 1, got 0"),
        ({"seed": "-1"}, "--seed must be an integer >= 0, got -1"),
        (
            {"likelihoods": "rna=gaussian,adt=gaussian"},
            "--likelihoods must give one of gaussian/bernoulli/poisson for each "
            "view (name=likelihood,...); missing for ['tf']",
        ),
        (
            {"likelihoods": "rna=gaussian,adt=normal,tf=gaussian"},
            "--likelihoods must give one of gaussian/bernoulli/poisson for each "
            "view (name=likelihood,...); unknown likelihood ['normal']",
        ),
        ({"metadata": "/no/such.tsv"}, "input file not found: /no/such.tsv"),
        ({"group_col": "batch"}, "column 'batch' not found in"),
        ({"factors": "60"}, "--factors (60) must be below the number of samples (60)"),
    ],
)
def test_argument_errors(tmp_path, changes, message):
    assert is_error(run_with(tmp_path, **changes), message)


@pytest.mark.parametrize(
    ("views", "message"),
    [
        (["rna"], "--view must be name=path, got 'rna'"),
        ([f"rna={VIEWS['rna']}"] * 2, "view 'rna' is given more than once"),
    ],
)
def test_view_flag_errors(tmp_path, views, message):
    argv = [x for k, v in ARGS.items() for x in (k, v)]
    argv += [x for view in views for x in ("--view", view)]
    assert is_error(run(*argv, "--outdir", str(tmp_path / "out")), message)


def test_too_many_views(tmp_path):
    views = {f"v{i}": VIEWS["rna"] for i in range(8)}
    lik = ",".join(f"v{i}=gaussian" for i in range(8))
    completed = run_with(tmp_path, views=views, likelihoods=lik)
    assert is_error(completed, "--view has 8 levels; the fixed palette holds 7")


@pytest.mark.parametrize(
    ("edit", "message"),
    [
        (
            lambda t: t.astype(object).assign(rna_001="x"),
            "view 'rna' has non-numeric values in 60 cells",
        ),
        (
            lambda t: pd.concat([t, t.iloc[:1]]),
            "duplicate sample IDs in view 'rna': S01",
        ),
        (
            lambda t: t.rename(columns={"rna_002": "rna_001"}),
            "duplicate feature names in view 'rna': rna_001",
        ),
        (
            lambda t: t.rename(columns={"rna_001": "adt_001"}),
            "1 feature names appear in more than one view, e.g. 'adt_001'",
        ),
        (
            lambda t: t.rename(index=lambda s: "X" + s),
            "inputs share no sample ids; check they come from one study and reference",
        ),
        (
            lambda t: t.rename(index={"S01": "S99"}),
            "1 samples in view 'rna' are not in the metadata, e.g. 'S99'",
        ),
    ],
)
def test_view_errors(tmp_path, edit, message):
    views = {**VIEWS, "rna": write_view(tmp_path, "rna", edit)}
    assert is_error(run_with(tmp_path, views=views), message)


@pytest.mark.parametrize(
    ("likelihood", "message"),
    [
        ("bernoulli", "view 'adt' is bernoulli but holds values other than 0 and 1"),
        ("poisson", "view 'adt' is poisson but holds negative or non-integer values"),
    ],
)
def test_likelihood_value_errors(tmp_path, likelihood, message):
    lik = f"rna=gaussian,adt={likelihood},tf=gaussian"
    assert is_error(run_with(tmp_path, likelihoods=lik), message)


@pytest.mark.parametrize(
    ("edit", "message"),
    [
        (lambda m: m.drop(columns="sample"), "column 'sample' not found in"),
        (lambda m: pd.concat([m, m.iloc[:1]]), "duplicate sample IDs in"),
        (
            lambda m: m.assign(group=m["group"].mask(m.index < 2)),
            "2 modelled samples have no value in 'group'",
        ),
        (
            lambda m: m.assign(group="healthy"),
            "'group' needs at least 2 levels, got ['healthy']",
        ),
        (
            lambda m: m.assign(group=[f"g{i % 8}" for i in m.index]),
            "group has 8 levels; the fixed palette holds 7",
        ),
    ],
)
def test_metadata_errors(tmp_path, edit, message):
    meta = write_meta(tmp_path, edit)
    assert is_error(run_with(tmp_path, metadata=meta), message)


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
