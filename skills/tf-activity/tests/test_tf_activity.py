"""CLI tests for scripts/tf_activity.py: demo goldens, output contract, fail strings.

Run from skills/tf-activity with pytest;
see references/maintainer.md "Running and testing".
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
from scipy import stats

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "tf_activity.py"
EXAMPLES = SKILL_DIR / "examples"
TEMPLATES = SKILL_DIR.parents[1] / "templates"

COMMON = {
    "--network": str(EXAMPLES / "demo_network.tsv"),
    "--method": "ulm",
    "--tmin": "5",
    "--n-top": "2",
    "--reference": "synthetic demo genes, no genome build",
    "--context": "synthetic demo",
    "--seed": "0",
}
ARGS = {
    "cells": {
        "--input": "cells",
        **COMMON,
        "--adata": str(EXAMPLES / "demo.h5ad"),
        "--groupby": "group",
        "--rank-method": "t-test_overestim_var",
        "--limit": "3",
    },
    "contrast": {
        "--input": "contrast",
        **COMMON,
        "--de-table": str(EXAMPLES / "demo_de.tsv"),
        "--stat-col": "stat",
        "--contrast": "g1_vs_g2",
    },
}


def run(*argv):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *argv],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )


def run_with(tmp_path, mode="cells", drop=(), **changes):
    """Mode arguments with `changes` ({"tmin": "0"}) and `drop` applied."""
    args = {**ARGS[mode], **{f"--{k.replace('_', '-')}": v for k, v in changes.items()}}
    argv = [x for k, v in args.items() if k not in drop for x in (k, v)]
    return run(*argv, "--outdir", str(tmp_path / "out"))


def result(out):
    return json.loads((out / "result.json").read_text())


def contract_paths(mode):
    """Output-tree paths for one mode; a `(cells)`/`(contrast)` note limits a file."""
    text = (SKILL_DIR / "SKILL.md").read_text()
    tree = text.split("## Output structure")[1].split("```")[1]
    paths, stack = [], []
    for line in tree.splitlines():
        match = re.search(r"[├└]── (\S+)", line)
        if not match:
            continue
        depth = match.start() // 4
        stack[depth:] = [match.group(1)]
        other = "(contrast)" if mode == "cells" else "(cells)"
        if not match.group(1).endswith("/") and other not in line:
            paths.append("".join(stack))
    return paths


def is_error(completed, message):
    return (
        completed.returncode == 2
        and f"ERROR: {message}" in completed.stderr
        and "Traceback" not in completed.stderr
    )


def written(tmp_path, name, frame):
    path = tmp_path / name
    frame.to_csv(path, sep="\t", index=False)
    return str(path)


@pytest.fixture(scope="module")
def demo_out(tmp_path_factory):
    out = tmp_path_factory.mktemp("demo") / "out"
    completed = run("--demo", "--outdir", str(out))
    assert completed.returncode == 0, completed.stderr
    return out


@pytest.fixture(scope="module")
def contrast_out(tmp_path_factory):
    out = tmp_path_factory.mktemp("contrast") / "out"
    completed = run("--demo", "--input", "contrast", "--outdir", str(out))
    assert completed.returncode == 0, completed.stderr
    return out


def test_demo_golden_values(demo_out):
    res = result(demo_out)
    assert (res["n_cells"], res["n_sources_scored"], res["n_sources_dropped"]) == (
        180,
        6,
        2,
    )
    assert [res["top_sources"][g][0] for g in ("g1", "g2")] == ["TF1", "TF2"]
    assert (res["n_matched"], res["n_in_resource"]) == (93, 99)
    assert res["reference"] == "synthetic demo genes, no genome build"
    table = pd.read_csv(demo_out / "tables" / "activity.tsv", sep="\t")
    table = table[table["status"] == "tested"]
    assert len(table) == 18 and set(table["n_targets"]) == {15}
    hits = table[(table["padj"] < 0.05) & (table["stat"] > 0)]
    assert sorted(zip(hits["source"], hits["group"], strict=True)) == [
        ("TF1", "g1"),
        ("TF2", "g2"),
    ]
    dropped = pd.read_csv(demo_out / "tables" / "sources_dropped.tsv", sep="\t")
    assert dropped.to_dict("list") == {"source": ["TF7", "TF8"], "n_targets": [3, 0]}


@pytest.mark.parametrize("mode", ["cells", "contrast"])
def test_status_rows_for_dropped_sources(demo_out, contrast_out, mode):
    out = demo_out if mode == "cells" else contrast_out
    table = pd.read_csv(out / "tables" / "activity.tsv", sep="\t")
    untested = table[table["status"] == "not_tested"]
    assert untested["source"].tolist() == ["TF7", "TF8"]
    assert untested[["score", "pvalue", "padj"]].isna().all().all()
    assert result(out)["status_counts"] == {"not_tested": 2, "tested": 6}


def test_context_leads_summary(demo_out):
    res = result(demo_out)
    assert res["context"].startswith("synthetic demo")
    assert res["chat_summary_lines"][0].startswith("Context: ")
    assert "Context: synthetic demo" in (demo_out / "report.md").read_text()


def test_activity_sorted_by_group_then_stat(demo_out):
    table = pd.read_csv(demo_out / "tables" / "activity.tsv", sep="\t")
    resorted = table.sort_values(
        ["group", "stat", "source"], ascending=[True, False, True]
    )
    assert table.index.equals(resorted.index)


def test_contrast_golden_values(contrast_out):
    res = result(contrast_out)
    assert set(res["top_sources"]) == {"TF1", "TF2"}
    assert (res["n_genes"], res["n_sources_scored"]) == (300, 6)
    assert (res["n_matched"], res["n_in_resource"]) == (93, 99)
    table = pd.read_csv(contrast_out / "tables" / "activity.tsv", sep="\t")
    table = table[table["status"] == "tested"]
    assert table["score"].is_monotonic_decreasing
    assert table["source"].iloc[[0, -1]].tolist() == ["TF1", "TF2"]


@pytest.mark.parametrize("method", ["ulm", "mlm"])
def test_contrast_padj_is_bh_of_pvalue(tmp_path, method):
    # ulm: raw p rebuilt in-script, padj from decoupler; mlm: the reverse.
    completed = run_with(tmp_path, "contrast", method=method)
    assert completed.returncode == 0, completed.stderr
    table = pd.read_csv(tmp_path / "out" / "tables" / "activity.tsv", sep="\t")
    table = table[table["status"] == "tested"]
    bh = stats.false_discovery_control(table["pvalue"], method="bh")
    np.testing.assert_allclose(table["padj"], bh, rtol=1e-5)
    assert set(result(tmp_path / "out")["top_sources"]) == {"TF1", "TF2"}


def test_mlm_cells_recovers_planted(tmp_path):
    completed = run_with(tmp_path, method="mlm", rank_method="wilcoxon")
    assert completed.returncode == 0, completed.stderr
    top = result(tmp_path / "out")["top_sources"]
    assert [top["g1"][0], top["g2"][0]] == ["TF1", "TF2"]


@pytest.mark.parametrize("mode", ["cells", "contrast"])
def test_output_contract(demo_out, contrast_out, mode):
    out = demo_out if mode == "cells" else contrast_out
    paths = contract_paths(mode)
    assert len(paths) == {"cells": 13, "contrast": 10}[mode]
    for rel in paths:
        assert (out / rel).stat().st_size > 0, rel
    check = subprocess.run(
        ["sha256sum", "-c", "reproducibility/checksums.sha256"],
        cwd=out,
        capture_output=True,
        text=True,
        check=False,
    )
    assert check.returncode == 0, check.stdout + check.stderr
    assert "--seed 0" in (out / "reproducibility" / "commands.sh").read_text()


@pytest.mark.parametrize("arg", ["--tmin", "--context"])
def test_missing_argument_fails_without_demo(tmp_path, arg):
    completed = run_with(tmp_path, drop=[arg])
    assert is_error(completed, f"missing required arguments: {arg}")


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--seed", "1", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --seed")


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


def test_mode_missing(tmp_path):
    completed = run_with(tmp_path, drop=["--groupby", "--limit"])
    assert is_error(completed, "--input cells needs: --groupby --limit")


def test_mode_unused(tmp_path):
    completed = run_with(tmp_path, "contrast", limit="3")
    assert is_error(completed, "--input contrast does not use: --limit")


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"tmin": "0"}, "--tmin must be an integer >= 1, got 0"),
        ({"limit": "0"}, "--limit must be > 0, got 0.0"),
        ({"adata": "/no/such.h5ad"}, "input file not found: /no/such.h5ad"),
        ({"groupby": "leiden"}, "obs column not found: leiden"),
        ({"tmin": "16"}, "no source has at least --tmin 16 targets in the data"),
    ],
)
def test_input_errors(tmp_path, changes, message):
    assert is_error(run_with(tmp_path, **changes), message)


def test_missing_stat_column(tmp_path):
    completed = run_with(tmp_path, "contrast", stat_col="t")
    assert is_error(completed, "DE table column not found: t")


def test_duplicate_genes(tmp_path):
    de = pd.read_csv(EXAMPLES / "demo_de.tsv", sep="\t")
    path = written(tmp_path, "de.tsv", pd.concat([de, de.head(1)]))
    completed = run_with(tmp_path, "contrast", de_table=path)
    assert is_error(completed, "duplicate gene names in the DE table: ")


@pytest.mark.parametrize(
    ("edit", "why"),
    [
        (lambda n: n.drop(columns="weight"), "missing weight"),
        (lambda n: n.assign(weight="up"), "non-numeric or empty weight"),
        (lambda n: pd.concat([n, n.head(1)]), "duplicate source-target rows"),
    ],
)
def test_bad_network(tmp_path, edit, why):
    net = pd.read_csv(EXAMPLES / "demo_network.tsv", sep="\t")
    path = written(tmp_path, "net.tsv", edit(net))
    assert is_error(run_with(tmp_path, network=path), f"network {path}: {why}")


def edited_demo(tmp_path, edit):
    adata = ad.read_h5ad(EXAMPLES / "demo.h5ad")
    edit(adata)
    path = tmp_path / "edited.h5ad"
    adata.write_h5ad(path)
    return str(path)


def test_not_lognorm(tmp_path):
    def counts(adata):
        adata.X = np.round(np.expm1(adata.X))

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, counts))
    assert is_error(completed, "X does not look log-normalised: all integers")


def test_all_zero_targets_count_as_absent(tmp_path):
    # TF3 owns G030-G044; 11 all-zero targets leave 4, under --tmin 5.
    def zero_tf3(adata):
        adata.X[:, 30:41] = 0

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, zero_tf3))
    assert completed.returncode == 0, completed.stderr
    dropped = pd.read_csv(tmp_path / "out" / "tables" / "sources_dropped.tsv", sep="\t")
    assert dropped.to_dict("list") == {
        "source": ["TF3", "TF7", "TF8"],
        "n_targets": [4, 3, 0],
    }


@pytest.mark.parametrize(("axis", "what"), [("var", "gene names"), ("obs", "cell ids")])
def test_duplicate_ids(tmp_path, axis, what):
    def duplicate(adata):
        names = list(getattr(adata, f"{axis}_names"))
        names[1] = names[0]
        setattr(adata, f"{axis}_names", names)

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, duplicate))
    assert is_error(completed, f"duplicate {what}: ")


def test_one_group(tmp_path):
    def one_group(adata):
        adata.obs["group"] = "g1"

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, one_group))
    assert is_error(completed, "--groupby group needs at least 2 groups, found 1")


def test_small_group(tmp_path):
    def singleton(adata):
        adata.obs["group"] = adata.obs["group"].astype(str)
        adata.obs.iloc[0, 0] = "g4"

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, singleton))
    assert is_error(completed, "groups with fewer than 2 cells in group: g4")


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
