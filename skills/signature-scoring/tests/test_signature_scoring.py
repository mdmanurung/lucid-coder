"""CLI tests for signature_scoring.py: demo goldens, output contract, fail strings.

Run from skills/signature-scoring with pytest;
see references/maintainer.md "Running and testing".
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
SCRIPT = SKILL_DIR / "scripts" / "signature_scoring.py"
DEMO = SKILL_DIR / "examples" / "demo.h5ad"
DEMO_SETS = SKILL_DIR / "examples" / "demo_gene_sets.tsv"
TEMPLATES = SKILL_DIR.parents[1] / "templates"

ARGS = {
    "--adata": str(DEMO),
    "--gene-sets": str(DEMO_SETS),
    "--sample-col": "sample",
    "--group-col": "group",
    "--method": "score_genes",
    "--ctrl-size": "50",
    "--n-bins": "25",
    "--min-genes": "5",
    "--n-random": "50",
    "--reference": "synthetic demo genes, no genome build",
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
    """Demo arguments with `changes` ({"min_genes": "0"}) and `drop` applied."""
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
    path = tmp_path / "edited.h5ad"
    adata.write_h5ad(path)
    return str(path)


def gene_sets(tmp_path, text):
    path = tmp_path / "sets.tsv"
    path.write_text(text)
    return str(path)


def read(out, name):
    return pd.read_csv(out / "tables" / name, sep="\t", keep_default_na=False)


@pytest.fixture(scope="module")
def demo_out(tmp_path_factory):
    out = tmp_path_factory.mktemp("demo") / "out"
    completed = run("--demo", "--outdir", str(out))
    assert completed.returncode == 0, completed.stderr
    return out


def test_demo_golden_values(demo_out):
    result = json.loads((demo_out / "result.json").read_text())
    assert (result["n_cells"], result["n_samples"], result["n_sets"]) == (400, 8, 2)
    assert result["n_found"] == {"sig_null": 10, "sig_up": 10}
    assert result["missing_genes"] == {"sig_null": ["GMISSING"]}
    assert result["beats_random"] == {"sig_null": False, "sig_up": True}
    assert result["p_empirical"]["sig_up"] == pytest.approx(1 / 51)
    assert result["seed"] == 0 and result["method"] == "score_genes"
    assert (result["n_matched"], result["n_in_resource"]) == (20, 21)
    assert result["reference"] == "synthetic demo genes, no genome build"
    # 4 samples per group is under FEW_POINTS: points and a median, no box.
    report = (demo_out / "report.md").read_text()
    assert "Dot plot of the per-sample mean sig_up score by group with a bar" in report
    found = read(demo_out, "genes_found.tsv")
    assert found.to_dict("list") == {
        "set": ["sig_null", "sig_up"],
        "n_listed": [11, 10],
        "n_found": [10, 10],
        "missing": ["GMISSING", ""],
    }
    scores = read(demo_out, "sample_scores.tsv")
    assert list(scores.columns) == ["feature", "subject", "group", "value", "n_cells"]
    assert len(scores) == 16 and set(scores["n_cells"]) == {50}
    up = scores[scores["feature"] == "sig_up"].groupby("group")["value"]
    assert up.min()["disease"] > up.max()["healthy"]
    assert read(demo_out, "cell_scores.tsv").shape == (400, 5)


def test_tables_sorted(demo_out):
    scores = read(demo_out, "sample_scores.tsv")
    assert scores.index.equals(scores.sort_values(["feature", "subject"]).index)
    cells = read(demo_out, "cell_scores.tsv")
    assert cells["cell"].is_monotonic_increasing
    assert read(demo_out, "random_control.tsv")["set"].is_monotonic_increasing


def test_aucell_engine(tmp_path):
    completed = run_with(
        tmp_path, drop=["--ctrl-size", "--n-bins"], method="aucell", n_up="30"
    )
    assert completed.returncode == 0, completed.stderr
    result = json.loads((tmp_path / "out" / "result.json").read_text())
    assert result["beats_random"] == {"sig_null": False, "sig_up": True}
    assert result["missing_genes"] == {"sig_null": ["GMISSING"]}


def test_output_contract(demo_out):
    paths = contract_paths()
    assert len(paths) == 25
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
    scored = ad.read_h5ad(demo_out / "data" / "scored.h5ad")
    assert {"score_genes_sig_up", "score_genes_sig_null"} <= set(scored.obs.columns)


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
    assert check.returncode == 0, check.stdout + check.stderr


def test_no_umap_skips_umap_figures(tmp_path):
    def drop_umap(adata):
        del adata.obsm["X_umap"]

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, drop_umap))
    assert completed.returncode == 0, completed.stderr
    figures = {p.name for p in (tmp_path / "out" / "figures").iterdir()}
    assert "random_control.pdf" in figures
    assert not any(name.startswith("score_umap") for name in figures)


def test_missing_argument_fails_without_demo(tmp_path):
    completed = run_with(tmp_path, drop=["--n-random"])
    assert is_error(completed, "missing required arguments: --n-random")


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--seed", "1", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --seed")


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


@pytest.mark.parametrize(
    ("changes", "drop", "message"),
    [
        ({"min_genes": "0"}, [], "--min-genes must be an integer >= 1, got 0"),
        ({"n_up": "30"}, [], "--method score_genes does not use --n-up"),
        (
            {"method": "aucell"},
            ["--ctrl-size", "--n-bins"],
            "--method aucell needs --n-up",
        ),
        (
            {"method": "aucell", "n_up": "301"},
            ["--ctrl-size", "--n-bins"],
            "--n-up must be <= 300, got 301",
        ),
        ({"adata": "/no/such.h5ad"}, [], "input file not found: /no/such.h5ad"),
        ({"group_col": "condition"}, [], "obs column not found: condition"),
        (
            {"min_genes": "11"},
            [],
            "sets with fewer than 11 genes in the data: sig_null (10), sig_up (10)",
        ),
        (
            {"group_col": "sample"},
            [],
            "--sample-col and --group-col must differ, both are sample",
        ),
    ],
)
def test_input_errors(tmp_path, changes, drop, message):
    assert is_error(run_with(tmp_path, drop=drop, **changes), message)


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("name\tgene\nA\tG000\n", "--gene-sets needs columns set, gene; found name"),
        ("set\tgene\nbad/name\tG000\n", "set names may use only letters"),
    ],
)
def test_gene_set_errors(tmp_path, text, message):
    completed = run_with(tmp_path, gene_sets=gene_sets(tmp_path, text))
    assert is_error(completed, message)


def test_raw_counts(tmp_path):
    def counts(adata):
        adata.X = adata.layers["counts"]

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, counts))
    assert is_error(completed, "X does not look log-normalised: all integers")


def test_negative_values(tmp_path):
    def scale(adata):
        adata.X = adata.X.toarray() - 1

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, scale))
    assert is_error(completed, "X does not look log-normalised: negative values")


def test_split_sample(tmp_path):
    def split(adata):
        adata.obs["group"] = adata.obs["group"].astype(str)
        adata.obs.iloc[0, adata.obs.columns.get_loc("group")] = "disease"

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, split))
    assert is_error(completed, "samples spanning several values of group: H1")


def test_many_groups(tmp_path):
    def per_sample(adata):
        adata.obs["group"] = adata.obs["sample"]

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, per_sample))
    assert is_error(completed, "--group-col group has 8 groups; box plot colours: 6")


def test_one_group(tmp_path):
    def one_group(adata):
        adata.obs["group"] = "healthy"

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, one_group))
    assert is_error(completed, "--group-col group needs at least 2 groups, found 1")


@pytest.mark.parametrize(("axis", "what"), [("var", "gene names"), ("obs", "cell ids")])
def test_duplicate_ids(tmp_path, axis, what):
    def duplicate(adata):
        names = list(getattr(adata, f"{axis}_names"))
        names[1] = names[0]
        setattr(adata, f"{axis}_names", names)

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, duplicate))
    assert is_error(completed, f"duplicate {what}: ")


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
