"""CLI tests for scripts/sc_embed.py: demo goldens, output contract, fail strings.

Run from skills/sc-embed with pytest;
see references/maintainer.md for the exact command.
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
SCRIPT = SKILL_DIR / "scripts" / "sc_embed.py"
DEMO = SKILL_DIR / "examples" / "demo.h5ad"
TEMPLATES = SKILL_DIR.parents[1] / "templates"

# Golden demo values, computed once with scanpy 1.12.4, harmonypy 2.1.0, igraph 1.0.0.
HARMONY_SIZES = {"0": 200, "1": 200, "2": 200}
NONE_SIZES = {"0": 94, "1": 106, "2": 95, "3": 105, "4": 94, "5": 106}

ARGS = {
    "--adata": str(DEMO),
    "--target-sum": "10000",
    "--n-top-genes": "100",
    "--flavor": "seurat",
    "--batch-key": "batch",
    "--n-comps": "20",
    "--integration": "harmony",
    "--n-neighbors": "15",
    "--n-pcs": "10",
    "--min-dist": "0.3",
    "--resolution": "0.3",
    "--seed": "0",
    "--reference": "synthetic demo genes, no genome build",
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
    """Demo arguments with `changes` ({"n_pcs": "30"}) and `drop` applied."""
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


def purity(outdir, truth):
    """Fraction of cells whose cluster's majority `truth` label matches their own."""
    obs = ad.read_h5ad(outdir / "data" / "embedded.h5ad").obs
    table = pd.crosstab(obs["leiden"], obs[truth])
    return table.max(axis=1).sum() / len(obs)


@pytest.fixture(scope="module")
def demo_out(tmp_path_factory):
    out = tmp_path_factory.mktemp("demo") / "out"
    completed = run("--demo", "--outdir", str(out))
    assert completed.returncode == 0, completed.stderr
    return out


@pytest.fixture(scope="module")
def none_out(tmp_path_factory):
    out = tmp_path_factory.mktemp("none") / "out"
    completed = run("--demo", "--integration", "none", "--outdir", str(out))
    assert completed.returncode == 0, completed.stderr
    return out


def test_demo_golden_values(demo_out):
    result = json.loads((demo_out / "result.json").read_text())
    assert result["n_cells"] == 600 and result["n_hvg"] == 100
    assert result["cluster_sizes"] == HARMONY_SIZES
    assert result["use_rep"] == "X_pca_harmony" and result["seed"] == 0
    assert 2 <= len(result["chat_summary_lines"]) <= 3
    assert result["reference"] == "synthetic demo genes, no genome build"
    report = (demo_out / "report.md").read_text()
    how = "10 PCs of X_pca_harmony, 15 neighbours, min_dist 0.3, seed 0"
    assert report.count(how) == 2  # both UMAP captions say how
    # Harmony recovers the 3 planted groups; each cluster mixes both batches.
    assert purity(demo_out, "true_group") == 1.0
    mix = pd.read_csv(demo_out / "tables" / "batch_by_cluster.tsv", sep="\t")
    assert mix["frac_of_cluster"].eq(0.5).all()


def test_no_integration_splits_by_batch(none_out):
    result = json.loads((none_out / "result.json").read_text())
    assert result["cluster_sizes"] == NONE_SIZES
    assert result["use_rep"] == "X_pca"
    assert purity(none_out, "batch") > 0.95


def test_rerun_is_identical(demo_out, tmp_path):
    out = tmp_path / "again"
    assert run("--demo", "--outdir", str(out)).returncode == 0
    for name in ["cluster_sizes.tsv", "batch_by_cluster.tsv"]:
        assert (out / "tables" / name).read_bytes() == (
            demo_out / "tables" / name
        ).read_bytes()
    first = ad.read_h5ad(demo_out / "data" / "embedded.h5ad")
    second = ad.read_h5ad(out / "data" / "embedded.h5ad")
    assert (first.obsm["X_umap"] == second.obsm["X_umap"]).all()


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
    embedded = ad.read_h5ad(demo_out / "data" / "embedded.h5ad")
    assert {"X_pca", "X_pca_harmony", "X_umap"} <= set(embedded.obsm)
    assert "counts" in embedded.layers and embedded.n_vars == 300


def test_missing_argument_fails_without_demo(tmp_path):
    assert is_error(
        run_with(tmp_path, drop=["--seed"]), "missing required arguments: --seed"
    )


def test_missing_reference_fails(tmp_path):
    completed = run_with(tmp_path, drop=["--reference"])
    assert is_error(completed, "missing required arguments: --reference")


@pytest.mark.parametrize("axis", ["var", "obs"])
def test_duplicate_ids(tmp_path, axis):
    def duplicate(adata):
        names = list(getattr(adata, f"{axis}_names"))
        names[1] = names[0]
        setattr(adata, f"{axis}_names", names)

    what = {"var": "gene names", "obs": "cell ids"}[axis]
    completed = run_with(tmp_path, adata=edited_demo(tmp_path, duplicate))
    assert is_error(completed, f"duplicate {what}: ")


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--seed", "1", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --seed")


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"n_pcs": "30"}, "--n-pcs must be in [2, --n-comps], got 30"),
        ({"n_top_genes": "400"}, "--n-top-genes must be <= the 300 genes, got 400"),
        ({"n_comps": "150"}, "--n-comps must be <= min(cells, --n-top-genes) - 1 = 99"),
        (
            {"batch_key": "none"},
            "--integration harmony needs a --batch-key, got none",
        ),
        ({"adata": "/no/such.h5ad"}, "input file not found: /no/such.h5ad"),
        ({"batch_key": "donor"}, "obs column not found: donor"),
    ],
)
def test_input_errors(tmp_path, changes, message):
    assert is_error(run_with(tmp_path, **changes), message)


def edited_demo(tmp_path, edit):
    adata = ad.read_h5ad(DEMO)
    edit(adata)
    adata.write_h5ad(tmp_path / "bad.h5ad")
    return str(tmp_path / "bad.h5ad")


def test_non_integer_counts(tmp_path):
    def halve(adata):
        adata.X = adata.X * 0.5

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, halve))
    assert is_error(completed, "counts in X must be non-negative integers")


def test_batch_na(tmp_path):
    def blank(adata):
        adata.obs["batch"] = adata.obs["batch"].astype(object)
        adata.obs.iloc[0, adata.obs.columns.get_loc("batch")] = None

    completed = run_with(tmp_path, adata=edited_demo(tmp_path, blank))
    assert is_error(completed, "obs column 'batch' has missing values")


def test_one_batch(tmp_path):
    completed = run_with(
        tmp_path,
        batch_key="true_group",
        adata=edited_demo(tmp_path, lambda a: a.obs.__setitem__("true_group", "g1")),
    )
    assert is_error(completed, "obs column 'true_group' has one level")


def test_no_batch_runs(tmp_path):
    completed = run_with(tmp_path, batch_key="none", integration="none")
    assert completed.returncode == 0, completed.stderr
    mix = pd.read_csv(tmp_path / "out" / "tables" / "batch_by_cluster.tsv", sep="\t")
    assert set(mix["batch"]) == {"none"}


@pytest.mark.parametrize("flavor", ["cell_ranger", "seurat_v3"])
def test_other_hvg_flavors(tmp_path, flavor):
    completed = run_with(tmp_path, flavor=flavor)
    assert completed.returncode == 0, completed.stderr
    result = json.loads((tmp_path / "out" / "result.json").read_text())
    assert result["n_hvg"] == 100


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
