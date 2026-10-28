"""CLI tests for scripts/velocity.py: demo goldens, output contract, fail strings.

Run from skills/velocity with pytest;
see references/maintainer.md for the uv line. The demo
run takes about a minute; the dynamical smoke test runs on a 150-cell subset.
"""

import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import pytest
from scipy import stats

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "velocity.py"
DEMO = SKILL_DIR / "examples" / "demo.h5ad"
TEMPLATES = SKILL_DIR.parents[1] / "templates"

# Golden demo values, computed once with scvelo 0.3.4 and cellrank 2.3.3. The
# Spearman rho of velocity pseudotime with the simulated time is 0.9963.
RHO_TRUE_TIME = 0.996
TERMINAL = ["s4"]
TOP_GENES = ["G08", "G20", "G28", "G44"]

ARGS = {
    "--adata": str(DEMO),
    "--gene-id-col": "gene_ids",
    "--mode": "stochastic",
    "--min-shared-counts": "10",
    "--n-top-genes": "30",
    "--n-pcs": "10",
    "--n-neighbors": "30",
    "--groupby": "stage",
    "--terminal-states": "s4",
    "--context": (
        "synthetic demo: scVelo simulation, 500 cells, 4 stages, no real tissue"
    ),
    "--reference": "synthetic demo, no genome build",
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
    """Demo arguments with `changes` ({"n_pcs": "1"}) and `drop` applied."""
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


def broken_demo(tmp_path, edit):
    adata = ad.read_h5ad(DEMO)
    out = edit(adata)
    adata = out if isinstance(out, ad.AnnData) else adata
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
    assert result["n_cells"] == 500 and result["n_genes"] == 30
    assert result["n_genes_input"] == 50
    assert result["reference"] == "synthetic demo, no genome build"
    assert result["terminal_states"] == TERMINAL
    assert result["top_genes"] == TOP_GENES
    medians = list(result["median_velocity_pseudotime"].values())
    assert list(result["median_velocity_pseudotime"]) == ["s1", "s2", "s3", "s4"]
    assert medians == sorted(medians)
    assert result["seed"] == 0
    assert len(result["chat_summary_lines"]) == 4
    assert result["chat_summary_lines"][0].startswith("Context: synthetic demo")
    assert result["context"].startswith("synthetic demo")
    assert result["status_counts"] == {"tested": 30, "not_tested": 20}
    assert result["gene_id_col"] == "gene_ids"
    genes = pd.read_csv(demo_out / "tables" / "velocity_genes.tsv", sep="\t")
    assert list(genes.columns[:2]) == ["gene_id", "gene_symbol"]
    assert len(genes) == 50 and genes["gene_symbol"].is_monotonic_increasing
    assert genes.loc[genes["gene_symbol"] == "G01", "gene_id"].item() == (
        "ENSG90000000001.1"
    )
    untested = genes[genes["status"] != "tested"]
    assert untested["velocity_r2"].isna().all()
    table = pd.read_csv(demo_out / "tables" / "velocity_pseudotime.tsv", sep="\t")
    truth = ad.read_h5ad(DEMO).obs["true_t"]
    rho = stats.spearmanr(
        table.set_index("cell").loc[truth.index, "velocity_pseudotime"], truth
    ).statistic
    assert round(rho, 3) == RHO_TRUE_TIME
    fates = pd.read_csv(demo_out / "tables" / "fate_probabilities.tsv", sep="\t")
    assert list(fates.columns) == ["cell", "stage", *TERMINAL]
    assert fates["cell"].is_monotonic_increasing


def test_output_contract(demo_out):
    paths = contract_paths()
    assert len(paths) == 22
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
    reread = ad.read_h5ad(demo_out / "data" / "velocity.h5ad")
    assert {"velocity", "Ms", "Mu"} <= set(reread.layers)


def test_commands_sh_reproduces_run(demo_out, tmp_path):
    """The resolved command in commands.sh, with a new outdir, gives the same tables."""
    argv = shlex.split(
        (demo_out / "reproducibility" / "commands.sh").read_text().splitlines()[-1]
    )
    argv[argv.index("--outdir") + 1] = str(tmp_path / "again")
    completed = subprocess.run(argv, capture_output=True, text=True, check=False)
    assert completed.returncode == 0, completed.stderr
    for table in sorted((demo_out / "tables").iterdir()):
        assert (tmp_path / "again" / "tables" / table.name).read_bytes() == (
            table.read_bytes()
        ), table.name


def test_auto_terminal_states(tmp_path):
    completed = run_with(tmp_path, terminal_states="auto", n_states="3")
    assert completed.returncode == 0, completed.stderr
    result = json.loads((tmp_path / "out" / "result.json").read_text())
    assert "s4" in result["terminal_states"]


def test_no_terminal_state(demo_out, monkeypatch):
    """predict_terminal_states raising (no stable macrostate) gives the fixed error."""
    sys.path.insert(0, str(SCRIPT.parent))
    import cellrank as cr
    import velocity
    from skillkit import SkillError

    def no_stable(self, **kwargs):
        raise ValueError("No macrostate has a stability >= 0.96")

    monkeypatch.setattr(cr.estimators.GPCCA, "predict_terminal_states", no_stable)
    adata = ad.read_h5ad(demo_out / "data" / "velocity.h5ad")
    p = {"n_states": 2, "groupby": "stage"}
    with pytest.raises(SkillError, match="no macrostate of 2 has stability >= 0.96"):
        velocity.fates(adata, None, p)


def test_dynamical_smoke(tmp_path):
    """Dynamical mode is slow; a 150-cell subset only checks that it runs."""
    path = broken_demo(tmp_path, lambda a: a[::3][:150].copy())
    completed = run_with(tmp_path, adata=path, mode="dynamical")
    assert completed.returncode == 0, completed.stderr
    genes = pd.read_csv(tmp_path / "out" / "tables" / "velocity_genes.tsv", sep="\t")
    assert "fit_likelihood" in genes.columns


@pytest.mark.parametrize("arg", ["--seed", "--reference", "--context", "--gene-id-col"])
def test_missing_argument_fails_without_demo(tmp_path, arg):
    assert is_error(
        run_with(tmp_path, drop=[arg]), f"missing required arguments: {arg}"
    )


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--seed", "1", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --seed")


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"mode": "steady"}, "--mode must be stochastic or dynamical, got steady"),
        ({"min_shared_counts": "-1"}, "--min-shared-counts must be an integer >= 0"),
        ({"n_top_genes": "1"}, "--n-top-genes must be an integer >= 2, got 1"),
        ({"n_pcs": "1"}, "--n-pcs must be an integer >= 2, got 1"),
        ({"n_neighbors": "1"}, "--n-neighbors must be an integer >= 2, got 1"),
        ({"seed": "-1"}, "--seed must be an integer >= 0, got -1"),
        (
            {"terminal_states": " , "},
            "--terminal-states must be a comma list of groups or auto, got ' , '",
        ),
        (
            {"terminal_states": "auto", "n_states": "1"},
            "--n-states must be an integer >= 2, got 1",
        ),
        ({"terminal_states": "auto"}, "--terminal-states auto needs --n-states"),
        ({"n_states": "3"}, "--n-states is used only with --terminal-states auto"),
        ({"adata": "/no/such.h5ad"}, "input file not found: /no/such.h5ad"),
        ({"groupby": "leiden"}, "obs column not found: leiden"),
        ({"gene_id_col": "ensembl"}, "var column not found: ensembl"),
        (
            {"terminal_states": "s4,s9"},
            "terminal states not in obs column 'stage': s9",
        ),
        (
            {"n_top_genes": "60"},
            "--n-top-genes 60 exceeds the 50 genes left after --min-shared-counts 10",
        ),
        ({"n_pcs": "30"}, "--n-pcs 30 must be below min(cells, genes) = 30"),
        (
            {"terminal_states": "auto", "n_states": "12"},
            "GPCCA could not form 12 macrostates",
        ),
    ],
)
def test_input_errors(tmp_path, changes, message):
    assert is_error(run_with(tmp_path, **changes), message)


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


@pytest.mark.parametrize("layer", ["spliced", "unspliced"])
def test_missing_layer(tmp_path, layer):
    path = broken_demo(tmp_path, lambda a: a.layers.pop(layer))
    assert is_error(run_with(tmp_path, adata=path), f"layer not found: {layer}")


def test_no_umap(tmp_path):
    path = broken_demo(tmp_path, lambda a: a.obsm.pop("X_umap"))
    assert is_error(run_with(tmp_path, adata=path), "obsm['X_umap'] not found")


def test_duplicate_gene_names(tmp_path):
    def dupe(a):
        a.var_names = ["G01", *a.var_names[1:-1], "G01"]

    path = broken_demo(tmp_path, dupe)
    assert is_error(run_with(tmp_path, adata=path), "duplicate gene names: G01")


def test_duplicate_cell_ids(tmp_path):
    def dupe(a):
        a.obs_names = [a.obs_names[1], *a.obs_names[1:]]

    path = broken_demo(tmp_path, dupe)
    assert is_error(run_with(tmp_path, adata=path), "duplicate cell ids: cell0001")


def test_one_group(tmp_path):
    path = broken_demo(tmp_path, lambda a: a.obs.__setitem__("one", "all"))
    completed = run_with(tmp_path, adata=path, groupby="one", terminal_states="all")
    assert is_error(completed, "obs column 'one' has 1 group(s); fates need at least 2")


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()


def test_duplicate_gene_id_fails(tmp_path):
    """One Ensembl id on two symbols is an ambiguous map: stop, never guess."""

    def dupe(a):
        a.var.loc["G02", "gene_ids"] = a.var.loc["G01", "gene_ids"]

    path = broken_demo(tmp_path, dupe)
    assert is_error(
        run_with(tmp_path, adata=path), "duplicate gene ids: ENSG90000000001.1"
    )


def test_gene_id_col_none_writes_na(tmp_path):
    completed = run_with(tmp_path, gene_id_col="none")
    assert completed.returncode == 0, completed.stderr
    result = json.loads((tmp_path / "out" / "result.json").read_text())
    assert result["gene_id_col"] == "none"
    genes = pd.read_csv(tmp_path / "out" / "tables" / "velocity_genes.tsv", sep="\t")
    assert genes["gene_id"].isna().all()
