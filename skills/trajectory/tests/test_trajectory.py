"""CLI tests for scripts/trajectory.py: demo goldens, output contract, fail strings.

Run from skills/trajectory with pytest; see references/maintainer.md for the uv line.
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
SCRIPT = SKILL_DIR / "scripts" / "trajectory.py"
DEMO = SKILL_DIR / "examples" / "demo.h5ad"
TEMPLATES = SKILL_DIR.parents[1] / "templates"

# Golden demo values, computed once with scanpy 1.12.4. The Spearman rho of DPT
# with the simulated time is 0.9791; rounded to 2 places it sits 0.004 from an edge.
RHO_TRUE_TIME = 0.98
ROOT_CELL = "cell0005"
EDGES = [("progenitor", "early"), ("early", "mid"), ("mid", "late")]

ARGS = {
    "--adata": str(DEMO),
    "--groups": "stage",
    "--root-group": "progenitor",
    "--n-dcs": "10",
    "--paga-threshold": "0.1",
    "--genes": "G00,G01,G02",
    "--reference": "synthetic demo, no genome build",
    "--context": "synthetic demo: 600 cells, 4 stages along one simulated path, "
    "no real tissue",
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
    """Demo arguments with `changes` ({"n_dcs": "1"}) and `drop` applied."""
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
    assert result["n_cells"] == 600 and result["n_groups"] == 4
    assert result["root_cell"] == ROOT_CELL and result["n_unreachable"] == 0
    edges = [(e["group_a"], e["group_b"]) for e in result["edges_above_threshold"]]
    assert edges == EDGES
    assert list(result["median_pseudotime"]) == ["progenitor", "early", "mid", "late"]
    assert result["seed"] == 0
    assert result["reference"] == "synthetic demo, no genome build"
    assert result["dropped_per_group"] == dict.fromkeys(
        ["progenitor", "early", "mid", "late"], 0
    )
    assert 3 <= len(result["chat_summary_lines"]) <= 4
    assert result["context"] == ARGS["--context"]
    assert result["chat_summary_lines"][0].startswith("Context: ")
    pseudotime = pd.read_csv(demo_out / "tables" / "pseudotime.tsv", sep="\t")
    truth = ad.read_h5ad(DEMO).obs["true_time"]
    rho = stats.spearmanr(
        pseudotime.set_index("cell").loc[truth.index, "dpt_pseudotime"], truth
    ).statistic
    assert round(rho, 2) == RHO_TRUE_TIME
    assert pseudotime.loc[pseudotime["is_root"], "cell"].tolist() == [ROOT_CELL]


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


def test_commands_sh_reproduces_run(demo_out, tmp_path):
    """The resolved command in commands.sh, with a new outdir, gives the same tables."""
    last = (demo_out / "reproducibility" / "commands.sh").read_text().splitlines()[-1]
    argv = shlex.split(last)
    argv[argv.index("--outdir") + 1] = str(tmp_path / "rerun")
    completed = subprocess.run(argv, capture_output=True, text=True, check=False)
    assert completed.returncode == 0, completed.stderr
    for name in ("pseudotime.tsv", "paga_connectivities.tsv"):
        rerun = (tmp_path / "rerun" / "tables" / name).read_bytes()
        assert rerun == (demo_out / "tables" / name).read_bytes(), name


def test_genes_from_file(tmp_path):
    genes = tmp_path / "genes.txt"
    genes.write_text("G03\nG04\n")
    completed = run_with(tmp_path, genes=str(genes))
    assert completed.returncode == 0, completed.stderr
    result = json.loads((tmp_path / "out" / "result.json").read_text())
    assert result["genes"] == ["G03", "G04"]


def test_missing_argument_fails_without_demo(tmp_path):
    assert is_error(
        run_with(tmp_path, drop=["--seed"]), "missing required arguments: --seed"
    )


def test_missing_context_fails(tmp_path):
    completed = run_with(tmp_path, drop=["--context"])
    assert is_error(completed, "missing required arguments: --context")


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--seed", "1", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --seed")


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"n_dcs": "1"}, "--n-dcs must be an integer >= 2, got 1"),
        ({"paga_threshold": "2"}, "--paga-threshold must be in [0, 1], got 2.0"),
        ({"seed": "-1"}, "--seed must be an integer >= 0, got -1"),
        ({"adata": "/no/such.h5ad"}, "input file not found: /no/such.h5ad"),
        ({"groups": "leiden"}, "obs column not found: leiden"),
        ({"root_group": "stem"}, "root group 'stem' not found in obs column 'stage'"),
        ({"n_dcs": "600"}, "--n-dcs 600 must be below the number of cells (600)"),
        ({"genes": "G00,NOPE"}, "genes not in var_names: NOPE"),
        ({"genes": " , "}, "--genes names no genes: ' , '"),
    ],
)
def test_input_errors(tmp_path, changes, message):
    assert is_error(run_with(tmp_path, **changes), message)


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


def test_one_group(tmp_path):
    path = broken_demo(tmp_path, lambda a: a.obs.__setitem__("one", "all"))
    completed = run_with(tmp_path, adata=path, groups="one", root_group="all")
    assert is_error(completed, "obs column 'one' has 1 group(s); PAGA needs at least 2")


def test_no_neighbors(tmp_path):
    path = broken_demo(tmp_path, lambda a: a.uns.pop("neighbors"))
    assert is_error(run_with(tmp_path, adata=path), "no neighbour graph")


@pytest.mark.parametrize(
    ("axis", "message"),
    [("var", "duplicate gene names: G00"), ("obs", "duplicate cell ids: cell0000")],
)
def test_duplicate_ids(tmp_path, axis, message):
    def edit(a):
        names = getattr(a, f"{axis}_names").tolist()
        names[1] = names[0]
        setattr(a, f"{axis}_names", names)

    assert is_error(run_with(tmp_path, adata=broken_demo(tmp_path, edit)), message)


def test_no_umap(tmp_path):
    path = broken_demo(tmp_path, lambda a: a.obsm.pop("X_umap"))
    assert is_error(run_with(tmp_path, adata=path), "obsm['X_umap'] not found")


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
