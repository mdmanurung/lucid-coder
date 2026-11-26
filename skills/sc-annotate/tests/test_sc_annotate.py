"""CLI tests for scripts/sc_annotate.py: demo goldens, output contract, fail strings.

Run from skills/sc-annotate with pytest;
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
SCRIPT = SKILL_DIR / "scripts" / "sc_annotate.py"
EXAMPLES = SKILL_DIR / "examples"
QUERY = EXAMPLES / "demo_query.h5ad"
REFERENCE = EXAMPLES / "demo_reference.h5ad"
MARKERS = EXAMPLES / "demo_markers.tsv"
TEMPLATES = SKILL_DIR.parents[1] / "templates"

# Golden demo values (celltypist 1.7.1, scanpy 1.12.4, scikit-learn 1.7.2):
# cells of 400 whose label matches obs["true_label"].
N_CORRECT = {"celltypist": 388, "ingest": 379, "model_no_voting": 385}
LABEL_COUNTS = {"B cell": 102, "Monocyte": 98, "NK cell": 102, "T cell": 98}

ARGS = {
    "--query": str(QUERY),
    "--markers": str(MARKERS),
    "--over-clustering": "cluster",
    "--engine": "celltypist",
    "--seed": "0",
    "--reference": "synthetic demo genes, no genome build",
    "--context": "synthetic demo, no real tissue",
    "--reference-adata": str(REFERENCE),
    "--reference-label": "cell_type",
    "--majority-voting": "yes",
}
INGEST = {"engine": "ingest", "majority_voting": None, "n_pcs": "10"}
INGEST["n_neighbors"] = "15"


def run(*argv):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *argv],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )


def run_with(tmp_path, drop=(), **changes):
    """ARGS with `changes` ({"seed": "1"}; None drops a flag) and `drop` applied."""
    args = {**ARGS, **{f"--{k.replace('_', '-')}": v for k, v in changes.items()}}
    argv = [
        x for k, v in args.items() if k not in drop and v is not None for x in (k, v)
    ]
    return run(*argv, "--outdir", str(tmp_path / "out"))


def is_error(completed, message):
    return (
        completed.returncode == 2
        and f"ERROR: {message}" in completed.stderr
        and "Traceback" not in completed.stderr
    )


def n_correct(outdir):
    truth = ad.read_h5ad(QUERY).obs["true_label"].astype(str)
    labels = pd.read_csv(outdir / "tables" / "labels.tsv", sep="\t", index_col="cell")
    return int((labels["label"].loc[truth.index] == truth).sum())


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


@pytest.fixture(scope="module")
def demo_out(tmp_path_factory):
    out = tmp_path_factory.mktemp("demo") / "out"
    completed = run("--demo", "--outdir", str(out))
    assert completed.returncode == 0, completed.stderr
    return out


def test_demo_golden_values(demo_out):
    result = json.loads((demo_out / "result.json").read_text())
    assert result["n_cells"] == 400 and result["n_labels"] == 4
    assert result["n_shared_genes"] == result["n_matched"] == 158
    assert result["n_in_resource"] == 158  # trained model features: shared genes
    assert result["reference"] == "synthetic demo genes, no genome build"
    assert result["context"].startswith("synthetic demo")
    assert result["chat_summary_lines"][0].startswith("Context: ")
    assert result["label_counts"] == LABEL_COUNTS
    assert result["seed"] == 0 and result["engine"] == "celltypist"
    assert 2 <= len(result["chat_summary_lines"]) <= 3
    assert n_correct(demo_out) == N_CORRECT["celltypist"]
    labels = pd.read_csv(demo_out / "tables" / "labels.tsv", sep="\t")
    assert list(labels.columns) == ["cell", "label", "confidence"]
    assert labels["cell"].is_monotonic_increasing
    assert labels["confidence"].between(0, 1).all()
    by_cluster = pd.read_csv(demo_out / "tables" / "label_by_cluster.tsv", sep="\t")
    assert by_cluster["n_cells"].sum() == 400
    assert by_cluster.groupby("cluster")["frac_of_cluster"].sum().round(9).eq(1).all()
    annotated = ad.read_h5ad(demo_out / "data" / "annotated.h5ad")
    assert (annotated.X.data == ad.read_h5ad(QUERY).X.data).all()  # raw counts kept
    assert "transferred_label" in annotated.obs


def test_ingest_demo(tmp_path):
    completed = run("--demo", "--engine", "ingest", "--outdir", str(tmp_path / "out"))
    assert completed.returncode == 0, completed.stderr
    assert n_correct(tmp_path / "out") == N_CORRECT["ingest"]


def test_trained_model_reused_without_voting(tmp_path, demo_out):
    model = demo_out / "data" / "celltypist_model.pkl"
    completed = run_with(
        tmp_path, reference_adata=None, reference_label=None, model=str(model),
        majority_voting="no",
    )  # fmt: skip
    assert completed.returncode == 0, completed.stderr
    assert n_correct(tmp_path / "out") == N_CORRECT["model_no_voting"]


def test_output_contract(demo_out):
    paths = contract_paths()
    assert len(paths) == 21
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


@pytest.mark.parametrize("flag", ["--seed", "--context"])
def test_missing_argument_fails_without_demo(tmp_path, flag):
    completed = run_with(tmp_path, drop=[flag])
    assert is_error(completed, f"missing required arguments: {flag}")


def test_missing_reference_fails(tmp_path):
    completed = run_with(tmp_path, drop=["--reference"])
    assert is_error(completed, "missing required arguments: --reference")


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--seed", "1", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --seed")


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


def write_h5ad(tmp_path, source, edit):
    adata = ad.read_h5ad(source)
    edit(adata)
    path = tmp_path / f"edited_{source.name}"
    adata.write_h5ad(path)
    return str(path)


def write_tsv(tmp_path, rows, columns=("label", "gene")):
    path = tmp_path / "markers.tsv"
    pd.DataFrame(rows, columns=list(columns)).to_csv(path, sep="\t", index=False)
    return str(path)


def scale_x(adata):
    adata.X = adata.X * 0.5


def drop_umap(adata):
    del adata.obsm["X_umap"]


def one_label(adata):
    adata.obs["cell_type"] = "T cell"


def duplicate_gene(adata):
    adata.var_names = [adata.var_names[0], *adata.var_names[1:-1], adata.var_names[0]]


def duplicate_cell(adata):
    adata.obs_names = [adata.obs_names[0], *adata.obs_names[1:-1], adata.obs_names[0]]


def rename_genes(adata):
    adata.var_names = [f"X{i}" for i in range(adata.n_vars)]


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"majority_voting": None}, "--engine celltypist needs: --majority-voting"),
        ({"n_pcs": "10"}, "--engine celltypist does not take: --n-pcs"),
        (
            {"model": str(MARKERS)},
            "--engine celltypist needs --model or --reference-adata, not both",
        ),
        ({"seed": "-1"}, "--seed must be an integer >= 0, got -1"),
        ({**INGEST, "n_neighbors": "1"}, "--n-neighbors must be an integer >= 2"),
        ({"query": "/no/such.h5ad"}, "input file not found: /no/such.h5ad"),
        ({"over_clustering": "batch"}, "obs column not found in the query: batch"),
        (
            {"reference_label": "label"},
            "obs column not found in the reference: label",
        ),
        (
            {"reference_adata": None, "reference_label": None, "model": str(MARKERS)},
            f"cannot load the CellTypist model: {MARKERS}",
        ),
        (
            {**INGEST, "n_pcs": "200"},
            "--n-pcs 200 needs more than 200 shared genes, got 158",
        ),
    ],
)
def test_input_errors(tmp_path, changes, message):
    assert is_error(run_with(tmp_path, **changes), message)


@pytest.mark.parametrize(
    ("which", "edit", "message"),
    [
        ("query", scale_x, "X of the query must hold raw counts"),
        ("query", drop_umap, "the query has no obsm['X_umap']"),
        ("reference", scale_x, "X of the reference must hold raw counts"),
        ("reference", one_label, "the reference has fewer than 2 labels"),
        ("reference", rename_genes, "the query shares no genes with the reference"),
        ("query", duplicate_gene, "duplicate gene names in the query: "),
        ("query", duplicate_cell, "duplicate cell ids in the query: "),
        ("reference", duplicate_gene, "duplicate gene names in the reference: "),
        ("reference", duplicate_cell, "duplicate cell ids in the reference: "),
    ],
)
def test_bad_h5ad(tmp_path, which, edit, message):
    source = QUERY if which == "query" else REFERENCE
    path = write_h5ad(tmp_path, source, edit)
    flag = "query" if which == "query" else "reference_adata"
    assert is_error(run_with(tmp_path, **{flag: path}), message)


@pytest.mark.parametrize(
    ("rows", "columns", "message"),
    [
        ([("T cell", "CD3E")], ("cell_type", "gene"), "--markers must be a TSV"),
        ([("Platelet", "CD3E")], ("label", "gene"), "marker labels not among"),
        (
            [("T cell", "PPBP")],
            ("label", "gene"),
            "marker genes not in the query: PPBP",
        ),
    ],
)
def test_bad_markers(tmp_path, rows, columns, message):
    path = write_tsv(tmp_path, rows, columns)
    assert is_error(run_with(tmp_path, markers=path), message)


def test_counts_layer_errors(tmp_path):
    completed = run_with(tmp_path, counts_layer="counts")
    assert is_error(completed, "layer 'counts' not found in the query")

    def layer_not_counts(adata):
        adata.layers["counts"] = adata.X * 0.5

    path = write_h5ad(tmp_path, QUERY, layer_not_counts)
    completed = run_with(tmp_path, query=path, counts_layer="counts")
    assert is_error(completed, "layer 'counts' of the query must hold raw counts")


def test_log_normalised_x_needs_counts_layer(tmp_path):
    def embed_like(adata):
        adata.layers["counts"] = adata.X.copy()
        adata.X = adata.X.astype("float32") * 0.5

    path = write_h5ad(tmp_path, QUERY, embed_like)
    completed = run_with(tmp_path, query=path)
    assert is_error(completed, "X of the query must hold raw counts")
    assert "pass --counts-layer <layer>" in completed.stderr


def test_embed_output_to_annotate_to_ccc(tmp_path):
    """sc-embed's embedded.h5ad (log X, raw counts layer) goes through annotate
    with --counts-layer, keeps its X, and is then accepted by ccc."""
    embed = SKILL_DIR.parent / "sc-embed" / "scripts" / "sc_embed.py"
    ccc = SKILL_DIR.parent / "ccc" / "scripts" / "ccc.py"
    if not (embed.is_file() and ccc.is_file()):
        pytest.skip("sc-embed and ccc are not beside sc-annotate")

    def go(script, *argv):
        return subprocess.run(
            [sys.executable, str(script), *argv],
            capture_output=True,
            text=True,
            check=False,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        )

    done = go(embed, "--demo", "--outdir", str(tmp_path / "embed"))
    assert done.returncode == 0, done.stderr
    embedded = tmp_path / "embed" / "data" / "embedded.h5ad"
    before = ad.read_h5ad(embedded)
    labels = sorted(before.obs["true_group"].astype(str).unique())
    genes = list(before.var_names[:6])
    markers = write_tsv(
        tmp_path, [(labels[0], genes[0]), (labels[1], genes[1])], ("label", "gene")
    )
    done = run(
        "--query", str(embedded), "--counts-layer", "counts",
        "--markers", markers, "--over-clustering", "leiden",
        "--engine", "celltypist", "--reference-adata", str(embedded),
        "--reference-counts-layer", "counts", "--reference-label", "true_group",
        "--majority-voting", "yes", "--seed", "0",
        "--reference", "synthetic demo genes, no genome build",
        "--context", "embed to annotate to ccc",
        "--outdir", str(tmp_path / "annotate"),
    )  # fmt: skip
    assert done.returncode == 0, done.stderr
    annotated = ad.read_h5ad(tmp_path / "annotate" / "data" / "annotated.h5ad")
    # X is exactly what sc-embed wrote: log-normalised, not counts; layer kept.
    assert (annotated.X != before.X).nnz == 0
    assert (annotated.layers["counts"] != before.layers["counts"]).nnz == 0
    assert "transferred_label" in annotated.obs

    resource = tmp_path / "resource.tsv"
    resource.write_text(f"ligand\treceptor\n{genes[2]}\t{genes[3]}\n")
    done = go(
        ccc,
        "--adata", str(tmp_path / "annotate" / "data" / "annotated.h5ad"),
        "--context", "embed to annotate to ccc", "--gene-id-col", "none",
        "--resource", str(resource), "--groupby", "transferred_label",
        "--methods", "rank_aggregate", "--expr-prop", "0.1", "--min-cells", "5",
        "--n-perms", "10", "--magnitude-max", "1", "--specificity-max", "1",
        "--sources", "all", "--targets", "all", "--n-top", "3",
        "--reference", "synthetic demo genes, no genome build", "--seed", "0",
        "--outdir", str(tmp_path / "ccc"),
    )  # fmt: skip
    assert "does not look log-normalised" not in done.stderr
    assert done.returncode == 0, done.stderr


def test_every_error_has_a_test():
    sys.path.insert(0, str(SCRIPT.parent))
    from sc_annotate import ERRORS

    tested = {
        "engine_missing", "engine_extra", "model_xor_reference", "bad_value",
        "missing_file", "missing_column", "no_umap", "bad_counts", "few_labels",
        "bad_model", "no_shared_genes", "too_few_genes", "bad_markers",
        "marker_label_absent", "marker_gene_absent", "bad_layer_counts", "no_layer",
    }  # fmt: skip
    assert set(ERRORS) == tested


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
