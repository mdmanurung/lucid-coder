"""Purpose: transfer cell-type labels from a reference to a query h5ad.
Inputs:  --query <file.h5ad> (raw counts in X, obsm X_umap), --markers <TSV>,
         and a CellTypist --model or a --reference-adata h5ad with --reference-label;
         --reference names the genome build and annotation release
Outputs: <outdir>/report.md, result.json, data/, tables/, figures/, reproducibility/
Run:     python scripts/sc_annotate.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""

import argparse
import logging
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skillkit as kit
from palettes import categorical_palette
from plotnine import scale_fill_gradient
from plots import (
    FEW_POINTS,
    ROTATE_X,
    plot_box,
    plot_dotplot,
    plot_embedding,
    plot_heatmap,
)
from scipy import sparse
from theme import save_panel

SKILL_DIR = Path(__file__).resolve().parents[1]
EXAMPLES = SKILL_DIR / "examples"

REQUIRED = [
    "query",
    "markers",
    "over_clustering",
    "engine",
    "seed",
    "reference",
    "context",
]
ENGINE_ARGS = {
    "celltypist": ["model", "reference_adata", "reference_label", "majority_voting"],
    "ingest": ["reference_adata", "reference_label", "n_pcs", "n_neighbors"],
}
OPTIONAL = ["model", "reference_adata", "reference_label", "majority_voting"]
OPTIONAL += ["n_pcs", "n_neighbors"]
KEYS = REQUIRED + OPTIONAL
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "query": str(EXAMPLES / "demo_query.h5ad"),
    "markers": str(EXAMPLES / "demo_markers.tsv"),
    "over_clustering": "cluster",
    "engine": "celltypist",
    "seed": 0,
    "reference": "synthetic demo genes, no genome build",
    "context": "synthetic demo: 4 PBMC-like cell types, reference and half-depth "
    "query, no real tissue",
    "outdir": "sc_annotate_demo",
}
DEMO_ENGINE = {
    "celltypist": {
        "reference_adata": str(EXAMPLES / "demo_reference.h5ad"),
        "reference_label": "cell_type",
        "majority_voting": "yes",
    },
    "ingest": {
        "reference_adata": str(EXAMPLES / "demo_reference.h5ad"),
        "reference_label": "cell_type",
        "n_pcs": 10,
        "n_neighbors": 15,
    },
}

ERRORS = {
    "engine_missing": "--engine {engine} needs: {args}",
    "engine_extra": "--engine {engine} does not take: {args}",
    "model_xor_reference": "--engine celltypist needs --model or --reference-adata, "
    "not both",
    "bad_value": "{arg} must be {rule}, got {value}",
    "missing_file": "input file not found: {path}",
    "missing_column": "obs column not found in the {which}: {col}",
    "no_umap": "the query has no obsm['X_umap']; compute the embedding first",
    "bad_counts": "X of the {which} must hold raw counts (non-negative integers)",
    "few_labels": "the reference has fewer than 2 labels in {col!r}",
    "bad_model": "cannot load the CellTypist model: {path}",
    "no_shared_genes": "the query shares no genes with the {which}",
    "too_few_genes": "--n-pcs {n_pcs} needs more than {n_pcs} shared genes, got {n}",
    "bad_markers": "--markers must be a TSV with columns label and gene",
    "marker_label_absent": "marker labels not among the reference labels: {labels}",
    "marker_gene_absent": "marker genes not in the query: {genes}",
}

PACKAGES = ["anndata", "numpy", "pandas", "scipy", "matplotlib", "plotnine", "scanpy"]
PACKAGES += ["celltypist", "scikit-learn", "pynndescent"]
# Fixed plot spec, not analysis thresholds.
# Hex strings: plotnine legends reject RGBA tuples.
DOT_MAX_SIZE = 4  # mm, dot diameter at frac = 1
JITTER_WIDTH = 0.3  # box widths each side
HEATMAP_LIMIT = 1  # fractions run 0 to 1


# fontTools prints "GSUB pruned" while subsetting PDF fonts.
logging.getLogger("fontTools").setLevel(logging.WARNING)


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--query")
    parser.add_argument("--markers")
    parser.add_argument("--over-clustering")
    parser.add_argument("--engine", choices=list(ENGINE_ARGS))
    parser.add_argument("--seed", type=int)
    parser.add_argument("--model")
    parser.add_argument("--reference")
    parser.add_argument("--context")
    parser.add_argument("--reference-adata")
    parser.add_argument("--reference-label")
    parser.add_argument("--majority-voting", choices=["yes", "no"])
    parser.add_argument("--n-pcs", type=int)
    parser.add_argument("--n-neighbors", type=int)
    p = kit.parse_required(
        parser, REQUIRED, DEMO, argv, OPTIONAL, demo_extra=("outdir", "engine")
    )
    p = {**dict.fromkeys(OPTIONAL), **p}
    return p | DEMO_ENGINE[p["engine"]] if "--demo" in argv else p


def check_args(p):
    engine = p["engine"]
    allowed = ENGINE_ARGS[engine]
    extra = [k for k in OPTIONAL if k not in allowed and p[k] is not None]
    if extra:
        kit.fail(ERRORS, "engine_extra", engine=engine, args=flags(extra))
    if engine == "celltypist":
        if p["model"] is not None and p["reference_adata"] is not None:
            kit.fail(ERRORS, "model_xor_reference")
        source = ["model"] if p["reference_adata"] is None else ["reference_label"]
        needed = ["majority_voting", *source]
        if p["reference_adata"] is not None:
            needed = ["reference_adata", *needed]
    else:
        needed = allowed
    missing = [k for k in needed if p[k] is None]
    if missing:
        kit.fail(ERRORS, "engine_missing", engine=engine, args=flags(missing))
    rules = [
        ("seed", 0, "an integer >= 0"),
        ("n_pcs", 1, "an integer >= 1"),
        ("n_neighbors", 2, "an integer >= 2"),
    ]
    for key, low, rule in rules:
        if p[key] is not None and p[key] < low:
            kit.fail(ERRORS, "bad_value", arg=kit.flag(key), rule=rule, value=p[key])
    for key in ["query", "markers", "model", "reference_adata"]:
        if p[key] is not None and not Path(p[key]).is_file():
            kit.fail(ERRORS, "missing_file", path=p[key])


def flags(keys):
    return " ".join(kit.flag(k) for k in keys)


def is_raw_counts(matrix):
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    return bool(np.all(values >= 0) and np.all(values == np.round(values)))


def read_counts(path, which, cols):
    adata = ad.read_h5ad(path)
    kit.check_unique(adata.var_names, f"gene names in the {which}")
    kit.check_unique(adata.obs_names, f"cell ids in the {which}")
    for col in cols:
        if col not in adata.obs.columns:
            kit.fail(ERRORS, "missing_column", which=which, col=col)
    if not is_raw_counts(adata.X):
        kit.fail(ERRORS, "bad_counts", which=which)
    return adata


def log_normalise(adata):
    """normalize_total(1e4) + log1p on a copy, as every corpus row does."""
    import scanpy as sc

    out = adata.copy()
    out.X = sparse.csr_matrix(out.X, dtype=np.float32)
    sc.pp.normalize_total(out, target_sum=1e4)
    sc.pp.log1p(out)
    return out


def read_markers(path, labels, genes):
    try:
        markers = pd.read_csv(path, sep="\t", dtype=str)
    except pd.errors.ParserError:
        kit.fail(ERRORS, "bad_markers")
    if list(markers.columns) != ["label", "gene"] or markers.isna().any().any():
        kit.fail(ERRORS, "bad_markers")
    absent = sorted(set(markers["label"]) - set(labels))
    if absent:
        kit.fail(ERRORS, "marker_label_absent", labels=", ".join(absent))
    absent = sorted(set(markers["gene"]) - set(genes))
    if absent:
        kit.fail(ERRORS, "marker_gene_absent", genes=", ".join(absent))
    return markers


def load_model(p):
    """A CellTypist model from --model, or trained on --reference-adata (seeded)."""
    import celltypist
    from celltypist import models

    if p["model"] is not None:
        try:
            return models.Model.load(str(Path(p["model"]).resolve())), "model"
        except Exception:  # noqa: BLE001 - celltypist raises a bare Exception
            kit.fail(ERRORS, "bad_model", path=p["model"])
    reference = read_counts(p["reference_adata"], "reference", [p["reference_label"]])
    if reference.obs[p["reference_label"]].nunique() < 2:
        kit.fail(ERRORS, "few_labels", col=p["reference_label"])
    # lbfgs, CellTypist's default solver, is deterministic; the seed reaches
    # the sklearn classifier for any solver that draws random numbers.
    model = celltypist.train(
        log_normalise(reference),
        labels=p["reference_label"],
        n_jobs=1,
        random_state=p["seed"],
    )
    return model, "reference"


def run_celltypist(query, p, model):
    """Label and confidence (CellTypist probability of the assigned label)."""
    import celltypist

    voting = p["majority_voting"] == "yes"
    result = celltypist.annotate(
        log_normalise(query),
        model=model,
        majority_voting=voting,
        over_clustering=p["over_clustering"] if voting else None,
    )
    pred = result.predicted_labels
    label = pred["majority_voting" if voting else "predicted_labels"].astype(str)
    prob = result.probability_matrix
    confidence = prob.to_numpy()[np.arange(len(prob)), prob.columns.get_indexer(label)]
    return label.to_numpy(), confidence


def run_ingest(query, reference, shared, p):
    """kNN vote of scanpy Ingest; confidence = fraction of the k neighbours agreeing."""
    import scanpy as sc

    col = p["reference_label"]
    ref = log_normalise(reference)[:, shared].copy()
    new = log_normalise(query)[:, shared].copy()
    sc.pp.pca(ref, n_comps=p["n_pcs"], random_state=p["seed"])
    sc.pp.neighbors(
        ref, n_neighbors=p["n_neighbors"], n_pcs=p["n_pcs"], random_state=p["seed"]
    )
    ing = sc.tl.Ingest(ref)
    ing.fit(new)
    ing.neighbors(random_state=p["seed"])
    ing.map_labels(col, "knn")
    label = ing._obs[col].astype(str).to_numpy()  # ties: first label in sorted order
    ref_labels = ref.obs[col].astype(str).to_numpy()
    confidence = (ref_labels[ing._indices] == label[:, None]).mean(axis=1)
    return label, confidence


def dotplot_frame(query, markers, labels):
    """Mean log-normalised expression (scaled 0-1 per gene) and fraction expressing."""
    genes = list(dict.fromkeys(markers["gene"]))
    expr = log_normalise(query)[:, genes].to_df()
    detected = (expr > 0).groupby(labels).mean()
    mean = expr.groupby(labels).mean()
    span = (mean.max() - mean.min()).replace(0, 1)
    scaled = (mean - mean.min()) / span
    frame = pd.DataFrame(
        {
            "gene": np.tile(genes, len(mean)),
            "label": np.repeat(mean.index, len(genes)),
            "scaled mean expression": scaled.to_numpy().ravel(),
            "fraction expressing": detected.to_numpy().ravel(),
        }
    )
    frame["gene"] = pd.Categorical(frame["gene"], categories=genes)
    frame["label"] = pd.Categorical(frame["label"], categories=sorted(mean.index))
    return frame


def draw_figures(figures, cells, dots, by_cluster, pal, seed, value_label):
    plot = plot_embedding(cells, "umap_1", "umap_2", "label", pal, axis_name="UMAP")
    save_panel(plot, cells, figures, "umap_label", 120, 89)

    plot = plot_dotplot(
        dots,
        "gene",
        "label",
        "scaled mean expression",
        "fraction expressing",
        max_size=DOT_MAX_SIZE,
        mean_label="Mean expression (scaled 0-1 per gene)",
    )
    save_panel(plot, dots, figures, "marker_dotplot", 140, 89)

    plot = plot_box(
        cells,
        "label",
        "confidence",
        pal,
        jitter_width=JITTER_WIDTH,
        seed=seed,
        value_label=value_label,
    )
    save_panel(plot + ROTATE_X, cells, figures, "confidence_box", 89, 89)

    plot = plot_heatmap(
        by_cluster, "label", "cluster", "fraction of cluster", limit=HEATMAP_LIMIT
    )
    # Fractions: white to red over 0-1, the upper half of the diverging scale.
    plot += scale_fill_gradient(low="white", high="#B2182B", limits=(0, HEATMAP_LIMIT))
    save_panel(plot, by_cluster, figures, "label_cluster_heatmap", 120, 89)


def analyse(argv, argv_full):
    """Steps 1-7 of SKILL.md: validate, transfer, check markers, write outputs."""
    p = parse_args(argv)
    check_args(p)
    query = read_counts(p["query"], "query", [p["over_clustering"]])
    if "X_umap" not in query.obsm:
        kit.fail(ERRORS, "no_umap")
    kit.seed_everything(p["seed"])

    from threadpoolctl import threadpool_limits

    with threadpool_limits(1):
        if p["engine"] == "celltypist":
            model, which = load_model(p)
            labels, genes = list(model.cell_types), model.features
        else:
            col = p["reference_label"]
            reference = read_counts(p["reference_adata"], "reference", [col])
            if reference.obs[col].nunique() < 2:
                kit.fail(ERRORS, "few_labels", col=col)
            which, labels = "reference", sorted(reference.obs[col].astype(str).unique())
            genes = reference.var_names
        shared = query.var_names.intersection(pd.Index(genes))
        if len(shared) == 0:
            kit.fail(ERRORS, "no_shared_genes", which=which)
        if p["engine"] == "ingest" and len(shared) <= p["n_pcs"]:
            kit.fail(ERRORS, "too_few_genes", n_pcs=p["n_pcs"], n=len(shared))
        markers = read_markers(p["markers"], labels, query.var_names)
        outdir = kit.make_outdir(p["outdir"])
        (outdir / "data").mkdir()
        if p["engine"] == "celltypist":
            label, confidence = run_celltypist(query, p, model)
            if which == "reference":
                model.write(str(outdir / "data" / "celltypist_model.pkl"))
        else:
            label, confidence = run_ingest(query, reference, shared, p)

    query.obs["transferred_label"] = pd.Categorical(label)
    query.obs["transferred_confidence"] = confidence
    query.write_h5ad(outdir / "data" / "annotated.h5ad")

    tables, figures = outdir / "tables", outdir / "figures"
    found = sorted(set(label))
    clusters = query.obs[p["over_clustering"]]
    order = (
        [str(c) for c in clusters.cat.categories]
        if isinstance(clusters.dtype, pd.CategoricalDtype)
        else sorted(clusters.astype(str).unique())
    )
    cells = pd.DataFrame(
        {
            "cell": query.obs_names,
            "label": pd.Categorical(label, categories=found),
            "confidence": confidence,
            "cluster": pd.Categorical(clusters.astype(str), categories=order),
            "umap_1": query.obsm["X_umap"][:, 0],
            "umap_2": query.obsm["X_umap"][:, 1],
        }
    )
    dots = dotplot_frame(query, markers, label)
    cells = cells.sort_values("cell", ignore_index=True)
    cells[["cell", "label", "confidence"]].to_csv(
        tables / "labels.tsv", sep="\t", index=False
    )
    by_cluster = (
        pd.crosstab(cells["cluster"], cells["label"], dropna=False)
        .stack()
        .rename("n_cells")
        .reset_index()
    )
    by_cluster["frac_of_cluster"] = by_cluster["n_cells"] / by_cluster.groupby(
        "cluster", observed=True
    )["n_cells"].transform("sum")
    by_cluster = by_cluster.sort_values(["cluster", "label"], ignore_index=True)
    by_cluster.to_csv(tables / "label_by_cluster.tsv", sep="\t", index=False)

    pal = categorical_palette(found)
    plt.switch_backend("Agg")
    heat = by_cluster.rename(columns={"frac_of_cluster": "fraction of cluster"})
    rule = (
        "CellTypist probability of the assigned label"
        if p["engine"] == "celltypist"
        else f"fraction of the {p['n_neighbors']} reference neighbours with that label"
    )
    value_label = (
        "CellTypist probability (0-1)"
        if p["engine"] == "celltypist"
        else "neighbour agreement (fraction)"
    )
    draw_figures(figures, cells, dots, heat, pal, p["seed"], value_label)

    counts = cells["label"].value_counts().reindex(found)
    method = (
        f"majority voting {p['majority_voting']}"
        if p["engine"] == "celltypist"
        else f"k = {p['n_neighbors']}, {p['n_pcs']} PCs"
    )
    summary = [
        f"{len(cells)} cells annotated with {len(found)} labels by {p['engine']} "
        f"from the {which} ({len(shared)} shared genes, {method}).",
        f"Median confidence {np.median(confidence):.3f}; largest label "
        f"{counts.idxmax()} ({counts.max()} cells).",
    ]
    captions = {
        "umap_label": f"UMAP of {len(cells)} query cells coloured by the transferred "
        "label; coordinates are the query's obsm['X_umap'], made upstream, so its "
        "PCs, neighbours, min_dist and seed are in that run's report.",
        "marker_dotplot": f"Dot plot of {dots['gene'].nunique()} marker genes from "
        "--markers by transferred label; colour is mean log-normalised expression "
        "scaled 0-1 per gene, dot size the fraction of cells expressing.",
        "confidence_box": f"Per-cell confidence ({rule}) by transferred label, one "
        f"point per cell, over a box where a label has {FEW_POINTS} or more cells "
        "and a median line otherwise.",
        "label_cluster_heatmap": "Heatmap of the fraction of each "
        f"{p['over_clustering']} cluster assigned to each label; columns sum to 1. "
        "Order fixed by label: rows sorted by label name, columns in cluster "
        "order; not clustered.",
    }
    kit.write_report(
        outdir, f"Cell-type annotation: {p['engine']}", p, KEYS, summary, captions
    )
    kit.write_result(
        outdir,
        {
            "n_cells": len(cells),
            "n_labels": len(found),
            "n_shared_genes": len(shared),
            "n_matched": len(shared),
            "n_in_resource": len(genes),
            "label_counts": {k: int(v) for k, v in counts.items()},
            "median_confidence": float(np.median(confidence)),
            "params": p,
            "engine": p["engine"],
            "seed": p["seed"],
        },
        summary,
        ["report.md", "tables/labels.tsv", "figures/marker_dotplot.png"],
        context=p["context"],
        reference=p["reference"],
    )
    kit.write_reproducibility(outdir, __file__, p, KEYS, argv_full, PACKAGES)
    return outdir


if __name__ == "__main__":
    raise SystemExit(kit.run("sc-annotate", analyse))
