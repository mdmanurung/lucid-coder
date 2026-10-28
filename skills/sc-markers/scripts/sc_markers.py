"""Purpose: marker genes per cluster or label from a log-normalised h5ad.
Inputs:  --adata <file.h5ad> with log-normalised values in X or --layer;
         obs column --groupby (e.g. sc-embed's data/embedded.h5ad)
Outputs: <outdir>/report.md, result.json, tables/, figures/, reproducibility/
Run:     python scripts/sc_markers.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""

import argparse
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from plots import plot_dotplot, plot_heatmap
from scipy import sparse
from skillkit import (
    check_gene_ids,
    check_unique,
    fail,
    make_outdir,
    parse_required,
    run,
    seed_everything,
    with_status,
    write_report,
    write_reproducibility,
    write_result,
)
from theme import save_panel

SKILL_DIR = Path(__file__).resolve().parents[1]

REQUIRED = [
    "adata",
    "groupby",
    "gene_id_col",
    "method",
    "padj_max",
    "lfc_min",
    "min_pct",
    "n_top",
    "limit",
    "reference",
    "context",
    "seed",
]
KEYS = [*REQUIRED[:3], "layer", *REQUIRED[3:]]
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "adata": str(SKILL_DIR / "examples" / "demo.h5ad"),
    "groupby": "cluster",
    "gene_id_col": "gene_ids",
    "layer": None,
    "method": "wilcoxon",
    "padj_max": 0.05,
    "lfc_min": 1.0,
    "min_pct": 0.25,
    "n_top": 5,
    "limit": 2.0,
    "reference": "synthetic demo genes, no genome build",
    "context": "synthetic demo: 3 clusters x 50 cells, no real tissue",
    "seed": 0,
    "outdir": "sc_markers_demo",
}

ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "missing_file": "input file not found: {path}",
    "missing_layer": "layer not found in the h5ad: {layer}",
    "missing_column": "obs column not found: {col}",
    "missing_gene_id_col": "gene id column not found in var: {col}",
    "not_lognorm": "{source} does not look log-normalised: {why}",
    "few_groups": "--groupby {col} needs at least 2 groups, found {n}",
    "small_group": "groups with fewer than 2 cells in {col}: {groups}",
    "no_markers": "zero genes pass {rule}; nothing to plot",
}

MARKER_COLUMNS = {
    "group": "group",
    "names": "gene_symbol",
    "scores": "score",
    "logfoldchanges": "logfc",
    "pvals": "pval",
    "pvals_adj": "padj",
    "pct_nz_group": "pct_in",
    "pct_nz_reference": "pct_out",
}
# Fixed plot spec, not a grill question: the dot plot widens 4 mm per gene, so a
# 4 mm dot at 100% expressing fills its column without touching the next one.
DOT_MAX_SIZE = 4
MARKER_COLUMNS_STATS = list(MARKER_COLUMNS.values())[2:]
MEAN_LABEL = "mean log-normalised expression"
FRAC_LABEL = "fraction expressing"
Z_LABEL = "z-score of group mean"


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--adata")
    parser.add_argument("--groupby")
    parser.add_argument("--gene-id-col")
    parser.add_argument("--layer")
    parser.add_argument("--method", choices=["wilcoxon", "t-test"])
    parser.add_argument("--padj-max", type=float)
    parser.add_argument("--lfc-min", type=float)
    parser.add_argument("--min-pct", type=float)
    parser.add_argument("--n-top", type=int)
    parser.add_argument("--limit", type=float)
    parser.add_argument("--reference")
    parser.add_argument("--context")
    parser.add_argument("--seed", type=int)
    p = parse_required(parser, REQUIRED, DEMO, argv, optional=["layer"])
    rules = [
        ("--padj-max", p["padj_max"], 0 < p["padj_max"] <= 1, "in (0, 1]"),
        ("--lfc-min", p["lfc_min"], p["lfc_min"] >= 0, ">= 0"),
        ("--min-pct", p["min_pct"], 0 <= p["min_pct"] <= 1, "in [0, 1]"),
        ("--n-top", p["n_top"], p["n_top"] >= 1, "an integer >= 1"),
        ("--limit", p["limit"], p["limit"] > 0, "> 0"),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ]
    for arg, value, ok, rule in rules:
        if not ok:
            fail(ERRORS, "bad_value", arg=arg, rule=rule, value=value)
    return p


def load(p):
    """Step 1: read the h5ad, check ids, the log-normalised matrix and the groups."""
    path = Path(p["adata"])
    if not path.is_file():
        fail(ERRORS, "missing_file", path=path)
    adata = ad.read_h5ad(path)
    # var_names are the symbols; --gene-id-col none: the source has no gene ids.
    ids = None
    if p["gene_id_col"] != "none":
        if p["gene_id_col"] not in adata.var.columns:
            fail(ERRORS, "missing_gene_id_col", col=p["gene_id_col"])
        ids = adata.var[p["gene_id_col"]].astype(str).tolist()
    gene_ids = check_gene_ids(ids, adata.var_names.tolist())
    check_unique(adata.obs_names, "cell ids")
    if p["layer"] is not None and p["layer"] not in adata.layers:
        fail(ERRORS, "missing_layer", layer=p["layer"])
    col = p["groupby"]
    if col not in adata.obs.columns:
        fail(ERRORS, "missing_column", col=col)

    source = f"layer {p['layer']}" if p["layer"] else "X"
    matrix = adata.layers[p["layer"]] if p["layer"] else adata.X
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if (values < 0).any():
        fail(ERRORS, "not_lognorm", source=source, why="negative values (scaled?)")
    if np.array_equal(values, np.round(values)):
        fail(ERRORS, "not_lognorm", source=source, why="all integers (raw counts?)")

    # Categorical order is the group order of every table and figure.
    groups = adata.obs[col]
    if not isinstance(groups.dtype, pd.CategoricalDtype):
        groups = pd.Categorical(groups.astype(str), categories=sorted(set(groups)))
    adata.obs[col] = pd.Categorical(groups).remove_unused_categories()
    sizes = adata.obs[col].value_counts(sort=False)
    if len(sizes) < 2:
        fail(ERRORS, "few_groups", col=col, n=len(sizes))
    if (sizes < 2).any():
        fail(ERRORS, "small_group", col=col, groups=", ".join(sizes.index[sizes < 2]))
    return adata, gene_ids


def rank_markers(adata, p, gene_ids):
    """Steps 2-3: each group against the rest, then the pass rule."""
    sc.tl.rank_genes_groups(
        adata,
        p["groupby"],
        method=p["method"],
        layer=p["layer"],
        # use_raw=False: an h5ad carrying .raw would otherwise be tested instead.
        use_raw=False,
        reference="rest",
        corr_method="benjamini-hochberg",
        pts=True,
    )
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers = markers.rename(columns=MARKER_COLUMNS)[list(MARKER_COLUMNS.values())]
    markers["group"] = pd.Categorical(
        markers["group"].astype(str), categories=adata.obs[p["groupby"]].cat.categories
    )
    markers.insert(2, "gene_id", markers["gene_symbol"].map(gene_ids))
    markers["passes"] = (
        (markers["padj"] < p["padj_max"])
        & (markers["logfc"] > p["lfc_min"])
        & (markers["pct_in"] >= p["min_pct"])
    )
    # Every gene is tested in every group: no filter, so every status is "tested".
    markers = with_status(markers, adata.var_names, "gene_symbol", {})
    return markers.sort_values(
        ["group", "score", "gene_symbol"], ascending=[True, False, True]
    ).reset_index(drop=True)


def group_profiles(adata, p, genes):
    """Mean and fraction expressing per group x gene, long, ordered for plotting."""
    col, groups = p["groupby"], adata.obs[p["groupby"]].cat.categories
    matrix = adata[:, genes].layers[p["layer"]] if p["layer"] else adata[:, genes].X
    matrix = sparse.csr_matrix(matrix)
    rows = []
    for group in groups:
        in_group = matrix[(adata.obs[col] == group).to_numpy()]
        rows.append(
            pd.DataFrame(
                {
                    "group": group,
                    "gene": genes,
                    MEAN_LABEL: np.asarray(in_group.mean(axis=0)).ravel(),
                    FRAC_LABEL: np.asarray((in_group > 0).mean(axis=0)).ravel(),
                }
            )
        )
    long = pd.concat(rows, ignore_index=True)
    means = long.groupby("gene")[MEAN_LABEL]
    long[Z_LABEL] = (long[MEAN_LABEL] - means.transform("mean")) / means.transform(
        "std"
    )
    long["group"] = pd.Categorical(long["group"], categories=groups)
    long["gene"] = pd.Categorical(long["gene"], categories=genes)
    return long


def draw_figures(figures, profiles, p):
    plot = plot_dotplot(
        profiles,
        "gene",
        "group",
        MEAN_LABEL,
        FRAC_LABEL,
        max_size=DOT_MAX_SIZE,
        mean_label="Mean expression (log-normalised)",
    )
    # mm: 89 single column up to 183 double; rows and columns grow with the data.
    width = min(183, max(89, 60 + 4 * profiles["gene"].nunique()))
    # 85 mm floor: the colourbar and size legend stack on the right.
    height = max(85, 25 + 8 * profiles["group"].nunique())
    save_panel(plot, profiles, figures, "marker_dotplot", width, height)

    plot = plot_heatmap(profiles, "group", "gene", Z_LABEL, limit=p["limit"])
    save_panel(plot, profiles, figures, "marker_heatmap", width, height)


def analyse(argv, full_argv):
    p = parse_args(argv)
    seed_everything(p["seed"])
    # Neither test draws random numbers; one job keeps the run single-threaded.
    sc.settings.n_jobs = 1
    adata, gene_ids = load(p)
    markers = rank_markers(adata, p, gene_ids)
    rule = f"padj < {p['padj_max']}, logfc > {p['lfc_min']}, pct_in >= {p['min_pct']}"
    top = markers[markers["passes"]].groupby("group", observed=False).head(p["n_top"])
    if top.empty:
        fail(ERRORS, "no_markers", rule=rule)
    top = top.assign(rank=top.groupby("group", observed=True).cumcount() + 1)
    top = top[["group", "rank", "gene_symbol", "gene_id", *MARKER_COLUMNS_STATS]]

    outdir = make_outdir(p["outdir"])
    tables, figures = outdir / "tables", outdir / "figures"
    markers.to_csv(tables / "markers.tsv", sep="\t", index=False)
    top.to_csv(tables / "top_markers.tsv", sep="\t", index=False)

    plt.switch_backend("Agg")
    # A gene topping two groups is plotted once, in its first group's block.
    genes = list(dict.fromkeys(top["gene_symbol"]))
    profiles = group_profiles(adata, p, genes)
    draw_figures(figures, profiles, p)

    col, n_groups = p["groupby"], adata.obs[p["groupby"]].nunique()
    n_pass = markers[markers["passes"]].groupby("group", observed=False).size()
    captions = {
        "marker_dotplot": f"Dot plot of the {len(genes)} top markers (up to "
        f"{p['n_top']} per {col} group, {rule}) across {n_groups} groups; colour is "
        "mean log-normalised expression, dot area the fraction of cells expressing.",
        "marker_heatmap": f"Heatmap of group-mean log-normalised expression for the "
        f"{len(genes)} top markers, z-scored per gene across {n_groups} {col} groups; "
        f"not clustered: rows in {col} order, columns in top-marker order; colours "
        f"saturate at ±{p['limit']:g}.",
    }
    empty = n_pass.index[n_pass == 0].astype(str).tolist()
    summary = [
        f"{adata.n_obs} cells in {n_groups} {col} groups; {adata.n_vars} genes "
        f"tested per group ({p['method']}, each group vs rest, BH).",
        f"{int(n_pass.sum())} group-gene pairs pass {rule}; "
        f"groups with none: {', '.join(empty) or 'none'}.",
    ]
    write_report(outdir, f"Marker genes by {col}", p, KEYS, summary, captions)
    write_result(
        outdir,
        {
            "n_cells": adata.n_obs,
            "n_groups": n_groups,
            "n_genes_tested": adata.n_vars,
            "status_counts": {
                k: int(v) for k, v in markers["status"].value_counts().items()
            },
            "gene_id_col": p["gene_id_col"],
            "n_pass": {str(g): int(n) for g, n in n_pass.items()},
            "top_markers": {
                str(g): d["gene_symbol"].tolist()
                for g, d in top.groupby("group", observed=True)
            },
            "method": p["method"],
            "seed": p["seed"],
            "params": p,
        },
        summary,
        ["report.md", "tables/top_markers.tsv", "figures/marker_dotplot.png"],
        context=p["context"],
        reference=p["reference"],
    )
    write_reproducibility(
        outdir,
        __file__,
        p,
        KEYS,
        full_argv,
        ["anndata", "scanpy", "numpy", "pandas", "scipy", "matplotlib", "plotnine"],
    )
    return outdir


if __name__ == "__main__":
    raise SystemExit(run("sc-markers", analyse))
