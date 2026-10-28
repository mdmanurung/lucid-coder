"""Purpose: PAGA graph and diffusion pseudotime (DPT) from an embedded h5ad.
Inputs:  --adata <file.h5ad> with neighbours, X_umap and a group column (sc-embed)
Outputs: <outdir>/report.md, result.json, data/, tables/, figures/, reproducibility/
Run:     python scripts/trajectory.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""

import os

# One thread for BLAS and numba before numpy loads: same numbers on any host.
for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(var, "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import argparse  # noqa: E402
from pathlib import Path  # noqa: E402

import anndata as ad  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
from matplotlib.colors import to_hex  # noqa: E402
from plotnine import element_text, theme  # noqa: E402
from plots import FEW_POINTS, plot_box  # noqa: E402
from scipy import sparse, stats  # noqa: E402
from skillkit import (  # noqa: E402
    check_unique,
    fail,
    make_outdir,
    parse_required,
    run,
    seed_everything,
    write_report,
    write_reproducibility,
    write_result,
)
from theme import save_panel  # noqa: E402

SKILL_DIR = Path(__file__).resolve().parents[1]
REQUIRED = [
    "adata",
    "groups",
    "root_group",
    "n_dcs",
    "paga_threshold",
    "genes",
    "reference",
    "context",
    "seed",
]
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "adata": str(SKILL_DIR / "examples" / "demo.h5ad"),
    "groups": "stage",
    "root_group": "progenitor",
    "n_dcs": 10,
    "paga_threshold": 0.1,
    "genes": "G00,G01,G02",
    "reference": "synthetic demo, no genome build",
    "context": "synthetic demo: 600 cells, 4 stages along one simulated path, "
    "no real tissue",
    "seed": 0,
    "outdir": "trajectory_demo",
}
ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "missing_file": "input file not found: {path}",
    "missing_column": "obs column not found: {col}",
    "one_group": "obs column {col!r} has {n} group(s); PAGA needs at least 2",
    "root_absent": "root group {root!r} not found in obs column {col!r}",
    "no_neighbors": "no neighbour graph (uns['neighbors'], obsp['connectivities']); "
    "run sc-embed or sc.pp.neighbors first",
    "no_umap": "obsm['X_umap'] not found; run sc-embed or sc.tl.umap first",
    "too_many_dcs": "--n-dcs {n_dcs} must be below the number of cells ({n_cells})",
    "missing_genes": "genes not in var_names: {genes}",
    "empty_genes": "--genes names no genes: {value!r}",
}
PACKAGES = ["anndata", "scanpy", "igraph", "numpy", "pandas", "scipy", "matplotlib"]
# Fixed plot spec, not analysis thresholds. save_panel saves without a tight
# bbox, so every figure sets a layout engine to keep its labels on the page.
N_TREND_BINS = 20
JITTER_WIDTH = 0.3
MAX_EDGE_WIDTH = 2
PSEUDOTIME = "diffusion pseudotime"
# DPT is rescaled to [0, 1] from the root: a unitless ordering, not time.
PSEUDOTIME_LABEL = "diffusion pseudotime (0-1, unitless)"


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--adata")
    parser.add_argument("--groups")
    parser.add_argument("--root-group")
    parser.add_argument("--n-dcs", type=int)
    parser.add_argument("--paga-threshold", type=float)
    parser.add_argument("--genes")
    parser.add_argument("--reference")
    parser.add_argument("--context")
    parser.add_argument("--seed", type=int)
    return parse_required(parser, REQUIRED, DEMO, argv)


def check_values(p):
    rules = [
        ("--n-dcs", p["n_dcs"], p["n_dcs"] >= 2, "an integer >= 2"),
        (
            "--paga-threshold",
            p["paga_threshold"],
            0 <= p["paga_threshold"] <= 1,
            "in [0, 1]",
        ),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ]
    for arg, value, ok, rule in rules:
        if not ok:
            fail(ERRORS, "bad_value", arg=arg, rule=rule, value=value)
    # A path to a file holds one gene per line; otherwise a comma list.
    genes_path = Path(p["genes"])
    raw = (
        genes_path.read_text().splitlines()
        if genes_path.is_file()
        else p["genes"].split(",")
    )
    genes = [g.strip() for g in raw if g.strip()]
    if not genes:
        fail(ERRORS, "empty_genes", value=p["genes"])
    return genes


def load(p, genes):
    """Step 1: read the h5ad and check every field the later steps use."""
    path = Path(p["adata"])
    if not path.is_file():
        fail(ERRORS, "missing_file", path=path)
    adata = ad.read_h5ad(path)
    check_unique(adata.var_names, "gene names")
    check_unique(adata.obs_names, "cell ids")
    col = p["groups"]
    if col not in adata.obs.columns:
        fail(ERRORS, "missing_column", col=col)
    if "neighbors" not in adata.uns or "connectivities" not in adata.obsp:
        fail(ERRORS, "no_neighbors")
    if "X_umap" not in adata.obsm:
        fail(ERRORS, "no_umap")
    if p["n_dcs"] >= adata.n_obs:
        fail(ERRORS, "too_many_dcs", n_dcs=p["n_dcs"], n_cells=adata.n_obs)
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        fail(ERRORS, "missing_genes", genes=", ".join(missing))
    # Step 2: PAGA needs a categorical; an existing category order is kept.
    if not isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
        adata.obs[col] = adata.obs[col].astype(str).astype("category")
    adata.obs[col] = adata.obs[col].cat.remove_unused_categories()
    levels = list(adata.obs[col].cat.categories)
    if len(levels) < 2:
        fail(ERRORS, "one_group", col=col, n=len(levels))
    if p["root_group"] not in levels:
        fail(ERRORS, "root_absent", root=p["root_group"], col=col)
    return adata, levels


def gene_trends(adata, genes, pseudotime):
    """Mean expression per equal-width pseudotime bin, with a t-based 95% CI."""
    reachable = np.isfinite(pseudotime)
    edges = np.linspace(0, pseudotime[reachable].max(), N_TREND_BINS + 1)
    bins = np.clip(np.digitize(pseudotime[reachable], edges) - 1, 0, N_TREND_BINS - 1)
    expr = adata[reachable, genes].X
    expr = expr.toarray() if sparse.issparse(expr) else np.asarray(expr)
    long = (
        pd.DataFrame(expr, columns=genes)
        .assign(bin=bins)
        .melt(id_vars="bin", var_name="gene")
    )
    trend = long.groupby(["gene", "bin"]).agg(
        mean=("value", "mean"), sd=("value", "std"), n=("value", "size")
    )
    half = stats.t.ppf(0.975, trend["n"] - 1) * trend["sd"] / np.sqrt(trend["n"])
    trend = trend.assign(lo=trend["mean"] - half, hi=trend["mean"] + half).reset_index()
    trend["pseudotime_mid"] = (edges[trend["bin"]] + edges[trend["bin"] + 1]) / 2
    trend["gene"] = pd.Categorical(trend["gene"], categories=genes)
    trend = trend.sort_values(["gene", "bin"]).reset_index(drop=True)
    return trend[["gene", "bin", "pseudotime_mid", "mean", "lo", "hi", "n"]]


def embedding_how(adata):
    """Caption text for the stored UMAP: neighbour parameters as sc-embed left them."""
    nn = adata.uns["neighbors"].get("params", {})
    pcs = nn.get("n_pcs", "all stored")
    return (
        f"UMAP as stored in the input ({nn.get('n_neighbors', '?')} neighbours, "
        f"{nn.get('metric', '?')} metric on {pcs} PCs, seed "
        f"{nn.get('random_state', '?')}; min_dist not recorded)"
    )


def draw_figures(figures, adata, p, levels, pal, cells, edges_long, trend, genes):
    col = p["groups"]
    umap = pd.DataFrame(adata.obsm["X_umap"][:, :2], columns=["UMAP 1", "UMAP 2"])
    umap[col] = adata.obs[col].to_numpy()

    # PAGA graph at group-median UMAP positions: no layout RNG, and it reads
    # against the UMAP panel (scv.pl.paga(basis="umap"), corpus A-pl-18).
    nodes = umap.groupby(col, observed=True)[["UMAP 1", "UMAP 2"]].median().loc[levels]
    sizes = adata.obs[col].value_counts().loc[levels]
    kept = edges_long[edges_long["above_threshold"]]
    fig, ax = plt.subplots(layout="constrained")
    for a, b, conn in kept[["group_a", "group_b", "connectivity"]].itertuples(
        index=False
    ):
        ax.plot(
            nodes.loc[[a, b], "UMAP 1"],
            nodes.loc[[a, b], "UMAP 2"],
            c="0.3",
            lw=MAX_EDGE_WIDTH * conn,
            zorder=1,
        )
    ax.scatter(
        nodes["UMAP 1"],
        nodes["UMAP 2"],
        s=20 + 80 * sizes / sizes.max(),
        c=[pal[lv] for lv in levels],
        zorder=2,
    )
    for lv, (x, y) in nodes.iterrows():
        ax.annotate(lv, (x, y), xytext=(4, 4), textcoords="offset points")
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    save_panel(fig, edges_long, figures, "paga_graph", 89, 89)
    plt.close(fig)

    umap[PSEUDOTIME] = cells[PSEUDOTIME].reindex(adata.obs_names).to_numpy()
    reachable = umap[np.isfinite(umap[PSEUDOTIME])]
    fig, ax = plt.subplots(layout="constrained")
    ax.scatter(umap["UMAP 1"], umap["UMAP 2"], s=1, c="0.85", rasterized=True)
    dots = ax.scatter(
        reachable["UMAP 1"],
        reachable["UMAP 2"],
        s=1,
        c=reachable[PSEUDOTIME],
        cmap="viridis",
        rasterized=True,
    )
    fig.colorbar(dots, ax=ax, label=PSEUDOTIME_LABEL)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    save_panel(fig, umap, figures, "umap_pseudotime", 89, 75)
    plt.close(fig)

    box = cells[np.isfinite(cells[PSEUDOTIME])][[col, PSEUDOTIME]]
    plot = plot_box(
        box,
        col,
        PSEUDOTIME,
        pal,
        jitter_width=JITTER_WIDTH,
        seed=p["seed"],
        value_label=PSEUDOTIME_LABEL,
    )
    plot += theme(axis_text_x=element_text(rotation=45, ha="right"))
    save_panel(plot, box, figures, "pseudotime_by_group", 89, 70)

    ncols = min(3, len(genes))
    nrows = -(-len(genes) // ncols)
    fig, axes = plt.subplots(
        nrows, ncols, sharex=True, layout="constrained", squeeze=False
    )
    for ax, gene in zip(axes.flat, genes, strict=False):
        g = trend[trend["gene"] == gene]
        ax.fill_between(
            g["pseudotime_mid"], g["lo"], g["hi"], color="0.6", alpha=0.4, lw=0
        )
        ax.plot(g["pseudotime_mid"], g["mean"], c="black", lw=0.8)
        ax.annotate(gene, (0.02, 0.9), xycoords="axes fraction")
        ax.set_ylabel("mean expression (X)")
    for ax in axes.flat[len(genes) :]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel(PSEUDOTIME_LABEL)
    save_panel(fig, trend, figures, "gene_trends", 60 * ncols, 45 * nrows)
    plt.close(fig)


def analyse(argv, full_argv):
    p = parse_args(argv)
    genes = check_values(p)
    adata, levels = load(p, genes)
    col = p["groups"]
    outdir = make_outdir(p["outdir"])
    (outdir / "data").mkdir()
    seed_everything(p["seed"])
    sc.settings.n_jobs = 1

    # Steps 3-5: PAGA on the groups, diffusion map, DPT from the root cell.
    sc.tl.paga(adata, groups=col)
    sc.tl.diffmap(adata, n_comps=p["n_dcs"], random_state=p["seed"])
    # Root rule: the first cell of --root-group in obs order (scanpy tutorial).
    root = int(np.flatnonzero(adata.obs[col] == p["root_group"])[0])
    adata.uns["iroot"] = root
    sc.tl.dpt(adata, n_dcs=p["n_dcs"])
    root_cell = adata.obs_names[root]

    cells = pd.DataFrame(
        {
            "cell": adata.obs_names,
            col: adata.obs[col].astype(str).to_numpy(),
            PSEUDOTIME: adata.obs["dpt_pseudotime"].to_numpy(),
            "is_root": adata.obs_names == root_cell,
        }
    ).sort_values("cell")
    cells.rename(columns={PSEUDOTIME: "dpt_pseudotime"}).to_csv(
        outdir / "tables" / "pseudotime.tsv", sep="\t", index=False
    )
    cells = cells.set_index("cell", drop=False)

    conn = adata.uns["paga"]["connectivities"].toarray()
    i, j = np.triu_indices(len(levels), k=1)
    edges_long = pd.DataFrame(
        {
            "group_a": np.asarray(levels)[i],
            "group_b": np.asarray(levels)[j],
            "connectivity": conn[i, j],
        }
    )
    edges_long["above_threshold"] = edges_long["connectivity"] >= p["paga_threshold"]
    edges_long.to_csv(
        outdir / "tables" / "paga_connectivities.tsv", sep="\t", index=False
    )
    adata.write_h5ad(outdir / "data" / "trajectory.h5ad")

    pseudotime = adata.obs["dpt_pseudotime"].to_numpy()
    n_unreachable = int((~np.isfinite(pseudotime)).sum())
    # Unreachable (inf) cells leave the box plot and trends; count them per group.
    dropped_per_group = {
        lv: int(n)
        for lv, n in (
            adata.obs.loc[~np.isfinite(pseudotime), col]
            .value_counts()
            .reindex(levels, fill_value=0)
            .items()
        )
    }
    trend = gene_trends(adata, genes, pseudotime)
    # ponytail: tab20 repeats colours past 20 groups; pass a palette when needed.
    tab20 = plt.get_cmap("tab20").colors
    pal = {lv: to_hex(tab20[k % 20]) for k, lv in enumerate(levels)}
    plt.switch_backend("Agg")
    draw_figures(
        outdir / "figures", adata, p, levels, pal, cells, edges_long, trend, genes
    )

    kept = edges_long[edges_long["above_threshold"]]
    medians = (
        cells[np.isfinite(cells[PSEUDOTIME])]
        .groupby(col)[PSEUDOTIME]
        .median()
        .loc[levels]
    )
    how = embedding_how(adata)
    reached = adata.obs.loc[np.isfinite(pseudotime), col].value_counts()
    # plots.plot_box draws a median tick, no box, under FEW_POINTS in any group.
    summary_mark = (
        "the bar is the group median"
        if reached[reached > 0].min() < FEW_POINTS
        else "boxes show median and quartiles, whiskers 1.5 IQR"
    )
    dpt_how = f"DPT on {p['n_dcs']} diffusion components, seed {p['seed']}"
    captions = {
        "paga_graph": f"PAGA graph of {len(levels)} {col} groups at their median "
        f"UMAP position; lines join groups with connectivity >= "
        f"{p['paga_threshold']}, width by connectivity, node size by group size. "
        f"{how}.",
        "umap_pseudotime": f"UMAP of {adata.n_obs} cells coloured by diffusion "
        f"pseudotime from root cell {root_cell} in {p['root_group']}; grey cells are "
        f"unreachable ({n_unreachable}). {how}; {dpt_how}.",
        "pseudotime_by_group": f"Diffusion pseudotime per {col} group; each point "
        f"is one cell, {summary_mark}; {n_unreachable} unreachable cells left out. "
        f"{dpt_how}.",
        "gene_trends": f"Mean expression (X) of {len(genes)} genes in {N_TREND_BINS} "
        "equal-width pseudotime bins; the band is the t-based 95% confidence interval "
        f"of each bin mean. {dpt_how}.",
    }
    summary = [
        f"{adata.n_obs} cells in {len(levels)} {col} groups; root cell {root_cell} "
        f"(first {p['root_group']} cell); {n_unreachable} cells unreachable.",
        f"{len(kept)} of {len(edges_long)} PAGA group pairs have connectivity >= "
        f"{p['paga_threshold']}.",
        "Median pseudotime by group: "
        + ", ".join(f"{lv} {v:.3f}" for lv, v in medians.items())
        + ".",
    ]
    keys = REQUIRED
    write_report(
        outdir, f"Trajectory: DPT from {p['root_group']}", p, keys, summary, captions
    )
    result = {
        "n_cells": adata.n_obs,
        "n_groups": len(levels),
        "root_cell": root_cell,
        "n_unreachable": n_unreachable,
        "dropped_per_group": dropped_per_group,
        "n_edges_above_threshold": len(kept),
        "edges_above_threshold": kept[["group_a", "group_b", "connectivity"]].to_dict(
            orient="records"
        ),
        "median_pseudotime": medians.round(6).to_dict(),
        "genes": genes,
        "params": p,
        "seed": p["seed"],
    }
    artifacts = ["report.md", "tables/pseudotime.tsv", "figures/umap_pseudotime.png"]
    write_result(
        outdir,
        result,
        summary,
        artifacts,
        context=p["context"],
        reference=p["reference"],
    )
    write_reproducibility(outdir, __file__, p, keys, full_argv, PACKAGES)
    return outdir


if __name__ == "__main__":
    raise SystemExit(run("trajectory", analyse))
