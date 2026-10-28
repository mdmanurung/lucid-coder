"""Purpose: normalise, HVG, PCA, integrate, neighbours, UMAP and Leiden for one h5ad.
Inputs:  --adata <file.h5ad> with QC'd raw integer counts in X (e.g. sc-qc output)
Outputs: <outdir>/report.md, result.json, data/embedded.h5ad, tables/, figures/,
         reproducibility/
Run:     python scripts/sc_embed.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""

import argparse
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator
from plots import plot_composition, plot_embedding
from scipy import sparse
from skillkit import (
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
from theme import save_panel
from threadpoolctl import threadpool_limits

SKILL_DIR = Path(__file__).resolve().parents[1]
ANALYSIS_ARGS = [
    "adata",
    "target_sum",
    "n_top_genes",
    "flavor",
    "batch_key",
    "n_comps",
    "integration",
    "n_neighbors",
    "n_pcs",
    "min_dist",
    "resolution",
    "seed",
    "reference",
    "context",
]
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "adata": str(SKILL_DIR / "examples" / "demo.h5ad"),
    "target_sum": 1e4,
    "n_top_genes": 100,
    "flavor": "seurat",
    "batch_key": "batch",
    "n_comps": 20,
    "integration": "harmony",
    "n_neighbors": 15,
    "n_pcs": 10,
    "min_dist": 0.3,
    "resolution": 0.3,
    "seed": 0,
    "reference": "synthetic demo genes, no genome build",
    "context": "synthetic demo: 2 batches x 3 groups, 300 genes, no real tissue",
    "outdir": "sc_embed_demo",
}

ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "harmony_needs_batch": "--integration harmony needs a --batch-key, got none",
    "missing_file": "input file not found: {path}",
    "missing_column": "obs column not found: {col}",
    "batch_na": "obs column {col!r} has missing values",
    "one_batch": "obs column {col!r} has one level; harmony needs at least 2",
    "bad_counts": "counts in X must be non-negative integers (raw counts)",
}

# Fixed steps, not analysis parameters: igraph Leiden with 2 iterations
# (C-an-05, C-an-10); Harmony stops at 20 rounds (B-an-06); `none` stands for a
# batch column when --batch-key none.
LEIDEN_ITERATIONS = 2
HARMONY_MAX_ITER = 20
NO_BATCH = "none"
PACKAGES = [
    "scanpy",
    "anndata",
    "harmonypy",
    "igraph",
    "umap-learn",
    "scikit-misc",
    "numpy",
    "pandas",
    "scipy",
    "matplotlib",
    "plotnine",
    "scikit-learn",
]
# Muted qualitative colours (Tol), assigned to labels in sorted order.
# ponytail: colours repeat past 10 labels; extend the list if runs need more.
COLOURS = [
    "#4477AA",
    "#EE6677",
    "#228833",
    "#CCBB44",
    "#66CCEE",
    "#AA3377",
    "#BBBBBB",
    "#332288",
    "#999933",
    "#882255",
]


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--adata")
    parser.add_argument("--target-sum", type=float)
    parser.add_argument("--n-top-genes", type=int)
    parser.add_argument("--flavor", choices=["seurat", "cell_ranger", "seurat_v3"])
    parser.add_argument("--batch-key")
    parser.add_argument("--n-comps", type=int)
    parser.add_argument("--integration", choices=["harmony", "none"])
    parser.add_argument("--n-neighbors", type=int)
    parser.add_argument("--n-pcs", type=int)
    parser.add_argument("--min-dist", type=float)
    parser.add_argument("--resolution", type=float)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--reference")
    parser.add_argument("--context")
    return parse_required(
        parser, ANALYSIS_ARGS, DEMO, argv, demo_extra=("outdir", "integration")
    )


def check_values(p):
    rules = [
        ("--target-sum", p["target_sum"], p["target_sum"] > 0, "> 0"),
        ("--n-top-genes", p["n_top_genes"], p["n_top_genes"] >= 2, "an integer >= 2"),
        ("--n-comps", p["n_comps"], p["n_comps"] >= 2, "an integer >= 2"),
        ("--n-neighbors", p["n_neighbors"], p["n_neighbors"] >= 2, "an integer >= 2"),
        ("--n-pcs", p["n_pcs"], 2 <= p["n_pcs"] <= p["n_comps"], "in [2, --n-comps]"),
        ("--min-dist", p["min_dist"], p["min_dist"] >= 0, ">= 0"),
        ("--resolution", p["resolution"], p["resolution"] > 0, "> 0"),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ]
    for arg, value, ok, rule in rules:
        if not ok:
            fail(ERRORS, "bad_value", arg=arg, rule=rule, value=value)
    if p["integration"] == "harmony" and p["batch_key"] == NO_BATCH:
        fail(ERRORS, "harmony_needs_batch")


def read_input(p):
    path = Path(p["adata"])
    if not path.is_file():
        fail(ERRORS, "missing_file", path=path)
    adata = ad.read_h5ad(path)
    check_unique(adata.var_names, "gene names")
    check_unique(adata.obs_names, "cell ids")
    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    if not (np.all(values >= 0) and np.all(values == np.round(values))):
        fail(ERRORS, "bad_counts")
    # Checked here, before the outdir exists, so a refused run leaves nothing.
    if p["n_top_genes"] > adata.n_vars:
        rule = f"<= the {adata.n_vars} genes"
        fail(
            ERRORS, "bad_value", arg="--n-top-genes", rule=rule, value=p["n_top_genes"]
        )
    n_max = min(adata.n_obs, p["n_top_genes"]) - 1
    if p["n_comps"] > n_max:
        rule = f"<= min(cells, --n-top-genes) - 1 = {n_max}"
        fail(ERRORS, "bad_value", arg="--n-comps", rule=rule, value=p["n_comps"])
    key = p["batch_key"]
    if key == NO_BATCH:
        adata.obs[NO_BATCH] = pd.Categorical([NO_BATCH] * adata.n_obs)
        return adata
    if key not in adata.obs.columns:
        fail(ERRORS, "missing_column", col=key)
    if adata.obs[key].isna().any():
        fail(ERRORS, "batch_na", col=key)
    adata.obs[key] = adata.obs[key].astype(str).astype("category")
    if p["integration"] == "harmony" and adata.obs[key].nunique() < 2:
        fail(ERRORS, "one_batch", col=key)
    return adata


def embed(adata, p):
    """Steps 2-8 on `adata` in place; PCA runs on a scaled HVG copy."""
    # Imported here: scanpy pulls in umap and numba (~20 s), which a refused run
    # does not need.
    import harmonypy
    import scanpy as sc

    seed, batch = p["seed"], p["batch_key"]
    sc.settings.n_jobs = 1
    hvg_batch = None if batch == NO_BATCH else batch
    adata.layers["counts"] = adata.X.copy()
    # seurat_v3 models raw counts, so it runs before normalisation; the others
    # rank genes on log-normalised data.
    if p["flavor"] == "seurat_v3":
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=p["n_top_genes"],
            flavor="seurat_v3",
            batch_key=hvg_batch,
            layer="counts",
        )
    sc.pp.normalize_total(adata, target_sum=p["target_sum"])
    sc.pp.log1p(adata)
    if p["flavor"] != "seurat_v3":
        sc.pp.highly_variable_genes(
            adata, n_top_genes=p["n_top_genes"], flavor=p["flavor"], batch_key=hvg_batch
        )

    hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hvg)
    sc.tl.pca(hvg, n_comps=p["n_comps"], svd_solver="arpack", random_state=seed)
    adata.obsm["X_pca"] = hvg.obsm["X_pca"]
    adata.uns["pca"] = hvg.uns["pca"]

    use_rep = "X_pca"
    if p["integration"] == "harmony":
        # Called directly: harmonypy 2.x returns Z_corr as cells x PCs, which
        # scanpy.external's harmony_integrate transposes the wrong way.
        harmony = harmonypy.run_harmony(
            adata.obsm["X_pca"],
            adata.obs,
            batch,
            max_iter_harmony=HARMONY_MAX_ITER,
            random_state=seed,
            ncores=1,
            verbose=False,
        )
        adata.obsm["X_pca_harmony"] = np.asarray(harmony.Z_corr)
        use_rep = "X_pca_harmony"

    sc.pp.neighbors(
        adata,
        n_neighbors=p["n_neighbors"],
        n_pcs=p["n_pcs"],
        use_rep=use_rep,
        random_state=seed,
    )
    sc.tl.umap(adata, min_dist=p["min_dist"], random_state=seed)
    sc.tl.leiden(
        adata,
        resolution=p["resolution"],
        flavor="igraph",
        n_iterations=LEIDEN_ITERATIONS,
        directed=False,
        random_state=seed,
        key_added="leiden",
    )
    clusters = sorted(adata.obs["leiden"].unique(), key=int)
    adata.obs["leiden"] = pd.Categorical(adata.obs["leiden"], categories=clusters)
    return use_rep


def palette(levels):
    return {lv: COLOURS[i % len(COLOURS)] for i, lv in enumerate(levels)}


def draw_figures(figures, adata, p, by_cluster):
    batch = p["batch_key"]
    ratio = adata.uns["pca"]["variance_ratio"]
    elbow = pd.DataFrame({"pc": np.arange(1, len(ratio) + 1), "variance_pct": ratio})
    elbow["variance_pct"] *= 100
    fig, ax = plt.subplots(layout="constrained")
    ax.plot(elbow["pc"], elbow["variance_pct"], c="black", lw=0.5, marker="o", ms=2)
    ax.axvline(p["n_pcs"], ls="--", lw=0.5, c="black")
    ax.set_ylim(bottom=0)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("principal component")
    ax.set_ylabel("variance explained (%)")
    save_panel(fig, elbow, figures, "pca_variance", 89, 60)
    plt.close(fig)

    umap = pd.DataFrame(adata.obsm["X_umap"], columns=["umap_1", "umap_2"])
    umap.insert(0, "cell", adata.obs_names)
    umap["batch"] = adata.obs[batch].to_numpy()
    umap["cluster"] = adata.obs["leiden"].to_numpy()
    pals = {
        "batch": palette(adata.obs[batch].cat.categories),
        "cluster": palette(adata.obs["leiden"].cat.categories),
    }
    for col, pal in pals.items():
        name = f"umap_{col}"
        plot = plot_embedding(umap, "umap_1", "umap_2", col, pal, axis_name="UMAP")
        save_panel(plot, umap, figures, name, 89, 75)

    plot = plot_composition(by_cluster, "cluster", "batch", "n_cells", pals["batch"])
    save_panel(plot, by_cluster, figures, "batch_composition", 89, 60)


def analyse(argv, argv_run):
    p = parse_args(argv)
    check_values(p)
    seed_everything(p["seed"])
    adata = read_input(p)
    batch = p["batch_key"]
    outdir = make_outdir(p["outdir"])
    # One thread for BLAS and numba-backed steps (n_jobs set in embed), so
    # floating-point sums run in one order on any host.
    with threadpool_limits(1):
        use_rep = embed(adata, p)

    data = outdir / "data"
    data.mkdir()
    adata.write_h5ad(data / "embedded.h5ad")

    clusters = adata.obs["leiden"].cat.categories
    sizes = adata.obs["leiden"].value_counts().reindex(clusters)
    cluster_sizes = pd.DataFrame(
        {"cluster": clusters, "n_cells": sizes.to_numpy(), "frac": sizes / adata.n_obs}
    )
    cluster_sizes.to_csv(outdir / "tables/cluster_sizes.tsv", sep="\t", index=False)
    by_cluster = (
        adata.obs.groupby(["leiden", batch], observed=False)
        .size()
        .rename("n_cells")
        .reset_index()
        .rename(columns={"leiden": "cluster", batch: "batch"})
    )
    by_cluster["frac_of_cluster"] = by_cluster["n_cells"] / by_cluster.groupby(
        "cluster", observed=False
    )["n_cells"].transform("sum")
    by_cluster.to_csv(outdir / "tables/batch_by_cluster.tsv", sep="\t", index=False)

    plt.switch_backend("Agg")
    draw_figures(outdir / "figures", adata, p, by_cluster)

    n_batches = adata.obs[batch].nunique()
    # fig.text: a UMAP caption says how it was made.
    how = (
        f"{p['n_pcs']} PCs of {use_rep}, {p['n_neighbors']} neighbours, "
        f"min_dist {p['min_dist']}, seed {p['seed']}, integration {p['integration']}"
    )
    captions = {
        "pca_variance": f"Line plot of variance explained by each of {p['n_comps']} "
        f"principal components of {p['n_top_genes']} scaled HVGs; the dashed line "
        f"marks the {p['n_pcs']} PCs used for the neighbour graph.",
        "umap_batch": f"UMAP of {adata.n_obs} cells coloured by {batch} "
        f"({n_batches} levels; {how}); well-mixed colours mean the embedding does "
        "not separate batches.",
        "umap_cluster": f"UMAP of {adata.n_obs} cells coloured by {len(clusters)} "
        f"Leiden clusters (resolution {p['resolution']}; {how}).",
        "batch_composition": f"Stacked bars of the percentage of each cluster's cells "
        f"from each {batch} level; the denominator is the cluster's cell count.",
    }
    summary = [
        f"{adata.n_obs} cells, {adata.n_vars} genes; {p['n_top_genes']} HVGs "
        f"({p['flavor']}), {p['n_comps']} PCs, {p['n_pcs']} used.",
        f"{len(clusters)} Leiden clusters at resolution {p['resolution']}, "
        f"integration {p['integration']} on {batch}.",
    ]
    write_report(outdir, "Embed and cluster", p, ANALYSIS_ARGS, summary, captions)
    write_result(
        outdir,
        {
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars,
            "n_hvg": int(adata.var["highly_variable"].sum()),
            "n_clusters": len(clusters),
            "cluster_sizes": dict(zip(clusters, sizes.tolist(), strict=True)),
            "use_rep": use_rep,
            "integration": p["integration"],
            "seed": p["seed"],
            "params": p,
        },
        summary,
        ["report.md", "figures/umap_cluster.png", "tables/cluster_sizes.tsv"],
        context=p["context"],
        reference=p["reference"],
    )
    write_reproducibility(outdir, __file__, p, ANALYSIS_ARGS, argv_run, PACKAGES)
    return outdir


if __name__ == "__main__":
    raise SystemExit(run("sc-embed", analyse))
