"""Purpose: Milo differential abundance of kNN neighbourhoods with pertpy (pydeseq2).
Inputs:  --adata <file.h5ad> with an embedding, UMAP and kNN graph (sc-embed output)
Outputs: <outdir>/report.md, result.json, data/milo.h5mu, tables/, figures/,
         reproducibility/
Run:     python scripts/milo_da.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""

import argparse
import re
import shutil
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from plots import diverging
from scipy import sparse
from skillkit import (
    SkillError,
    check_design_rank,
    check_replication,
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
    "context",
    "sample_col",
    "donor_col",
    "design",
    "model_contrasts",
    "prop",
    "use_rep",
    "n_neighbors",
    "spatial_fdr_max",
    "annotation_col",
    "min_annotation_frac",
    "seed",
]
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "adata": str(SKILL_DIR / "examples" / "demo.h5ad"),
    "context": "synthetic demo: 6 samples, condition B vs A, "
    "4 Gaussian-blob cell types, no real tissue",
    "sample_col": "sample",
    "donor_col": "donor",
    "design": "~condition",
    "model_contrasts": "conditionB-conditionA",
    "prop": 0.1,
    "use_rep": "X_pca",
    "n_neighbors": 30,
    "spatial_fdr_max": 0.1,
    "annotation_col": "cell_type",
    "min_annotation_frac": 0.6,
    "seed": 0,
    "outdir": "milo_da_demo",
}

ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "missing_file": "input file not found: {path}",
    "missing_column": "obs column not found: {col}",
    "column_na": "obs column {col!r} has missing values",
    "missing_rep": "obsm key not found: {name}",
    "not_nested": "design column {col!r} takes more than one value within a "
    "{sample} level; each sample needs one value",
    "model_refused": "pertpy da_nhoods refused the model: {message}",
}

# Fixed, not analysis parameters: the label for nhoods under --min-annotation-frac
# (milopy annotate_nhoods), the layout basis of the nhood graph, the solver.
MIXED = "Mixed"
BASIS = "X_umap"
SOLVER = "pydeseq2"
NOT_SIGNIFICANT = "#D9D9D9"
PACKAGES = [
    "pertpy",
    "pydeseq2",
    "scanpy",
    "anndata",
    "mudata",
    "numpy",
    "pandas",
    "scipy",
    "matplotlib",
    "threadpoolctl",
]


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--adata")
    parser.add_argument("--context")
    parser.add_argument("--sample-col")
    parser.add_argument("--donor-col")
    parser.add_argument("--design")
    parser.add_argument("--model-contrasts")
    parser.add_argument("--prop", type=float)
    parser.add_argument("--use-rep")
    parser.add_argument("--n-neighbors", type=int)
    parser.add_argument("--spatial-fdr-max", type=float)
    parser.add_argument("--annotation-col")
    parser.add_argument("--min-annotation-frac", type=float)
    parser.add_argument("--seed", type=int)
    return parse_required(
        parser, ANALYSIS_ARGS, DEMO, argv, demo_extra=("outdir", "n_neighbors")
    )


def design_columns(design):
    """obs columns named by a `~a + b` or `~a * b` formula, as pertpy splits it."""
    return [x.strip() for x in re.split(r"[+*:]", design.lstrip("~ ")) if x.strip()]


def check_values(p):
    rules = [
        ("--design", p["design"], p["design"].startswith("~"), "a formula like ~x"),
        ("--prop", p["prop"], 0 < p["prop"] <= 1, "in (0, 1]"),
        ("--n-neighbors", p["n_neighbors"], p["n_neighbors"] >= 2, "an integer >= 2"),
        (
            "--spatial-fdr-max",
            p["spatial_fdr_max"],
            0 < p["spatial_fdr_max"] < 1,
            "in (0, 1)",
        ),
        (
            "--min-annotation-frac",
            p["min_annotation_frac"],
            0 < p["min_annotation_frac"] <= 1,
            "in (0, 1]",
        ),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ]
    for arg, value, ok, rule in rules:
        if not ok:
            fail(ERRORS, "bad_value", arg=arg, rule=rule, value=value)


def read_input(p):
    path = Path(p["adata"])
    if not path.is_file():
        fail(ERRORS, "missing_file", path=path)
    adata = ad.read_h5ad(path)
    for key in [p["use_rep"], BASIS]:
        if key not in adata.obsm:
            fail(ERRORS, "missing_rep", name=key)
    sample = p["sample_col"]
    covariates = design_columns(p["design"])
    check_unique(adata.obs_names, "cell ids")
    for col in [sample, p["donor_col"], p["annotation_col"], *covariates]:
        if col not in adata.obs.columns:
            fail(ERRORS, "missing_column", col=col)
        if adata.obs[col].isna().any():
            fail(ERRORS, "column_na", col=col)
    # Each sample is one replicate: a covariate that varies inside a sample would
    # be silently collapsed by count_nhoods.
    per_sample = adata.obs.groupby(sample, observed=True)[covariates].nunique()
    for col in covariates:
        if (per_sample[col] > 1).any():
            fail(ERRORS, "not_nested", col=col, sample=sample)
    # The per-sample rows are what the NB GLM fits: check rank and replication there.
    per_sample_meta = adata.obs.drop_duplicates(sample)
    check_design_rank(per_sample_meta, covariates)
    n_donors, _ = check_replication(adata.obs, sample, p["donor_col"], covariates)
    if p["n_neighbors"] >= adata.n_obs:
        rule = f"< the {adata.n_obs} cells"
        fail(
            ERRORS, "bad_value", arg="--n-neighbors", rule=rule, value=p["n_neighbors"]
        )
    for col in [sample, p["annotation_col"]]:
        adata.obs[col] = adata.obs[col].astype(str).astype("category")
    return adata, n_donors


def graph_matches(adata, p):
    params = adata.uns.get("neighbors", {}).get("params", {})
    return (
        "connectivities" in adata.obsp
        and params.get("use_rep") == p["use_rep"]
        and params.get("n_neighbors") == p["n_neighbors"]
    )


def run_milo(adata, p):
    """Steps 2-6; returns the MuData and whether the kNN graph was rebuilt."""
    # Imported here: pertpy and scanpy take ~20 s, which a refused run does not need.
    import pertpy as pt
    import scanpy as sc

    rebuilt = not graph_matches(adata, p)
    if rebuilt:
        sc.pp.neighbors(
            adata,
            n_neighbors=p["n_neighbors"],
            use_rep=p["use_rep"],
            random_state=p["seed"],
        )
    milo = pt.tl.Milo()
    mdata = milo.load(adata)
    milo.make_nhoods(mdata["rna"], prop=p["prop"], seed=p["seed"])
    mdata = milo.count_nhoods(mdata, sample_col=p["sample_col"])
    try:
        milo.da_nhoods(
            mdata,
            design=p["design"],
            model_contrasts=p["model_contrasts"],
            solver=SOLVER,
        )
    except ValueError as exc:
        fail(ERRORS, "model_refused", message=str(exc).splitlines()[0])
    milo.annotate_nhoods(mdata, anno_col=p["annotation_col"])
    milo.build_nhood_graph(mdata, basis=BASIS)
    var = mdata["milo"].var
    var["nhood_annotation"] = pd.Categorical(var["nhood_annotation"].astype(str))
    var["nhood_annotation"] = var["nhood_annotation"].cat.add_categories(MIXED)
    var.loc[
        var["nhood_annotation_frac"] < p["min_annotation_frac"], "nhood_annotation"
    ] = MIXED
    # count_nhoods names the sample index after --sample-col, which h5mu cannot
    # write next to the column of the same name.
    mdata["milo"].obs.index.name = None
    return mdata, rebuilt


def draw_figures(figures, nhoods, mdata, p):
    fdr = p["spatial_fdr_max"]
    sig = nhoods["SpatialFDR"] < fdr
    limit = float(nhoods.loc[sig, "logFC"].abs().max()) if sig.any() else 1.0
    colour = np.where(sig, nhoods["logFC"], np.nan)

    graph = sparse.triu(mdata["milo"].varp["nhood_connectivities"]).tocoo()
    xy = nhoods[["umap_1", "umap_2"]].to_numpy()
    # nhoods rows are sorted by significance; the graph is indexed by nhood id.
    by_id = pd.DataFrame(xy, index=nhoods["nhood"]).sort_index().to_numpy()
    fig, ax = plt.subplots(layout="constrained")
    edges = np.stack([by_id[graph.row], by_id[graph.col]], axis=1)
    ax.add_collection(LineCollection(edges, colors="0.7", linewidths=0.2, zorder=1))
    area = 40 * nhoods["size"] / nhoods["size"].max()
    ax.scatter(*xy[~sig].T, s=area[~sig], c=NOT_SIGNIFICANT, lw=0, zorder=2)
    dots = ax.scatter(
        *xy[sig].T,
        s=area[sig],
        c=colour[sig],
        cmap=diverging,
        vmin=-limit,
        vmax=limit,
        lw=0,
        zorder=3,
    )
    fig.colorbar(
        dots, ax=ax, shrink=0.7, label=f"log2 fold change (SpatialFDR < {fdr})"
    )
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    save_panel(fig, nhoods, figures, "nhood_graph", 120, 90)
    plt.close(fig)

    # Annotations ordered by median logFC, as milopy's beeswarm orders them.
    order = nhoods.groupby("annotation", observed=True)["logFC"].median().sort_values()
    levels = list(order.index)
    y = pd.Categorical(nhoods["annotation"], categories=levels).codes
    jitter = np.random.default_rng(p["seed"]).uniform(-0.3, 0.3, len(nhoods))
    fig, ax = plt.subplots(layout="constrained")
    ax.axvline(0, ls="--", lw=0.5, c="black")
    ax.scatter(nhoods.loc[~sig, "logFC"], (y + jitter)[~sig], s=4, c=NOT_SIGNIFICANT)
    dots = ax.scatter(
        nhoods.loc[sig, "logFC"],
        (y + jitter)[sig],
        s=4,
        c=colour[sig],
        cmap=diverging,
        vmin=-limit,
        vmax=limit,
    )
    fig.colorbar(dots, ax=ax, label=f"log2 fold change (SpatialFDR < {fdr})")
    counts = nhoods["annotation"].value_counts()
    ax.set_yticks(range(len(levels)), [f"{lv} (n = {counts[lv]})" for lv in levels])
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("neighbourhood label")
    save_panel(fig, nhoods, figures, "da_beeswarm", 120, 75)
    plt.close(fig)

    fig, ax = plt.subplots(layout="constrained")
    ax.hist(nhoods["size"], bins="auto", color="0.4")
    ax.set_xlabel("cells per neighbourhood")
    ax.set_ylabel("neighbourhoods")
    save_panel(fig, nhoods[["nhood", "size"]], figures, "nhood_size", 89, 60)
    plt.close(fig)


def analyse(argv, argv_run):
    p = parse_args(argv)
    check_values(p)
    seed_everything(p["seed"])
    adata, n_donors = read_input(p)
    outdir = make_outdir(p["outdir"])
    # One BLAS thread, so sums run in one order. pydeseq2's loky workers still use
    # every CPU (pertpy passes no n_cpus); each nhood is fitted alone, so the
    # result does not depend on the worker count.
    with threadpool_limits(1):
        try:
            mdata, rebuilt = run_milo(adata, p)
        except SkillError:
            # A refused model leaves no half-written outdir behind.
            shutil.rmtree(outdir)
            raise

    data = outdir / "data"
    data.mkdir()
    mdata.write_h5mu(data / "milo.h5mu")

    var = mdata["milo"].var
    nhoods = (
        pd.DataFrame(
            {
                "nhood": var.index.astype(int).to_numpy(),
                "index_cell": var["index_cell"].to_numpy(),
                "size": var["Nhood_size"].astype(int).to_numpy(),
                "logFC": var["logFC"].to_numpy(),
                "PValue": var["PValue"].to_numpy(),
                "SpatialFDR": var["SpatialFDR"].to_numpy(),
                "annotation": var["nhood_annotation"].astype(str).to_numpy(),
                "annotation_frac": var["nhood_annotation_frac"].to_numpy(),
                "umap_1": mdata["milo"].varm["X_milo_graph"][:, 0],
                "umap_2": mdata["milo"].varm["X_milo_graph"][:, 1],
            }
        )
        .sort_values(["SpatialFDR", "PValue", "nhood"])
        .reset_index(drop=True)
    )
    nhoods.drop(columns=["umap_1", "umap_2"]).to_csv(
        outdir / "tables/nhoods.tsv", sep="\t", index=False
    )
    fdr = p["spatial_fdr_max"]
    sig = nhoods["SpatialFDR"] < fdr
    by_annotation = (
        nhoods.assign(up=sig & (nhoods["logFC"] > 0), down=sig & (nhoods["logFC"] < 0))
        .groupby("annotation")
        .agg(n_nhoods=("nhood", "size"), n_up=("up", "sum"), n_down=("down", "sum"))
        .reset_index()
        .sort_values("annotation")
    )
    by_annotation.to_csv(outdir / "tables/da_by_annotation.tsv", sep="\t", index=False)

    plt.switch_backend("Agg")
    draw_figures(outdir / "figures", nhoods, mdata, p)

    n_sig, n_up = int(sig.sum()), int((sig & (nhoods["logFC"] > 0)).sum())
    sig_counts = nhoods.loc[sig, "annotation"].value_counts()
    # Ties go to the alphabetically first label, so the summary is stable.
    dominant = (
        sorted(sig_counts[sig_counts == sig_counts.max()].index)[0] if n_sig else None
    )
    n_samples = adata.obs[p["sample_col"]].nunique()
    graph = (
        f"rebuilt with k = {p['n_neighbors']} on all columns of {p['use_rep']}"
        if rebuilt
        else f"kept from the input (k = {p['n_neighbors']}, {p['use_rep']})"
    )
    captions = {
        "nhood_graph": f"Neighbourhood graph on the input UMAP (obsm {BASIS}): one "
        "dot per neighbourhood at its index cell, area by cell count, edges for "
        "shared cells. Colour is "
        f"log2 fold change ({p['model_contrasts']}) where SpatialFDR < {fdr}; "
        "grey dots are not significant.",
        "da_beeswarm": f"Jittered strip of neighbourhood log2 fold change "
        f"({p['model_contrasts']}) per {p['annotation_col']} label, ordered by "
        f"median; coloured where SpatialFDR < {fdr}. {MIXED} holds neighbourhoods "
        f"whose top label is under {p['min_annotation_frac']:.0%} of cells.",
        "nhood_size": f"Histogram of cells per neighbourhood for {len(nhoods)} "
        f"neighbourhoods (prop {p['prop']}, k = {p['n_neighbors']}).",
    }
    summary = [
        f"{adata.n_obs} cells in {n_samples} samples from {n_donors} donors; "
        f"{len(nhoods)} neighbourhoods "
        f"(prop {p['prop']}, seed {p['seed']}), median size "
        f"{int(nhoods['size'].median())} cells.",
        f"kNN graph {graph}.",
        f"{n_sig} neighbourhoods at SpatialFDR < {fdr} ({n_up} up, "
        f"{n_sig - n_up} down) for {p['model_contrasts']}"
        + (f"; most are {dominant}." if n_sig else "."),
    ]
    write_report(
        outdir, "Milo differential abundance", p, ANALYSIS_ARGS, summary, captions
    )
    write_result(
        outdir,
        {
            "n_cells": adata.n_obs,
            "n_samples": n_samples,
            "n_donors": n_donors,
            "n_nhoods": len(nhoods),
            "n_significant": n_sig,
            "n_up": n_up,
            "n_down": n_sig - n_up,
            "dominant_annotation": dominant,
            "significant_by_annotation": sig_counts.sort_index().to_dict(),
            "graph_rebuilt": rebuilt,
            "solver": SOLVER,
            "seed": p["seed"],
            "params": p,
        },
        summary,
        ["report.md", "figures/nhood_graph.png", "tables/da_by_annotation.tsv"],
        context=p["context"],
    )
    write_reproducibility(outdir, __file__, p, ANALYSIS_ARGS, argv_run, PACKAGES)
    return outdir


if __name__ == "__main__":
    raise SystemExit(run("milo-da", analyse))
