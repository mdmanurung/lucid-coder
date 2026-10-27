"""Purpose: gene-program discovery by consensus NMF (cNMF-style) on raw counts.
Inputs:  --adata <file.h5ad> with raw counts in X or --layer; obs column --groupby
         (e.g. sc-qc's filtered h5ad, or sc-embed's data/embedded.h5ad with counts)
Outputs: <outdir>/report.md, result.json, tables/, figures/, data/programs.h5ad,
         reproducibility/
Run:     python scripts/gene_programs.py --demo --outdir demo_out
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
from plots import FEW_POINTS, plot_box, plot_heatmap
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
from sklearn.cluster import KMeans
from sklearn.decomposition import NMF, non_negative_factorization
from sklearn.metrics import pairwise_distances, roc_auc_score, silhouette_samples
from theme import save_panel
from threadpoolctl import threadpool_limits

SKILL_DIR = Path(__file__).resolve().parents[1]

REQUIRED = [
    "adata",
    "groupby",
    "reference",
    "k",
    "n_hvg",
    "n_replicates",
    "density_threshold",
    "active_threshold",
    "n_top_genes",
    "seed",
]
KEYS = [*REQUIRED[:2], "layer", *REQUIRED[2:]]
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "adata": str(SKILL_DIR / "examples" / "demo.h5ad"),
    "groupby": "group",
    "reference": "synthetic demo, no genome build",
    "layer": None,
    "k": 3,
    "n_hvg": 200,
    "n_replicates": 20,
    "density_threshold": 0.5,
    "active_threshold": 0.1,
    "n_top_genes": 10,
    "seed": 0,
    "outdir": "gene_programs_demo",
}

ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "missing_file": "input file not found: {path}",
    "missing_layer": "layer not found in the h5ad: {layer}",
    "missing_column": "obs column not found: {col}",
    "not_counts": "{source} does not look like raw counts: {why}",
    "few_groups": "--groupby {col} needs at least 2 groups, found {n}",
    "too_many_hvg": "--n-hvg {n_hvg} exceeds the {n_expressed} expressed genes",
    "zero_cells": "{n} cells have zero counts over the {n_hvg} HVGs; filter them first",
    "too_few_kept": "{kept} of {total} replicate components pass "
    "--density-threshold {threshold}; need at least --k {k}",
}

# Fixed method constants (cNMF), not analysis parameters.
NEIGHBOUR_FRACTION = 0.3  # cNMF --local-neighborhood-size
NMF_MAX_ITER = 1000
KMEANS_N_INIT = 10
TPM_SCALE = 1e6
# Fixed plot spec: heatmap saturation (z units), box jitter (box widths).
HEATMAP_LIMIT = 2.0
JITTER_WIDTH = 0.2
AUC_LABEL = "AUC, one vs rest (0 to 1)"
# Muted qualitative colours (Tol), assigned to labels in order.
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
    parser.add_argument("--groupby")
    parser.add_argument("--reference")
    parser.add_argument("--layer")
    parser.add_argument("--k", type=int)
    parser.add_argument("--n-hvg", type=int)
    parser.add_argument("--n-replicates", type=int)
    parser.add_argument("--density-threshold", type=float)
    parser.add_argument("--active-threshold", type=float)
    parser.add_argument("--n-top-genes", type=int)
    parser.add_argument("--seed", type=int)
    p = parse_required(parser, REQUIRED, DEMO, argv, optional=["layer"])
    rules = [
        ("--k", p["k"], p["k"] >= 2, "an integer >= 2"),
        ("--n-hvg", p["n_hvg"], p["n_hvg"] > p["k"], "an integer > --k"),
        # int(0.3 * n) neighbours must be at least 1.
        ("--n-replicates", p["n_replicates"], p["n_replicates"] >= 4, ">= 4"),
        (
            "--density-threshold",
            p["density_threshold"],
            p["density_threshold"] > 0,
            "> 0",
        ),
        (
            "--active-threshold",
            p["active_threshold"],
            0 < p["active_threshold"] < 1,
            "in (0, 1)",
        ),
        ("--n-top-genes", p["n_top_genes"], p["n_top_genes"] >= 1, "an integer >= 1"),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ]
    for arg, value, ok, rule in rules:
        if not ok:
            fail(ERRORS, "bad_value", arg=arg, rule=rule, value=value)
    return p


def load(p):
    """Step 1: read the h5ad, check raw counts and the groups."""
    path = Path(p["adata"])
    if not path.is_file():
        fail(ERRORS, "missing_file", path=path)
    adata = ad.read_h5ad(path)
    check_unique(adata.var_names, "gene names")
    check_unique(adata.obs_names, "cell ids")
    if p["layer"] is not None and p["layer"] not in adata.layers:
        fail(ERRORS, "missing_layer", layer=p["layer"])
    col = p["groupby"]
    if col not in adata.obs.columns:
        fail(ERRORS, "missing_column", col=col)

    source = f"layer {p['layer']}" if p["layer"] else "X"
    counts = sparse.csr_matrix(adata.layers[p["layer"]] if p["layer"] else adata.X)
    counts = counts.astype(np.float64)
    if (counts.data < 0).any():
        fail(ERRORS, "not_counts", source=source, why="negative values")
    if not np.array_equal(counts.data, np.round(counts.data)):
        fail(
            ERRORS, "not_counts", source=source, why="non-integer values (normalised?)"
        )

    # Categorical order is the group order of every table and figure.
    groups = adata.obs[col]
    if not isinstance(groups.dtype, pd.CategoricalDtype):
        groups = pd.Categorical(
            groups.astype(str), categories=sorted(set(groups.astype(str)))
        )
    adata.obs[col] = pd.Categorical(groups).remove_unused_categories()
    n_groups = len(adata.obs[col].cat.categories)
    if n_groups < 2:
        fail(ERRORS, "few_groups", col=col, n=n_groups)
    return adata, counts


def select_hvg(adata, counts, p):
    """Step 2: top --n-hvg expressed genes by Seurat dispersion (log-normalised)."""
    expressed = np.asarray(counts.sum(axis=0)).ravel() > 0
    if p["n_hvg"] > expressed.sum():
        fail(ERRORS, "too_many_hvg", n_hvg=p["n_hvg"], n_expressed=int(expressed.sum()))
    lognorm = ad.AnnData(counts[:, expressed].copy())
    sc.pp.normalize_total(lognorm, target_sum=1e4)
    sc.pp.log1p(lognorm)
    sc.pp.highly_variable_genes(lognorm, flavor="seurat", n_top_genes=p["n_hvg"])
    hvg = np.flatnonzero(expressed)[lognorm.var["highly_variable"].to_numpy()]
    hvg_counts = counts[:, hvg]
    zero = np.asarray(hvg_counts.sum(axis=1)).ravel() == 0
    if zero.any():
        fail(ERRORS, "zero_cells", n=int(zero.sum()), n_hvg=len(hvg))
    # cNMF input: HVG counts scaled to unit variance per gene, not centred.
    mean = np.asarray(hvg_counts.mean(axis=0)).ravel()
    sd = np.sqrt(np.asarray(hvg_counts.power(2).mean(axis=0)).ravel() - mean**2)
    return hvg, sparse.csr_matrix(hvg_counts.multiply(1 / sd))


def replicate_spectra(x, p):
    """Step 3: --n-replicates NMF fits, one seed each, drawn from --seed."""
    seeds = np.random.default_rng(p["seed"]).integers(1, 2**31 - 1, p["n_replicates"])
    spectra = [
        NMF(
            n_components=p["k"],
            init="random",
            solver="cd",
            max_iter=NMF_MAX_ITER,
            random_state=int(s),
        )
        .fit(x)
        .components_
        for s in seeds
    ]
    return np.vstack(spectra), np.repeat(np.arange(p["n_replicates"]), p["k"])


def consensus(spectra, replicate, p):
    """Step 4: local-density filter, k-means of components, median per cluster."""
    unit = spectra / np.linalg.norm(spectra, axis=1, keepdims=True)
    n_neighbours = int(NEIGHBOUR_FRACTION * p["n_replicates"])
    dist = np.sort(pairwise_distances(unit), axis=1)
    # Column 0 is the distance to itself.
    density = dist[:, 1 : n_neighbours + 1].mean(axis=1)
    kept = density < p["density_threshold"]
    if kept.sum() < p["k"]:
        fail(
            ERRORS,
            "too_few_kept",
            kept=int(kept.sum()),
            total=len(kept),
            threshold=p["density_threshold"],
            k=p["k"],
        )
    labels = KMeans(
        n_clusters=p["k"], n_init=KMEANS_N_INIT, random_state=p["seed"]
    ).fit_predict(unit[kept])
    centres = np.vstack(
        [np.median(unit[kept][labels == c], axis=0) for c in range(p["k"])]
    )
    centres /= centres.sum(axis=1, keepdims=True)
    silhouette = silhouette_samples(unit[kept], labels)
    stability = pd.DataFrame(
        {
            "cluster": labels,
            "replicate": replicate[kept],
            "silhouette": silhouette,
        }
    )
    return centres, stability, density, kept


def fit_usage(x, centres, p):
    """Step 5: usage with the consensus spectra fixed, then per-cell sum to 1."""
    usage, _, _ = non_negative_factorization(
        x,
        H=centres,
        n_components=p["k"],
        update_H=False,
        solver="cd",
        max_iter=NMF_MAX_ITER,
    )
    return usage / usage.sum(axis=1, keepdims=True)


def gene_scores(counts, usage):
    """Step 7: cNMF gene_spectra_score, OLS of z-scored TPM on usage, all genes.

    Z = (T - mean) / sd, scores = pinv(U) @ Z, computed without densifying T.
    Unexpressed genes (sd 0) get no score.
    """
    tpm = sparse.csr_matrix(counts.multiply(TPM_SCALE / counts.sum(axis=1)))
    mean = np.asarray(tpm.mean(axis=0)).ravel()
    sd = np.sqrt(np.asarray(tpm.power(2).mean(axis=0)).ravel() - mean**2)
    pinv = np.linalg.pinv(usage)
    projected = (tpm.T @ pinv.T).T
    with np.errstate(divide="ignore", invalid="ignore"):
        scores = (projected - np.outer(pinv.sum(axis=1), mean)) / sd
    scores[:, sd == 0] = np.nan
    return scores


def summarise_groups(usage_long, col, programs, groups):
    """Mean usage, fraction active and one-vs-rest AUC per program x group."""
    rows = []
    for program in programs:
        one = usage_long[usage_long["program"] == program]
        for group in groups:
            in_group = (one[col] == group).to_numpy()
            rows.append(
                {
                    "program": program,
                    col: group,
                    "n_cells": int(in_group.sum()),
                    "mean_usage": one.loc[in_group, "usage"].mean(),
                    "frac_active": one.loc[in_group, "active"].mean(),
                    "auc": roc_auc_score(in_group, one["usage"]),
                }
            )
    by_group = pd.DataFrame(rows)
    by_group["program"] = pd.Categorical(by_group["program"], categories=programs)
    by_group[col] = pd.Categorical(by_group[col], categories=groups)
    return by_group.sort_values(["program", col]).reset_index(drop=True)


def draw_figures(
    figures, usage_long, by_group, top, loadings, col, programs, groups, p
):
    pal_program = dict(zip(programs, COLOURS * len(programs), strict=False))
    pal_group = {g: COLOURS[i % len(COLOURS)] for i, g in enumerate(groups)}

    # Usage heatmap: group means centred and scaled per program (D-pl-34).
    heat = by_group[["program", col, "mean_usage"]].copy()
    means = heat.groupby("program", observed=True)["mean_usage"]
    heat["z-score of mean usage"] = (
        heat["mean_usage"] - means.transform("mean")
    ) / means.transform("std")
    fig = plot_heatmap(
        heat, "program", col, "z-score of mean usage", limit=HEATMAP_LIMIT
    )
    fig.set_layout_engine("tight")
    width = min(183, max(89, 50 + 6 * len(groups)))
    height = max(50, 25 + 6 * len(programs))
    save_panel(fig, heat, figures, "usage_heatmap", width, height)
    plt.close(fig)

    # Structure plot: cells by group, then dominant program, then its usage.
    wide = usage_long.pivot(index="cell", columns="program", values="usage")[programs]
    meta = usage_long.drop_duplicates("cell").set_index("cell")[[col]]
    meta["dominant"] = pd.Categorical(wide.idxmax(axis=1), categories=programs)
    meta["top_usage"] = wide.max(axis=1)
    order = meta.sort_values(
        [col, "dominant", "top_usage"], ascending=[True, True, False]
    )
    wide = wide.loc[order.index]
    fig, ax = plt.subplots()
    ax.stackplot(
        np.arange(len(wide)),
        wide.T.to_numpy(),
        labels=programs,
        colors=[pal_program[g] for g in programs],
        lw=0,
    )
    sizes = order[col].value_counts(sort=False).reindex(groups).to_numpy()
    edges = np.concatenate([[0], np.cumsum(sizes)])
    for edge in edges[1:-1]:
        ax.axvline(edge - 0.5, c="black", lw=0.5)
    ax.set_xticks((edges[:-1] + edges[1:]) / 2 - 0.5, groups, rotation=90)
    ax.set_xlim(-0.5, len(wide) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_xlabel(f"cells by {col}")
    ax.set_ylabel("program usage (fraction of cell)")
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    fig.set_layout_engine("tight")
    structure = wide.reset_index().melt(
        id_vars="cell", var_name="program", value_name="usage"
    )
    save_panel(fig, structure, figures, "usage_structure", 183, 60)
    plt.close(fig)

    # Top-gene heatmap: scores z-scored per gene across programs.
    genes = list(dict.fromkeys(top["gene"]))
    loads = (
        top[["gene"]]
        .drop_duplicates()
        .merge(loadings, on="gene", validate="1:m")
        .pivot(index="program", columns="gene", values="score")
    )
    loads = loads.reindex(index=programs, columns=genes)
    z = (loads - loads.mean()) / loads.std()
    tg = z.reset_index().melt(
        id_vars="program", var_name="gene", value_name="z-score of gene score"
    )
    tg["program"] = pd.Categorical(tg["program"], categories=programs)
    tg["gene"] = pd.Categorical(tg["gene"], categories=genes)
    fig = plot_heatmap(
        tg, "program", "gene", "z-score of gene score", limit=HEATMAP_LIMIT
    )
    fig.set_layout_engine("tight")
    save_panel(
        fig,
        tg,
        figures,
        "top_gene_heatmap",
        min(183, max(89, 40 + 3 * len(genes))),
        height,
    )
    plt.close(fig)

    # AUC per program by group (D-pl-30): one point per program.
    auc = by_group[["program", col, "auc"]].copy()
    auc[col] = auc[col].astype(str)
    fig = plot_box(
        auc,
        col,
        "auc",
        pal_group,
        jitter_width=JITTER_WIDTH,
        seed=p["seed"],
        value_label=AUC_LABEL,
    )
    fig.axes[0].axhline(0.5, ls="--", lw=0.5, c="black")
    fig.set_layout_engine("tight")
    save_panel(
        fig, auc, figures, "auc_by_group", min(183, max(89, 40 + 15 * len(groups))), 60
    )
    plt.close(fig)


def analyse(argv, full_argv):
    p = parse_args(argv)
    seed_everything(p["seed"])
    sc.settings.n_jobs = 1
    adata, counts = load(p)
    col = p["groupby"]
    groups = list(adata.obs[col].cat.categories)

    with threadpool_limits(1):
        hvg, x = select_hvg(adata, counts, p)
        spectra, replicate = replicate_spectra(x, p)
        centres, members, density, kept = consensus(spectra, replicate, p)
        usage = fit_usage(x, centres, p)
        # Step 6: name programs GP1.. by total usage, largest first.
        rank = np.argsort(-usage.sum(axis=0), kind="stable")
        usage, programs = usage[:, rank], [f"GP{i + 1}" for i in range(p["k"])]
        cluster_name = {int(c): programs[i] for i, c in enumerate(rank)}
        scores = gene_scores(counts, usage)

    # Active: usage over the program's max across cells (immgenT Figure2.R).
    normalised = usage / usage.max(axis=0, keepdims=True)
    usage_long = pd.DataFrame(
        {
            "cell": np.repeat(adata.obs_names.to_numpy(), p["k"]),
            col: np.repeat(adata.obs[col].astype(str).to_numpy(), p["k"]),
            "program": np.tile(programs, adata.n_obs),
            "usage": usage.ravel(),
            "usage_max_scaled": normalised.ravel(),
        }
    )
    usage_long["active"] = usage_long["usage_max_scaled"] > p["active_threshold"]
    usage_long["program"] = pd.Categorical(usage_long["program"], categories=programs)
    usage_long = usage_long.sort_values(["program", "cell"]).reset_index(drop=True)

    is_hvg = np.zeros(adata.n_vars, dtype=bool)
    is_hvg[hvg] = True
    loadings = pd.DataFrame(scores.T, index=adata.var_names, columns=programs)
    loadings = loadings.assign(hvg=is_hvg).reset_index(names="gene")
    loadings = loadings.melt(
        id_vars=["gene", "hvg"], var_name="program", value_name="score"
    )
    loadings = loadings.dropna(subset="score")
    loadings["program"] = pd.Categorical(loadings["program"], categories=programs)
    loadings = loadings[["program", "gene", "score", "hvg"]].sort_values(
        ["program", "score", "gene"], ascending=[True, False, True]
    )
    loadings = loadings.reset_index(drop=True)
    top = loadings.groupby("program", observed=True).head(p["n_top_genes"]).copy()
    top.insert(1, "rank", top.groupby("program", observed=True).cumcount() + 1)
    top = top.drop(columns="hvg")

    by_group = summarise_groups(usage_long, col, programs, groups)

    members["program"] = pd.Categorical(
        members["cluster"].map(cluster_name), categories=programs
    )
    stability = (
        members.groupby("program", observed=False)
        .agg(
            n_components=("replicate", "size"),
            n_replicates=("replicate", "nunique"),
            silhouette=("silhouette", "mean"),
        )
        .reset_index()
    )

    outdir = make_outdir(p["outdir"])
    tables, figures = outdir / "tables", outdir / "figures"
    usage_long.to_csv(tables / "usage.tsv", sep="\t", index=False)
    loadings.to_csv(tables / "gene_loadings.tsv", sep="\t", index=False)
    top.to_csv(tables / "top_genes.tsv", sep="\t", index=False)
    by_group.to_csv(tables / "program_by_group.tsv", sep="\t", index=False)
    stability.to_csv(tables / "stability.tsv", sep="\t", index=False)

    plt.switch_backend("Agg")
    draw_figures(figures, usage_long, by_group, top, loadings, col, programs, groups, p)

    (outdir / "data").mkdir()
    adata.obsm["usage"] = pd.DataFrame(usage, index=adata.obs_names, columns=programs)
    adata.varm["gene_scores"] = pd.DataFrame(
        scores.T, index=adata.var_names, columns=programs
    )
    adata.var["highly_variable"] = is_hvg
    adata.write_h5ad(outdir / "data" / "programs.h5ad")

    n_kept, n_total = int(kept.sum()), len(kept)
    mean_silhouette = float(members["silhouette"].mean())
    top_genes = {
        g: d["gene"].tolist() for g, d in top.groupby("program", observed=True)
    }
    best = by_group.loc[by_group.groupby("program", observed=True)["auc"].idxmax()]
    how = (
        f"NMF k = {p['k']} on {len(hvg)} HVGs, consensus of {p['n_replicates']} "
        f"replicates, seed {p['seed']}"
    )
    summary_form = (
        "Points with the median" if p["k"] < FEW_POINTS else "Box plot with points"
    )
    captions = {
        "usage_heatmap": f"Heatmap of mean usage of {p['k']} consensus-NMF "
        f"programs in each {col} level, centred and scaled per program across "
        f"{len(groups)} groups ({how}); rows fixed by program name (GP1 has the "
        f"largest total usage), columns fixed by {col} order, not clustered; "
        f"colours saturate at ±{HEATMAP_LIMIT:g}.",
        "usage_structure": f"Structure plot of program usage in {adata.n_obs} cells "
        f"({how}), each a column summing to 1, ordered by {col}, then dominant "
        "program, then its usage.",
        "top_gene_heatmap": f"Heatmap of gene scores (cNMF z-score units) for the top "
        f"{p['n_top_genes']} genes per program ({how}), scaled per gene across "
        "programs; rows fixed by program name, columns by program then score rank, "
        f"not clustered; colours saturate at ±{HEATMAP_LIMIT:g}.",
        "auc_by_group": f"{summary_form} of one-vs-rest AUC of program usage per "
        f"{col} group ({how}), one point per program (n = {p['k']}); dashed line "
        "at 0.5, no separation.",
    }
    summary = [
        f"{adata.n_obs} cells, {len(hvg)} HVGs; {p['k']} programs from "
        f"{p['n_replicates']} NMF replicates.",
        f"{loadings['gene'].nunique()} of {adata.n_vars} genes scored; the other "
        f"{adata.n_vars - loadings['gene'].nunique()} have no counts and are not "
        "tested.",
        f"{n_kept} of {n_total} replicate components pass density < "
        f"{p['density_threshold']}; mean silhouette {mean_silhouette:.3f}.",
        "Top group per program by AUC: "
        + ", ".join(
            f"{gp} {g} ({a:.2f})"
            for gp, g, a in zip(best["program"], best[col], best["auc"], strict=True)
        )
        + ".",
    ]
    write_report(
        outdir,
        f"Gene programs by consensus NMF (k = {p['k']})",
        p,
        KEYS,
        summary,
        captions,
    )
    write_result(
        outdir,
        {
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars,
            "n_hvg": len(hvg),
            "n_genes_scored": int(loadings["gene"].nunique()),
            "k": p["k"],
            "n_components_kept": n_kept,
            "n_components_total": n_total,
            "mean_silhouette": round(mean_silhouette, 4),
            "top_genes": top_genes,
            "seed": p["seed"],
            "params": p,
        },
        summary,
        ["report.md", "tables/program_by_group.tsv", "figures/usage_heatmap.png"],
        reference=p["reference"],
    )
    write_reproducibility(
        outdir,
        __file__,
        p,
        KEYS,
        full_argv,
        ["anndata", "scanpy", "numpy", "pandas", "scipy", "scikit-learn", "matplotlib"],
    )
    return outdir


if __name__ == "__main__":
    raise SystemExit(run("gene-programs", analyse))
