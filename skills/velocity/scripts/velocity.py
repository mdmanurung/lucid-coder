"""Purpose: RNA velocity (scVelo) and fate probabilities (CellRank 2) from an h5ad.
Inputs:  --adata <file.h5ad> with spliced/unspliced layers, X_umap and a group column
Outputs: <outdir>/report.md, result.json, data/, tables/, figures/, reproducibility/
Run:     python scripts/velocity.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""

import os

# One thread for BLAS and numba before numpy loads: same numbers on any host.
for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(var, "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import argparse  # noqa: E402
import warnings  # noqa: E402
from pathlib import Path  # noqa: E402

import anndata as ad  # noqa: E402
import cellrank as cr  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import scvelo as scv  # noqa: E402
from matplotlib.colors import to_hex  # noqa: E402
from plots import plot_composition  # noqa: E402
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
    "mode",
    "min_shared_counts",
    "n_top_genes",
    "n_pcs",
    "n_neighbors",
    "groupby",
    "terminal_states",
    "reference",
    "seed",
]
OPTIONAL = ["n_states"]
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "adata": str(SKILL_DIR / "examples" / "demo.h5ad"),
    "mode": "stochastic",
    "min_shared_counts": 10,
    "n_top_genes": 30,
    "n_pcs": 10,
    "n_neighbors": 30,
    "groupby": "stage",
    "terminal_states": "s4",
    "n_states": None,
    "reference": "synthetic demo, no genome build",
    "seed": 0,
    "outdir": "velocity_demo",
}
ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "auto_needs_n_states": "--terminal-states auto needs --n-states",
    "n_states_unused": "--n-states is used only with --terminal-states auto",
    "missing_file": "input file not found: {path}",
    "missing_layer": "layer not found: {layer}; velocity needs spliced and "
    "unspliced counts (velocyto, STARsolo or alevin-fry)",
    "no_umap": "obsm['X_umap'] not found; run sc-embed or sc.tl.umap first",
    "missing_column": "obs column not found: {col}",
    "one_group": "obs column {col!r} has {n} group(s); fates need at least 2",
    "terminal_absent": "terminal states not in obs column {col!r}: {groups}",
    "too_few_genes": "--n-top-genes {n_top} exceeds the {n_kept} genes left after "
    "--min-shared-counts {min_shared}",
    "too_many_pcs": "--n-pcs {n_pcs} must be below min(cells, genes) = {limit}",
    "macrostates_failed": "GPCCA could not form {n} macrostates ({reason}); "
    "ask the user for another --n-states",
    "no_terminal": "no macrostate of {n} has stability >= {t}; "
    "ask the user for another --n-states or name the terminal groups",
}
PACKAGES = [
    "scvelo",
    "cellrank",
    "anndata",
    "scanpy",
    "numpy",
    "pandas",
    "scipy",
    "matplotlib",
]
# Fixed method choices, not user thresholds; SKILL.md "Steps" names each one.
# Auto terminal states: CellRank's tutorial rule, macrostates with
# coarse-grained self-transition >= 0.96 (predict_terminal_states default).
STABILITY_THRESHOLD = 0.96
# Gene rank for the phase portraits: scVelo's own fit score per mode.
SCORE = {"stochastic": "velocity_r2", "dynamical": "fit_likelihood"}
N_PHASE_GENES = 4
PSEUDOTIME = "velocity pseudotime"


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--adata")
    parser.add_argument("--mode")
    parser.add_argument("--min-shared-counts", type=int)
    parser.add_argument("--n-top-genes", type=int)
    parser.add_argument("--n-pcs", type=int)
    parser.add_argument("--n-neighbors", type=int)
    parser.add_argument("--groupby")
    parser.add_argument("--terminal-states")
    parser.add_argument("--n-states", type=int)
    parser.add_argument("--reference")
    parser.add_argument("--seed", type=int)
    return parse_required(parser, REQUIRED, DEMO, argv, optional=OPTIONAL)


def check_values(p):
    """Return the terminal groups as a list, or None for `auto`."""
    rules = [
        ("--mode", p["mode"], p["mode"] in SCORE, "stochastic or dynamical"),
        (
            "--min-shared-counts",
            p["min_shared_counts"],
            p["min_shared_counts"] >= 0,
            "an integer >= 0",
        ),
        ("--n-top-genes", p["n_top_genes"], p["n_top_genes"] >= 2, "an integer >= 2"),
        ("--n-pcs", p["n_pcs"], p["n_pcs"] >= 2, "an integer >= 2"),
        ("--n-neighbors", p["n_neighbors"], p["n_neighbors"] >= 2, "an integer >= 2"),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ]
    for arg, value, ok, rule in rules:
        if not ok:
            fail(ERRORS, "bad_value", arg=arg, rule=rule, value=value)
    auto = p["terminal_states"] == "auto"
    if auto and p["n_states"] is None:
        fail(ERRORS, "auto_needs_n_states")
    if not auto and p["n_states"] is not None:
        fail(ERRORS, "n_states_unused")
    if auto:
        if p["n_states"] < 2:
            fail(
                ERRORS,
                "bad_value",
                arg="--n-states",
                rule="an integer >= 2",
                value=p["n_states"],
            )
        return None
    terminal = [g.strip() for g in p["terminal_states"].split(",") if g.strip()]
    if not terminal:
        fail(
            ERRORS,
            "bad_value",
            arg="--terminal-states",
            rule="a comma list of groups or auto",
            value=repr(p["terminal_states"]),
        )
    return terminal


def load(p, terminal):
    """Step 1: read the h5ad and check every field the later steps use."""
    path = Path(p["adata"])
    if not path.is_file():
        fail(ERRORS, "missing_file", path=path)
    adata = ad.read_h5ad(path)
    check_unique(adata.var_names, "gene names")
    check_unique(adata.obs_names, "cell ids")
    for layer in ("spliced", "unspliced"):
        if layer not in adata.layers:
            fail(ERRORS, "missing_layer", layer=layer)
    if "X_umap" not in adata.obsm:
        fail(ERRORS, "no_umap")
    col = p["groupby"]
    if col not in adata.obs.columns:
        fail(ERRORS, "missing_column", col=col)
    # Step 2: groups become categorical; an existing category order is kept.
    if not isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
        adata.obs[col] = adata.obs[col].astype(str).astype("category")
    adata.obs[col] = adata.obs[col].cat.remove_unused_categories()
    levels = list(adata.obs[col].cat.categories)
    if len(levels) < 2:
        fail(ERRORS, "one_group", col=col, n=len(levels))
    absent = [g for g in terminal or [] if g not in levels]
    if absent:
        fail(ERRORS, "terminal_absent", col=col, groups=", ".join(absent))
    return adata, levels


def preprocess(adata, p):
    """Steps 3-4: gene filter, normalise, HVG subset, PCA, neighbours, moments.

    Returns the input gene count, so the report can say how many were not tested.
    """
    n_genes_input = adata.n_vars
    scv.pp.filter_genes(adata, min_shared_counts=p["min_shared_counts"])
    if p["n_top_genes"] > adata.n_vars:
        fail(
            ERRORS,
            "too_few_genes",
            n_top=p["n_top_genes"],
            n_kept=adata.n_vars,
            min_shared=p["min_shared_counts"],
        )
    limit = min(adata.n_obs, p["n_top_genes"])
    if p["n_pcs"] >= limit:
        fail(ERRORS, "too_many_pcs", n_pcs=p["n_pcs"], limit=limit)
    scv.pp.normalize_per_cell(adata)
    # sc-embed output is already log1p; logging twice would flatten X.
    if "log1p" not in adata.uns:
        sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=p["n_top_genes"], subset=True)
    sc.pp.pca(adata, n_comps=p["n_pcs"], random_state=p["seed"])
    sc.pp.neighbors(
        adata, n_pcs=p["n_pcs"], n_neighbors=p["n_neighbors"], random_state=p["seed"]
    )
    scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
    return n_genes_input


def fates(adata, terminal, p):
    """Steps 7-8: VelocityKernel, GPCCA terminal states, fate probabilities."""
    vk = cr.kernels.VelocityKernel(adata).compute_transition_matrix(
        model="deterministic", n_jobs=1, show_progress_bar=False
    )
    g = cr.estimators.GPCCA(vk)
    if terminal is None:
        try:
            g.compute_macrostates(n_states=p["n_states"], cluster_key=p["groupby"])
        except ValueError as exc:
            fail(ERRORS, "macrostates_failed", n=p["n_states"], reason=exc)
        try:
            g.predict_terminal_states(
                method="stability", stability_threshold=STABILITY_THRESHOLD
            )
        except ValueError:
            fail(ERRORS, "no_terminal", n=p["n_states"], t=STABILITY_THRESHOLD)
    else:
        # Every cell of a named group is absorbing (the dict form of the API).
        col = adata.obs[p["groupby"]]
        g.set_terminal_states(
            states={s: adata.obs_names[col == s].tolist() for s in terminal},
            cluster_key=p["groupby"],
        )
    g.compute_fate_probabilities(
        solver="gmres", use_petsc=False, n_jobs=1, show_progress_bar=False
    )
    lineages = g.fate_probabilities
    probs = pd.DataFrame(
        np.asarray(lineages), index=adata.obs_names, columns=list(lineages.names)
    )
    # Lineage objects do not round-trip through h5ad; store a plain array.
    adata.obsm["lineages_fwd"] = probs.to_numpy()
    adata.uns["lineages_fwd_names"] = list(probs.columns)
    adata.uns.pop("eigendecomposition_fwd", None)
    adata.obs["term_states_fwd"] = adata.obs["term_states_fwd"].astype(str)
    return probs


def frame_umap(ax):
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")


def phase_slope(var, mode):
    """Steady-state ratio Mu / Ms of each gene's fit, drawn on its phase portrait."""
    if mode == "stochastic":
        return var["velocity_gamma"]
    # scVelo's dynamical model fits u / scaling; steady state is u = gamma/beta s.
    return var["fit_scaling"] * var["fit_gamma"] / var["fit_beta"]


def draw_figures(figures, adata, p, levels, pal, cells, by_group, top):
    col = p["groupby"]
    umap = pd.DataFrame(adata.obsm["X_umap"][:, :2], columns=["UMAP 1", "UMAP 2"])
    umap[col] = adata.obs[col].to_numpy()
    umap[["velocity UMAP 1", "velocity UMAP 2"]] = adata.obsm["velocity_umap"][:, :2]

    adata.uns[f"{col}_colors"] = [pal[lv] for lv in levels]
    fig, ax = plt.subplots(layout="constrained")
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color=col,
        ax=ax,
        show=False,
        title="",
        legend_loc="right margin",
        size=8,
        alpha=0.6,
        linewidth=0.6,
        arrow_size=0.6,
        frameon=True,
    )
    frame_umap(ax)
    save_panel(fig, umap, figures, "velocity_stream", 120, 89)
    plt.close(fig)

    umap[PSEUDOTIME] = cells.loc[adata.obs_names, PSEUDOTIME].to_numpy()
    fig, ax = plt.subplots(layout="constrained")
    dots = ax.scatter(
        umap["UMAP 1"],
        umap["UMAP 2"],
        s=2,
        c=umap[PSEUDOTIME],
        cmap="viridis",
        rasterized=True,
    )
    fig.colorbar(dots, ax=ax, label=PSEUDOTIME)
    frame_umap(ax)
    save_panel(fig, umap, figures, "umap_velocity_pseudotime", 89, 75)
    plt.close(fig)

    fate_pal = dict(
        zip(
            by_group["terminal_state"].cat.categories,
            [to_hex(c) for c in plt.get_cmap("Set2").colors],
            strict=False,
        )
    )
    fig = plot_composition(
        by_group, col, "terminal_state", "fate probability", fate_pal
    )
    fig.axes[0].tick_params(axis="x", labelrotation=45)
    fig.set_layout_engine("constrained")
    save_panel(fig, by_group, figures, "fate_by_group", 89, 70)
    plt.close(fig)

    slope = phase_slope(adata.var, p["mode"])
    phase = pd.concat(
        [
            pd.DataFrame(
                {
                    "gene": gene,
                    "cell": adata.obs_names,
                    col: adata.obs[col].to_numpy(),
                    "Ms": adata[:, gene].layers["Ms"].ravel(),
                    "Mu": adata[:, gene].layers["Mu"].ravel(),
                }
            )
            for gene in top
        ]
    )
    fig, axes = plt.subplots(1, len(top), layout="constrained", squeeze=False)
    for ax, gene in zip(axes.flat, top, strict=True):
        g = phase[phase["gene"] == gene]
        ax.scatter(g["Ms"], g["Mu"], s=2, c=[pal[lv] for lv in g[col]], rasterized=True)
        xs = np.array([0, g["Ms"].max()])
        ax.plot(xs, slope[gene] * xs, ls="--", lw=0.6, c="black")
        ax.annotate(gene, (0.02, 0.9), xycoords="axes fraction")
        ax.set_xlabel("spliced (Ms)")
    axes.flat[0].set_ylabel("unspliced (Mu)")
    save_panel(fig, phase, figures, "phase_portraits", 45 * len(top), 45)
    plt.close(fig)


def analyse(argv, full_argv):
    p = parse_args(argv)
    terminal = check_values(p)
    adata, levels = load(p, terminal)
    col = p["groupby"]
    seed_everything(p["seed"])
    sc.settings.n_jobs = 1
    scv.settings.verbosity = 0
    cr.settings.logging_level = "ERROR"
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    # Preprocess before make_outdir: its gene and PC checks leave no outdir.
    n_genes_input = preprocess(adata, p)
    outdir = make_outdir(p["outdir"])
    (outdir / "data").mkdir()
    # Steps 5-6: velocity, velocity graph, velocity pseudotime, UMAP projection.
    if p["mode"] == "dynamical":
        scv.tl.recover_dynamics(adata, n_jobs=1, show_progress_bar=False)
    scv.tl.velocity(adata, mode=p["mode"])
    scv.tl.velocity_graph(adata, n_jobs=1, show_progress_bar=False)
    scv.tl.velocity_pseudotime(adata)
    scv.tl.velocity_embedding(adata, basis="umap")
    probs = fates(adata, terminal, p)
    states = list(probs.columns)

    groups = adata.obs[col].astype(str).to_numpy()
    cells = pd.DataFrame(
        {
            "cell": adata.obs_names,
            col: groups,
            PSEUDOTIME: adata.obs["velocity_pseudotime"].to_numpy(),
        }
    ).sort_values("cell")
    cells.rename(columns={PSEUDOTIME: "velocity_pseudotime"}).to_csv(
        outdir / "tables" / "velocity_pseudotime.tsv", sep="\t", index=False
    )
    cells = cells.set_index("cell", drop=False)

    fate_table = probs.assign(cell=probs.index, **{col: groups})[
        ["cell", col, *states]
    ].sort_values("cell")
    fate_table.to_csv(
        outdir / "tables" / "fate_probabilities.tsv", sep="\t", index=False
    )
    by_group = (
        probs.assign(**{col: groups})
        .groupby(col)[states]
        .mean()
        .reindex(levels)
        .rename_axis(index=col, columns="terminal_state")
        .stack()
        .rename("fate probability")
        .reset_index()
    )
    by_group[col] = pd.Categorical(by_group[col], categories=levels)
    by_group["terminal_state"] = pd.Categorical(
        by_group["terminal_state"], categories=states
    )
    n_cells = adata.obs[col].value_counts()
    by_group["n_cells"] = n_cells.loc[by_group[col]].to_numpy()
    by_group = by_group.sort_values([col, "terminal_state"]).reset_index(drop=True)
    by_group.rename(columns={"fate probability": "mean_fate_probability"}).to_csv(
        outdir / "tables" / "fate_by_group.tsv", sep="\t", index=False
    )

    score = SCORE[p["mode"]]
    genes = (
        adata.var.filter(regex="^(velocity|fit)_")
        .rename_axis("gene")
        .reset_index()
        .sort_values("gene")
    )
    genes.to_csv(outdir / "tables" / "velocity_genes.tsv", sep="\t", index=False)
    top = (
        genes.sort_values([score, "gene"], ascending=[False, True])["gene"]
        .head(N_PHASE_GENES)
        .tolist()
    )
    adata.write_h5ad(outdir / "data" / "velocity.h5ad")

    # ponytail: tab20 repeats colours past 20 groups; pass a palette when needed.
    tab20 = plt.get_cmap("tab20").colors
    pal = {lv: to_hex(tab20[k % 20]) for k, lv in enumerate(levels)}
    plt.switch_backend("Agg")
    draw_figures(outdir / "figures", adata, p, levels, pal, cells, by_group, top)

    medians = cells.groupby(col)[PSEUDOTIME].median().loc[levels]
    mode_note = "steady-state model" if p["mode"] == "stochastic" else "dynamical model"
    graph = (
        f"velocity graph on {p['n_pcs']} PCs of {adata.n_vars} highly variable genes, "
        f"{p['n_neighbors']} neighbours, seed {p['seed']}; UMAP from the input "
        "X_umap, not recomputed"
    )
    how = (
        f"every cell of {', '.join(terminal)} absorbing"
        if terminal
        else f"GPCCA on {p['n_states']} macrostates, stability >= {STABILITY_THRESHOLD}"
    )
    captions = {
        "velocity_stream": f"UMAP of {adata.n_obs} cells coloured by {col}, with "
        f"scVelo velocity ({mode_note}) projected as streamlines ({graph}); arrows "
        "are a model projection, not measured movement.",
        "umap_velocity_pseudotime": f"UMAP of {adata.n_obs} cells coloured by "
        "scVelo velocity pseudotime (0 at the inferred root cells, 1 at the end "
        f"points; {graph}).",
        "fate_by_group": f"Mean CellRank fate probability per {col} group towards "
        f"each of {len(states)} terminal states ({how}); bars sum to 100%.",
        "phase_portraits": f"Phase portraits of the {len(top)} genes with the "
        f"highest {score}: smoothed unspliced (Mu) against spliced (Ms) per cell, "
        f"moments over {p['n_neighbors']} neighbours, coloured by group; the dashed "
        "line is the fitted steady-state ratio.",
    }
    summary = [
        f"{adata.n_obs} cells, {adata.n_vars} of {n_genes_input} genes tested "
        f"(the rest filtered, not tested), {p['mode']} velocity; top "
        f"genes by {score}: {', '.join(top)}.",
        "Median velocity pseudotime by group: "
        + ", ".join(f"{lv} {v:.3f}" for lv, v in medians.items())
        + ".",
        f"Terminal states: {', '.join(states)} ({how}).",
    ]
    keys = [*REQUIRED, *OPTIONAL]
    write_report(outdir, f"Velocity: {p['mode']} model", p, keys, summary, captions)
    result = {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "n_genes_input": n_genes_input,
        "mode": p["mode"],
        "terminal_states": states,
        "top_genes": top,
        "median_velocity_pseudotime": medians.round(6).to_dict(),
        "params": p,
        "seed": p["seed"],
    }
    artifacts = [
        "report.md",
        "tables/fate_by_group.tsv",
        "figures/velocity_stream.png",
    ]
    write_result(outdir, result, summary, artifacts, reference=p["reference"])
    write_reproducibility(outdir, __file__, p, keys, full_argv, PACKAGES)
    return outdir


if __name__ == "__main__":
    raise SystemExit(run("velocity", analyse))
