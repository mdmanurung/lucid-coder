"""Purpose: ligand-receptor cell-cell communication with LIANA+ rank_aggregate.
Inputs:  --adata <file.h5ad> with log-normalised X and an obs cell-type column
         (e.g. sc-annotate's output); --resource TSV with ligand, receptor columns
Outputs: <outdir>/report.md, result.json, tables/, figures/, reproducibility/
Run:     python scripts/ccc.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""

import argparse
import itertools
from pathlib import Path

import anndata as ad
import liana as li
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from liana.method import AggregateClass, aggregate_meta
from matplotlib.transforms import offset_copy
from plots import grey_to_red, plot_dotplot
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

REQUIRED = [
    "adata",
    "resource",
    "groupby",
    "methods",
    "expr_prop",
    "min_cells",
    "n_perms",
    "magnitude_max",
    "specificity_max",
    "sources",
    "targets",
    "n_top",
    "reference",
    "seed",
]
KEYS = [*REQUIRED[:3], "sample", *REQUIRED[3:]]
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "adata": str(SKILL_DIR / "examples" / "demo.h5ad"),
    "resource": str(SKILL_DIR / "examples" / "demo_resource.tsv"),
    "groupby": "cell_type",
    "sample": None,
    "methods": "rank_aggregate",
    "expr_prop": 0.1,
    "min_cells": 5,
    "n_perms": 100,
    "magnitude_max": 0.002,
    "specificity_max": 0.01,
    "sources": "all",
    "targets": "all",
    "n_top": 10,
    "reference": "synthetic demo genes, no genome build",
    "seed": 0,
    "outdir": "ccc_demo",
}
# The five methods liana's rank_aggregate combines; `rank_aggregate` means all.
METHODS = ["cellphonedb", "connectome", "logfc", "natmi", "singlecellsignalr"]

ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "bad_method": "--methods: unknown {names}; use rank_aggregate or a comma list "
    "of " + ", ".join(METHODS),
    "no_score": "--methods {methods} has no {kind} score; add a method that has one",
    "missing_file": "input file not found: {path}",
    "missing_column": "{what} column not found: {col}",
    "not_lognorm": "X does not look log-normalised: {why}",
    "few_groups": "--groupby {col} needs at least 2 groups, found {n}",
    "unknown_label": "{arg}: labels not in {col}: {labels}",
    "no_overlap": "no resource pair has its ligand and receptor in var_names; "
    "check species and gene identifiers",
    "no_pass": "zero interactions pass {rule}; nothing to plot",
}

COLUMNS = ["source", "target", "ligand", "receptor", "magnitude_rank"]
COLUMNS += ["specificity_rank", "passes"]
# Fixed plot spec, not a grill question: the dot plot widens 6 mm per pair, so a
# 4 mm dot at size 1 fills its column without touching the next one.
DOT_MAX_SIZE = 4
MAGNITUDE_LABEL = "-log10 magnitude rank"
SPECIFICITY_LABEL = "1 - specificity rank"
COUNT_LABEL = "passing interactions (n)"


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--adata")
    parser.add_argument("--resource")
    parser.add_argument("--groupby")
    parser.add_argument("--sample")
    parser.add_argument("--methods")
    parser.add_argument("--expr-prop", type=float)
    parser.add_argument("--min-cells", type=int)
    parser.add_argument("--n-perms", type=int)
    parser.add_argument("--magnitude-max", type=float)
    parser.add_argument("--specificity-max", type=float)
    parser.add_argument("--sources")
    parser.add_argument("--targets")
    parser.add_argument("--n-top", type=int)
    parser.add_argument("--reference")
    parser.add_argument("--seed", type=int)
    p = parse_required(parser, REQUIRED, DEMO, argv, optional=["sample"])
    rules = [
        ("--expr-prop", p["expr_prop"], 0 <= p["expr_prop"] <= 1, "in [0, 1]"),
        ("--min-cells", p["min_cells"], p["min_cells"] >= 1, "an integer >= 1"),
        ("--n-perms", p["n_perms"], p["n_perms"] >= 1, "an integer >= 1"),
        (
            "--magnitude-max",
            p["magnitude_max"],
            0 < p["magnitude_max"] <= 1,
            "in (0, 1]",
        ),
        (
            "--specificity-max",
            p["specificity_max"],
            0 < p["specificity_max"] <= 1,
            "in (0, 1]",
        ),
        ("--n-top", p["n_top"], p["n_top"] >= 1, "an integer >= 1"),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ]
    for arg, value, ok, rule in rules:
        if not ok:
            fail(ERRORS, "bad_value", arg=arg, rule=rule, value=value)
    return p


def pick_methods(spec):
    """--methods to liana method objects; each score kind needs one method."""
    names = METHODS if spec == "rank_aggregate" else spec.split(",")
    unknown = [name for name in names if name not in METHODS]
    if unknown:
        fail(ERRORS, "bad_method", names=", ".join(unknown))
    methods = [getattr(li.mt, name) for name in names]
    for kind in ("magnitude", "specificity"):
        if all(getattr(method, kind) is None for method in methods):
            fail(ERRORS, "no_score", methods=spec, kind=kind)
    return methods


def labels_of(spec, arg, groups, col):
    if spec == "all":
        return list(groups)
    labels = spec.split(",")
    unknown = [label for label in labels if label not in groups]
    if unknown:
        fail(ERRORS, "unknown_label", arg=arg, col=col, labels=", ".join(unknown))
    return labels


def load(p):
    """Step 1: read the h5ad and resource; check log-normalised X and labels."""
    for key in ("adata", "resource"):
        if not Path(p[key]).is_file():
            fail(ERRORS, "missing_file", path=p[key])
    adata = ad.read_h5ad(p["adata"])
    check_unique(adata.var_names, "gene names")
    check_unique(adata.obs_names, "cell ids")
    resource = pd.read_csv(p["resource"], sep="\t")
    for col in ("ligand", "receptor"):
        if col not in resource.columns:
            fail(ERRORS, "missing_column", what="resource", col=col)
    for col in (p["groupby"], p["sample"]):
        if col is not None and col not in adata.obs.columns:
            fail(ERRORS, "missing_column", what="obs", col=col)

    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    if (values < 0).any():
        fail(ERRORS, "not_lognorm", why="negative values (scaled?)")
    if np.array_equal(values, np.round(values)):
        fail(ERRORS, "not_lognorm", why="all integers (raw counts?)")

    # Complex subunits are joined by "_"; a pair counts if every subunit is a gene.
    genes = set(adata.var_names)
    present = resource[["ligand", "receptor"]].map(
        lambda x: set(str(x).split("_")) <= genes
    )
    if not present.all(axis=1).any():
        fail(ERRORS, "no_overlap")
    # Resource overlap: distinct ligand and receptor genes (subunits) in var_names.
    listed = {
        g
        for x in resource[["ligand", "receptor"]].to_numpy().ravel()
        for g in str(x).split("_")
    }
    overlap = {"n_matched": len(listed & genes), "n_in_resource": len(listed)}

    # Categorical order is the group order of every table and figure.
    col = p["groupby"]
    groups = adata.obs[col]
    if not isinstance(groups.dtype, pd.CategoricalDtype):
        groups = pd.Categorical(groups.astype(str), categories=sorted(set(groups)))
    adata.obs[col] = pd.Categorical(groups).remove_unused_categories()
    categories = list(adata.obs[col].cat.categories)
    if len(categories) < 2:
        fail(ERRORS, "few_groups", col=col, n=len(categories))
    sources = labels_of(p["sources"], "--sources", categories, col)
    targets = labels_of(p["targets"], "--targets", categories, col)
    return adata, resource, overlap, categories, sources, targets


def infer(adata, resource, methods, sources, targets, categories, p):
    """Step 2: liana rank_aggregate over the chosen methods, one thread, seeded."""
    aggregate = AggregateClass(aggregate_meta, methods=methods)
    kwargs = {
        "groupby": p["groupby"],
        "resource": resource,
        "expr_prop": p["expr_prop"],
        "min_cells": p["min_cells"],
        "n_perms": p["n_perms"],
        "seed": p["seed"],
        "n_jobs": 1,
        # X holds the log-normalised values; .raw is never read.
        "use_raw": False,
        "inplace": False,
    }
    # Restricting pairs changes the ranked universe; only pass it when asked.
    if sources != categories or targets != categories:
        kwargs["groupby_pairs"] = pd.DataFrame(
            itertools.product(sources, targets), columns=["source", "target"]
        )
    if p["sample"] is None:
        return aggregate(adata, **kwargs)
    return aggregate.by_sample(adata, sample_key=p["sample"], **kwargs)


def tidy(res, p):
    """Step 3: rename, call passes, sort by sample, ranks, then labels."""
    res = res.rename(
        columns={"ligand_complex": "ligand", "receptor_complex": "receptor"}
    )
    res["passes"] = (res["magnitude_rank"] <= p["magnitude_max"]) & (
        res["specificity_rank"] <= p["specificity_max"]
    )
    lead = ["sample"] if p["sample"] else []
    order = [*lead, "magnitude_rank", "specificity_rank", *COLUMNS[:4]]
    return res[lead + COLUMNS].sort_values(order).reset_index(drop=True)


def count_pairs(interactions, categories, p):
    """Step 4: passing interactions per (sample,) source x target, zeros kept."""
    by = (["sample"] if p["sample"] else []) + ["source", "target"]
    levels = [categories, categories]
    if p["sample"]:
        levels.insert(0, sorted(interactions["sample"].unique()))
    grid = pd.MultiIndex.from_product(levels, names=by)
    counts = interactions[interactions["passes"]].groupby(by).size()
    return counts.reindex(grid, fill_value=0).rename("n_pass").reset_index()


def draw_counts(figures, pair_counts, categories):
    """Sequential heatmap: counts are unsigned, so not the diverging plot_heatmap."""
    data = pair_counts.groupby(["source", "target"], as_index=False)["n_pass"].sum()
    matrix = data.pivot(index="source", columns="target", values="n_pass")
    matrix = matrix.reindex(index=categories, columns=categories)
    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(matrix, cmap=grey_to_red, vmin=0)
    for (i, j), n in np.ndenumerate(matrix.to_numpy()):
        ax.text(j + 0.5, i + 0.5, n, ha="center", va="center")
    ticks = np.arange(len(categories)) + 0.5
    ax.set_xticks(ticks, categories, rotation=90)
    ax.set_yticks(ticks, categories)
    ax.invert_yaxis()
    ax.set_xlabel("target")
    ax.set_ylabel("source")
    fig.colorbar(mesh, ax=ax, label=COUNT_LABEL)
    fig.set_layout_engine("tight")
    side = max(70, 30 + 8 * len(categories))
    save_panel(fig, data, figures, "pair_counts_heatmap", side + 20, side)
    plt.close(fig)


def draw_dotplot(figures, top, p):
    """Top passing interactions: colour is magnitude, dot area specificity."""
    prefix = top["sample"] + ": " if p["sample"] else ""
    pair = prefix + top["source"] + " → " + top["target"]
    lr = top["ligand"] + " → " + top["receptor"]
    data = top.assign(
        pair=pd.Categorical(pair, categories=list(dict.fromkeys(pair))),
        interaction=pd.Categorical(lr, categories=list(dict.fromkeys(lr))),
        **{
            MAGNITUDE_LABEL: -np.log10(top["magnitude_rank"]),
            SPECIFICITY_LABEL: 1 - top["specificity_rank"],
        },
    )
    fig = plot_dotplot(
        data,
        "pair",
        "interaction",
        MAGNITUDE_LABEL,
        SPECIFICITY_LABEL,
        max_size=DOT_MAX_SIZE,
    )
    ax = fig.axes[0]
    ax.invert_yaxis()
    # plot_dotplot labels sizes as percents; 1 - rank is a plain fraction.
    for text, value in zip(
        ax.get_legend().get_texts(), ("0.25", "0.5", "1"), strict=True
    ):
        text.set_text(value)
    # 36 pt right of the colorbar clears its ticks and label at 8 pt type.
    cbar_right = offset_copy(fig.axes[1].transAxes, fig=fig, x=36, units="points")
    ax.get_legend().set_bbox_to_anchor((1, 1), transform=cbar_right)
    fig.set_layout_engine("tight")
    width = min(183, max(110, 70 + 6 * data["pair"].nunique()))
    height = max(60, 25 + 6 * data["interaction"].nunique())
    save_panel(fig, data, figures, "lr_dotplot", width, height)
    plt.close(fig)


def analyse(argv, full_argv):
    p = parse_args(argv)
    methods = pick_methods(p["methods"])
    seed_everything(p["seed"])
    adata, resource, overlap, categories, sources, targets = load(p)
    # One BLAS thread and n_jobs=1 keep the permutations single-threaded.
    with threadpool_limits(1):
        res = infer(adata, resource, methods, sources, targets, categories, p)
    interactions = tidy(res, p)
    rule = (
        f"magnitude_rank <= {p['magnitude_max']}, "
        f"specificity_rank <= {p['specificity_max']}"
    )
    passing = interactions[interactions["passes"]]
    if passing.empty:
        fail(ERRORS, "no_pass", rule=rule)
    top = passing.head(p["n_top"])
    pair_counts = count_pairs(interactions, categories, p)

    outdir = make_outdir(p["outdir"])
    tables, figures = outdir / "tables", outdir / "figures"
    interactions.to_csv(tables / "interactions.tsv", sep="\t", index=False)
    pair_counts.to_csv(tables / "pair_counts.tsv", sep="\t", index=False)

    plt.switch_backend("Agg")
    draw_counts(figures, pair_counts, categories)
    draw_dotplot(figures, top, p)

    col = p["groupby"]
    per = f" per sample ({p['sample']}), summed" if p["sample"] else ""
    captions = {
        "pair_counts_heatmap": f"Heatmap of the number of interactions passing {rule} "
        f"for each source (rows) and target (columns) {col}{per}; not clustered: "
        f"rows and columns in {col} order.",
        "lr_dotplot": f"Dot plot of the {len(top)} top passing ligand-receptor "
        "interactions by magnitude rank; colour is -log10 magnitude rank, dot area "
        "is 1 - specificity rank.",
    }
    method_names = ", ".join(m.method_name for m in methods)
    summary = [
        f"{adata.n_obs} cells in {len(categories)} {col} groups; "
        f"{len(resource)} resource pairs; methods {method_names}.",
        f"{len(passing)} of {len(interactions)} interactions pass {rule}.",
        "Top: "
        + "; ".join(
            f"{r.source} → {r.target} {r.ligand}-{r.receptor}"
            for r in top.head(3).itertuples()
        ),
    ]
    write_report(
        outdir, f"Cell-cell communication by {col}", p, KEYS, summary, captions
    )
    write_result(
        outdir,
        {
            "n_cells": adata.n_obs,
            "n_groups": len(categories),
            "n_resource_pairs": len(resource),
            **overlap,
            "n_interactions": len(interactions),
            "n_pass": len(passing),
            "top_interactions": top[COLUMNS[:4]].to_dict("records"),
            "methods": [m.method_name for m in methods],
            "seed": p["seed"],
            "params": p,
        },
        summary,
        ["report.md", "tables/interactions.tsv", "figures/lr_dotplot.png"],
        reference=p["reference"],
    )
    write_reproducibility(
        outdir,
        __file__,
        p,
        KEYS,
        full_argv,
        ["liana", "anndata", "numpy", "pandas", "scipy", "matplotlib"],
    )
    return outdir


if __name__ == "__main__":
    raise SystemExit(run("ccc", analyse))
