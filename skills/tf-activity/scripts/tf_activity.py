"""Purpose: TF or pathway activity from a prior network with decoupler (ULM or MLM).
Inputs:  --input cells:    --adata <file.h5ad>, log-normalised X, obs column --groupby
         --input contrast: --de-table <file.tsv> with gene and --stat-col columns
         --network <file.tsv> with source, target, weight (e.g. CollecTRI, PROGENy)
Outputs: <outdir>/report.md, result.json, tables/, figures/, reproducibility/
Run:     python scripts/tf_activity.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""
# ruff: noqa: E402

import os

# One thread for BLAS and numba: set before numpy is imported anywhere.
for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(var, "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import argparse
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plots import plot_heatmap
from scipy import sparse, stats
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

SKILL_DIR = Path(__file__).resolve().parents[1]

REQUIRED = ["input", "network", "method", "tmin", "n_top", "reference", "seed"]
MODE_ARGS = {
    "cells": ["adata", "groupby", "rank_method", "limit"],
    "contrast": ["de_table", "stat_col", "contrast"],
}
OPTIONAL = [a for args in MODE_ARGS.values() for a in args]
EXAMPLES = SKILL_DIR / "examples"
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "cells": {
        "input": "cells",
        "network": str(EXAMPLES / "demo_network.tsv"),
        "method": "ulm",
        "tmin": 5,
        "n_top": 2,
        "reference": "synthetic demo genes, no genome build",
        "seed": 0,
        "adata": str(EXAMPLES / "demo.h5ad"),
        "groupby": "group",
        "rank_method": "t-test_overestim_var",
        "limit": 3.0,
        "outdir": "tf_activity_demo",
    },
    "contrast": {
        "input": "contrast",
        "network": str(EXAMPLES / "demo_network.tsv"),
        "method": "ulm",
        "tmin": 5,
        "n_top": 2,
        "reference": "synthetic demo genes, no genome build",
        "seed": 0,
        "de_table": str(EXAMPLES / "demo_de.tsv"),
        "stat_col": "stat",
        "contrast": "g1_vs_g2",
        "outdir": "tf_activity_demo",
    },
}

ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "mode_missing": "--input {mode} needs: {args}; ask the user, or use --demo",
    "mode_unused": "--input {mode} does not use: {args}",
    "missing_file": "input file not found: {path}",
    "bad_network": "network {path}: {why}",
    "missing_column": "{what} column not found: {col}",
    "not_lognorm": "X does not look log-normalised: {why}",
    "few_groups": "--groupby {col} needs at least 2 groups, found {n}",
    "small_group": "groups with fewer than 2 cells in {col}: {groups}",
    "bad_de": "DE table {path}: {why}",
    "no_sources": "no source has at least --tmin {tmin} targets in the data",
}
# decoupler's ulm and mlm scores are t-values of the regression slope.
UNIT = "t-value"


def parse_args(argv):
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--input", choices=list(MODE_ARGS))
    mode = pre.parse_known_args(argv)[0].input or "cells"
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", choices=list(MODE_ARGS))
    parser.add_argument("--network")
    parser.add_argument("--method", choices=["ulm", "mlm"])
    parser.add_argument("--tmin", type=int)
    parser.add_argument("--n-top", type=int)
    parser.add_argument("--reference")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--adata")
    parser.add_argument("--groupby")
    parser.add_argument(
        "--rank-method", choices=["t-test_overestim_var", "t-test", "wilcoxon"]
    )
    parser.add_argument("--limit", type=float)
    parser.add_argument("--de-table")
    parser.add_argument("--stat-col")
    parser.add_argument("--contrast")
    p = parse_required(
        parser,
        REQUIRED,
        DEMO[mode],
        argv,
        optional=OPTIONAL,
        demo_extra=("outdir", "input"),
    )
    mode = p["input"]
    missing = [a for a in MODE_ARGS[mode] if p.get(a) is None]
    if missing:
        fail(ERRORS, "mode_missing", mode=mode, args=" ".join(flags(missing)))
    unused = [a for a in OPTIONAL if a not in MODE_ARGS[mode] and p.get(a) is not None]
    if unused:
        fail(ERRORS, "mode_unused", mode=mode, args=" ".join(flags(unused)))
    rules = [
        ("--tmin", p["tmin"], p["tmin"] >= 1, "an integer >= 1"),
        ("--n-top", p["n_top"], p["n_top"] >= 1, "an integer >= 1"),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ]
    if mode == "cells":
        rules.append(("--limit", p["limit"], p["limit"] > 0, "> 0"))
    for arg, value, ok, rule in rules:
        if not ok:
            fail(ERRORS, "bad_value", arg=arg, rule=rule, value=value)
    return p


def flags(keys):
    return [f"--{k.replace('_', '-')}" for k in keys]


def existing(path):
    path = Path(path)
    if not path.is_file():
        fail(ERRORS, "missing_file", path=path)
    return path


def load_network(p):
    """Step 1: read the network; decoupler asserts on these, so check them here."""
    path = existing(p["network"])
    net = pd.read_csv(path, sep="\t")
    missing = [c for c in ("source", "target", "weight") if c not in net.columns]
    if missing:
        fail(ERRORS, "bad_network", path=path, why="missing " + ", ".join(missing))
    net = net[["source", "target", "weight"]].astype({"source": str, "target": str})
    net["weight"] = pd.to_numeric(net["weight"], errors="coerce")
    if net["weight"].isna().any():
        fail(ERRORS, "bad_network", path=path, why="non-numeric or empty weight")
    if net.duplicated(["source", "target"]).any():
        fail(ERRORS, "bad_network", path=path, why="duplicate source-target rows")
    return net


def load_cells(p):
    """Step 2 (cells): log-normalised h5ad with at least 2 groups of 2 cells."""
    adata = ad.read_h5ad(existing(p["adata"]))
    check_unique(adata.var_names, "gene names")
    check_unique(adata.obs_names, "cell ids")
    col = p["groupby"]
    if col not in adata.obs.columns:
        fail(ERRORS, "missing_column", what="obs", col=col)
    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    if (values < 0).any():
        fail(ERRORS, "not_lognorm", why="negative values (scaled?)")
    if np.array_equal(values, np.round(values)):
        fail(ERRORS, "not_lognorm", why="all integers (raw counts?)")
    groups = adata.obs[col]
    if not isinstance(groups.dtype, pd.CategoricalDtype):
        groups = pd.Categorical(groups.astype(str), categories=sorted(set(groups)))
    adata.obs[col] = pd.Categorical(groups).remove_unused_categories()
    sizes = adata.obs[col].value_counts(sort=False)
    if len(sizes) < 2:
        fail(ERRORS, "few_groups", col=col, n=len(sizes))
    if (sizes < 2).any():
        fail(ERRORS, "small_group", col=col, groups=", ".join(sizes.index[sizes < 2]))
    return adata


def load_contrast(p):
    """Step 2 (contrast): one row of gene statistics; NA rows are dropped, counted."""
    path = existing(p["de_table"])
    de = pd.read_csv(path, sep="\t")
    for col in ("gene", p["stat_col"]):
        if col not in de.columns:
            fail(ERRORS, "missing_column", what="DE table", col=col)
    check_unique(de["gene"], "gene names in the DE table")
    stat = pd.to_numeric(de[p["stat_col"]], errors="coerce")
    keep = np.isfinite(stat)
    if keep.sum() < 3:
        fail(ERRORS, "bad_de", path=path, why="fewer than 3 genes with a finite stat")
    mat = pd.DataFrame(
        [stat[keep].to_numpy()], index=[p["contrast"]], columns=de["gene"][keep]
    )
    return mat, int((~keep).sum())


def prune(net, features, tmin):
    """Step 3: targets per source present in the data; drop sources under tmin."""
    shared = net[net["target"].isin(set(features))]
    counts = (
        shared.groupby("source").size().reindex(net["source"].unique(), fill_value=0)
    )
    n_targets = counts.rename("n_targets").rename_axis("source").reset_index()
    kept = n_targets[n_targets["n_targets"] >= tmin]
    if kept.empty:
        fail(ERRORS, "no_sources", tmin=tmin)
    dropped = n_targets[n_targets["n_targets"] < tmin].sort_values("source")
    return shared[shared["source"].isin(kept["source"])], n_targets, dropped


def score_cells(adata, net, p):
    """Step 4 (cells): per-cell scores, then each group against the rest."""
    import decoupler as dc  # deferred: a slow import, not needed to reject bad input

    method = getattr(dc.mt, p["method"])
    repaired = method(adata, net, tmin=p["tmin"])
    # decoupler returns a copy only when it dropped empty cells.
    adata = adata if repaired is None else repaired
    acts = dc.pp.get_obsm(adata, f"score_{p['method']}")
    ranked = dc.tl.rankby_group(
        acts, groupby=p["groupby"], reference="rest", method=p["rank_method"]
    )
    groups = acts.obs[p["groupby"]].cat.categories
    means = (
        pd.DataFrame(acts.X, columns=acts.var_names)
        .assign(group=acts.obs[p["groupby"]].to_numpy())
        .groupby("group", observed=True)
        .mean()
        .stack()
        .rename("score")
    )
    table = ranked.astype({"group": str, "name": str}).rename(
        columns={"name": "source", "pval": "pvalue"}
    )
    table["score"] = means.loc[
        list(zip(table["group"], table["source"], strict=True))
    ].to_numpy()
    table["group"] = pd.Categorical(table["group"], categories=groups)
    table = table.sort_values(
        ["group", "stat", "source"], ascending=[True, False, True]
    )
    return acts, table[["source", "group", "score", "stat", "pvalue", "padj"]]


def score_contrast(mat, net, p):
    """Step 4 (contrast): one score per source for the contrast statistic."""
    import decoupler as dc  # deferred: a slow import, not needed to reject bad input

    method = getattr(dc.mt, p["method"])
    # empty=False keeps zero-stat genes, so the gene count below is decoupler's.
    es, pv = method(mat, net, tmin=p["tmin"], empty=False)
    score, pv = es.iloc[0], pv.iloc[0]
    if p["method"] == "ulm":
        # decoupler returns BH-adjusted p for ulm; rebuild raw p from the t value
        # with its own df (genes - 2). A test checks BH(pvalue) == decoupler's padj.
        pvalue = 2 * stats.t.sf(np.abs(score), mat.shape[1] - 2)
        padj = pv.to_numpy()
    else:
        # mlm p is raw; BH across sources, the same family decoupler uses for ulm.
        pvalue = pv.to_numpy()
        padj = stats.false_discovery_control(pvalue, method="bh")
    table = pd.DataFrame(
        {
            "source": score.index,
            "contrast": p["contrast"],
            "score": score.to_numpy(),
            "pvalue": pvalue,
            "padj": padj,
        }
    )
    return table.sort_values(["score", "source"], ascending=[False, True])


def draw_cells(figures, acts, table, top, p):
    sources = list(dict.fromkeys(s for srcs in top.values() for s in srcs))
    heat = table[table["source"].isin(sources)].copy()
    heat["source"] = pd.Categorical(heat["source"], categories=sources)
    label = f"mean {p['method'].upper()} activity ({UNIT})"
    heat = heat.rename(columns={"score": label})
    fig = plot_heatmap(heat, "source", "group", label, limit=p["limit"])
    fig.set_layout_engine("tight")
    width = min(183, max(89, 50 + 8 * heat["group"].nunique()))
    height = max(50, 20 + 5 * len(sources))
    save_panel(fig, heat, figures, "activity_heatmap", width, height)
    plt.close(fig)

    groups = list(acts.obs[p["groupby"]].cat.categories)
    long = (
        pd.DataFrame(acts[:, sources].X, columns=sources, index=acts.obs_names)
        .assign(group=acts.obs[p["groupby"]].to_numpy())
        .rename_axis("cell")
        .reset_index()
        .melt(id_vars=["cell", "group"], var_name="source", value_name="score")
    )
    fig, axes = plt.subplots(len(sources), 1, sharex=True, squeeze=False)
    for ax, source in zip(axes[:, 0], sources, strict=True):
        sub = long[long["source"] == source]
        parts = ax.violinplot(
            [sub.loc[sub["group"] == g, "score"] for g in groups],
            positions=range(len(groups)),
            showmedians=True,
        )
        for body in parts["bodies"]:
            body.set_facecolor("grey")
        for key in ("cbars", "cmins", "cmaxes", "cmedians"):
            parts[key].set_color("black")
        ax.axhline(0, c="grey", lw=0.5)
        ax.set_ylabel(f"{source}\n{p['method'].upper()} ({UNIT})")
    counts = acts.obs[p["groupby"]].value_counts()
    axes[-1, 0].set_xticks(
        range(len(groups)), [f"{g} (n = {counts[g]})" for g in groups], rotation=90
    )
    axes[-1, 0].set_xlabel(p["groupby"])
    fig.set_layout_engine("constrained")
    save_panel(fig, long, figures, "activity_violin", width, 15 + 20 * len(sources))
    plt.close(fig)
    return sources


def draw_contrast(figures, table, top, p):
    bars = table[table["source"].isin(top)].sort_values("score")
    fig, ax = plt.subplots()
    ax.barh(bars["source"], bars["score"], color="grey")
    ax.axvline(0, c="black", lw=0.5)
    ax.set_xlabel(f"{p['method'].upper()} activity ({UNIT})")
    ax.set_ylabel("source")
    fig.set_layout_engine("tight")
    save_panel(fig, bars, figures, "activity_bar", 89, max(40, 15 + 5 * len(top)))
    plt.close(fig)


def analyse(argv, full_argv):
    p = parse_args(argv)
    seed_everything(p["seed"])
    net = load_network(p)
    cells = p["input"] == "cells"
    if cells:
        adata = load_cells(p)
        # decoupler drops genes that are zero in every cell before its own prune.
        detected = np.asarray((adata.X != 0).sum(axis=0)).ravel() > 0
        features, n_dropped_genes = adata.var_names[detected], 0
    else:
        mat, n_dropped_genes = load_contrast(p)
        features = mat.columns
    # Resource overlap: distinct network targets present in the data.
    targets = net["target"].drop_duplicates()
    n_in_resource, n_matched = len(targets), int(targets.isin(set(features)).sum())
    net, n_targets, dropped = prune(net, features, p["tmin"])

    if cells:
        acts, table = score_cells(adata, net, p)
        top = {
            str(g): d["source"].head(p["n_top"]).tolist()
            for g, d in table.groupby("group", observed=True)
        }
    else:
        table = score_contrast(mat, net, p)
        top = table.reindex(table["score"].abs().sort_values(ascending=False).index)
        top = top["source"].head(p["n_top"]).tolist()
    table = table.merge(n_targets, on="source", how="left")
    if set(table["source"]) != set(net["source"]):
        raise RuntimeError("decoupler kept different sources than the in-script prune")

    outdir = make_outdir(p["outdir"])
    tables, figures = outdir / "tables", outdir / "figures"
    table.to_csv(tables / "activity.tsv", sep="\t", index=False)
    dropped.to_csv(tables / "sources_dropped.tsv", sep="\t", index=False)

    plt.switch_backend("Agg")
    n_src = table["source"].nunique()
    if cells:
        col = p["groupby"]
        shown = draw_cells(figures, acts, table, top, p)
        captions = {
            "activity_heatmap": f"Heatmap of mean {p['method'].upper()} activity per "
            f"group of {col} for the {len(shown)} top sources (up to {p['n_top']} per "
            f"group by {p['rank_method']} statistic against the rest); colours "
            f"saturate at ±{p['limit']:g}; not clustered: rows in top-source order, "
            f"columns in {col} order.",
            "activity_violin": f"Violin plot of per-cell {p['method'].upper()} "
            f"activity of the same {len(shown)} sources by {col}, one panel "
            "per source; y axes differ between panels and the middle bar marks "
            "each median.",
        }
        summary = [
            f"{acts.n_obs} cells in {len(top)} groups of {col} scored with "
            f"{p['method']} on {n_src} sources (tmin {p['tmin']}); "
            f"{len(dropped)} sources dropped under tmin.",
            "Top source per group: "
            + "; ".join(f"{g}: {s[0]}" for g, s in top.items())
            + f" ({p['rank_method']}, group vs rest, BH within group).",
        ]
        n_units, artifacts = acts.n_obs, "figures/activity_heatmap.png"
    else:
        draw_contrast(figures, table, top, p)
        captions = {
            "activity_bar": f"Bar plot of {p['method'].upper()} activity scores for "
            f"the {len(top)} sources with the largest absolute score in contrast "
            f"{p['contrast']} ({p['stat_col']} column); positive means targets move "
            "with the source's signed weights.",
        }
        summary = [
            f"{mat.shape[1]} genes of {p['contrast']} ({p['stat_col']}) scored with "
            f"{p['method']} on {n_src} sources (tmin {p['tmin']}); "
            f"{len(dropped)} sources dropped under tmin; {n_dropped_genes} genes "
            "dropped for a missing stat.",
            f"Top {len(top)} sources by |score|: {', '.join(top)} (BH across sources).",
        ]
        n_units, artifacts = mat.shape[1], "figures/activity_bar.png"

    keys = [*REQUIRED, *MODE_ARGS[p["input"]]]
    write_report(outdir, "TF and pathway activity", p, keys, summary, captions)
    write_result(
        outdir,
        {
            "input": p["input"],
            "n_cells" if cells else "n_genes": n_units,
            "n_sources_scored": n_src,
            "n_sources_dropped": len(dropped),
            "n_matched": n_matched,
            "n_in_resource": n_in_resource,
            "top_sources": top,
            "method": p["method"],
            "seed": p["seed"],
            "params": p,
        },
        summary,
        ["report.md", "tables/activity.tsv", artifacts],
        reference=p["reference"],
    )
    write_reproducibility(
        outdir,
        __file__,
        p,
        keys,
        full_argv,
        ["decoupler", "anndata", "numpy", "pandas", "scipy", "matplotlib"],
    )
    return outdir


if __name__ == "__main__":
    raise SystemExit(run("tf-activity", analyse))
