"""Purpose: cell-type composition / differential abundance with crumblr + dream.
Inputs:  --adata <cells.h5ad> (obs only) with --sample-col, --label-col; or
         --counts <samples x labels .tsv> with --meta <samples .tsv>
Outputs: <outdir>/report.md, result.json, tables/, figures/, reproducibility/
Run:     python scripts/composition.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""

import argparse
import os
import re
import shutil
import subprocess
import tempfile
from itertools import cycle
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from palettes import categorical_palette
from plots import FEW_POINTS, ROTATE_X, plot_box, plot_composition
from skillkit import (
    check_design_rank,
    check_replication,
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
R_ENGINE = Path(__file__).resolve().parent / "crumblr.R"

ANALYSIS_ARGS = [
    "context",
    "formula",
    "coef",
    "group_col",
    "tau",
    "min_count",
    "padj_max",
    "seed",
]
INPUT_ARGS = ["adata", "sample_col", "label_col", "counts", "meta"]
KEYS = [*INPUT_ARGS, "donor_col", *ANALYSIS_ARGS]
NAME = r"[A-Za-z_][\w.]*"
RANDOM_TERM = re.compile(rf"\(\s*1\s*\|\s*({NAME})\s*\)")
TERM = rf"(?:{NAME}|{RANDOM_TERM.pattern})"
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "adata": str(SKILL_DIR / "examples" / "demo.h5ad"),
    "sample_col": "sample",
    "label_col": "cell_type",
    "counts": None,
    "meta": None,
    "donor_col": "donor",
    "context": "synthetic demo: 8 donors, ctrl vs stim, "
    "7 PBMC-like labels, no real tissue",
    "formula": "~ condition + (1|donor)",
    "coef": "conditionstim",
    "group_col": "condition",
    "tau": 1.0,
    "min_count": 50,
    "padj_max": 0.05,
    "seed": 0,
    "outdir": "composition_demo",
}

ERRORS = {
    "input_mode": "give --adata with --sample-col and --label-col, "
    "or --counts with --meta, not both",
    "bad_value": "{arg} must be {rule}, got {value}",
    "bad_formula": "--formula must be '~ a + b' with sample-level column names, "
    "or '~ a + (1|b)' with a random effect, got {value!r}",
    "missing_file": "input file not found: {path}",
    "missing_column": "column not found: {col}",
    "varies_within_sample": "column {col!r} varies within a sample; "
    "formula terms and --group-col must be sample-level",
    "bad_counts": "--counts must hold non-negative integers, one column per label",
    "sample_mismatch": "--counts and --meta disagree on samples: {samples}",
    "too_few_labels": "{n} label(s) keep >= --min-count {min_count} cells; "
    "crumblr needs at least 2",
    "bad_coef": "--coef {coef!r} is not a coefficient of the fit; choose one of: "
    "{choices}",
    "no_rscript": "Rscript not found; set RSCRIPT or put Rscript on PATH",
    "r_failed": "{script} failed with exit {code}: {tail}",
}

DA_COLUMNS = ["label", "logFC", "se", "t", "ci_low", "ci_high", "pvalue", "padj"]
# Fixed plot spec, not analysis thresholds. Groups take Paul Tol's bright set
# (reference blue, then red); labels take categorical_palette (colour-blind safe).
GROUP_COLOURS = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377"]
SIG, OTHER = "#B2182B", "0.5"
JITTER_WIDTH = 0.2


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--adata")
    parser.add_argument("--sample-col")
    parser.add_argument("--label-col")
    parser.add_argument("--counts")
    parser.add_argument("--meta")
    parser.add_argument("--donor-col")
    parser.add_argument("--context")
    parser.add_argument("--formula")
    parser.add_argument("--coef")
    parser.add_argument("--group-col")
    parser.add_argument("--tau", type=float)
    parser.add_argument("--min-count", type=int)
    parser.add_argument("--padj-max", type=float)
    parser.add_argument("--seed", type=int)
    p = parse_required(
        parser, ["donor_col", *ANALYSIS_ARGS], DEMO, argv, optional=INPUT_ARGS
    )
    h5ad_mode = [p[k] is not None for k in ("adata", "sample_col", "label_col")]
    tsv_mode = [p[k] is not None for k in ("counts", "meta")]
    if not (all(h5ad_mode) and not any(tsv_mode)) and not (
        all(tsv_mode) and not any(h5ad_mode)
    ):
        fail(ERRORS, "input_mode")
    return p


def check_values(p):
    rules = [
        ("--tau", p["tau"], p["tau"] > 0, "> 0"),
        ("--min-count", p["min_count"], p["min_count"] >= 0, "an integer >= 0"),
        ("--padj-max", p["padj_max"], 0 < p["padj_max"] <= 1, "in (0, 1]"),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ]
    for arg, value, ok, rule in rules:
        if not ok:
            fail(ERRORS, "bad_value", arg=arg, rule=rule, value=value)
    if not re.fullmatch(rf"~\s*{TERM}(\s*\+\s*{TERM})*\s*", p["formula"]):
        fail(ERRORS, "bad_formula", value=p["formula"])
    # Every term, fixed or random, is a column; the "1" of (1|x) is not matched.
    terms = re.findall(NAME, p["formula"])
    fixed = [t for t in terms if t not in RANDOM_TERM.findall(p["formula"])]
    columns = list(dict.fromkeys([*terms, p["group_col"], p["donor_col"]]))
    return terms, fixed, columns


def read_file(path, reader):
    if not Path(path).is_file():
        fail(ERRORS, "missing_file", path=path)
    return reader(path)


def sample_level(table, sample_col, columns):
    """One row per sample; each column must be constant within a sample."""
    for col in [sample_col, *columns]:
        if col not in table.columns:
            fail(ERRORS, "missing_column", col=col)
    varies = table.groupby(sample_col, observed=True)[columns].nunique() > 1
    for col in columns:
        if varies[col].any():
            fail(ERRORS, "varies_within_sample", col=col)
    meta = table.drop_duplicates(sample_col).set_index(sample_col)[columns]
    return meta.rename_axis("sample")


def read_counts(p, columns):
    """Step 1-2: samples x labels counts and sample metadata, both sorted by sample."""
    if p["adata"]:
        obs = read_file(p["adata"], lambda f: ad.read_h5ad(f, backed="r")).obs
        # Categories become strings; numeric columns stay covariates.
        obs = obs.apply(
            lambda c: c.astype(str) if isinstance(c.dtype, pd.CategoricalDtype) else c
        )
        check_unique(obs.index, "cell ids in --adata")
        for col in (p["sample_col"], p["label_col"]):
            if col not in obs.columns:
                fail(ERRORS, "missing_column", col=col)
        ids = [p["sample_col"], p["label_col"]]
        obs[ids] = obs[ids].astype(str)
        meta = sample_level(obs, p["sample_col"], columns)
        counts = pd.crosstab(obs[p["sample_col"]], obs[p["label_col"]])
        return counts.sort_index().rename_axis("sample"), meta.sort_index()

    counts = read_file(p["counts"], lambda f: pd.read_csv(f, sep="\t", index_col=0))
    values = counts.to_numpy()
    if (
        counts.shape[1] == 0
        or not np.issubdtype(values.dtype, np.number)
        or (values < 0).any()
        or (values != np.round(values)).any()
    ):
        fail(ERRORS, "bad_counts")
    check_unique(counts.index.astype(str), "sample ids in --counts")
    # Numeric columns enter the formula as covariates, as in R's read.delim.
    table = read_file(p["meta"], lambda f: pd.read_csv(f, sep="\t"))
    table[table.columns[0]] = table[table.columns[0]].astype(str)
    check_unique(table[table.columns[0]], "sample ids in --meta")
    meta = sample_level(table, table.columns[0], columns)
    counts.index, counts.columns = counts.index.astype(str), counts.columns.astype(str)
    differ = sorted(set(counts.index) ^ set(meta.index))
    if differ:
        fail(ERRORS, "sample_mismatch", samples=", ".join(differ[:5]))
    counts = counts.astype(np.int64).sort_index().rename_axis("sample")
    return counts, meta.loc[counts.index]


def run_crumblr(tmp, counts, meta, p):
    """Step 4: crumblr + dream + eBayes + topTable via crumblr.R."""
    rscript = os.environ.get("RSCRIPT") or shutil.which("Rscript")
    if not rscript or not Path(rscript).is_file():
        fail(ERRORS, "no_rscript")
    files = {k: Path(tmp) / f"{k}.tsv" for k in ("counts", "meta", "da", "clr")}
    counts.to_csv(files["counts"], sep="\t")
    meta.to_csv(files["meta"], sep="\t")
    cmd = [
        rscript,
        str(R_ENGINE),
        str(files["counts"]),
        str(files["meta"]),
        p["formula"],
        p["coef"],
        str(p["tau"]),
        str(p["seed"]),
        str(files["da"]),
        str(files["clr"]),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        choices = re.search(r"not in the fit; choose one of: (.*)", proc.stderr)
        if choices:
            fail(ERRORS, "bad_coef", coef=p["coef"], choices=choices.group(1))
        tail = " | ".join(proc.stderr.strip().splitlines()[-5:])
        fail(ERRORS, "r_failed", script=R_ENGINE.name, code=proc.returncode, tail=tail)
    da = pd.read_csv(files["da"], sep="\t", dtype={"label": str})[DA_COLUMNS]
    clr = pd.read_csv(files["clr"], sep="\t", index_col=0)
    return da, clr, [f"Rscript {rscript}", *proc.stdout.strip().splitlines()]


def paired_by(p, meta, levels):
    """The random-effect term, if each of its values has one sample per group level."""
    for term in RANDOM_TERM.findall(p["formula"]):
        per = meta.groupby(term)[p["group_col"]].agg(["size", "nunique"])
        if ((per["size"] == len(levels)) & (per["nunique"] == len(levels))).all():
            return term
    return None


def draw_figures(figures, counts, props, meta, da, p):
    group = p["group_col"]
    levels = sorted(meta[group].astype(str).unique())
    group_pal = dict(zip(levels, cycle(GROUP_COLOURS)))
    label_pal = categorical_palette(counts.columns)
    order = meta.assign(g=meta[group].astype(str)).sort_values(["g", "sample"]).index
    captions = {}

    long = (
        counts.loc[order]
        .reset_index()
        .melt("sample", var_name="label", value_name="cells")
    )
    long["sample"] = pd.Categorical(long["sample"], categories=order)
    long[group] = meta.loc[long["sample"].astype(str), group].to_numpy()
    plot = plot_composition(long, "sample", "label", "cells", label_pal) + ROTATE_X
    save_panel(plot, long, figures, "composition", 120, 70)
    captions["composition"] = (
        f"Stacked bars of the share of cells per label in each of {len(order)} "
        f"samples, ordered by {group} ({', '.join(levels)}); the denominator is "
        f"the sample's cells in the {counts.shape[1]} labels tested."
    )

    pair = paired_by(p, meta, levels)
    joined = "lines join one " + pair + "'s samples" if pair else "points jittered"
    per_level = meta[group].astype(str).value_counts()
    # plots.plot_box draws points + median under FEW_POINTS per group, else a box.
    shape = (
        "Points with the group median" if per_level.min() < FEW_POINTS else "Box plot"
    )
    for label in da["label"]:
        value = f"{label} proportion"
        rows = meta[[group]].astype(str).assign(**{value: props[label]})
        if pair:
            rows[pair] = meta[pair]
        rows = rows.loc[order].reset_index()
        plot = plot_box(
            rows,
            group,
            value,
            group_pal,
            jitter_width=JITTER_WIDTH,
            seed=p["seed"],
            value_label=f"{label} (fraction of sample's cells)",
            pair=pair,
        )
        name = "proportion_" + re.sub(r"\W+", "_", label).strip("_")
        save_panel(plot, rows, figures, name, 60, 60)
        captions[name] = (
            f"{shape} of the {label} share of each sample's cells per {group}; "
            f"each point is one sample, {joined}. crumblr {p['coef']} logFC "
            f"{da.set_index('label').loc[label, 'logFC']:.2f}, padj "
            f"{da.set_index('label').loc[label, 'padj']:.2g}."
        )

    forest = da.sort_values("logFC").reset_index(drop=True)
    fig, ax = plt.subplots(layout="constrained")
    colours = np.where(forest["significant"], SIG, OTHER)
    y = np.arange(len(forest))
    ax.errorbar(
        forest["logFC"],
        y,
        xerr=[forest["logFC"] - forest["ci_low"], forest["ci_high"] - forest["logFC"]],
        fmt="none",
        ecolor=colours,
        lw=0.8,
    )
    ax.scatter(forest["logFC"], y, s=10, c=colours, zorder=3)
    ax.axvline(0, lw=0.5, ls="--", c="black")
    ax.set_yticks(y, forest["label"])
    ax.set_xlabel(f"{p['coef']} logFC (CLR scale, 95% CI)")
    ax.set_ylabel("cell-type label")
    save_panel(fig, forest, figures, "effect_forest", 89, 70)
    plt.close(fig)
    captions["effect_forest"] = (
        f"Forest plot of the crumblr {p['coef']} effect per label: logFC on the "
        f"CLR scale with its 95% confidence interval; red labels have padj < "
        f"{p['padj_max']} (BH across {len(da)} labels), the dashed line marks zero."
    )
    return captions, pair


def analyse(argv, full_argv):
    p = parse_args(argv)
    terms, fixed, columns = check_values(p)
    seed_everything(p["seed"])
    counts, meta = read_counts(p, columns)
    check_design_rank(meta, fixed)
    n_donors, n_samples = check_replication(
        meta.reset_index(), "sample", p["donor_col"], terms
    )

    # Step 3: drop labels with fewer than --min-count cells over all samples.
    totals = counts.sum()
    dropped = sorted(totals.index[totals < p["min_count"]])
    # Cells removed with the dropped labels, per --group-col level (no sample drops).
    by_level = counts[dropped].sum(axis=1).groupby(meta[p["group_col"]].astype(str))
    dropped_cells = {k: int(v) for k, v in by_level.sum().items()}
    counts = counts.drop(columns=dropped)
    if counts.shape[1] < 2:
        fail(ERRORS, "too_few_labels", n=counts.shape[1], min_count=p["min_count"])

    outdir = make_outdir(p["outdir"])
    tables, figures = outdir / "tables", outdir / "figures"
    with tempfile.TemporaryDirectory() as tmp:
        da, clr, r_versions = run_crumblr(tmp, counts, meta, p)

    props = counts.div(counts.sum(axis=1), axis=0)
    counts.to_csv(tables / "counts.tsv", sep="\t")
    props.to_csv(tables / "proportions.tsv", sep="\t", float_format="%.6g")
    clr.rename_axis("sample").to_csv(tables / "clr.tsv", sep="\t", float_format="%.6g")
    da["significant"] = da["padj"] < p["padj_max"]
    # Labels under --min-count stay in the table with NA statistics (data.status).
    da = with_status(da, totals.index, "label", dict.fromkeys(dropped, "not_tested"))
    da = da.sort_values(["padj", "label"]).reset_index(drop=True)
    status_counts = {k: int(v) for k, v in da["status"].value_counts().items()}
    da.to_csv(tables / "da_results.tsv", sep="\t", index=False, float_format="%.6g")

    plt.switch_backend("Agg")
    tested = da[da["status"] == "tested"]
    captions, pair = draw_figures(figures, counts, props, meta, tested, p)

    sig = da[da["significant"].eq(True)]
    signs = {r.label: "up" if r.logFC > 0 else "down" for r in sig.itertuples()}
    summary = [
        f"{n_samples} samples from {n_donors} donors, {counts.shape[1]} labels "
        f"tested ({', '.join(dropped) or 'none'} dropped by --min-count "
        f"{p['min_count']}; cells removed per {p['group_col']}: "
        + ", ".join(f"{k} {v}" for k, v in dropped_cells.items())
        + ").",
        f"{len(sig)} labels change ({p['coef']}, padj < {p['padj_max']}, BH): "
        + (", ".join(f"{k} {v}" for k, v in signs.items()) or "none")
        + ".",
    ]
    # A bare "|" in (1|x) would split the markdown table cell.
    shown = {**p, "formula": p["formula"].replace("|", "\\|")}
    write_report(outdir, f"Composition: {p['coef']}", shown, KEYS, summary, captions)
    result = {
        "n_samples": n_samples,
        "n_donors": n_donors,
        "dropped_cells_by_level": dropped_cells,
        "n_labels_tested": counts.shape[1],
        "dropped_labels": dropped,
        "status_counts": status_counts,
        "n_significant": len(sig),
        "significant": signs,
        "paired_by": pair,
        "engine": "crumblr",
        "params": p,
        "seed": p["seed"],
    }
    artifacts = ["report.md", "tables/da_results.tsv", "figures/effect_forest.png"]
    write_result(outdir, result, summary, artifacts, context=p["context"])
    packages = ["anndata", "numpy", "pandas", "matplotlib", "plotnine"]
    write_reproducibility(
        outdir, __file__, p, KEYS, full_argv, packages, extra=r_versions
    )
    return outdir


if __name__ == "__main__":
    raise SystemExit(run("composition", analyse))
