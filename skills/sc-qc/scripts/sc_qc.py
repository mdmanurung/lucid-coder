"""Purpose: scRNA cell QC and Scrublet doublet removal, per sample, from a raw h5ad.
Inputs:  --adata <file.h5ad> with raw integer counts in X or --layer; a sample column
Outputs: <outdir>/report.md, result.json, data/qc_filtered.h5ad, tables/, figures/,
         reproducibility/
Run:     python scripts/sc_qc.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""

import argparse
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plots import grey_to_red
from scipy import sparse
from scipy.stats import median_abs_deviation
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
    "sample_col",
    "mito_prefix",
    "mode",
    "expected_doublet_rate",
    "doublet_threshold",
    "seed",
    "reference",
    "context",
]
MODE_ARGS = {
    "fixed": ["min_genes", "max_genes", "min_counts", "max_pct_mito"],
    "mad": ["mads_counts", "mads_genes", "mads_mito"],
}
OPTIONAL = ["layer", *MODE_ARGS["fixed"], *MODE_ARGS["mad"]]
KEYS = [*REQUIRED, *OPTIONAL]
# Demo-only values, not defaults: SKILL.md lists them under the same label.
# MAD counts follow cns-haematopoiesis (5, 5, 3 MADs).
DEMO = {
    "adata": str(SKILL_DIR / "examples" / "demo.h5ad"),
    "sample_col": "sample",
    "mito_prefix": "MT-",
    "mode": "mad",
    "expected_doublet_rate": 0.08,
    "doublet_threshold": "auto",
    "seed": 0,
    "reference": "synthetic demo genes, no genome build",
    "context": "synthetic demo: 3 samples of 5 populations, planted "
    "low-quality cells and doublets, no real tissue",
    "layer": None,
    **dict.fromkeys(MODE_ARGS["fixed"]),
    "mads_counts": 5.0,
    "mads_genes": 5.0,
    "mads_mito": 3.0,
    "outdir": "sc_qc_demo",
}

ERRORS = {
    "mode_args_missing": "--mode {mode} needs: {args}; ask the user, or use --demo",
    "mode_args_unused": "--mode {mode} does not use: {args}; drop them",
    "bad_value": "{arg} must be {rule}, got {value}",
    "missing_file": "input file not found: {path}",
    "missing_layer": "layer not found in the h5ad: {layer}",
    "missing_column": "obs column not found: {col}",
    "bad_counts": "counts in {source} must be non-negative integers (raw counts)",
    "no_mito_genes": "no gene name starts with --mito-prefix {prefix!r}; "
    "check the species and case (MT- human, mt- mouse)",
    "too_few_cells": "sample {sample!r} has {n} cells; Scrublet needs more than "
    "{pcs} (one per principal component)",
    "scrublet_failed": "Scrublet failed: {message}",
    "no_threshold": "Scrublet found no doublet-score threshold for sample {sample!r}; "
    "ask the user for a numeric --doublet-threshold",
    "no_cells_left": "zero cells pass QC and doublet removal",
}

# Ribosomal genes are flagged for the report only, never thresholded.
RIBO_PREFIXES = ("RPS", "RPL")
# Fixed plot spec, not analysis thresholds.
N_SCORE_BINS = 40
# Scrublet's own default PC count, written out so the cell guard matches it.
SCRUBLET_PCS = 30
FACET_COLS = 4
STATUS_COLOURS = {"kept": "0.3", "failed QC": "#B2182B", "doublet": "#2166AC"}
METRICS = {
    "n_genes_by_counts": ("genes detected", "min_genes", "max_genes"),
    "total_counts": ("UMI counts", "min_counts", "max_counts"),
    "pct_counts_mito": ("% mitochondrial counts", None, "max_pct_mito"),
}


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--adata")
    parser.add_argument("--layer")
    parser.add_argument("--sample-col")
    parser.add_argument("--mito-prefix")
    parser.add_argument("--mode", choices=list(MODE_ARGS))
    parser.add_argument("--min-genes", type=int)
    parser.add_argument("--max-genes", type=int)
    parser.add_argument("--min-counts", type=int)
    parser.add_argument("--max-pct-mito", type=float)
    parser.add_argument("--mads-counts", type=float)
    parser.add_argument("--mads-genes", type=float)
    parser.add_argument("--mads-mito", type=float)
    parser.add_argument("--expected-doublet-rate", type=float)
    parser.add_argument("--doublet-threshold")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--reference")
    parser.add_argument("--context")
    p = parse_required(parser, REQUIRED, DEMO, argv, optional=OPTIONAL)
    flags = " ".join

    missing = [f"--{k.replace('_', '-')}" for k in MODE_ARGS[p["mode"]] if p[k] is None]
    if missing:
        fail(ERRORS, "mode_args_missing", mode=p["mode"], args=flags(missing))
    unused = [
        f"--{k.replace('_', '-')}"
        for mode, keys in MODE_ARGS.items()
        if mode != p["mode"]
        for k in keys
        if p[k] is not None
    ]
    if unused:
        fail(ERRORS, "mode_args_unused", mode=p["mode"], args=flags(unused))
    return p


def check_values(p):
    rate = p["expected_doublet_rate"]
    rules = [
        ("--expected-doublet-rate", rate, 0 < rate < 1, "in (0, 1)"),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ]
    if p["mode"] == "fixed":
        rules += [
            ("--min-genes", p["min_genes"], p["min_genes"] >= 0, ">= 0"),
            ("--min-counts", p["min_counts"], p["min_counts"] >= 0, ">= 0"),
            (
                "--max-genes",
                p["max_genes"],
                p["max_genes"] > p["min_genes"],
                "> --min-genes",
            ),
            (
                "--max-pct-mito",
                p["max_pct_mito"],
                0 < p["max_pct_mito"] <= 100,
                "in (0, 100]",
            ),
        ]
    else:
        rules += [
            (f"--{k.replace('_', '-')}", p[k], p[k] > 0, "> 0")
            for k in MODE_ARGS["mad"]
        ]
    threshold = p["doublet_threshold"]
    try:
        ok = threshold == "auto" or 0 < float(threshold) < 1
    except ValueError:
        ok = False
    rules.append(("--doublet-threshold", threshold, ok, "'auto' or a number in (0, 1)"))
    for arg, value, ok, rule in rules:
        if not ok:
            fail(ERRORS, "bad_value", arg=arg, rule=rule, value=value)


def load(p):
    """Step 1 input: raw counts and the sample column, validated."""
    path = Path(p["adata"])
    if not path.is_file():
        fail(ERRORS, "missing_file", path=path)
    adata = ad.read_h5ad(path)
    check_unique(adata.var_names, "gene names")
    check_unique(adata.obs_names, "cell ids")
    if p["layer"] is not None and p["layer"] not in adata.layers:
        fail(ERRORS, "missing_layer", layer=p["layer"])
    if p["sample_col"] not in adata.obs.columns:
        fail(ERRORS, "missing_column", col=p["sample_col"])
    counts = sparse.csr_matrix(adata.layers[p["layer"]] if p["layer"] else adata.X)
    if not (np.all(counts.data >= 0) and np.all(counts.data == np.round(counts.data))):
        fail(ERRORS, "bad_counts", source=f"layer {p['layer']}" if p["layer"] else "X")
    if not adata.var_names.str.startswith(p["mito_prefix"]).any():
        fail(ERRORS, "no_mito_genes", prefix=p["mito_prefix"])
    return adata, counts


def qc_metrics(adata, counts, p):
    """Step 1: per-cell QC metrics on raw counts, mito and ribo flags."""
    # Imported here: scanpy's import is slow, and bad input should fail fast.
    import scanpy as sc

    cells = ad.AnnData(
        X=counts.astype(np.float32),
        obs=pd.DataFrame(
            {"sample": adata.obs[p["sample_col"]].astype(str).to_numpy()},
            index=adata.obs_names,
        ),
        var=pd.DataFrame(index=adata.var_names),
    )
    cells.var["mito"] = cells.var_names.str.startswith(p["mito_prefix"])
    cells.var["ribo"] = cells.var_names.str.upper().str.startswith(RIBO_PREFIXES)
    # percent_top=None: the top-gene shares are unused and fail on small panels.
    sc.pp.calculate_qc_metrics(
        cells, qc_vars=["mito", "ribo"], percent_top=None, inplace=True, log1p=True
    )
    return cells


def scrublet(cells, p):
    """Step 2: Scrublet per sample, seeded; returns scores, calls, thresholds."""
    import scanpy as sc

    sc.settings.n_jobs = 1
    threshold = (
        None if p["doublet_threshold"] == "auto" else float(p["doublet_threshold"])
    )
    sizes = cells.obs["sample"].value_counts()
    if (sizes <= SCRUBLET_PCS).any():
        small = sizes.idxmin()
        fail(ERRORS, "too_few_cells", sample=small, n=sizes[small], pcs=SCRUBLET_PCS)
    scored = cells.copy()
    try:
        sc.pp.scrublet(
            scored,
            batch_key="sample",
            expected_doublet_rate=p["expected_doublet_rate"],
            threshold=threshold,
            n_prin_comps=SCRUBLET_PCS,
            # Exact kNN: annoy's approximate search is not reproducible.
            use_approx_neighbors=False,
            random_state=p["seed"],
            verbose=False,
        )
    except ValueError as exc:
        fail(ERRORS, "scrublet_failed", message=str(exc).splitlines()[0])
    thresholds = {}
    for sample, info in scored.uns["scrublet"]["batches"].items():
        if "threshold" not in info:
            fail(ERRORS, "no_threshold", sample=sample)
        thresholds[sample] = float(info["threshold"])
    calls = scored.obs[["doublet_score"]].copy()
    # Scrublet's own rule (score > threshold); a NaN score (empty cell) is no call.
    cutoff = scored.obs["sample"].map(thresholds).astype(float)
    calls["predicted_doublet"] = calls["doublet_score"] > cutoff
    return calls, pd.Series(thresholds, name="doublet_threshold")


def limits_per_sample(obs, p):
    """Step 3: QC limits per sample, on the raw scale, for both modes."""
    samples = sorted(obs["sample"].unique())
    if p["mode"] == "fixed":
        return pd.DataFrame(
            {
                "min_counts": float(p["min_counts"]),
                "max_counts": np.inf,
                "min_genes": float(p["min_genes"]),
                "max_genes": float(p["max_genes"]),
                "max_pct_mito": p["max_pct_mito"],
            },
            index=pd.Index(samples, name="sample"),
        )
    # Unscaled MAD (scipy default scale=1.0), as cns-haematopoiesis; two-sided on
    # log1p counts and genes, upper only on % mito.
    by = obs.groupby("sample")
    cols = ["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mito"]
    median, mad = by[cols].median(), by[cols].agg(median_abs_deviation)

    def log_limit(col, n, sign):
        return np.expm1(median[f"log1p_{col}"] + sign * n * mad[f"log1p_{col}"])

    return pd.DataFrame(
        {
            "min_counts": log_limit("total_counts", p["mads_counts"], -1),
            "max_counts": log_limit("total_counts", p["mads_counts"], 1),
            "min_genes": log_limit("n_genes_by_counts", p["mads_genes"], -1),
            "max_genes": log_limit("n_genes_by_counts", p["mads_genes"], 1),
            "max_pct_mito": median["pct_counts_mito"]
            + p["mads_mito"] * mad["pct_counts_mito"],
        }
    ).loc[samples]


def qc_status(obs, limits):
    """First failing rule names the cell's status; a cell passes inside every limit."""
    lim = limits.loc[obs["sample"]].set_axis(obs.index)
    rules = {
        "low_counts": obs["total_counts"] < lim["min_counts"],
        "high_counts": obs["total_counts"] > lim["max_counts"],
        "low_genes": obs["n_genes_by_counts"] < lim["min_genes"],
        "high_genes": obs["n_genes_by_counts"] > lim["max_genes"],
        "high_mito": obs["pct_counts_mito"] > lim["max_pct_mito"],
    }
    return np.select(list(rules.values()), list(rules), default="pass")


def facet_grid(n):
    ncols = min(n, FACET_COLS)
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(
        nrows, ncols, sharex=True, sharey=True, squeeze=False, layout="constrained"
    )
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    return fig, list(axes.flat)


def draw_figures(figures, per_cell, limits, thresholds):
    samples = list(limits.index)
    n = len(samples)
    width = min(183, 50 * min(n, FACET_COLS) + 20)
    height = 45 * -(-n // FACET_COLS) + 10

    # QC violins: one panel per metric, one violin per sample, limits as segments.
    long = per_cell.melt(
        id_vars=["cell", "sample"], value_vars=list(METRICS), var_name="metric"
    )
    fig, axes = plt.subplots(1, len(METRICS), layout="constrained")
    for ax, (metric, (label, lo, hi)) in zip(axes, METRICS.items(), strict=True):
        values = [per_cell.loc[per_cell["sample"] == s, metric] for s in samples]
        parts = ax.violinplot(values, positions=range(n), showextrema=False)
        for body in parts["bodies"]:
            body.set_facecolor("0.6")
            body.set_alpha(0.8)
        for i, s in enumerate(samples):
            for col in (lo, hi):
                if col and np.isfinite(limits.loc[s, col]):
                    ax.hlines(
                        limits.loc[s, col],
                        i - 0.4,
                        i + 0.4,
                        ls="--",
                        lw=0.5,
                        colors="black",
                    )
        if metric != "pct_counts_mito":
            ax.set_yscale("log")
        counts = per_cell["sample"].value_counts()
        ax.set_xticks(
            range(n), [f"{s} (n = {counts[s]})" for s in samples], rotation=90
        )
        ax.set_ylabel(label)
    save_panel(fig, long, figures, "qc_violins", 183, 70)
    plt.close(fig)

    # Counts vs genes, coloured by % mito, one panel per sample with its limits.
    fig, axes = facet_grid(n)
    for ax, s in zip(axes, samples, strict=False):
        cells = per_cell[per_cell["sample"] == s]
        dots = ax.scatter(
            cells["total_counts"],
            cells["n_genes_by_counts"],
            c=cells["pct_counts_mito"],
            cmap=grey_to_red,
            vmin=0,
            vmax=per_cell["pct_counts_mito"].max(),
            s=2,
            rasterized=True,
        )
        for col in ("min_counts", "max_counts"):
            if np.isfinite(limits.loc[s, col]):
                ax.axvline(limits.loc[s, col], ls="--", lw=0.5, c="black")
        for col in ("min_genes", "max_genes"):
            ax.axhline(limits.loc[s, col], ls="--", lw=0.5, c="black")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.annotate(s, (0.03, 0.97), xycoords="axes fraction", va="top")
        ax.set_xlabel("UMI counts")
        ax.set_ylabel("genes detected")
    fig.colorbar(dots, ax=axes[:n], label="% mitochondrial counts")
    save_panel(
        fig,
        per_cell[
            ["cell", "sample", "total_counts", "n_genes_by_counts", "pct_counts_mito"]
        ],
        figures,
        "counts_vs_genes",
        width + 20,
        height,
    )
    plt.close(fig)

    # Doublet-score histogram per sample, with that sample's Scrublet cutoff.
    bins = np.linspace(0, per_cell["doublet_score"].max(), N_SCORE_BINS + 1)
    fig, axes = facet_grid(n)
    for ax, s in zip(axes, samples, strict=False):
        scores = per_cell.loc[per_cell["sample"] == s, "doublet_score"]
        ax.hist(scores, bins=bins, color="0.6", edgecolor="white", lw=0.3)
        ax.axvline(thresholds[s], ls="--", lw=0.5, c="black")
        ax.annotate(s, (0.03, 0.97), xycoords="axes fraction", va="top")
        ax.set_xlabel("Scrublet doublet score")
        ax.set_ylabel("cells")
    save_panel(
        fig,
        per_cell[["cell", "sample", "doublet_score", "predicted_doublet"]],
        figures,
        "doublet_scores",
        width,
        height,
    )
    plt.close(fig)

    # Cells per sample, stacked by fate; the bar height is the count before QC.
    fate = pd.crosstab(per_cell["sample"], per_cell["fate"]).reindex(
        index=samples, columns=list(STATUS_COLOURS), fill_value=0
    )
    fig, ax = plt.subplots(layout="constrained")
    bottom = np.zeros(n)
    for status, colour in STATUS_COLOURS.items():
        ax.bar(samples, fate[status], bottom=bottom, color=colour, label=status)
        bottom += fate[status].to_numpy()
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_xlabel("sample")
    ax.set_ylabel("cells")
    source = fate.reset_index().melt(
        id_vars="sample", var_name="fate", value_name="cells"
    )
    save_panel(fig, source, figures, "cells_per_sample", 89, 70)
    plt.close(fig)


def analyse(argv, argv_full):
    p = parse_args(argv)
    check_values(p)
    seed_everything(p["seed"])
    adata, counts = load(p)
    with threadpool_limits(1):
        cells = qc_metrics(adata, counts, p)
        calls, thresholds = scrublet(cells, p)
    obs = cells.obs.join(calls, validate="1:1")
    limits = limits_per_sample(obs, p)
    obs["qc_status"] = qc_status(obs, limits)
    obs["kept"] = (obs["qc_status"] == "pass") & ~obs["predicted_doublet"]
    obs["fate"] = np.select(
        [obs["kept"], obs["qc_status"] != "pass"], ["kept", "failed QC"], "doublet"
    )
    if not obs["kept"].any():
        fail(ERRORS, "no_cells_left")
    # Created only after every input check, so a refused run leaves no outdir.
    outdir = make_outdir(p["outdir"])
    (outdir / "data").mkdir()

    columns = [
        "sample",
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mito",
        "pct_counts_ribo",
        "doublet_score",
        "predicted_doublet",
        "qc_status",
        "fate",
        "kept",
    ]
    per_cell = (
        obs[columns].rename_axis("cell").reset_index().sort_values(["sample", "cell"])
    )
    per_cell.to_csv(outdir / "tables" / "qc_per_cell.tsv", sep="\t", index=False)
    per_sample = (
        per_cell.groupby("sample")
        .agg(
            n_cells_before=("cell", "size"),
            n_failed_qc=("qc_status", lambda s: (s != "pass").sum()),
            n_doublets=("predicted_doublet", "sum"),
            n_cells_after=("kept", "sum"),
        )
        .join(limits)
        .join(thresholds)
        .rename_axis("sample")
    )
    per_sample.to_csv(outdir / "tables" / "cells_per_sample.tsv", sep="\t")

    filtered = adata[obs.index[obs["kept"]]].copy()
    qc_cols = [c for c in columns if c != "sample"]
    filtered.obs[qc_cols] = obs.loc[filtered.obs_names, qc_cols]
    filtered.write_h5ad(outdir / "data" / "qc_filtered.h5ad")

    plt.switch_backend("Agg")
    draw_figures(outdir / "figures", per_cell, limits, thresholds)

    n_before, n_after = len(per_cell), int(per_cell["kept"].sum())
    n_failed = int((per_cell["qc_status"] != "pass").sum())
    n_doublets = int(per_cell["predicted_doublet"].sum())
    rule = (
        f"{p['mads_counts']}/{p['mads_genes']}/{p['mads_mito']} MADs per sample "
        "(counts/genes/% mito)"
        if p["mode"] == "mad"
        else f"genes {p['min_genes']}-{p['max_genes']}, counts >= {p['min_counts']}, "
        f"% mito <= {p['max_pct_mito']}"
    )
    summary = [
        f"{n_after} of {n_before} cells kept across {len(per_sample)} samples.",
        f"{n_failed} cells failed QC ({p['mode']} mode: {rule}); "
        f"{n_doublets} Scrublet doublets (threshold {p['doublet_threshold']}).",
    ]
    captions = {
        "qc_violins": "Violin plots of genes detected, UMI counts (log scale) and "
        "% mitochondrial counts per sample, before filtering; dashed segments mark "
        "each sample's QC limits.",
        "counts_vs_genes": "Scatter of UMI counts against genes detected (log axes), "
        "one panel per sample, coloured by % mitochondrial counts; dashed lines mark "
        "that sample's count and gene limits.",
        "doublet_scores": f"Histograms of Scrublet doublet scores per sample in "
        f"{N_SCORE_BINS} bins; the dashed line marks that sample's doublet-score "
        "threshold.",
        "cells_per_sample": "Stacked bars of cells per sample by fate: kept, failed "
        "QC, or passed QC but called a doublet; bar height is the count before QC.",
    }
    write_report(outdir, "scRNA QC and doublet removal", p, KEYS, summary, captions)
    write_result(
        outdir,
        {
            "n_cells_before": n_before,
            "n_cells_after": n_after,
            "n_failed_qc": n_failed,
            "n_doublets": n_doublets,
            "cells_after_per_sample": per_sample["n_cells_after"].astype(int).to_dict(),
            "failed_qc_per_sample": per_sample["n_failed_qc"].astype(int).to_dict(),
            "doublets_per_sample": per_sample["n_doublets"].astype(int).to_dict(),
            "doublet_thresholds": thresholds.to_dict(),
            "mode": p["mode"],
            "seed": p["seed"],
            "params": p,
        },
        summary,
        ["report.md", "tables/cells_per_sample.tsv", "figures/counts_vs_genes.png"],
        context=p["context"],
        reference=p["reference"],
    )
    packages = [
        "anndata",
        "scanpy",
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "plotnine",
        "scikit-image",
    ]
    write_reproducibility(outdir, __file__, p, KEYS, argv_full, packages)
    return outdir


if __name__ == "__main__":
    raise SystemExit(run("sc-qc", analyse))
