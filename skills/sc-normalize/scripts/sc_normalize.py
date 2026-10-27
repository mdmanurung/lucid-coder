"""Purpose: normalise one count h5ad by log1p, analytic Pearson residuals or ADT dsb.
Inputs:  --adata <file.h5ad> with raw integer counts in X (RNA, or ADT for dsb);
         --empty-droplets <file.h5ad> with ADT counts of empty drops (dsb only)
Outputs: <outdir>/report.md, result.json, data/normalized.h5ad, tables/, figures/,
         reproducibility/
Run:     python scripts/sc_normalize.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""

import argparse
import warnings
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plots import plot_scatter_fit
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
from threadpoolctl import threadpool_limits

SKILL_DIR = Path(__file__).resolve().parents[1]
EXAMPLES = SKILL_DIR / "examples"
REQUIRED = ["adata", "method", "layer_out", "seed", "reference"]
METHOD_ARGS = {
    "log1p": ["target_sum"],
    "pearson": ["theta", "clip"],
    "dsb": ["empty_droplets", "isotype_controls"],
}
OPTIONAL = [k for keys in METHOD_ARGS.values() for k in keys]
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "method": "pearson",
    "layer_out": "normalized",
    "seed": 0,
    "reference": "synthetic demo features, no genome build",
    "outdir": "sc_normalize_demo",
}
DEMO_BY_METHOD = {
    "log1p": {"adata": str(EXAMPLES / "demo.h5ad"), "target_sum": 1e4},
    "pearson": {"adata": str(EXAMPLES / "demo.h5ad"), "theta": 100.0, "clip": 25.0},
    "dsb": {
        "adata": str(EXAMPLES / "demo_adt.h5ad"),
        "empty_droplets": str(EXAMPLES / "demo_empty.h5ad"),
        "isotype_controls": "IgG1,IgG2a,IgG2b",
    },
}

ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "method_needs": "--method {method} needs: {flags}",
    "method_unused": "--method {method} does not use: {flags}",
    "missing_file": "input file not found: {path}",
    "bad_counts": "counts in X of {path} must be non-negative integers (raw counts)",
    "empty_cells": "{n} cells have zero total counts; filter them first (sc-qc)",
    "isotype_missing": "isotype controls not in the ADT features: {names}",
    "panel_mismatch": "--adata and --empty-droplets must have the same features "
    "in the same order",
    "barcode_overlap": "{n} barcodes are in both --adata and --empty-droplets",
    "flat_empty": "features with zero variance in the empty droplets: {names}",
}

# Fixed steps, not analysis parameters. dsb: log(x + 10) and standardise by the
# empty droplets, then regress out the technical component (Mule et al. 2022,
# muon's defaults). Diagnostics: PC1 from 2 seeded arpack PCs; the 6 most
# variable features in the value-distribution figure; 20 quantile bins for the
# mean-variance trend.
DSB_PSEUDOCOUNT = 10
N_DIAG_PCS = 2
N_TOP_SHOWN = 6
TREND_BINS = 20
PACKAGES = [
    "scanpy",
    "anndata",
    "muon",
    "mudata",
    "scikit-learn",
    "numpy",
    "pandas",
    "scipy",
    "matplotlib",
    "threadpoolctl",
]


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--adata")
    parser.add_argument("--method", choices=list(METHOD_ARGS))
    parser.add_argument("--layer-out")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--reference")
    parser.add_argument("--target-sum", type=float)
    parser.add_argument("--theta", type=float)
    parser.add_argument("--clip", type=float)
    parser.add_argument("--empty-droplets")
    parser.add_argument("--isotype-controls")
    p = parse_required(
        parser, REQUIRED, DEMO, argv, optional=OPTIONAL, demo_extra=("outdir", "method")
    )
    p = {k: p.get(k) for k in [*REQUIRED, *OPTIONAL, "outdir"]}
    if "--demo" in argv:
        p.update(DEMO_BY_METHOD[p["method"]])
    return p


def check_values(p):
    method = p["method"]
    needs = [k for k in METHOD_ARGS[method] if p[k] is None]
    if needs:
        flags = " ".join(f"--{k.replace('_', '-')}" for k in needs)
        fail(ERRORS, "method_needs", method=method, flags=flags)
    unused = [k for k in OPTIONAL if k not in METHOD_ARGS[method] and p[k] is not None]
    if unused:
        flags = " ".join(f"--{k.replace('_', '-')}" for k in unused)
        fail(ERRORS, "method_unused", method=method, flags=flags)
    rules = [
        ("--target-sum", p["target_sum"], "> 0"),
        ("--theta", p["theta"], "> 0"),
        ("--clip", p["clip"], "> 0"),
    ]
    for arg, value, rule in rules:
        if value is not None and value <= 0:
            fail(ERRORS, "bad_value", arg=arg, rule=rule, value=value)
    if p["seed"] < 0:
        fail(ERRORS, "bad_value", arg="--seed", rule="an integer >= 0", value=p["seed"])
    if p["layer_out"] in ("counts", ""):
        rule = "a layer name other than 'counts', which keeps the input"
        fail(ERRORS, "bad_value", arg="--layer-out", rule=rule, value=p["layer_out"])


def read_counts(path):
    path = Path(path)
    if not path.is_file():
        fail(ERRORS, "missing_file", path=path)
    adata = ad.read_h5ad(path)
    check_unique(adata.var_names, f"feature names in {path}")
    check_unique(adata.obs_names, f"cell ids in {path}")
    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    if not (np.all(values >= 0) and np.all(values == np.round(values))):
        fail(ERRORS, "bad_counts", path=path)
    return adata


def read_inputs(p):
    """Cells (and, for dsb, empty droplets); refused before the outdir exists."""
    adata = read_counts(p["adata"])
    n_empty_cells = int((np.asarray(adata.X.sum(axis=1)).ravel() == 0).sum())
    if n_empty_cells:
        fail(ERRORS, "empty_cells", n=n_empty_cells)
    if p["method"] != "dsb":
        return adata, None
    empty = read_counts(p["empty_droplets"])
    if list(empty.var_names) != list(adata.var_names):
        fail(ERRORS, "panel_mismatch")
    isotypes = p["isotype_controls"].split(",")
    missing = [name for name in isotypes if name not in adata.var_names]
    if missing:
        fail(ERRORS, "isotype_missing", names=", ".join(missing))
    # muon silently drops empty droplets whose barcode is also a cell.
    overlap = int(empty.obs_names.isin(adata.obs_names).sum())
    if overlap:
        fail(ERRORS, "barcode_overlap", n=overlap)
    # A constant count gives SD 0 after log, so dsb would divide by zero.
    empty_counts = dense(empty.X)
    flat = empty.var_names[(empty_counts == empty_counts[0]).all(axis=0)]
    if len(flat):
        fail(ERRORS, "flat_empty", names=", ".join(flat))
    return adata, empty


def dense(x):
    return x.toarray() if sparse.issparse(x) else np.asarray(x)


def normalise(adata, empty, p):
    """Steps 3-5 in place: drop all-zero features, keep counts, normalise X."""
    # Imported here: scanpy and muon pull in numba (~20 s), which a refused run
    # does not need.
    import scanpy as sc

    keep = np.asarray(adata.X.sum(axis=0)).ravel() > 0
    n_dropped = int((~keep).sum())
    adata._inplace_subset_var(keep)
    if empty is not None:
        empty._inplace_subset_var(keep)
    adata.layers["counts"] = adata.X.copy()
    method = p["method"]
    if method == "log1p":
        sc.pp.normalize_total(adata, target_sum=p["target_sum"])
        sc.pp.log1p(adata)
    if method == "pearson":
        sc.experimental.pp.normalize_pearson_residuals(
            adata, theta=p["theta"], clip=p["clip"]
        )
    if method == "dsb":
        import muon as mu

        with warnings.catch_warnings():
            # No empty_counts_range: muon warns that every non-cell is empty,
            # which is what --empty-droplets holds.
            warnings.simplefilter("ignore", UserWarning)
            mu.prot.pp.dsb(
                adata,
                empty,
                pseudocount=DSB_PSEUDOCOUNT,
                denoise_counts=True,
                isotype_controls=p["isotype_controls"].split(","),
                scale_factor="standardize",
                random_state=p["seed"],
            )
    adata.X = dense(adata.X)
    adata.layers[p["layer_out"]] = adata.X.copy()
    return n_dropped


def tables(adata, seed):
    """Per-cell and per-feature summaries, plus PC1 of the normalised values."""
    from sklearn.decomposition import PCA

    counts = adata.layers["counts"]
    per_cell = pd.DataFrame(
        {
            "cell": adata.obs_names,
            "library_size": np.asarray(counts.sum(axis=1)).ravel(),
            "normalised_total": adata.X.sum(axis=1),
            "n_detected": np.asarray((counts > 0).sum(axis=1)).ravel(),
        }
    )
    pca = PCA(n_components=N_DIAG_PCS, svd_solver="arpack", random_state=seed)
    per_cell["pc1"] = pca.fit_transform(adata.X)[:, 0]
    per_cell = per_cell.sort_values("cell", ignore_index=True)

    counts = dense(counts).astype(float)
    per_feature = pd.DataFrame(
        {
            "feature": adata.var_names,
            "mean_counts": counts.mean(axis=0),
            "var_counts": counts.var(axis=0, ddof=1),
            "mean_normalised": adata.X.mean(axis=0),
            "var_normalised": adata.X.var(axis=0, ddof=1),
        }
    ).sort_values("feature", ignore_index=True)
    return per_cell, per_feature


def draw_figures(figures, adata, per_cell, per_feature, method, isotypes):
    stages = {"counts": "var_counts", "normalised": "var_normalised"}
    long = pd.concat(
        [
            pd.DataFrame(
                {
                    "feature": per_feature["feature"],
                    "stage": stage,
                    "log10_mean_counts": np.log10(per_feature["mean_counts"]),
                    "log10_variance": np.log10(per_feature[col]),
                }
            )
            for stage, col in stages.items()
        ],
        ignore_index=True,
    )
    long = long[np.isfinite(long["log10_variance"])]
    fig, axes = plt.subplots(1, 2, sharex=True, layout="constrained")
    for ax, stage in zip(axes, stages, strict=True):
        part = long[long["stage"] == stage]
        ax.scatter(
            part["log10_mean_counts"], part["log10_variance"], s=2, c="0.5", lw=0
        )
        bins = pd.qcut(part["log10_mean_counts"], TREND_BINS, duplicates="drop")
        trend = part.groupby(bins, observed=True)[
            ["log10_mean_counts", "log10_variance"]
        ].median()
        ax.plot(trend["log10_mean_counts"], trend["log10_variance"], c="#B2182B")
        ax.set_xlabel("log10 mean count")
        ax.set_ylabel(f"log10 variance ({stage})")
    lims = axes[0].get_xlim()
    axes[0].plot(lims, lims, ls="--", lw=0.5, c="black")
    save_panel(fig, long, figures, "mean_variance", 183, 70)
    plt.close(fig)

    depth = pd.DataFrame(
        {
            "log10_library_size": np.log10(per_cell["library_size"]),
            "pc1": per_cell["pc1"],
        }
    )
    fig = plot_scatter_fit(depth, "log10_library_size", "pc1")
    fig.set_layout_engine("constrained")
    fig.axes[0].set_xlabel("log10 library size (counts)")
    fig.axes[0].set_ylabel(f"PC1 ({method})")
    save_panel(fig, depth, figures, "depth_pc1", 89, 70)
    plt.close(fig)

    top = per_feature.nlargest(N_TOP_SHOWN, "var_normalised")["feature"].tolist()
    values = pd.DataFrame(adata[:, top].X, columns=top).melt(
        var_name="feature", value_name="value"
    )
    fig, ax = plt.subplots(layout="constrained")
    ax.boxplot(
        [values.loc[values["feature"] == f, "value"] for f in top],
        tick_labels=top,
        showfliers=False,
        medianprops={"color": "black"},
    )
    ax.set_xlabel("feature (by normalised variance)")
    ax.set_ylabel(f"normalised value ({method})")
    save_panel(fig, values, figures, "top_variable_values", 89, 70)
    plt.close(fig)

    if method != "dsb":
        return
    counts = dense(adata.layers["counts"]).astype(float)
    step_one = np.log(counts + DSB_PSEUDOCOUNT)
    is_iso = adata.var_names.isin(isotypes)
    # ponytail: the median non-isotype log count stands in for dsb's per-cell
    # GMM background mean, which muon does not return; fine for panels where
    # most proteins are background in most cells.
    iso = pd.DataFrame(
        {
            "isotype_mean_log": step_one[:, is_iso].mean(axis=1),
            "background_median_log": np.median(step_one[:, ~is_iso], axis=1),
        }
    )
    fig = plot_scatter_fit(iso, "isotype_mean_log", "background_median_log")
    fig.set_layout_engine("constrained")
    fig.axes[0].set_xlabel(f"mean isotype log(count + {DSB_PSEUDOCOUNT})")
    fig.axes[0].set_ylabel(f"median protein log(count + {DSB_PSEUDOCOUNT})")
    save_panel(fig, iso, figures, "isotype_vs_background", 89, 70)
    plt.close(fig)


def analyse(argv, argv_run):
    p = parse_args(argv)
    check_values(p)
    seed_everything(p["seed"])
    adata, empty = read_inputs(p)
    outdir = make_outdir(p["outdir"])
    method = p["method"]
    keys = [*REQUIRED, *METHOD_ARGS[method]]
    # One BLAS thread: dsb fits a GMM per cell and the PCA sums in a fixed order.
    with threadpool_limits(1):
        n_dropped = normalise(adata, empty, p)
        per_cell, per_feature = tables(adata, p["seed"])

    data = outdir / "data"
    data.mkdir()
    adata.write_h5ad(data / "normalized.h5ad")
    per_cell.to_csv(outdir / "tables/per_cell_summary.tsv", sep="\t", index=False)
    per_feature.to_csv(outdir / "tables/gene_mean_variance.tsv", sep="\t", index=False)

    rho = abs(stats.spearmanr(per_cell["library_size"], per_cell["pc1"])[0])
    isotypes = (p["isotype_controls"] or "").split(",")
    plt.switch_backend("Agg")
    draw_figures(outdir / "figures", adata, per_cell, per_feature, method, isotypes)

    n_cells, n_features = adata.n_obs, adata.n_vars
    captions = {
        "mean_variance": f"Scatter of log10 variance against log10 mean count for "
        f"{n_features} features, raw counts (left, dashed line: Poisson, variance "
        f"= mean) and {method} values (right); red lines are binned medians.",
        "depth_pc1": f"Scatter of PC1 of the {method} values against log10 library "
        f"size for {n_cells} cells, with a linear fit and 95% CI; |Spearman rho| "
        f"= {rho:.3f}. A strong trend means depth still drives the first axis.",
        "top_variable_values": f"Box plots of {method} values in {n_cells} cells for "
        f"the {N_TOP_SHOWN} features with the highest normalised variance; "
        "whiskers 1.5 IQR, outliers hidden.",
    }
    if method == "dsb":
        captions["isotype_vs_background"] = (
            f"Scatter of each cell's median protein log(count + {DSB_PSEUDOCOUNT}) "
            f"against its mean isotype log count, {n_cells} cells, with a linear "
            "fit and 95% CI; a positive trend is the technical component dsb "
            "regresses out."
        )
    summary = [
        f"{n_cells} cells, {n_features} features normalised by {method} into "
        f"X and layers['{p['layer_out']}']; {n_dropped} all-zero features dropped.",
        f"|Spearman rho| of PC1 with library size: {rho:.3f}.",
    ]
    write_report(outdir, f"Normalise ({method})", p, keys, summary, captions)
    write_result(
        outdir,
        {
            "method": method,
            "n_cells": n_cells,
            "n_features": n_features,
            "n_features_dropped": n_dropped,
            "pc1_depth_abs_spearman": round(float(rho), 4),
            "layer_out": p["layer_out"],
            "seed": p["seed"],
            "params": {k: p[k] for k in [*keys, "outdir"]},
        },
        summary,
        ["report.md", "figures/depth_pc1.png", "figures/mean_variance.png"],
        reference=p["reference"],
    )
    write_reproducibility(outdir, __file__, p, keys, argv_run, PACKAGES)
    return outdir


if __name__ == "__main__":
    raise SystemExit(run("sc-normalize", analyse))
