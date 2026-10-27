"""Purpose: pseudobulk DESeq2 for one cell type, from a cell-level h5ad.
Inputs:  --adata <file.h5ad> with raw integer counts in X or --layer; cell-level obs
Outputs: <outdir>/report.md, result.json, tables/, figures/, reproducibility/
Run:     python scripts/pseudobulk_de.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""

import argparse
import hashlib
import json
import os
import random
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plots import plot_heatmap, plot_volcano
from scipy import sparse
from theme import save_panel

SKILL_DIR = Path(__file__).resolve().parents[1]
R_ENGINE = Path(__file__).resolve().parent / "deseq2.R"

ANALYSIS_ARGS = [
    "adata",
    "cell_type_col",
    "cell_type",
    "sample_col",
    "design",
    "contrast",
    "min_cells",
    "min_frac",
    "padj_max",
    "lfc_min",
    "seed",
    "engine",
]
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "adata": str(SKILL_DIR / "examples" / "demo.h5ad"),
    "cell_type_col": "cell_type",
    "cell_type": "T",
    "sample_col": "sample",
    "design": "~ patient + condition",
    "contrast": "condition,stim,ctrl",
    "min_cells": 10,
    "min_frac": 0.1,
    "padj_max": 0.05,
    "lfc_min": 1.0,
    "seed": 0,
    "engine": "python",
    "layer": None,
    "outdir": "pseudobulk_de_demo",
}

ERRORS = {
    "missing_args": "missing required arguments: {args}; ask the user, or use --demo",
    "demo_conflict": "--demo fixes the analysis arguments; drop: {args}",
    "bad_value": "{arg} must be {rule}, got {value}",
    "bad_design": "--design must be '~ a + b' with obs column names, got {value!r}",
    "bad_contrast": "--contrast must be factor,numerator,denominator, got {value!r}",
    "outdir_exists": "outdir already exists: {path}; pick a new one",
    "missing_file": "input file not found: {path}",
    "missing_layer": "layer not found in the h5ad: {layer}",
    "missing_column": "obs column not found: {col}",
    "cell_type_absent": "cell type {cell_type!r} not found in obs column {col!r}",
    "bad_counts": "counts in {source} must be non-negative integers (raw counts)",
    "contrast_not_in_design": "contrast factor {factor!r} is not a term of --design",
    "level_absent": "contrast level {level!r} not found in obs column {factor!r}",
    "varies_within_sample": "obs column {col!r} varies within a sample; "
    "design terms must be sample-level",
    "too_few_samples": "fewer than 2 samples for {factor}={level!r} "
    "after --min-cells {min_cells}",
    "no_genes": "zero genes pass --min-frac {min_frac}",
    "no_rscript": "Rscript not found; set RSCRIPT or put Rscript on PATH",
    "r_failed": "deseq2.R failed with exit {code}: {tail}",
}

DE_COLUMNS = ["gene", "base_mean", "log2_fc", "lfc_se", "stat", "pvalue", "padj"]
PYDESEQ2_NAMES = {
    "baseMean": "base_mean",
    "log2FoldChange": "log2_fc",
    "lfcSE": "lfc_se",
}
# Fixed plot spec, not analysis thresholds. save_panel saves without a tight
# bbox, so every figure sets a layout engine to keep its labels on the page.
N_PCA_GENES = 500
N_HEATMAP_GENES = 20
HEATMAP_LIMIT = 2
UP, DOWN, OTHER = "#B2182B", "#2166AC", "0.6"


def fail(key, **fields):
    raise ValueError(ERRORS[key].format(**fields))


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--adata")
    parser.add_argument("--layer")
    parser.add_argument("--cell-type-col")
    parser.add_argument("--cell-type")
    parser.add_argument("--sample-col")
    parser.add_argument("--design")
    parser.add_argument("--contrast")
    parser.add_argument("--min-cells", type=int)
    parser.add_argument("--min-frac", type=float)
    parser.add_argument("--padj-max", type=float)
    parser.add_argument("--lfc-min", type=float)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--engine", choices=["python", "r"])
    parser.add_argument("--outdir")
    parser.add_argument("--demo", action="store_true")
    args = vars(parser.parse_args(argv))
    given = [k for k in [*ANALYSIS_ARGS, "layer"] if args[k] is not None]
    flag = "--{}".format
    if args.pop("demo"):
        # --engine and --outdir may be set beside --demo; the analysis may not.
        conflicts = [k for k in given if k != "engine"]
        if conflicts:
            fail(
                "demo_conflict",
                args=" ".join(flag(k.replace("_", "-")) for k in conflicts),
            )
        return {**DEMO, **{k: v for k, v in args.items() if v is not None}}
    missing = [k for k in [*ANALYSIS_ARGS, "outdir"] if args[k] is None]
    if missing:
        fail("missing_args", args=" ".join(flag(k.replace("_", "-")) for k in missing))
    return args


def check_values(p):
    rules = [
        ("--min-cells", p["min_cells"], p["min_cells"] >= 1, "an integer >= 1"),
        ("--min-frac", p["min_frac"], 0 <= p["min_frac"] <= 1, "in [0, 1]"),
        ("--padj-max", p["padj_max"], 0 < p["padj_max"] <= 1, "in (0, 1]"),
        ("--lfc-min", p["lfc_min"], p["lfc_min"] >= 0, ">= 0"),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ]
    for arg, value, ok, rule in rules:
        if not ok:
            fail("bad_value", arg=arg, rule=rule, value=value)
    if not re.fullmatch(
        r"~\s*[A-Za-z_][\w.]*(\s*\+\s*[A-Za-z_][\w.]*)*\s*", p["design"]
    ):
        fail("bad_design", value=p["design"])
    contrast = [part.strip() for part in p["contrast"].split(",")]
    if len(contrast) != 3 or not all(contrast):
        fail("bad_contrast", value=p["contrast"])
    return re.findall(r"[A-Za-z_][\w.]*", p["design"]), contrast


def is_raw_counts(matrix):
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    return bool(np.all(values >= 0) and np.all(values == np.round(values)))


def build_pseudobulk(p, terms, factor, levels):
    """Steps 1-5: validate, subset, sum per sample, drop small samples, filter genes."""
    path = Path(p["adata"])
    if not path.is_file():
        fail("missing_file", path=path)
    adata = ad.read_h5ad(path)
    if p["layer"] is not None and p["layer"] not in adata.layers:
        fail("missing_layer", layer=p["layer"])
    for col in [p["cell_type_col"], p["sample_col"], *terms, factor]:
        if col not in adata.obs.columns:
            fail("missing_column", col=col)
    if factor not in terms:
        fail("contrast_not_in_design", factor=factor)

    in_type = (adata.obs[p["cell_type_col"]].astype(str) == p["cell_type"]).to_numpy()
    if not in_type.any():
        fail("cell_type_absent", cell_type=p["cell_type"], col=p["cell_type_col"])
    obs = adata.obs.loc[in_type, [p["sample_col"], *terms]].astype(str)
    counts = adata.layers[p["layer"]] if p["layer"] else adata.X
    counts = sparse.csr_matrix(counts[in_type])
    if not is_raw_counts(counts):
        fail("bad_counts", source=f"layer {p['layer']}" if p["layer"] else "X")
    for level in levels:
        if level not in set(obs[factor]):
            fail("level_absent", level=level, factor=factor)

    per_sample = obs.groupby(p["sample_col"])[terms].nunique()
    for col in terms:
        if (per_sample[col] > 1).any():
            fail("varies_within_sample", col=col)
    meta = obs.drop_duplicates(p["sample_col"]).set_index(p["sample_col"]).sort_index()
    meta.index.name = "sample"
    meta["n_cells"] = obs[p["sample_col"]].value_counts()

    samples = pd.Categorical(obs[p["sample_col"]], categories=meta.index)
    indicator = sparse.csr_matrix(
        (np.ones(len(obs)), (samples.codes, np.arange(len(obs)))),
        shape=(len(meta), len(obs)),
    )
    summed = pd.DataFrame(
        (indicator @ counts).toarray().astype(np.int64).T,
        index=adata.var_names,
        columns=meta.index,
    )

    kept = meta["n_cells"] >= p["min_cells"]
    dropped = meta.index[~kept].tolist()
    meta = meta[kept]
    for level in levels:
        if (meta[factor] == level).sum() < 2:
            fail(
                "too_few_samples", factor=factor, level=level, min_cells=p["min_cells"]
            )

    # Detection is measured over the cells of the samples kept in step 4.
    kept_cells = samples.isin(meta.index)
    detected = np.asarray((counts[kept_cells] > 0).mean(axis=0)).ravel()
    genes = detected >= p["min_frac"]
    if not genes.any():
        fail("no_genes", min_frac=p["min_frac"])
    return summed.loc[genes, meta.index], meta, dropped, int(adata.n_vars)


def run_python(counts, meta, p, contrast):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    # PyDESeq2's default path draws no random numbers; the seed is set for the record.
    # One CPU: joblib otherwise spawns a worker per host core, slow on shared nodes.
    inference = DefaultInference(n_cpus=1)
    dds = DeseqDataSet(
        counts=counts.T,
        metadata=meta.drop(columns="n_cells"),
        design=p["design"],
        inference=inference,
        quiet=True,
    )
    dds.deseq2()
    # alpha = --padj-max, so independent filtering targets the user's cutoff.
    stats = DeseqStats(
        dds, contrast=contrast, alpha=p["padj_max"], inference=inference, quiet=True
    )
    stats.summary()
    res = stats.results_df.rename(columns=PYDESEQ2_NAMES)
    return res.rename_axis("gene").reset_index()[DE_COLUMNS], []


def run_r(tables, p, contrast):
    rscript = os.environ.get("RSCRIPT") or shutil.which("Rscript")
    if not rscript or not Path(rscript).is_file():
        fail("no_rscript")
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "de.tsv"
        cmd = [
            rscript,
            str(R_ENGINE),
            str(tables / "pseudobulk_counts.tsv"),
            str(tables / "sample_meta.tsv"),
            p["design"],
            *contrast,
            str(p["padj_max"]),
            str(p["seed"]),
            str(out),
        ]
        run = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if run.returncode != 0:
            tail = " | ".join(run.stderr.strip().splitlines()[-5:])
            fail("r_failed", code=run.returncode, tail=tail)
        res = pd.read_csv(out, sep="\t")
    return res[DE_COLUMNS], [f"Rscript {rscript}", *run.stdout.strip().splitlines()]


def log_cpm(counts):
    return np.log2(counts / counts.sum(axis=0) * 1e6 + 1)


def draw_figures(figures, counts, meta, de, p, factor, num, den):
    logcpm = log_cpm(counts)
    level_colour = {lv: OTHER for lv in meta[factor]} | {num: UP, den: DOWN}

    top_var = logcpm.var(axis=1).nlargest(min(N_PCA_GENES, len(logcpm))).index
    centred = logcpm.loc[top_var].T
    centred = centred - centred.mean()
    u, s, _ = np.linalg.svd(centred.to_numpy(), full_matrices=False)
    explained = 100 * s**2 / np.sum(s**2)
    pca = pd.DataFrame(
        {
            "sample": meta.index,
            factor: meta[factor],
            "pc1": u[:, 0] * s[0],
            "pc2": u[:, 1] * s[1],
        }
    )
    fig, ax = plt.subplots(layout="constrained")
    ax.scatter(pca["pc1"], pca["pc2"], s=12, c=pca[factor].map(level_colour))
    for name, x, y in zip(pca["sample"], pca["pc1"], pca["pc2"], strict=True):
        ax.annotate(name, (x, y), xytext=(2, 2), textcoords="offset points", fontsize=5)
    handles = [
        plt.Line2D([], [], marker="o", ls="", c=colour, label=lv)
        for lv, colour in level_colour.items()
    ]
    ax.legend(handles=handles, title=factor, loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_xlabel(f"PC1 ({explained[0]:.0f}% variance)")
    ax.set_ylabel(f"PC2 ({explained[1]:.0f}% variance)")
    save_panel(fig, pca, figures, "sample_pca", 89, 70)
    plt.close(fig)

    tested = de.dropna(subset=["pvalue"])
    fig, ax = plt.subplots(layout="constrained")
    ax.hist(
        tested["pvalue"], bins=np.linspace(0, 1, 21), color=OTHER, edgecolor="white"
    )
    ax.set_xlabel("raw p-value")
    ax.set_ylabel("genes")
    save_panel(fig, tested[["gene", "pvalue"]], figures, "pvalue_hist", 89, 70)
    plt.close(fig)

    fig, ax = plt.subplots(layout="constrained")
    ax.scatter(
        de["base_mean"],
        de["log2_fc"],
        s=4,
        c=np.where(de["significant"], "darkred", OTHER),
    )
    ax.set_xscale("log")
    ax.axhline(0, lw=0.5, c="black")
    for y in (-p["lfc_min"], p["lfc_min"]):
        ax.axhline(y, ls="--", lw=0.5, c="black")
    ax.set_xlabel("mean normalised count")
    ax.set_ylabel(f"log2 fold change ({num} / {den})")
    save_panel(fig, de, figures, "ma", 89, 70)
    plt.close(fig)

    # padj of 0 would be an infinite -log10; floor it at the smallest double.
    volcano = de.dropna(subset=["padj"]).assign(
        padj=lambda d: d["padj"].clip(lower=1e-300)
    )
    fig = plot_volcano(
        volcano, "log2_fc", "padj", padj_max=p["padj_max"], lfc_min=p["lfc_min"]
    )
    fig.set_layout_engine("constrained")
    save_panel(fig, volcano, figures, "volcano", 89, 70)
    plt.close(fig)

    top = de.dropna(subset=["padj"])["gene"].head(N_HEATMAP_GENES).tolist()
    order = meta.assign(rank=meta[factor].map({den: 0, num: 1}).fillna(2))
    order = order.sort_values(["rank", factor]).index.tolist()
    rows = logcpm.loc[top, order]
    z = rows.sub(rows.mean(axis=1), axis=0).div(rows.std(axis=1), axis=0)
    value = "row z-score of log2 CPM"
    heat = (
        z.rename_axis("gene")
        .reset_index()
        .melt(id_vars="gene", var_name="sample", value_name=value)
    )
    heat["gene"] = pd.Categorical(heat["gene"], categories=top)
    heat["sample"] = pd.Categorical(heat["sample"], categories=order)
    fig = plot_heatmap(heat, "gene", "sample", value, limit=HEATMAP_LIMIT)
    # Constrained layout cannot adopt a colorbar added before it; tight can.
    fig.set_layout_engine("tight")
    save_panel(fig, heat, figures, "top_genes_heatmap", 89, 89)
    plt.close(fig)
    return len(top_var), len(top)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_reproducibility(outdir, p, argv, r_versions):
    repro = outdir / "reproducibility"
    repro.mkdir()
    resolved = [sys.executable, str(Path(__file__).resolve())]
    for key in [*ANALYSIS_ARGS, "layer", "outdir"]:
        if p[key] is not None:
            resolved += [f"--{key.replace('_', '-')}", str(p[key])]
    (repro / "commands.sh").write_text(
        "#!/usr/bin/env bash\n"
        f"# argv as run: {shlex.join(argv)}\n"
        f"{shlex.join(resolved)}\n"
    )
    packages = ["anndata", "numpy", "pandas", "scipy", "matplotlib", "pydeseq2"]
    env = [f"python {sys.version.split()[0]}"] + [
        f"{pkg} {version(pkg)}" for pkg in packages
    ]
    (repro / "environment.txt").write_text("\n".join(env + r_versions) + "\n")
    hashed = sorted(outdir.glob("tables/*.tsv")) + sorted(outdir.glob("figures/*.pdf"))
    lines = [f"{sha256(path)}  {path.relative_to(outdir)}" for path in hashed]
    (repro / "checksums.sha256").write_text("\n".join(lines) + "\n")


def write_report(outdir, p, summary, captions):
    rows = "\n".join(
        f"| `--{k.replace('_', '-')}` | `{p[k] if p[k] is not None else 'X (none)'}` |"
        for k in [*ANALYSIS_ARGS, "layer"]
    )
    figures = "\n\n".join(
        f"![{name}](figures/{name}.png)\n\n{caption}"
        for name, caption in captions.items()
    )
    (outdir / "report.md").write_text(
        f"# Pseudobulk DE: {p['cell_type']}, {p['contrast']}\n\n"
        f"| parameter | value |\n|---|---|\n{rows}\n\n"
        + "\n".join(f"- {line}" for line in summary)
        + f"\n\n## Figures\n\n{figures}\n"
    )


def analyse(p, argv):
    terms, contrast = check_values(p)
    factor, num, den = contrast
    outdir = Path(p["outdir"])
    if outdir.exists():
        fail("outdir_exists", path=outdir)
    random.seed(p["seed"])
    np.random.seed(p["seed"])
    counts, meta, dropped, n_genes_in = build_pseudobulk(p, terms, factor, [num, den])

    outdir.mkdir(parents=True)
    tables, figures = outdir / "tables", outdir / "figures"
    tables.mkdir()
    figures.mkdir()
    counts.rename_axis("gene").to_csv(tables / "pseudobulk_counts.tsv", sep="\t")
    meta.drop(columns="n_cells").to_csv(tables / "sample_meta.tsv", sep="\t")

    if p["engine"] == "r":
        de, r_versions = run_r(tables, p, contrast)
    else:
        de, r_versions = run_python(counts, meta, p, contrast)
    de["significant"] = (de["padj"] < p["padj_max"]) & (
        de["log2_fc"].abs() > p["lfc_min"]
    )
    de = de.sort_values(["padj", "gene"], na_position="last").reset_index(drop=True)
    de.to_csv(tables / "de_results.tsv", sep="\t", index=False)

    plt.switch_backend("Agg")
    n_pca, n_heat = draw_figures(figures, counts, meta, de, p, factor, num, den)
    rule = f"padj < {p['padj_max']} and |log2FC| > {p['lfc_min']}"
    captions = {
        "sample_pca": f"PCA of {len(meta)} pseudobulk samples on log2 CPM of the "
        f"{n_pca} "
        f"most variable genes, coloured by {factor}; each point is one sample.",
        "pvalue_hist": "Histogram of raw p-values for "
        f"{de['pvalue'].notna().sum()} tested "
        "genes in 20 bins of width 0.05; a flat floor with a peak near 0 is expected.",
        "ma": f"MA plot of log2 fold change ({num} / {den}) against mean normalised "
        f"count; red genes pass {rule}, dashed lines mark the fold-change cutoff.",
        "volcano": f"Volcano plot of log2 fold change ({num} / {den}) against -log10 "
        f"adjusted p; red genes pass {rule}, dashed lines mark both cutoffs.",
        "top_genes_heatmap": f"Heatmap of row z-scored log2 CPM for the {n_heat} genes "
        f"with the lowest adjusted p, samples ordered by {factor}; colours saturate at "
        f"±{HEATMAP_LIMIT}.",
    }

    n_sig = int(de["significant"].sum())
    top = de.dropna(subset=["padj"]).head(10)
    summary = [
        f"{len(meta)} samples kept ({', '.join(dropped) or 'none'} dropped by "
        f"--min-cells {p['min_cells']}); {len(de)} of {n_genes_in} genes tested.",
        f"{n_sig} genes significant ({num} vs {den}, {rule}), engine {p['engine']}.",
    ]
    write_report(outdir, p, summary, captions)
    result = {
        "n_samples": len(meta),
        "n_genes_tested": len(de),
        "n_significant": n_sig,
        "dropped_samples": dropped,
        "top_genes": top[["gene", "log2_fc", "padj"]].to_dict(orient="records"),
        "params": p,
        "engine": p["engine"],
        "seed": p["seed"],
        "chat_summary_lines": summary,
        "preferred_artifacts": [
            "report.md",
            "tables/de_results.tsv",
            "figures/volcano.png",
        ],
    }
    (outdir / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    write_reproducibility(outdir, p, argv, r_versions)
    return outdir


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    # Fixed PDF CreationDate so checksums match across reruns.
    os.environ.setdefault("SOURCE_DATE_EPOCH", "0")
    try:
        outdir = analyse(parse_args(argv), [sys.argv[0], *argv])
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"pseudobulk-de wrote {outdir / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
