"""Purpose: pseudobulk DE (DESeq2 or dreamlet) for one cell type, from cell-level h5ad.
Inputs:  --adata <file.h5ad> with raw integer counts in X or --layer; cell-level obs
Outputs: <outdir>/report.md, result.json, tables/, figures/, reproducibility/
Run:     python scripts/pseudobulk_de.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""

import argparse
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plots import plot_heatmap, plot_volcano
from scipy import sparse
from skillkit import (
    check_design_rank,
    check_gene_ids,
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
R_ENGINES = {
    "r": Path(__file__).resolve().parent / "deseq2.R",
    "dreamlet": Path(__file__).resolve().parent / "dreamlet.R",
}

ANALYSIS_ARGS = [
    "adata",
    "gene_id_col",
    "cell_type_col",
    "cell_type",
    "sample_col",
    "donor_col",
    "design",
    "contrast",
    "min_cells",
    "min_frac",
    "padj_max",
    "lfc_min",
    "seed",
    "engine",
    "reference",
    "context",
]
KEYS = [*ANALYSIS_ARGS, "layer"]
NAME = r"[A-Za-z_][\w.]*"
RANDOM_TERM = re.compile(rf"\(\s*1\s*\|\s*({NAME})\s*\)")
TERM = rf"(?:{NAME}|{RANDOM_TERM.pattern})"
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "adata": str(SKILL_DIR / "examples" / "demo.h5ad"),
    "gene_id_col": "gene_ids",
    "cell_type_col": "cell_type",
    "cell_type": "T",
    "sample_col": "sample",
    "donor_col": "patient",
    "design": "~ patient + condition",
    "contrast": "condition,stim,ctrl",
    "min_cells": 10,
    "min_frac": 0.1,
    "padj_max": 0.05,
    "lfc_min": 1.0,
    "seed": 0,
    "engine": "python",
    "reference": "synthetic demo, no genome",
    "context": "synthetic demo: 7 patients, ctrl vs stim T cells, no real tissue",
    "layer": None,
    "outdir": "pseudobulk_de_demo",
}
# `--demo --engine dreamlet` swaps in this demo-only design: same pairing, as a
# random effect (dreamlet's repeated-measures form).
DEMO_DREAMLET_DESIGN = "~ condition + (1|patient)"

ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "bad_design": "--design must be '~ a + b' with obs column names, "
    "or '~ a + (1|b)' with a random effect, got {value!r}",
    "random_needs_dreamlet": "--design has the random-effect term (1|{term}): "
    "repeated measures need --engine dreamlet, not --engine {engine}",
    "bad_contrast": "--contrast must be factor,numerator,denominator, got {value!r}",
    "missing_file": "input file not found: {path}",
    "missing_layer": "layer not found in the h5ad: {layer}",
    "missing_column": "obs column not found: {col}",
    "missing_gene_id_col": "gene id column not found in var: {col}",
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
    "r_failed": "{script} failed with exit {code}: {tail}",
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
MA_X = ("mean normalised count", "mean log2 CPM (dreamlet AveExpr)")


def parse_args(argv):
    """Parse argv. A random-effect design selects dreamlet when --engine is absent."""
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--design")
    random_terms = RANDOM_TERM.findall(pre.parse_known_args(argv)[0].design or "")
    required = [k for k in ANALYSIS_ARGS if k != "engine" or not random_terms]
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--adata")
    parser.add_argument("--gene-id-col")
    parser.add_argument("--layer")
    parser.add_argument("--cell-type-col")
    parser.add_argument("--cell-type")
    parser.add_argument("--sample-col")
    parser.add_argument("--donor-col")
    parser.add_argument("--design")
    parser.add_argument("--contrast")
    parser.add_argument("--min-cells", type=int)
    parser.add_argument("--min-frac", type=float)
    parser.add_argument("--padj-max", type=float)
    parser.add_argument("--lfc-min", type=float)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--engine", choices=["python", "r", "dreamlet"])
    parser.add_argument("--reference")
    parser.add_argument("--context")
    # --engine and --outdir may be set beside --demo; the analysis may not.
    p = parse_required(
        parser,
        required,
        DEMO,
        argv,
        optional=["layer", "engine"],
        demo_extra=("outdir", "engine"),
    )
    if "--demo" in argv and p["engine"] == "dreamlet":
        p["design"] = DEMO_DREAMLET_DESIGN
    random_terms = RANDOM_TERM.findall(p["design"])
    if random_terms and p["engine"] in ("python", "r"):
        fail(ERRORS, "random_needs_dreamlet", term=random_terms[0], engine=p["engine"])
    # Only a random-effect design reaches here without --engine (see `required`).
    p["engine"] = p["engine"] or "dreamlet"
    return p


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
            fail(ERRORS, "bad_value", arg=arg, rule=rule, value=value)
    if not re.fullmatch(rf"~\s*{TERM}(\s*\+\s*{TERM})*\s*", p["design"]):
        fail(ERRORS, "bad_design", value=p["design"])
    contrast = [part.strip() for part in p["contrast"].split(",")]
    if len(contrast) != 3 or not all(contrast):
        fail(ERRORS, "bad_contrast", value=p["contrast"])
    # Every term, fixed or random, is an obs column; the "1" of (1|x) is not matched.
    return re.findall(NAME, p["design"]), contrast


def is_raw_counts(matrix):
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    return bool(np.all(values >= 0) and np.all(values == np.round(values)))


def build_pseudobulk(p, terms, factor, levels):
    """Steps 1-5: validate, subset, sum per sample, drop small samples, filter genes."""
    path = Path(p["adata"])
    if not path.is_file():
        fail(ERRORS, "missing_file", path=path)
    adata = ad.read_h5ad(path)
    # var_names are the symbols; --gene-id-col none: the source has no gene ids.
    ids = None
    if p["gene_id_col"] != "none":
        if p["gene_id_col"] not in adata.var.columns:
            fail(ERRORS, "missing_gene_id_col", col=p["gene_id_col"])
        ids = adata.var[p["gene_id_col"]].astype(str).tolist()
    gene_ids = check_gene_ids(ids, adata.var_names.tolist())
    check_unique(adata.obs_names, "cell ids")
    if p["layer"] is not None and p["layer"] not in adata.layers:
        fail(ERRORS, "missing_layer", layer=p["layer"])
    # `--cell-type-col none`: the input is one population; test all cells, no subset.
    all_cells = p["cell_type_col"] == "none"
    type_col = [] if all_cells else [p["cell_type_col"]]
    for col in [*type_col, p["sample_col"], p["donor_col"], *terms, factor]:
        if col not in adata.obs.columns:
            fail(ERRORS, "missing_column", col=col)
    if factor not in terms or factor in RANDOM_TERM.findall(p["design"]):
        fail(ERRORS, "contrast_not_in_design", factor=factor)

    if all_cells:
        in_type = np.ones(adata.n_obs, dtype=bool)
    else:
        in_type = (
            adata.obs[p["cell_type_col"]].astype(str) == p["cell_type"]
        ).to_numpy()
    if not in_type.any():
        fail(
            ERRORS, "cell_type_absent", cell_type=p["cell_type"], col=p["cell_type_col"]
        )
    cols = list(dict.fromkeys([p["sample_col"], p["donor_col"], *terms]))
    obs = adata.obs.loc[in_type, cols].astype(str)
    counts = adata.layers[p["layer"]] if p["layer"] else adata.X
    counts = sparse.csr_matrix(counts[in_type])
    if not is_raw_counts(counts):
        fail(ERRORS, "bad_counts", source=f"layer {p['layer']}" if p["layer"] else "X")
    for level in levels:
        if level not in set(obs[factor]):
            fail(ERRORS, "level_absent", level=level, factor=factor)

    sample_level = [c for c in cols if c != p["sample_col"]]
    per_sample = obs.groupby(p["sample_col"])[sample_level].nunique()
    for col in sample_level:
        if (per_sample[col] > 1).any():
            fail(ERRORS, "varies_within_sample", col=col)
    meta = obs.drop_duplicates(p["sample_col"]).set_index(p["sample_col"]).sort_index()
    meta = meta[list(dict.fromkeys(terms))]
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
    dropped_by_level = {
        lv: sorted(meta.index[~kept & (meta[factor] == lv)]) for lv in levels
    }
    meta = meta[kept]
    for level in levels:
        if (meta[factor] == level).sum() < 2:
            fail(
                ERRORS,
                "too_few_samples",
                factor=factor,
                level=level,
                min_cells=p["min_cells"],
            )

    # Design checks run on the kept samples, the rows the engine fits.
    kept_cells = samples.isin(meta.index)
    random = RANDOM_TERM.findall(p["design"])
    fixed = list(dict.fromkeys(t for t in terms if t not in random))
    # obs was cast to str; a numeric covariate is one column, not one dummy per value.
    numeric = [c for c in fixed if pd.api.types.is_numeric_dtype(adata.obs[c])]
    check_design_rank(meta[fixed].astype(dict.fromkeys(numeric, float)), fixed)
    n_donors, _ = check_replication(
        obs[kept_cells], p["sample_col"], p["donor_col"], terms
    )

    # Detection is measured over the cells of the samples kept in step 4.
    detected = np.asarray((counts[kept_cells] > 0).mean(axis=0)).ravel()
    genes = detected >= p["min_frac"]
    if not genes.any():
        fail(ERRORS, "no_genes", min_frac=p["min_frac"])
    # Why an untested gene has no result: never detected, or under --min-frac.
    reasons = {
        g: "not_detected" if d == 0 else "not_tested"
        for g, d in zip(adata.var_names, detected, strict=True)
    }
    drops = {
        "dropped": dropped,
        "by_level": dropped_by_level,
        "n_donors": n_donors,
        "reasons": reasons,
        "gene_ids": gene_ids,
    }
    return summed.loc[genes, meta.index], meta, drops, int(adata.n_vars)


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


def run_r(counts, tables, meta, p, contrast):
    """deseq2.R or dreamlet.R; dreamlet also gets n_cells per sample."""
    rscript = os.environ.get("RSCRIPT") or shutil.which("Rscript")
    if not rscript or not Path(rscript).is_file():
        fail(ERRORS, "no_rscript")
    engine = R_ENGINES[p["engine"]]
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "de.tsv"
        # The engines read a symbol-keyed matrix; the table copy also has gene_id.
        counts_file = Path(tmp) / "pseudobulk_counts.tsv"
        counts.rename_axis("gene").to_csv(counts_file, sep="\t")
        meta_file = tables / "sample_meta.tsv"
        if p["engine"] == "dreamlet":
            # n_cells stays out of the checksummed sample_meta.tsv.
            meta_file = Path(tmp) / "sample_meta.tsv"
            meta.to_csv(meta_file, sep="\t")
        cmd = [
            rscript,
            str(engine),
            str(counts_file),
            str(meta_file),
            p["design"],
            *contrast,
            str(p["padj_max"]),
            str(p["seed"]),
            str(out),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            tail = " | ".join(proc.stderr.strip().splitlines()[-5:])
            fail(
                ERRORS, "r_failed", script=engine.name, code=proc.returncode, tail=tail
            )
        res = pd.read_csv(out, sep="\t")
    return res[DE_COLUMNS], [f"Rscript {rscript}", *proc.stdout.strip().splitlines()]


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
    # dreamlet's base_mean is AveExpr, mean log2 CPM: already a log scale.
    if p["engine"] != "dreamlet":
        ax.set_xscale("log")
    ax.axhline(0, lw=0.5, c="black")
    for y in (-p["lfc_min"], p["lfc_min"]):
        ax.axhline(y, ls="--", lw=0.5, c="black")
    ax.set_xlabel(MA_X[p["engine"] == "dreamlet"])
    ax.set_ylabel(f"log2 fold change ({num} / {den})")
    save_panel(fig, de, figures, "ma", 89, 70)
    plt.close(fig)

    # padj of 0 would be an infinite -log10; floor it at the smallest double.
    volcano = de.dropna(subset=["padj"]).assign(
        padj=lambda d: d["padj"].clip(lower=1e-300)
    )
    plot = plot_volcano(
        volcano, "log2_fc", "padj", padj_max=p["padj_max"], lfc_min=p["lfc_min"]
    )
    save_panel(plot, volcano, figures, "volcano", 89, 70)

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
    plot = plot_heatmap(heat, "gene", "sample", value, limit=HEATMAP_LIMIT)
    save_panel(plot, heat, figures, "top_genes_heatmap", 89, 89)
    return len(top_var), len(top)


def analyse(argv, full_argv):
    p = parse_args(argv)
    terms, contrast = check_values(p)
    factor, num, den = contrast
    seed_everything(p["seed"])
    counts, meta, drops, n_genes_in = build_pseudobulk(p, terms, factor, [num, den])
    dropped = drops["dropped"]

    outdir = make_outdir(p["outdir"])
    tables, figures = outdir / "tables", outdir / "figures"
    # gene_id, gene_symbol first: variance-partition reads this table as --counts.
    counts.set_index(
        pd.MultiIndex.from_arrays(
            [counts.index.map(drops["gene_ids"]), counts.index],
            names=["gene_id", "gene_symbol"],
        )
    ).to_csv(tables / "pseudobulk_counts.tsv", sep="\t")
    meta.drop(columns="n_cells").to_csv(tables / "sample_meta.tsv", sep="\t")

    if p["engine"] == "python":
        de, r_versions = run_python(counts, meta, p, contrast)
    else:
        de, r_versions = run_r(counts, tables, meta, p, contrast)
    # A gene the engine drops (dreamlet: all-zero or failed fit) stays as an NA row,
    # tested; a gene filtered before the engine keeps NA with its reason.
    de = de.set_index("gene").reindex(counts.index).rename_axis("gene").reset_index()
    reasons = drops["reasons"]
    de = with_status(de, list(reasons), "gene", reasons)
    de["significant"] = (de["padj"] < p["padj_max"]) & (
        de["log2_fc"].abs() > p["lfc_min"]
    )
    de = de.sort_values(["padj", "gene"], na_position="last").reset_index(drop=True)
    de.insert(1, "gene_id", de["gene"].map(drops["gene_ids"]))
    de.rename(columns={"gene": "gene_symbol"}).to_csv(
        tables / "de_results.tsv", sep="\t", index=False
    )
    n_tested = int((de["status"] == "tested").sum())
    tested = de[de["status"] == "tested"]

    plt.switch_backend("Agg")
    n_pca, n_heat = draw_figures(figures, counts, meta, tested, p, factor, num, den)
    rule = f"padj < {p['padj_max']} and |log2FC| > {p['lfc_min']}"
    captions = {
        "sample_pca": f"PCA of {len(meta)} pseudobulk samples on log2 CPM of the "
        f"{n_pca} "
        f"most variable genes, coloured by {factor}; each point is one sample.",
        "pvalue_hist": "Histogram of raw p-values for "
        f"{tested['pvalue'].notna().sum()} tested "
        "genes in 20 bins of width 0.05; a flat floor with a peak near 0 is expected.",
        "ma": f"MA plot of log2 fold change ({num} / {den}) against "
        f"{MA_X[p['engine'] == 'dreamlet']}; red genes pass {rule}, dashed lines "
        "mark the fold-change cutoff.",
        "volcano": f"Volcano plot of log2 fold change ({num} / {den}) against -log10 "
        f"adjusted p; red genes pass {rule} and go up, blue pass it and go down, grey "
        "do not; dashed lines mark both cutoffs.",
        "top_genes_heatmap": f"Heatmap of row z-scored log2 CPM for the {n_heat} genes "
        "with the lowest adjusted p; not clustered: rows fixed by adjusted-p rank, "
        f"columns fixed by {factor} label; colours saturate at ±{HEATMAP_LIMIT}.",
    }

    n_sig = int(de["significant"].sum())
    top = de.dropna(subset=["padj"]).head(10)
    by_level = "; ".join(
        f"{lv}: {', '.join(s) or 'none'}" for lv, s in drops["by_level"].items()
    )
    summary = [
        f"{len(meta)} samples from {drops['n_donors']} donors kept; dropped by "
        f"--min-cells {p['min_cells']}: {by_level}. {n_tested} of {n_genes_in} genes "
        "tested; filtered genes are not tested, not absent.",
        f"{n_sig} genes significant ({num} vs {den}, {rule}), engine {p['engine']}.",
    ]
    title = f"Pseudobulk DE: {p['cell_type']}, {p['contrast']}"
    # A bare "|" in (1|x) would split the markdown table cell.
    shown = {
        **p,
        "layer": p["layer"] or "X (none)",
        "design": p["design"].replace("|", "\\|"),
    }
    write_report(outdir, title, shown, KEYS, summary, captions)
    result = {
        "n_samples": len(meta),
        "n_donors": drops["n_donors"],
        "n_genes_tested": n_tested,
        "status_counts": {k: int(v) for k, v in de["status"].value_counts().items()},
        "gene_id_col": p["gene_id_col"],
        "n_significant": n_sig,
        "dropped_samples": dropped,
        "dropped_by_level": drops["by_level"],
        "top_genes": top[["gene", "log2_fc", "padj"]].to_dict(orient="records"),
        "params": p,
        "engine": p["engine"],
        "seed": p["seed"],
    }
    artifacts = ["report.md", "tables/de_results.tsv", "figures/volcano.png"]
    write_result(
        outdir,
        result,
        summary,
        artifacts,
        context=p["context"],
        reference=p["reference"],
    )
    packages = [
        "anndata",
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "plotnine",
        "pydeseq2",
    ]
    write_reproducibility(
        outdir, __file__, p, KEYS, full_argv, packages, extra=r_versions
    )
    return outdir


if __name__ == "__main__":
    raise SystemExit(run("pseudobulk-de", analyse))
