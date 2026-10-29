"""Purpose: variance partitioning of pseudobulk expression across covariates.
Inputs:  --counts/--meta pseudobulk TSVs, or --adata cells summed per --sample-col
Outputs: <outdir>/report.md, result.json, tables/, figures/, reproducibility/
Run:     python scripts/variance_partition.py --demo --outdir demo_out
Notes:   variancePartition runs in R (varpart.R via Rscript). Every analysis argument
         is required; --demo alone pins the demo values. Steps run in a fixed order;
         see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""

import argparse
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from palettes import categorical_palette, palette_for
from plots import plot_heatmap
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
ENGINE = Path(__file__).resolve().parent / "varpart.R"

MODEL_ARGS = [
    "formula",
    "min_frac",
    "n_top",
    "max_cancor",
    "reference",
    "context",
    "seed",
]
COUNTS_ARGS = ["counts", "meta"]
ADATA_ARGS = ["adata", "gene_id_col", "sample_col", "donor_col", "min_cells"]
KEYS = [*COUNTS_ARGS, *ADATA_ARGS, "layer", *MODEL_ARGS]
NAME = r"[A-Za-z_][\w.]*"
RANDOM_TERM = re.compile(rf"\(\s*1\s*\|\s*({NAME})\s*\)")
TERM = rf"(?:{NAME}|{RANDOM_TERM.pattern})"
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "counts": str(SKILL_DIR / "examples" / "demo_counts.tsv"),
    "meta": str(SKILL_DIR / "examples" / "demo_meta.tsv"),
    "adata": None,
    "gene_id_col": None,
    "sample_col": None,
    "donor_col": None,
    "min_cells": None,
    "layer": None,
    "formula": "~ (1|donor) + (1|condition) + rin",
    "min_frac": 0.1,
    "n_top": 5,
    "max_cancor": 0.9,
    "reference": "synthetic demo genes, no genome build",
    "context": "synthetic demo: 12 donors, ctrl vs stim pseudobulk, no real tissue",
    "seed": 0,
    "outdir": "variance_partition_demo",
}

ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "bad_formula": "--formula must be '~ (1|a) + b' with column names, got {value!r}",
    "too_many_terms": "--formula has {n} terms; the figures colour at most 23",
    "input_mode": "give --counts with --meta, or --adata with --sample-col, not both",
    "min_cells_needs_adata": "--min-cells applies only to --adata input; "
    "--counts samples are already summed",
    "missing_file": "input file not found: {path}",
    "missing_layer": "layer not found in the h5ad: {layer}",
    "missing_gene_id_col": "gene id column not found in var: {col}",
    "missing_id_columns": "--counts must start with gene_id and gene_symbol columns, "
    "then one column per sample",
    "missing_column": "column not found in the sample metadata: {col}",
    "bad_counts": "counts in {source} must be non-negative integers (raw counts)",
    "samples_mismatch": "--counts columns and --meta rows name different samples: "
    "{only}",
    "varies_within_sample": "obs column {col!r} varies within a sample; "
    "formula terms must be sample-level",
    "categorical_fixed": "covariate {col!r} is categorical; write it as (1|{col})",
    "too_few_samples": "{n} samples left after --min-cells {min_cells}; "
    "need at least 3",
    "too_few_count_samples": "{n} samples in --counts; need at least 3",
    "no_genes": "zero genes pass --min-frac {min_frac}",
    "collinear": "canonical correlation {value:.3f} between {a} and {b} exceeds "
    "--max-cancor {max_cancor}; all pairs: {pairs}",
    "no_rscript": "Rscript not found; set RSCRIPT or put Rscript on PATH",
    "r_failed": "{script} failed with exit {code}: {tail}",
}

# Fixed plot spec, not analysis thresholds. Terms take categorical_palette
# (colour-blind safe; the corpus used RColorBrewer Set1, E-pl-01), grey85 for
# Residuals. 23 is the largest colour-blind-safe set.
MAX_TERMS = 23
RESIDUAL_COLOUR = "#D9D9D9"


def parse_args(argv):
    """Parse argv. --adata switches the required input arguments to the h5ad path."""
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--adata")
    pre.add_argument("--counts")
    pre.add_argument("--min-cells")
    given = pre.parse_known_args(argv)[0]
    if given.adata and given.counts:
        fail(ERRORS, "input_mode")
    if given.counts and given.min_cells:
        fail(ERRORS, "min_cells_needs_adata")
    inputs = ADATA_ARGS if given.adata else COUNTS_ARGS
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--counts")
    parser.add_argument("--meta")
    parser.add_argument("--adata")
    parser.add_argument("--gene-id-col")
    parser.add_argument("--sample-col")
    parser.add_argument("--donor-col")
    parser.add_argument("--layer")
    parser.add_argument("--min-cells", type=int)
    parser.add_argument("--formula")
    parser.add_argument("--min-frac", type=float)
    parser.add_argument("--n-top", type=int)
    parser.add_argument("--max-cancor", type=float)
    parser.add_argument("--reference")
    parser.add_argument("--context")
    parser.add_argument("--seed", type=int)
    optional = [k for k in [*COUNTS_ARGS, *ADATA_ARGS, "layer"] if k not in inputs]
    return parse_required(parser, [*inputs, *MODEL_ARGS], DEMO, argv, optional)


def check_values(p):
    rules = [
        ("--min-frac", p["min_frac"], 0 <= p["min_frac"] <= 1, "in [0, 1]"),
        ("--n-top", p["n_top"], p["n_top"] >= 1, "an integer >= 1"),
        ("--max-cancor", p["max_cancor"], 0 < p["max_cancor"] <= 1, "in (0, 1]"),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ]
    if p["min_cells"] is not None:
        rules.append(
            ("--min-cells", p["min_cells"], p["min_cells"] >= 1, "an integer >= 1")
        )
    for arg, value, ok, rule in rules:
        if not ok:
            fail(ERRORS, "bad_value", arg=arg, rule=rule, value=value)
    if not re.fullmatch(rf"~\s*{TERM}(\s*\+\s*{TERM})*\s*", p["formula"]):
        fail(ERRORS, "bad_formula", value=p["formula"])
    terms = re.findall(NAME, p["formula"])
    if len(terms) > MAX_TERMS:
        fail(ERRORS, "too_many_terms", n=len(terms))
    return terms, set(RANDOM_TERM.findall(p["formula"]))


def check_counts(values, source):
    if not (np.all(values >= 0) and np.all(values == np.round(values))):
        fail(ERRORS, "bad_counts", source=source)


def read_tsv_input(p, terms):
    """--counts/--meta: summed counts already; detection is over samples."""
    for path in (p["counts"], p["meta"]):
        if not Path(path).is_file():
            fail(ERRORS, "missing_file", path=path)
    # pandas renames a repeated header ("S1.1"), so read the header as written.
    with open(p["counts"]) as handle:
        header = handle.readline().rstrip("\n").split("\t")
    if header[:2] != ["gene_id", "gene_symbol"]:
        fail(ERRORS, "missing_id_columns")
    check_unique(header[2:], "sample ids")
    counts = pd.read_csv(p["counts"], sep="\t", dtype={"gene_id": str})
    # An empty gene_id column (pseudobulk-de --gene-id-col none): symbols only.
    ids = None if counts["gene_id"].isna().all() else counts["gene_id"].tolist()
    gene_ids = check_gene_ids(ids, counts["gene_symbol"].astype(str).tolist())
    counts = counts.drop(columns="gene_id").set_index("gene_symbol").rename_axis("gene")
    counts.index = counts.index.astype(str)
    meta = pd.read_csv(p["meta"], sep="\t", index_col=0)
    check_unique(meta.index.astype(str), "sample ids")
    check_counts(counts.to_numpy(), p["counts"])
    only = sorted(set(counts.columns) ^ set(meta.index.astype(str)))
    if only:
        fail(ERRORS, "samples_mismatch", only=", ".join(only))
    if counts.shape[1] < 3:
        fail(ERRORS, "too_few_count_samples", n=counts.shape[1])
    for col in terms:
        if col not in meta.columns:
            fail(ERRORS, "missing_column", col=col)
    meta.index = meta.index.astype(str)
    meta = meta.loc[sorted(counts.columns), terms].rename_axis("sample")
    counts = counts[meta.index].astype(np.int64)
    detected = (counts > 0).mean(axis=1)
    return counts, meta, [], detected, None, gene_ids


def read_adata_input(p, terms):
    """--adata: sum per sample as in pseudobulk-de; detection over kept cells."""
    import anndata as ad
    from scipy import sparse

    path = Path(p["adata"])
    if not path.is_file():
        fail(ERRORS, "missing_file", path=path)
    adata = ad.read_h5ad(path)
    if p["layer"] is not None and p["layer"] not in adata.layers:
        fail(ERRORS, "missing_layer", layer=p["layer"])
    for col in [p["sample_col"], p["donor_col"], *terms]:
        if col not in adata.obs.columns:
            fail(ERRORS, "missing_column", col=col)
    # var_names are the symbols; --gene-id-col none: the source has no gene ids.
    ids = None
    if p["gene_id_col"] != "none":
        if p["gene_id_col"] not in adata.var.columns:
            fail(ERRORS, "missing_gene_id_col", col=p["gene_id_col"])
        ids = adata.var[p["gene_id_col"]].astype(str).tolist()
    gene_ids = check_gene_ids(ids, adata.var_names.tolist())
    check_unique(adata.obs_names, "cell ids")
    matrix = sparse.csr_matrix(adata.layers[p["layer"]] if p["layer"] else adata.X)
    check_counts(matrix.data, f"layer {p['layer']}" if p["layer"] else "X")

    sample_level = list(dict.fromkeys([*terms, p["donor_col"]]))
    obs = adata.obs[[p["sample_col"], *sample_level]].copy()
    obs[p["sample_col"]] = obs[p["sample_col"]].astype(str)
    per_sample = obs.groupby(p["sample_col"], observed=True)[sample_level].nunique()
    for col in sample_level:
        if (per_sample[col] > 1).any():
            fail(ERRORS, "varies_within_sample", col=col)
    meta = obs.drop_duplicates(p["sample_col"]).set_index(p["sample_col"]).sort_index()
    meta.index.name = "sample"
    n_cells = obs[p["sample_col"]].value_counts()

    samples = pd.Categorical(obs[p["sample_col"]], categories=meta.index)
    indicator = sparse.csr_matrix(
        (np.ones(len(obs)), (samples.codes, np.arange(len(obs)))),
        shape=(len(meta), len(obs)),
    )
    counts = pd.DataFrame(
        (indicator @ matrix).toarray().astype(np.int64).T,
        index=adata.var_names.rename("gene"),
        columns=meta.index,
    )
    kept = n_cells.reindex(meta.index) >= p["min_cells"]
    dropped = meta.index[~kept].tolist()
    meta = meta[kept.to_numpy()]
    if len(meta) < 3:
        fail(ERRORS, "too_few_samples", n=len(meta), min_cells=p["min_cells"])
    # Donor versus sample: several samples per donor need a donor term.
    n_donors, _ = check_replication(meta.reset_index(), "sample", p["donor_col"], terms)
    meta = meta[terms]
    in_kept = samples.isin(meta.index)
    detected = pd.Series(
        np.asarray((matrix[in_kept] > 0).mean(axis=0)).ravel(), index=counts.index
    )
    return counts[meta.index], meta, dropped, detected, n_donors, gene_ids


def run_r(counts, meta, p):
    rscript = os.environ.get("RSCRIPT") or shutil.which("Rscript")
    if not rscript or not Path(rscript).is_file():
        fail(ERRORS, "no_rscript")
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        counts.to_csv(tmp / "counts.tsv", sep="\t")
        meta.to_csv(tmp / "meta.tsv", sep="\t")
        cmd = [
            rscript,
            str(ENGINE),
            str(tmp / "counts.tsv"),
            str(tmp / "meta.tsv"),
            p["formula"],
            str(p["max_cancor"]),
            str(p["seed"]),
            str(tmp),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            tail = " | ".join(proc.stderr.strip().splitlines()[-5:])
            fail(
                ERRORS, "r_failed", script=ENGINE.name, code=proc.returncode, tail=tail
            )
        cancor = pd.read_csv(tmp / "cancor.tsv", sep="\t", index_col=0)
        fractions_file = tmp / "variance_fractions.tsv"
        fractions = (
            pd.read_csv(fractions_file, sep="\t") if fractions_file.is_file() else None
        )
    return cancor, fractions, [f"Rscript {rscript}", *proc.stdout.strip().splitlines()]


def draw_figures(figures, fractions, order, top_genes, cancor, pal):
    long = fractions.melt(id_vars="gene", var_name="term", value_name="fraction")
    long["percent"] = 100 * long["fraction"]
    long["term"] = pd.Categorical(long["term"], categories=order)
    long = long.sort_values(["term", "gene"]).reset_index(drop=True)

    fig, ax = plt.subplots(layout="constrained")
    groups = [long.loc[long["term"] == t, "percent"].to_numpy() for t in order]
    violins = ax.violinplot(groups, positions=range(len(order)), showextrema=False)
    for body, colour in zip(violins["bodies"], palette_for(order, pal), strict=True):
        body.set_facecolor(colour)
        body.set_edgecolor("black")
        body.set_linewidth(0.3)
        body.set_alpha(1)
    ax.boxplot(
        groups,
        positions=range(len(order)),
        widths=0.1,
        showfliers=False,
        medianprops={"color": "black"},
    )
    ax.set_xticks(range(len(order)), order, rotation=45, ha="right")
    ax.set_ylim(0, 100)
    ax.set_ylabel("variance explained (%)")
    save_panel(fig, long, figures, "variance_violin", 89, 70)
    plt.close(fig)

    bars = long[long["gene"].isin(top_genes)].copy()
    bars["gene"] = pd.Categorical(bars["gene"], categories=top_genes)
    bars = bars.sort_values(["gene", "term"]).reset_index(drop=True)
    wide = bars.pivot(index="gene", columns="term", values="percent").loc[
        top_genes, order
    ]
    fig, ax = plt.subplots(layout="constrained")
    left = np.zeros(len(wide))
    for term, colour in zip(order, palette_for(order, pal), strict=True):
        ax.barh(range(len(wide)), wide[term], left=left, color=colour, label=term)
        left += wide[term].to_numpy()
    ax.set_yticks(range(len(wide)), top_genes)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("variance explained (%)")
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    height = max(70, 4 * len(top_genes) + 20)
    save_panel(fig, bars, figures, "percent_bars", 89, height)
    plt.close(fig)

    terms = list(cancor.index)
    value = "canonical correlation"
    heat = (
        cancor.rename_axis("term_row")
        .reset_index()
        .melt(id_vars="term_row", var_name="term_col", value_name=value)
    )
    heat["term_row"] = pd.Categorical(heat["term_row"], categories=terms)
    heat["term_col"] = pd.Categorical(heat["term_col"], categories=terms)
    heat = heat.sort_values(["term_row", "term_col"]).reset_index(drop=True)
    # canCor lies in [0, 1]: white at 0, red at 1 on the diverging scale.
    plot = plot_heatmap(heat, "term_row", "term_col", value, limit=1)
    save_panel(plot, heat, figures, "cancor_heatmap", 89, 70)


def analyse(argv, full_argv):
    p = parse_args(argv)
    terms, random_terms = check_values(p)
    seed_everything(p["seed"])
    read = read_adata_input if p["adata"] else read_tsv_input
    counts, meta, dropped, detected, n_donors, gene_ids = read(p, terms)
    for col in terms:
        if col not in random_terms and not pd.api.types.is_numeric_dtype(meta[col]):
            fail(ERRORS, "categorical_fixed", col=col)
    fixed = [t for t in terms if t not in random_terms]
    if fixed:
        check_design_rank(meta, fixed)
    genes = detected[detected >= p["min_frac"]].index
    if genes.empty:
        fail(ERRORS, "no_genes", min_frac=p["min_frac"])
    n_genes_in = len(counts)
    counts = counts.loc[sorted(genes)]

    if Path(p["outdir"]).exists():
        make_outdir(p["outdir"])  # fails fast, before the slow R run
    # Every refusal comes before make_outdir, so a refused run leaves no outdir.
    cancor, fractions, r_versions = run_r(counts, meta, p)
    cancor = cancor.loc[terms, terms]
    pairs = (
        cancor.where(np.triu(np.ones(cancor.shape, dtype=bool), k=1)).stack().dropna()
    )
    if len(pairs) and pairs.max() > p["max_cancor"]:
        a, b = pairs.idxmax()
        listed = ", ".join(f"{x}/{y} {v:.3f}" for (x, y), v in pairs.items())
        fail(
            ERRORS,
            "collinear",
            value=pairs.max(),
            a=a,
            b=b,
            max_cancor=p["max_cancor"],
            pairs=listed,
        )
    outdir = make_outdir(p["outdir"])
    tables, figures = outdir / "tables", outdir / "figures"
    cancor.rename_axis("term").to_csv(tables / "cancor.tsv", sep="\t")

    fractions = fractions[["gene", *terms, "Residuals"]].sort_values("gene")
    # Every input gene keeps a row; a gene under --min-frac has NA fractions.
    reasons = {
        g: "not_detected" if d == 0 else "not_tested" for g, d in detected.items()
    }
    table = with_status(fractions, list(reasons), "gene", reasons).sort_values("gene")
    table.insert(1, "gene_id", table["gene"].map(gene_ids))
    table.rename(columns={"gene": "gene_symbol"}).to_csv(
        tables / "variance_fractions.tsv", sep="\t", index=False
    )
    long = fractions.melt(id_vars="gene", var_name="term", value_name="fraction")
    summary_table = (
        long.groupby("term")["fraction"]
        .agg(median_fraction="median", mean_fraction="mean")
        .reset_index()
        .sort_values(["median_fraction", "term"], ascending=[False, True])
    )
    summary_table.to_csv(tables / "summary.tsv", sep="\t", index=False)

    # sortCols order (E-pl-01): covariates by median fraction, Residuals last.
    order = [t for t in summary_table["term"] if t != "Residuals"] + ["Residuals"]
    top_covariate = order[0]
    top_by_term = {
        term: fractions.sort_values([term, "gene"], ascending=[False, True])["gene"]
        .head(p["n_top"])
        .tolist()
        for term in order
    }
    top_genes = list(dict.fromkeys(g for genes in top_by_term.values() for g in genes))
    pal = categorical_palette(order[:-1]) | {"Residuals": RESIDUAL_COLOUR}

    plt.switch_backend("Agg")
    draw_figures(figures, fractions, order, top_genes, cancor, pal)
    medians = summary_table.set_index("term")["median_fraction"]
    captions = {
        "variance_violin": f"Violin plots of the percent of expression variance "
        f"explained by each formula term across {len(fractions)} genes, terms ordered "
        "by median, Residuals last; inner boxes mark the quartiles.",
        "percent_bars": f"Stacked bars of the percent variance explained for the "
        f"{p['n_top']} genes with the largest fraction for each term, "
        f"{len(top_genes)} genes in all, grouped by the term they top.",
        "cancor_heatmap": f"Heatmap of canonical correlations between the "
        f"{len(terms)} formula terms over {len(meta)} samples, rows and columns in "
        f"formula order (not clustered); the run stops above --max-cancor "
        f"{p['max_cancor']}.",
    }
    summary = [
        f"{len(meta)} samples ({', '.join(dropped) or 'none'} dropped by --min-cells); "
        f"{len(fractions)} of {n_genes_in} genes pass --min-frac {p['min_frac']}.",
        f"Largest median fraction: {top_covariate} ({medians[top_covariate]:.3f}); "
        f"Residuals {medians['Residuals']:.3f}.",
        f"Largest canonical correlation {pairs.max():.3f} "
        f"({' / '.join(pairs.idxmax())}), under --max-cancor {p['max_cancor']}."
        if len(pairs)
        else "One term: no canonical correlation to check.",
    ]
    write_report(
        outdir, f"Variance partition: {p['formula']}", p, KEYS, summary, captions
    )
    result = {
        "n_samples": len(meta),
        "n_donors": n_donors,
        "n_genes": len(fractions),
        "status_counts": {k: int(v) for k, v in table["status"].value_counts().items()},
        "dropped_samples": dropped,
        "top_covariate": top_covariate,
        "median_fractions": medians.round(6).to_dict(),
        "max_cancor": round(float(pairs.max()), 6) if len(pairs) else None,
        "top_genes": top_by_term,
        "params": p,
        "seed": p["seed"],
    }
    artifacts = [
        "report.md",
        "tables/variance_fractions.tsv",
        "figures/variance_violin.png",
    ]
    if p["adata"]:
        result["gene_id_col"] = p["gene_id_col"]
    write_result(
        outdir,
        result,
        summary,
        artifacts,
        context=p["context"],
        reference=p["reference"],
    )
    packages = ["numpy", "pandas", "matplotlib", "plotnine", "anndata", "scipy"]
    write_reproducibility(
        outdir, __file__, p, KEYS, full_argv, packages, extra=r_versions
    )
    return outdir


if __name__ == "__main__":
    raise SystemExit(run("variance-partition", analyse))
