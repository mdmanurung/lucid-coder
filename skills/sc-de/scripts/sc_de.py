"""Purpose: per-feature DE: MAST hurdle model on cells, or a linear model on samples.
Inputs:  --mode mast: --adata <log-normalised h5ad>; --mode lm: --input <long.tsv>
Outputs: <outdir>/report.md, result.json, tables/, figures/, reproducibility/
Run:     python scripts/sc_de.py --demo --outdir demo_out   (add --mode lm for NPX)
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Both modes fit in R (mast.R, lm.R via Rscript); BH runs here. Bad input:
         "ERROR: ..." and exit 2. Steps run in a fixed order; see SKILL.md.
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
import skillkit
from plotnine import labs
from plots import plot_box, plot_scatter_fit, plot_volcano
from scipy import sparse, stats
from theme import save_panel

SCRIPTS = Path(__file__).resolve().parent
SKILL_DIR = SCRIPTS.parent
REQUIRED = ["mode", "formula", "padj_max", "lfc_min", "seed", "reference", "context"]
MAST_ARGS = [
    "adata",
    "gene_id_col",
    "cell_type_col",
    "cell_type",
    "donor_col",
    "contrast_term",
    "n_cells_per_donor",
    "min_frac",
    "cdr",
    "log_base",
]
LM_ARGS = ["input", "term"]
OPTIONAL = [*MAST_ARGS, "layer", *LM_ARGS]
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO_COMMON = {
    "padj_max": 0.05,
    "lfc_min": 0.5,
    "seed": 0,
    "layer": None,
    "reference": "synthetic demo, no genome",
}
DEMO = {
    "mast": {
        **DEMO_COMMON,
        "mode": "mast",
        "adata": str(SKILL_DIR / "examples" / "demo.h5ad"),
        "gene_id_col": "gene_ids",
        "context": "synthetic demo: 6 donors, ctrl vs stim T cells, no real tissue",
        "cell_type_col": "cell_type",
        "cell_type": "T",
        "donor_col": "donor",
        "formula": "~ condition",
        "contrast_term": "condition,stim,ctrl",
        "n_cells_per_donor": "50",
        "min_frac": 0.1,
        "cdr": "yes",
        "log_base": "e",
        "input": None,
        "term": None,
    },
    "lm": {
        **DEMO_COMMON,
        "mode": "lm",
        "input": str(SKILL_DIR / "examples" / "demo_npx.tsv"),
        "formula": "~ status + age + sex",
        "term": "status,disease,healthy",
        "context": "synthetic demo: 40 samples, disease vs healthy NPX, no real cohort",
        **dict.fromkeys(MAST_ARGS),
    },
}

ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "needs_args": "--mode {mode} needs {args}",
    "unused_args": "--mode {mode} does not use {args}",
    "bad_formula": "--formula must start with ~ and name columns, got {value!r}",
    "bad_term": "{arg} must be column or column,level,reference, got {value!r}",
    "term_not_in_formula": "{arg} column {col!r} is not in --formula",
    "missing_file": "input file not found: {path}",
    "missing_layer": "layer not found in h5ad: {layer}",
    "missing_gene_id_col": "gene id column not found in var: {col}",
    "missing_column": "input column not found: {col}",
    "cell_type_absent": "cell type {cell_type!r} not found in column {col!r}",
    "level_absent": "level {level!r} not found in column {col!r}",
    "not_numeric": "column {col!r} must be numeric",
    "not_log_normalised": (
        "{source} must hold log-normalised expression (non-negative, not all "
        "integers); raw counts or scaled data cannot go into MAST"
    ),
    "duplicate_rows": "{n} rows repeat a (feature, sample) combination; one value each",
    "no_genes": "no gene is detected in at least {min_frac} of the kept cells",
    "no_rscript": "Rscript not found; set RSCRIPT or put Rscript on PATH",
    "r_failed": "{script} failed with exit {code}: {tail}",
}

# Fixed plot spec, not analysis thresholds: corpus jitter width (position_jitter(0.2)),
# reference blue and level red from palettes.pal_condition.
JITTER_WIDTH = 0.2
LEVEL_COLOURS = ("#4477AA", "#EE6677")
PANEL_MM = 60
DE_COLUMNS = {
    "mast": ["feature", "log2_fc", "se", "stat", "pvalue"],
    "lm": ["feature", "estimate", "se", "stat", "pvalue"],
}


def fail(key, **fields):
    skillkit.fail(ERRORS, key, **fields)


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--mode", choices=["mast", "lm"])
    parser.add_argument("--formula")
    parser.add_argument("--padj-max", type=float)
    parser.add_argument("--lfc-min", type=float)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--reference")
    parser.add_argument("--context")
    parser.add_argument("--adata")
    parser.add_argument("--gene-id-col")
    parser.add_argument("--layer")
    parser.add_argument("--cell-type-col")
    parser.add_argument("--cell-type")
    parser.add_argument("--donor-col")
    parser.add_argument("--contrast-term")
    parser.add_argument("--n-cells-per-donor")
    parser.add_argument("--min-frac", type=float)
    parser.add_argument("--cdr", choices=["yes", "no"])
    parser.add_argument("--log-base", choices=["2", "e"])
    parser.add_argument("--input")
    parser.add_argument("--term")
    # --mode picks which demo values --demo pins; it may stand beside --demo.
    mode = parser.parse_known_args(argv)[0].mode or "mast"
    return skillkit.parse_required(
        parser, REQUIRED, DEMO[mode], argv, OPTIONAL, ("outdir", "mode")
    )


def check_values(p):
    """Step 1: per-mode arguments, ranges, formula and term syntax."""
    needs = MAST_ARGS if p["mode"] == "mast" else LM_ARGS
    missing = [k for k in needs if p[k] is None]
    if missing:
        fail("needs_args", mode=p["mode"], args=" ".join(map(skillkit.flag, missing)))
    unused = [k for k in OPTIONAL if k not in [*needs, "layer"] and p[k] is not None]
    if p["mode"] == "lm" and p["layer"] is not None:
        unused.append("layer")
    if unused:
        fail("unused_args", mode=p["mode"], args=" ".join(map(skillkit.flag, unused)))
    rules = [
        ("--padj-max", p["padj_max"], 0 < p["padj_max"] < 1, "in (0, 1)"),
        ("--lfc-min", p["lfc_min"], p["lfc_min"] >= 0, ">= 0"),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ]
    if p["mode"] == "mast":
        n = p["n_cells_per_donor"]
        ok = n == "all" or (n.isdigit() and int(n) >= 1)
        rules += [
            ("--n-cells-per-donor", n, ok, "'all' or an integer >= 1"),
            ("--min-frac", p["min_frac"], 0 <= p["min_frac"] < 1, "in [0, 1)"),
        ]
    for arg, value, ok, rule in rules:
        if not ok:
            fail("bad_value", arg=arg, rule=rule, value=value)

    # Names followed by "(" are functions (log, poly), not columns.
    names = re.findall(r"([A-Za-z_.][\w.]*)(\s*\()?", p["formula"])
    columns = [name for name, call in names if not call]
    # Fixed-effect columns: the formula without its (1|x) random terms.
    fixed_part = re.sub(r"\(\s*1\s*\|[^)]*\)", "", p["formula"])
    found = re.findall(r"([A-Za-z_.][\w.]*)(\s*\()?", fixed_part)
    fixed = list(dict.fromkeys(name for name, call in found if not call))
    if not p["formula"].strip().startswith("~") or not columns:
        fail("bad_formula", value=p["formula"])
    key = "contrast_term" if p["mode"] == "mast" else "term"
    parts = [x.strip() for x in p[key].split(",")]
    if (
        len(parts) not in (1, 3)
        or not all(parts)
        or len(set(parts[1:])) < len(parts) - 1
    ):
        fail("bad_term", arg=skillkit.flag(key), value=p[key])
    if parts[0] not in columns:
        fail("term_not_in_formula", arg=skillkit.flag(key), col=parts[0])
    return list(dict.fromkeys(columns)), fixed, parts


def check_rank(frame, fixed):
    """Design rank on the fitted rows; an all-numeric column enters as a number."""
    meta = pd.DataFrame(index=frame.index)
    for col in fixed:
        number = pd.to_numeric(frame[col], errors="coerce")
        meta[col] = number if number.notna().all() else frame[col].astype(str)
    skillkit.check_design_rank(meta, fixed)


def check_term(frame, term):
    """The term column is numeric, or holds both level and reference."""
    col = term[0]
    if len(term) == 1:
        values = pd.to_numeric(frame[col], errors="coerce")
        if values.isna().any():
            fail("not_numeric", col=col)
        frame[col] = values
        return
    present = set(frame[col].astype(str))
    for level in term[1:]:
        if level not in present:
            fail("level_absent", level=level, col=col)


def read_cells(p, columns, term):
    """Steps 2-6: read, subset, check log scale, downsample, CDR, gene filter."""
    import anndata as ad

    path = Path(p["adata"])
    if not path.is_file():
        fail("missing_file", path=path)
    adata = ad.read_h5ad(path)
    # var_names are the symbols; --gene-id-col none: the source has no gene ids.
    ids = None
    if p["gene_id_col"] != "none":
        if p["gene_id_col"] not in adata.var.columns:
            fail("missing_gene_id_col", col=p["gene_id_col"])
        ids = adata.var[p["gene_id_col"]].astype(str).tolist()
    gene_ids = skillkit.check_gene_ids(ids, adata.var_names.tolist())
    skillkit.check_unique(adata.obs_names, "cell ids")
    if p["layer"] is not None and p["layer"] not in adata.layers:
        fail("missing_layer", layer=p["layer"])
    needed = [p["cell_type_col"], p["donor_col"], *columns]
    for col in dict.fromkeys(needed):
        if col not in adata.obs.columns:
            fail("missing_column", col=col)
    in_type = (adata.obs[p["cell_type_col"]].astype(str) == p["cell_type"]).to_numpy()
    if not in_type.any():
        fail("cell_type_absent", cell_type=p["cell_type"], col=p["cell_type_col"])
    adata = adata[in_type]
    x = adata.layers[p["layer"]] if p["layer"] else adata.X
    x = sparse.csr_matrix(x, dtype=np.float64)
    values = x.data
    if values.size and ((values < 0).any() or np.all(values == np.round(values))):
        fail("not_log_normalised", source=f"layer {p['layer']}" if p["layer"] else "X")
    if p["log_base"] == "e":
        x = x / np.log(2)  # log1p is natural log; MAST logFC is then log2
    obs = adata.obs[list(dict.fromkeys(needed))].astype(str)
    check_term(obs, term)

    # Seeded downsampling per donor (and per level of a categorical term).
    rng = np.random.default_rng(p["seed"])
    keys = [p["donor_col"], *([term[0]] if len(term) == 3 else [])]
    n = p["n_cells_per_donor"]
    keep = []
    for _, idx in sorted(obs.groupby(keys, observed=True).indices.items()):
        if n != "all" and len(idx) > int(n):
            idx = np.sort(rng.choice(idx, int(n), replace=False))
        keep.append(idx)
    keep = np.sort(np.concatenate(keep))
    x, obs = x[keep], obs.iloc[keep]

    detected = np.asarray((x > 0).sum(axis=1)).ravel()
    if p["cdr"] == "yes":
        # Cellular detection rate over all genes, centred and scaled (MAST cngeneson).
        obs = obs.assign(cngeneson=(detected - detected.mean()) / detected.std(ddof=1))
    frac = np.asarray((x > 0).mean(axis=0)).ravel()
    genes = frac >= p["min_frac"]
    if not genes.any():
        fail("no_genes", min_frac=p["min_frac"])
    expr = pd.DataFrame(
        x[:, genes].toarray().T,
        index=adata.var_names[genes],
        columns=obs.index,
    )
    # Why an untested gene has no result: never detected, or under --min-frac.
    reasons = {
        g: "not_detected" if f == 0 else "not_tested"
        for g, f in zip(adata.var_names, frac, strict=True)
    }
    return expr, obs, gene_ids, reasons


def read_long(p, columns, term):
    """Steps 2-3 (lm): read the long TSV, drop incomplete rows, refuse repeats."""
    path = Path(p["input"])
    if not path.is_file():
        fail("missing_file", path=path)
    data = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    cols = list(dict.fromkeys(["sample", "feature", "value", *columns]))
    for col in cols:
        if col not in data.columns:
            fail("missing_column", col=col)
    data = data[cols]
    # Rows lost to a missing value, counted per term level (all rows if numeric).
    lost = data[data.isna().any(axis=1)]
    by = lost[term[0]].fillna("missing") if len(term) == 3 else ["all"] * len(lost)
    dropped = pd.Series(1, index=lost.index).groupby(by).sum().sort_index()
    # A feature whose rows all drop, or whose term is not estimable, is not tested.
    reasons = dict.fromkeys(sorted(data["feature"].dropna().unique()), "not_tested")
    data = data.dropna()
    value = pd.to_numeric(data["value"], errors="coerce")
    if value.isna().any():
        fail("not_numeric", col="value")
    data["value"] = value
    check_term(data, term)
    repeated = int(data.duplicated(["feature", "sample"], keep=False).sum())
    if repeated:
        fail("duplicate_rows", n=repeated)
    return (
        data.sort_values(["feature", "sample"]),
        {str(k): int(v) for k, v in dropped.items()},
        reasons,
    )


def run_r(script, files, args, tmp):
    """Run an R engine; return its table and the version lines it prints."""
    rscript = os.environ.get("RSCRIPT") or shutil.which("Rscript")
    if not rscript or not Path(rscript).is_file():
        fail("no_rscript")
    out = Path(tmp) / "de.tsv"
    cmd = [rscript, str(SCRIPTS / script), *map(str, files), *args, str(out)]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        tail = " | ".join(proc.stderr.strip().splitlines()[-5:])
        fail("r_failed", script=script, code=proc.returncode, tail=tail)
    res = pd.read_csv(out, sep="\t", keep_default_na=False, na_values=["NA", "NaN"])
    return res, [f"Rscript {rscript}", *proc.stdout.strip().splitlines()]


def term_args(term):
    return [term[0], *(term[1:] if len(term) == 3 else ["", ""])]


def draw_figures(figures, table, data, p, term, effect):
    """Step 9: volcano, and the top feature by group (or against a numeric term)."""
    captions = {}
    tested = table.dropna(subset=["padj"])
    plot = plot_volcano(
        tested, effect, "padj", padj_max=p["padj_max"], lfc_min=p["lfc_min"]
    )
    save_panel(plot, tested, figures, "volcano", PANEL_MM, PANEL_MM)
    unit = "cells" if p["mode"] == "mast" else "samples"
    captions["volcano"] = (
        f"Volcano plot of {len(tested)} features: {effect} against -log10 BH-adjusted "
        f"p; red (up) and blue (down) mark padj < {p['padj_max']} and "
        f"|{effect}| > {p['lfc_min']} (dashed lines), grey the rest. "
        f"P-values treat {unit} as independent replicates."
    )
    if tested.empty:
        return captions
    top = tested.iloc[0]
    feature = top["feature"]
    rows = data(feature)
    unit_note = "one cell" if p["mode"] == "mast" else "one sample"
    value_label = (
        f"{feature} expression (log2 normalised)"
        if p["mode"] == "mast"
        else f"{feature} value (units of --input)"
    )
    if len(term) == 3:
        pal = {term[2]: LEVEL_COLOURS[0], term[1]: LEVEL_COLOURS[1]}
        plot = plot_box(
            rows,
            term[0],
            feature,
            pal,
            jitter_width=JITTER_WIDTH,
            seed=p["seed"],
            value_label=value_label,
        )
        shown = f"Box plot of {feature} per {term[0]} ({term[1]} vs {term[2]})"
    else:
        plot = plot_scatter_fit(rows, term[0], feature) + labs(
            x=f"{term[0]} (units of --input)", y=value_label
        )
        shown = f"Scatter plot of {feature} against {term[0]} with a least-squares line"
    save_panel(plot, rows, figures, "top_feature", PANEL_MM, PANEL_MM)
    captions["top_feature"] = (
        f"{shown}, the feature with the lowest padj; each point is {unit_note}. "
        f"{effect} {top[effect]:.3g}, padj {top['padj']:.2g} (BH across "
        f"{len(tested)} features)."
    )
    return captions


def analyse(argv, argv_as_run):
    p = parse_args(argv)
    columns, fixed, term = check_values(p)
    skillkit.seed_everything(p["seed"])
    mode = p["mode"]
    formula = p["formula"]
    with tempfile.TemporaryDirectory() as tmp:
        if mode == "mast":
            expr, obs, gene_ids, reasons = read_cells(p, columns, term)
            if p["cdr"] == "yes":
                formula = f"{formula} + cngeneson"
                fixed = [*fixed, "cngeneson"]
            check_rank(obs, fixed)
            # ponytail: dense TSV hand-off; switch to Matrix Market past ~1e8 values.
            expr.to_csv(Path(tmp) / "expr.tsv", sep="\t", float_format="%.6g")
            obs.to_csv(Path(tmp) / "meta.tsv", sep="\t")
            files = [Path(tmp) / "expr.tsv", Path(tmp) / "meta.tsv"]
            args = [formula, *term_args(term), str(p["seed"])]
            res, r_env = run_r("mast.R", files, args, tmp)
            n_units = {"n_cells": len(obs), "n_donors": obs[p["donor_col"]].nunique()}
            dropped = {}
        else:
            data, dropped, reasons = read_long(p, columns, term)
            gene_ids = None
            check_rank(data, fixed)
            data.to_csv(Path(tmp) / "long.tsv", sep="\t", index=False)
            args = [formula, *term_args(term)]
            res, r_env = run_r("lm.R", [Path(tmp) / "long.tsv"], args, tmp)
            n_units = {"n_samples": data["sample"].nunique()}

    # Steps 7-8: BH across tested features, call significance, sort.
    effect = DE_COLUMNS[mode][1]
    table = res[DE_COLUMNS[mode]].copy()
    table["feature"] = table["feature"].astype(str)
    table["padj"] = np.nan
    tested = table["pvalue"].notna()
    if tested.any():
        table.loc[tested, "padj"] = stats.false_discovery_control(
            table.loc[tested, "pvalue"], method="bh"
        )
    table["significant"] = (table["padj"] < p["padj_max"]) & (
        table[effect].abs() > p["lfc_min"]
    )
    # Every input feature keeps a row; untested ones have NA statistics.
    table = skillkit.with_status(table[tested], list(reasons), "feature", reasons)
    table["significant"] = table["significant"].eq(True)
    if gene_ids is not None:
        table.insert(1, "gene_id", table["feature"].map(gene_ids))
    table = table.sort_values(["padj", "feature"], na_position="last")

    outdir = skillkit.make_outdir(p["outdir"])
    # Mast features are genes: gene_symbol, gene_id (the names mashr-sharing reads).
    shown = (
        table.rename(columns={"feature": "gene_symbol"}) if mode == "mast" else table
    )
    shown.to_csv(
        outdir / "tables" / "de_results.tsv", sep="\t", index=False, float_format="%.6g"
    )
    plt.switch_backend("Agg")
    if mode == "mast":

        def rows_for(feature):
            return obs.assign(**{feature: expr.loc[feature].to_numpy()})

    else:

        def rows_for(feature):
            rows = data[data["feature"] == feature]
            return rows.rename(columns={"value": feature})

    captions = draw_figures(outdir / "figures", table, rows_for, p, term, effect)

    n_tested = int(tested.sum())
    n_sig = int(table["significant"].sum())
    sig = table.loc[table["significant"], "feature"].tolist()
    contrast = f"{term[1]} vs {term[2]}" if len(term) == 3 else f"per unit {term[0]}"
    unit = (
        f"{n_units['n_cells']} cells from {n_units['n_donors']} donors; the cell "
        "is the unit, so p-values overstate donor-level evidence"
        if mode == "mast"
        else f"{n_units['n_samples']} samples; the sample is the unit"
        f" ({sum(dropped.values())} rows with a missing value dropped"
        + (f": {dropped}" if dropped else "")
        + ")"
    )
    summary = [
        f"{mode}: {n_tested} of {len(table)} features tested for {term[0]} "
        f"({contrast}) with {formula}; {unit}.",
        f"{n_sig} features with padj < {p['padj_max']} (BH) and "
        f"|{effect}| > {p['lfc_min']}.",
        f"Significant: {', '.join(sig[:10]) or 'none'}"
        + (f" and {len(sig) - 10} more." if len(sig) > 10 else "."),
    ]
    keys = [*REQUIRED, *(MAST_ARGS + ["layer"] if mode == "mast" else LM_ARGS)]
    title = f"Differential expression ({mode}): {term[0]} {contrast}"
    skillkit.write_report(outdir, title, p, keys, summary, captions)
    result = {
        "mode": mode,
        "formula_fitted": formula,
        **n_units,
        "dropped_by_level": dropped,
        "n_features": len(table),
        "n_tested": n_tested,
        "status_counts": {k: int(v) for k, v in table["status"].value_counts().items()},
        "n_significant": n_sig,
        "significant_features": sig,
        "params": p,
        "seed": p["seed"],
    }
    artifacts = [
        "report.md",
        "tables/de_results.tsv",
        *[f"figures/{n}.png" for n in captions],
    ]
    if mode == "mast":
        result["gene_id_col"] = p["gene_id_col"]
    skillkit.write_result(
        outdir,
        result,
        summary,
        artifacts,
        context=p["context"],
        reference=p["reference"],
    )
    packages = ["anndata", "numpy", "pandas", "scipy", "matplotlib", "plotnine"]
    skillkit.write_reproducibility(
        outdir, __file__, p, keys, argv_as_run, packages, extra=r_env
    )
    return outdir


if __name__ == "__main__":
    raise SystemExit(skillkit.run("sc-de", analyse))
