"""Purpose: cross-cell-type effect sharing with mashr, from one DE table per cell type.
Inputs:  --de name=path (repeatable), each a TSV with gene, effect and SE columns
Outputs: <outdir>/report.md, result.json, tables/, figures/, reproducibility/
Run:     python scripts/mashr_sharing.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""

import argparse
import os
import shutil
import subprocess
import tempfile
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from skillkit import (
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
MASH_R = Path(__file__).resolve().parent / "mash.R"

ANALYSIS_ARGS = [
    "de",
    "gene_col",
    "effect_col",
    "se_col",
    "align",
    "covariances",
    "null_corr",
    "strong_lfsr",
    "lfsr_max",
    "sharing_factor",
    "reference",
    "seed",
]
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "de": [
        f"{ct}={SKILL_DIR / 'examples' / f'demo_{ct}.tsv'}"
        for ct in ["B", "CD4T", "CD8T", "Mono"]
    ],
    "gene_col": "gene",
    "effect_col": "log2_fc",
    "se_col": "lfc_se",
    "align": "intersect",
    "covariances": "canonical+ed",
    "null_corr": "identity",
    "strong_lfsr": 0.05,
    "lfsr_max": 0.05,
    "sharing_factor": 0.5,
    "reference": "synthetic demo genes, no genome build",
    "seed": 0,
    "outdir": "mashr_sharing_demo",
}

ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "bad_de": "--de must be name=path, got {value!r}",
    "duplicate_name": "--de name {name!r} is given twice",
    "too_few_tables": "need --de for at least 2 cell types, got {n}",
    "missing_file": "input file not found: {path}",
    "missing_column": "column {col!r} not found in the --de table for {name!r}",
    "duplicate_gene": "{n} genes repeat in the --de table for {name!r}, e.g. {gene!r}",
    "not_numeric": "column {col!r} in the --de table for {name!r} "
    "must be numeric with no missing values",
    "bad_se": "{n} rows of {col!r} in the --de table for {name!r} are <= 0",
    "no_genes": "zero genes are shared by all --de tables (--align intersect)",
    "no_rscript": "Rscript not found; set RSCRIPT or put Rscript on PATH",
    "missing_r_package": "{tail}; ask for it to be installed, do not install it here",
    "too_few_strong": "{n} genes pass --strong-lfsr {strong_lfsr} in any cell type, "
    "fewer than the {k} cell types that canonical+ed needs; ask the user about "
    "--strong-lfsr or --covariances canonical",
    "r_failed": "mash.R failed with exit {code}: {tail}",
}

# Fixed plot spec, not analysis thresholds. save_panel saves without a tight
# bbox, so every figure sets a layout engine to keep its labels on the page.
SEQUENTIAL = LinearSegmentedColormap.from_list("white_to_red", ["white", "#B2182B"])
SIG, OTHER = "#B2182B", "0.6"


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--de", action="append")
    parser.add_argument("--gene-col")
    parser.add_argument("--effect-col")
    parser.add_argument("--se-col")
    parser.add_argument("--align", choices=["intersect", "union"])
    parser.add_argument("--covariances", choices=["canonical", "canonical+ed"])
    parser.add_argument("--null-corr", choices=["identity", "em"])
    parser.add_argument("--strong-lfsr", type=float)
    parser.add_argument("--lfsr-max", type=float)
    parser.add_argument("--sharing-factor", type=float)
    parser.add_argument("--reference")
    parser.add_argument("--seed", type=int)
    return parse_required(parser, ANALYSIS_ARGS, DEMO, argv)


def check_values(p):
    rules = [
        ("--strong-lfsr", p["strong_lfsr"], 0 < p["strong_lfsr"] < 1, "in (0, 1)"),
        ("--lfsr-max", p["lfsr_max"], 0 < p["lfsr_max"] < 1, "in (0, 1)"),
        (
            "--sharing-factor",
            p["sharing_factor"],
            0 < p["sharing_factor"] < 1,
            "in (0, 1)",
        ),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ]
    for arg, value, ok, rule in rules:
        if not ok:
            fail(ERRORS, "bad_value", arg=arg, rule=rule, value=value)
    tables = {}
    for spec in p["de"]:
        name, sep, path = spec.partition("=")
        if not sep or not name or not path:
            fail(ERRORS, "bad_de", value=spec)
        if name in tables:
            fail(ERRORS, "duplicate_name", name=name)
        tables[name] = Path(path)
    if len(tables) < 2:
        fail(ERRORS, "too_few_tables", n=len(tables))
    return tables


def read_de(name, path, p):
    """One DE table -> gene-indexed frame with effect and se columns."""
    if not path.is_file():
        fail(ERRORS, "missing_file", path=path)
    de = pd.read_csv(path, sep="\t", dtype={p["gene_col"]: str})
    cols = [p["gene_col"], p["effect_col"], p["se_col"]]
    for col in cols:
        if col not in de.columns:
            fail(ERRORS, "missing_column", col=col, name=name)
    repeated = de[p["gene_col"]][de[p["gene_col"]].duplicated()]
    if len(repeated):
        fail(
            ERRORS,
            "duplicate_gene",
            n=repeated.nunique(),
            name=name,
            gene=repeated.iloc[0],
        )
    for col in cols[1:]:
        values = pd.to_numeric(de[col], errors="coerce")
        if not np.isfinite(values).all():
            fail(ERRORS, "not_numeric", col=col, name=name)
        de[col] = values
    if (de[p["se_col"]] <= 0).any():
        fail(
            ERRORS, "bad_se", n=(de[p["se_col"]] <= 0).sum(), col=p["se_col"], name=name
        )
    return de.set_index(p["gene_col"])[cols[1:]].set_axis(["effect", "se"], axis=1)


def align(tables, p):
    """Steps 1-2: read every table, align genes; NA marks a gene a table lacks."""
    frames = {name: read_de(name, path, p) for name, path in tables.items()}
    join = "inner" if p["align"] == "intersect" else "outer"
    bhat = pd.concat({n: f["effect"] for n, f in frames.items()}, axis=1, join=join)
    shat = pd.concat({n: f["se"] for n, f in frames.items()}, axis=1, join=join)
    if bhat.empty:
        fail(ERRORS, "no_genes")
    order = sorted(bhat.index)
    # Drops per group: genes each table loses to --align intersect (0 for union).
    dropped = {n: int((~f.index.isin(bhat.index)).sum()) for n, f in frames.items()}
    return bhat.loc[order], shat.loc[order], dropped


def run_mash(bhat, shat, p):
    """Steps 3-6 in mash.R; returns pm, psd, lfsr, sharing, n_strong, versions."""
    rscript = os.environ.get("RSCRIPT") or shutil.which("Rscript")
    if not rscript or not Path(rscript).is_file():
        fail(ERRORS, "no_rscript")
    # One BLAS thread so the posterior digits, and the checksums, repeat.
    threads = {k: "1" for k in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"]}
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        bhat.rename_axis("gene").to_csv(tmp / "bhat.tsv", sep="\t")
        shat.rename_axis("gene").to_csv(tmp / "shat.tsv", sep="\t")
        cmd = [
            rscript,
            str(MASH_R),
            str(tmp / "bhat.tsv"),
            str(tmp / "shat.tsv"),
            p["covariances"],
            p["null_corr"],
            str(p["strong_lfsr"]),
            str(p["lfsr_max"]),
            str(p["sharing_factor"]),
            str(p["seed"]),
            str(tmp),
        ]
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
            env={**os.environ, **threads},
        )
        tail = " | ".join(proc.stderr.strip().splitlines()[-5:])
        if proc.returncode == 3:
            fail(ERRORS, "missing_r_package", tail=tail)
        if proc.returncode == 4:
            fail(
                ERRORS,
                "too_few_strong",
                n=tail.split()[-1],
                strong_lfsr=p["strong_lfsr"],
                k=bhat.shape[1],
            )
        if proc.returncode != 0:
            fail(ERRORS, "r_failed", code=proc.returncode, tail=tail)
        out = {
            name: pd.read_csv(tmp / f"{name}.tsv", sep="\t", index_col=0)
            for name in ["pm", "psd", "lfsr", "sharing", "null_corr"]
        }
    lines = proc.stdout.strip().splitlines()
    n_strong = lines[0].split()[-1]
    out["n_strong"] = None if n_strong == "NA" else int(n_strong)
    return out, [f"Rscript {rscript}", *lines[1:]]


def pairs_table(matrix, value):
    """Upper triangle of a square cell-type matrix as long rows, sorted by pair."""
    return pd.DataFrame(
        [
            {"cell_type_a": a, "cell_type_b": b, value: matrix.loc[a, b]}
            for a, b in combinations(matrix.index, 2)
        ]
    )


def raw_spearman(bhat):
    """C-an-40: Spearman of raw effects per pair, over genes measured in both."""
    rows = []
    for a, b in combinations(bhat.columns, 2):
        both = bhat[[a, b]].dropna()
        rho, pvalue = stats.spearmanr(both[a], both[b])
        rows.append(
            {
                "cell_type_a": a,
                "cell_type_b": b,
                "n_genes": len(both),
                "rho": rho,
                "pvalue": pvalue,
            }
        )
    return pd.DataFrame(rows)


def draw_figures(figures, mash, posterior, counts, pair, p):
    names = list(mash["sharing"].index)
    sharing = mash["sharing"].loc[names, names]
    fig, ax = plt.subplots(layout="constrained")
    mesh = ax.pcolormesh(sharing, cmap=SEQUENTIAL, vmin=0, vmax=1)
    for i in range(len(names)):
        for j in range(len(names)):
            ax.text(
                j + 0.5, i + 0.5, f"{sharing.iloc[i, j]:.2f}", ha="center", va="center"
            )
    ax.set_xticks(np.arange(len(names)) + 0.5, names, rotation=90)
    ax.set_yticks(np.arange(len(names)) + 0.5, names)
    ax.invert_yaxis()
    fig.colorbar(mesh, ax=ax, label="fraction of effects shared")
    heat = sharing.rename_axis("cell_type_a").reset_index()
    save_panel(fig, heat, figures, "sharing_heatmap", 89, 75)
    plt.close(fig)

    fig, ax = plt.subplots(layout="constrained")
    ax.bar(counts["n_cell_types"], counts["n_genes"], color=OTHER)
    ax.set_xticks(counts["n_cell_types"])
    ax.set_xlabel(f"cell types with lfsr < {p['lfsr_max']}")
    ax.set_ylabel("genes")
    save_panel(fig, counts, figures, "significant_hist", 89, 70)
    plt.close(fig)

    a, b = pair
    wide = posterior.pivot(index="gene", columns="cell_type")
    data = pd.DataFrame(
        {
            f"raw_{a}": wide[("raw_effect", a)],
            f"raw_{b}": wide[("raw_effect", b)],
            f"pm_{a}": wide[("pm", a)],
            f"pm_{b}": wide[("pm", b)],
            "significant_either": wide[("significant", a)].eq(True)
            | wide[("significant", b)].eq(True),
        }
    ).dropna()
    colour = np.where(data["significant_either"], SIG, OTHER)
    fig, axes = plt.subplots(1, 2, layout="constrained", sharex=True, sharey=True)
    for ax, kind, label in zip(
        axes, ["raw", "pm"], ["raw effect", "posterior mean"], strict=True
    ):
        ax.scatter(data[f"{kind}_{a}"], data[f"{kind}_{b}"], s=4, c=colour)
        ax.axline((0, 0), slope=1, ls="--", lw=0.5, c="black")
        ax.set_xlabel(f"{label}, {a}")
        ax.set_ylabel(f"{label}, {b}")
    save_panel(fig, data.reset_index(), figures, "raw_vs_posterior", 160, 75)
    plt.close(fig)


def analyse(argv, full_argv):
    p = parse_args(argv)
    tables = check_values(p)
    seed_everything(p["seed"])
    bhat, shat, dropped = align(tables, p)
    if Path(p["outdir"]).exists():
        make_outdir(p["outdir"])  # fails fast, before the R run
    # Every refusal comes before make_outdir, so a refused run leaves no outdir.
    mash, r_versions = run_mash(bhat, shat, p)
    outdir = make_outdir(p["outdir"])

    tables_dir, figures = outdir / "tables", outdir / "figures"
    names = list(bhat.columns)
    long = {
        "raw_effect": bhat,
        "raw_se": shat,
        "pm": mash["pm"],
        "psd": mash["psd"],
        "lfsr": mash["lfsr"],
    }
    posterior = pd.concat(
        {k: v.loc[bhat.index, names].stack(future_stack=True) for k, v in long.items()},
        axis=1,
    ).rename_axis(["gene", "cell_type"])
    posterior = posterior.dropna(subset=["raw_effect"]).reset_index()
    posterior["significant"] = posterior["lfsr"] < p["lfsr_max"]
    posterior = posterior.sort_values(["gene", "cell_type"]).reset_index(drop=True)
    posterior.to_csv(tables_dir / "posterior.tsv", sep="\t", index=False)

    n_sig = posterior.groupby("gene")["significant"].sum()
    counts = (
        n_sig.value_counts()
        .reindex(range(1, len(names) + 1), fill_value=0)
        .rename_axis("n_cell_types")
        .reset_index(name="n_genes")
    )
    counts.to_csv(tables_dir / "significant_counts.tsv", sep="\t", index=False)

    sharing = pairs_table(mash["sharing"].loc[names, names], "sharing")
    sharing.to_csv(tables_dir / "pairwise_sharing.tsv", sep="\t", index=False)
    spearman = raw_spearman(bhat)
    spearman.to_csv(tables_dir / "raw_spearman.tsv", sep="\t", index=False)
    null_corr = pairs_table(mash["null_corr"].loc[names, names], "v")
    null_corr.to_csv(tables_dir / "null_corr.tsv", sep="\t", index=False)

    # Scatter pair: highest sharing; a tie goes to the first pair in name order.
    top = sharing.sort_values(
        ["sharing", "cell_type_a", "cell_type_b"], ascending=[False, True, True]
    ).iloc[0]
    pair = (top["cell_type_a"], top["cell_type_b"])
    plt.switch_backend("Agg")
    draw_figures(figures, mash, posterior, counts, pair, p)

    k = len(names)
    n_shared = int((n_sig == k).sum())
    n_specific = int((n_sig == 1).sum())
    n_any = int((n_sig >= 1).sum())
    rule = f"lfsr < {p['lfsr_max']}"
    captions = {
        "sharing_heatmap": "Heatmap of mashr pairwise sharing between cell types: "
        f"of genes with {rule} in either type, the fraction whose posterior means "
        f"share sign and lie within a factor {p['sharing_factor']} of each other; "
        "rows and columns in --de order (not clustered).",
        "significant_hist": f"Bar chart of genes by the number of cell types (of {k}) "
        f"in which they pass {rule}; genes passing in none are not shown.",
        "raw_vs_posterior": f"Scatter of {pair[0]} against {pair[1]}, the pair with "
        "the highest sharing: raw effects (left) and mashr posterior means (right); "
        f"red genes pass {rule} in either type, the dashed line is y = x.",
    }
    summary = [
        f"{len(bhat)} genes x {k} cell types ({', '.join(names)}), "
        f"--align {p['align']}, covariances {p['covariances']}; genes dropped by "
        "--align per cell type: "
        + ", ".join(f"{n} {v}" for n, v in dropped.items())
        + ".",
        f"{n_any} genes pass {rule} in at least one cell type: {n_shared} in all {k}, "
        f"{n_specific} in exactly one.",
        f"Pairwise sharing ranges {sharing['sharing'].min():.2f}-"
        f"{sharing['sharing'].max():.2f}; raw Spearman rho "
        f"{spearman['rho'].min():.2f}-{spearman['rho'].max():.2f}.",
    ]
    shown = {**p, "de": " ".join(p["de"])}
    write_report(
        outdir, "mashr cross-cell-type sharing", shown, ANALYSIS_ARGS, summary, captions
    )
    result = {
        "n_genes": len(bhat),
        "cell_types": names,
        "genes_dropped_by_align": dropped,
        "n_strong": mash["n_strong"],
        "null_corr": p["null_corr"],
        "n_significant_any": n_any,
        "n_shared_all": n_shared,
        "n_specific": n_specific,
        "significant_per_cell_type": posterior.groupby("cell_type")["significant"]
        .sum()
        .astype(int)
        .to_dict(),
        "scatter_pair": list(pair),
        "params": p,
        "seed": p["seed"],
    }
    artifacts = ["report.md", "tables/posterior.tsv", "figures/sharing_heatmap.png"]
    write_result(outdir, result, summary, artifacts, reference=p["reference"])
    write_reproducibility(
        outdir,
        __file__,
        p,
        ANALYSIS_ARGS,
        full_argv,
        ["numpy", "pandas", "scipy", "matplotlib"],
        extra=r_versions,
    )
    return outdir


if __name__ == "__main__":
    raise SystemExit(run("mashr-sharing", analyse))
