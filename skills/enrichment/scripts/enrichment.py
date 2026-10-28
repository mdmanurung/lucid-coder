"""Purpose: gene-set enrichment (fgsea or ORA) on DE tables, BH within each contrast.
Inputs:  repeated --de name=<DE table TSV> (pseudobulk-de tables/de_results.tsv);
         --gene-sets <GMT, or TSV with set and gene columns>
Outputs: <outdir>/report.md, result.json, tables/enrichment.tsv, figures/,
         reproducibility/
Run:     python scripts/enrichment.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         fgsea runs in R (scripts/fgsea.R); ORA is a one-sided Fisher test in scipy.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
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
from plotnine import aes, geom_text
from plots import diverging, plot_heatmap
from scipy import stats
from theme import save_panel

SKILL_DIR = Path(__file__).resolve().parents[1]
EXAMPLES = SKILL_DIR / "examples"
FGSEA_R = Path(__file__).resolve().parent / "fgsea.R"
REQUIRED = [
    "de",
    "gene_sets",
    "method",
    "min_size",
    "max_size",
    "padj_max_sets",
    "limit",
    "seed",
    "reference",
    "context",
]
METHOD_ARGS = {
    "fgsea": ["rank_col", "n_perm_simple"],
    "ora": ["universe", "padj_max", "lfc_min", "direction"],
}
OPTIONAL = [*METHOD_ARGS["fgsea"], *METHOD_ARGS["ora"]]
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "de": [f"{c}={EXAMPLES / f'demo_{c}.tsv'}" for c in ("T", "Mono")],
    "gene_sets": str(EXAMPLES / "demo.gmt"),
    "method": "fgsea",
    "rank_col": "stat",
    "n_perm_simple": 10000,
    "universe": None,
    "padj_max": None,
    "lfc_min": None,
    "direction": None,
    "min_size": 10,
    "max_size": 500,
    "padj_max_sets": 0.05,
    "limit": 2.5,
    "seed": 0,
    "reference": "synthetic demo, no genome build",
    "context": "synthetic demo: 1000 genes, planted T contrast and null Mono "
    "contrast, no real tissue",
}
# `--demo --method ora` swaps the fgsea-only values for these demo-only ones.
DEMO_ORA = {
    "rank_col": None,
    "n_perm_simple": None,
    "universe": "tested",
    "padj_max": 0.05,
    "lfc_min": 0.5,
    "direction": "both",
}

ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "needs_args": "--method {method} needs {args}",
    "unused_args": "--method {method} does not use {args}",
    "bad_de": "--de must be name=path, the name made of letters, digits, '.', "
    "'_' or '-', got {value!r}",
    "duplicate_de": "contrast {name!r} is given more than once",
    "missing_file": "input file not found: {path}",
    "missing_column": "column {col!r} not found in {path}",
    "not_numeric": "column {col!r} in {path} must be numeric and finite",
    "bad_gene_sets": "--gene-sets must be a .gmt file or a .tsv with set and gene "
    "columns, got {path}",
    "no_overlap": "no gene of --gene-sets is in contrast {contrast!r}; "
    "check the gene ID type (symbol or Ensembl)",
    "no_sets": "no gene set has {min_size} to {max_size} genes in contrast "
    "{contrast!r}",
    "no_foreground": "no contrast has a {direction} gene with "
    "padj < {padj_max} and |log2_fc| > {lfc_min}",
    "no_rscript": "Rscript not found; set RSCRIPT or put Rscript on PATH",
    "missing_r_package": "R package fgsea is not installed for {rscript}; "
    "ask for it to be installed",
    "r_failed": "fgsea.R failed with exit {code}: {tail}",
}

# Fixed plot spec, not analysis thresholds: the 20 sets with the lowest padj in
# any contrast, dot area scaled to 4 mm at the largest -log10 padj (ggplot2
# scale_size_area), running-ES curve in the volcano's significant red.
N_SETS = 20
MAX_DOT_MM = 4
CURVE_COLOUR = "darkred"
PACKAGES = ["numpy", "pandas", "scipy", "matplotlib", "plotnine"]


def fail(key, **fields):
    skillkit.fail(ERRORS, key, **fields)


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--de", action="append")
    parser.add_argument("--gene-sets")
    parser.add_argument("--method", choices=list(METHOD_ARGS))
    parser.add_argument("--rank-col")
    parser.add_argument("--n-perm-simple", type=int)
    parser.add_argument("--universe", choices=["tested", "tested-in-sets"])
    parser.add_argument("--padj-max", type=float)
    parser.add_argument("--lfc-min", type=float)
    parser.add_argument("--direction", choices=["up", "down", "both"])
    parser.add_argument("--min-size", type=int)
    parser.add_argument("--max-size", type=int)
    parser.add_argument("--padj-max-sets", type=float)
    parser.add_argument("--limit", type=float)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--reference")
    parser.add_argument("--context")
    p = skillkit.parse_required(
        parser,
        REQUIRED,
        DEMO,
        argv,
        optional=OPTIONAL,
        demo_extra=("outdir", "method"),
    )
    if "--demo" in argv and p["method"] == "ora":
        p |= DEMO_ORA
    return p


def check_values(p):
    """Step 1: method arguments, value ranges, --de as {name: path}."""
    needs = METHOD_ARGS[p["method"]]
    missing = [k for k in needs if p[k] is None]
    if missing:
        fail(
            "needs_args", method=p["method"], args=" ".join(map(skillkit.flag, missing))
        )
    unused = [k for k in OPTIONAL if k not in needs and p[k] is not None]
    if unused:
        fail(
            "unused_args", method=p["method"], args=" ".join(map(skillkit.flag, unused))
        )
    rules = [
        ("--min-size", p["min_size"], p["min_size"] >= 1, "an integer >= 1"),
        ("--max-size", p["max_size"], p["max_size"] >= p["min_size"], ">= --min-size"),
        (
            "--padj-max-sets",
            p["padj_max_sets"],
            0 < p["padj_max_sets"] <= 1,
            "in (0, 1]",
        ),
        ("--limit", p["limit"], p["limit"] > 0, "> 0"),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ]
    if p["method"] == "fgsea":
        rules.append(
            ("--n-perm-simple", p["n_perm_simple"], p["n_perm_simple"] >= 1, ">= 1")
        )
    else:
        rules += [
            ("--padj-max", p["padj_max"], 0 < p["padj_max"] <= 1, "in (0, 1]"),
            ("--lfc-min", p["lfc_min"], p["lfc_min"] >= 0, ">= 0"),
        ]
    for arg, value, ok, rule in rules:
        if not ok:
            fail("bad_value", arg=arg, rule=rule, value=value)
    contrasts = {}
    for item in p["de"]:
        name, sep, path = item.partition("=")
        if not sep or not path or not re.fullmatch(r"[\w.-]+", name):
            fail("bad_de", value=item)
        if name in contrasts:
            fail("duplicate_de", name=name)
        contrasts[name] = Path(path)
    return contrasts


def read_de(path, cols):
    """Step 2: one DE table and its count of rows dropped for an NA in `cols`."""
    if not path.is_file():
        fail("missing_file", path=path)
    de = pd.read_csv(path, sep="\t", dtype={"gene_symbol": str})
    for col in ["gene_symbol", *cols]:
        if col not in de.columns:
            fail("missing_column", col=col, path=path)
    # Gene sets name genes by symbol, so the symbol is the join key.
    de = de.rename(columns={"gene_symbol": "gene"})
    for col in cols:
        values = pd.to_numeric(de[col], errors="coerce")
        if (values.isna() & de[col].notna()).any() or np.isinf(values).any():
            fail("not_numeric", col=col, path=path)
        de[col] = values
    skillkit.check_unique(de["gene"], f"gene names in {path}")
    # Rows the engine cannot use (NA rank; NA padj or log2_fc) leave the universe.
    kept = de[["gene", *cols]].dropna()
    return kept, len(de) - len(kept)


def read_gene_sets(path):
    """Step 3: long (set, gene) table from a GMT or a TSV."""
    if not path.is_file():
        fail("missing_file", path=path)
    if path.suffix == ".gmt":
        rows = [line.split("\t") for line in path.read_text().splitlines() if line]
        skillkit.check_unique([row[0] for row in rows], f"gene set names in {path}")
        sets = pd.DataFrame(
            [(row[0], gene) for row in rows for gene in row[2:] if gene],
            columns=["set", "gene"],
        )
    elif path.suffix == ".tsv":
        sets = pd.read_csv(path, sep="\t", dtype=str)
        if not {"set", "gene"} <= set(sets.columns):
            fail("bad_gene_sets", path=path)
    else:
        fail("bad_gene_sets", path=path)
    return sets[["set", "gene"]].drop_duplicates().sort_values(["set", "gene"])


def sets_in(sets, genes, p, contrast):
    """Step 4: sets cut to `genes`, kept when --min-size <= size <= --max-size."""
    inside = sets[sets["gene"].isin(genes)]
    if inside.empty:
        fail("no_overlap", contrast=contrast)
    size = inside.groupby("set")["gene"].transform("size")
    kept = inside[size.between(p["min_size"], p["max_size"])]
    if kept.empty:
        fail(
            "no_sets", min_size=p["min_size"], max_size=p["max_size"], contrast=contrast
        )
    return kept


def run_ora(tested, sets, p, contrast):
    """Step 5 (ora): one-sided Fisher (hypergeometric upper tail) per set."""
    universe = tested["gene"]
    if p["universe"] == "tested-in-sets":
        universe = universe[universe.isin(sets["gene"])]
    kept = sets_in(sets, universe, p, contrast)
    lfc = tested["log2_fc"]
    changed = {
        "up": lfc > p["lfc_min"],
        "down": lfc < -p["lfc_min"],
        "both": lfc.abs() > p["lfc_min"],
    }[p["direction"]]
    foreground = set(tested.loc[(tested["padj"] < p["padj_max"]) & changed, "gene"])
    foreground &= set(universe)
    # An empty foreground gives overlap 0 and p = 1 for every set, not an error:
    # a null contrast beside others is a result. analyse() fails if all are empty.
    rows = []
    for name, genes in kept.groupby("set")["gene"]:
        hits = sorted(set(genes) & foreground)
        k, size, n_fg, n = len(hits), len(genes), len(foreground), len(universe)
        table = [[k, size - k], [n_fg - k, n - size - n_fg + k]]
        # scipy's statistic is the sample odds ratio; its p is phyper's upper tail.
        odds_ratio, pvalue = stats.fisher_exact(table, alternative="greater")
        rows.append(
            {
                "contrast": contrast,
                "set": name,
                "size": size,
                "overlap": k,
                "odds_ratio": odds_ratio,
                "pvalue": pvalue,
                "leading_edge": ",".join(hits),
            }
        )
    return pd.DataFrame(rows), {"n_universe": len(universe), "n_foreground": n_fg}


def run_fgsea(ranks, sets, p):
    """Step 5 (fgsea): fgsea.R on every contrast, seeded, BPPARAM = SerialParam()."""
    rscript = os.environ.get("RSCRIPT") or shutil.which("Rscript")
    if not rscript or not Path(rscript).is_file():
        fail("no_rscript")
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        ranks.to_csv(tmp / "ranks.tsv", sep="\t", index=False)
        sets.to_csv(tmp / "sets.tsv", sep="\t", index=False)
        cmd = [
            rscript,
            str(FGSEA_R),
            str(tmp / "ranks.tsv"),
            str(tmp / "sets.tsv"),
            *(str(p[k]) for k in ("min_size", "max_size", "n_perm_simple", "seed")),
            str(tmp / "out.tsv"),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode == 3:
            fail("missing_r_package", rscript=rscript)
        if proc.returncode != 0:
            tail = " | ".join(proc.stderr.strip().splitlines()[-5:])
            fail("r_failed", code=proc.returncode, tail=tail)
        res = pd.read_csv(
            tmp / "out.tsv", sep="\t", dtype={"contrast": str, "set": str}
        )
    return res, [f"Rscript {rscript}", *proc.stdout.strip().splitlines()]


def bh(pvalues):
    """Step 6: BH over the sets tested in one contrast; NA stays NA."""
    padj = pd.Series(np.nan, index=pvalues.index)
    tested = pvalues.dropna()
    padj[tested.index] = stats.false_discovery_control(tested, method="bh")
    return padj


def with_status(table, sets, tables, names):
    """Step 7b: every set in the file per contrast; untested sets keep NA stats.

    A set outside --min-size/--max-size in a contrast (too few or too many of its
    genes found there) is `not_tested`; `size` still says how many were found.
    """
    every = sets["set"].unique()
    parts = []
    for name in names:
        tested = table[table["contrast"] == name]
        reasons = dict.fromkeys(set(every) - set(tested["set"]), "not_tested")
        part = skillkit.with_status(tested, every, "set", reasons)
        found = sets[sets["gene"].isin(tables[name]["gene"])].groupby("set").size()
        part["contrast"] = name
        part["size"] = part["size"].fillna(part["set"].map(found).fillna(0))
        parts.append(part)
    full = pd.concat(parts, ignore_index=True)
    full["set_size"] = full["set"].map(sets.groupby("set").size())
    full["contrast"] = pd.Categorical(full["contrast"], categories=names)
    full = full.sort_values(["contrast", "padj", "set"], na_position="last")
    return full.assign(
        contrast=full["contrast"].astype(str), size=full["size"].astype(int)
    )


def running_es(ranks, members):
    """GSEA running score (weight 1, as fgsea::plotEnrichment) along the ranking."""
    ranked = ranks.sort_values(["rank", "gene"], ascending=[False, True])
    hit = ranked["gene"].isin(members).to_numpy()
    weight = np.where(hit, ranked["rank"].abs(), 0)
    es = np.cumsum(weight) / weight.sum() - np.cumsum(~hit) / (~hit).sum()
    return ranked.assign(position=np.arange(1, len(ranked) + 1), hit=hit, running_es=es)


def draw_figures(figures, table, p, contrasts, effect, label, ranks, sets):
    """Step 8: dotplot, heatmap and (fgsea) running-ES curves; returns captions."""
    best = table.groupby("set")["padj"].min().reset_index()
    top = best.sort_values(["padj", "set"], na_position="last")["set"].head(N_SETS)
    shown = table[table["set"].isin(top)].assign(
        neg_log10_padj=-np.log10(table["padj"].clip(lower=1e-300)),
    )
    shown["set"] = pd.Categorical(shown["set"], categories=top.tolist())
    shown["contrast"] = pd.Categorical(shown["contrast"], categories=contrasts)
    shown = shown.sort_values(["set", "contrast"])
    rule = f"padj < {p['padj_max_sets']}"
    order = "Order fixed, not clustered: sets by lowest padj, contrasts as given."
    height = 25 + 4 * len(top)

    # Area is proportional to -log10 padj (scale_size_area), with the cutoff as floor
    # of the scale so an all-null run does not blow up tiny dots.
    peak = max(shown["neg_log10_padj"].max(), -np.log10(p["padj_max_sets"]))
    max_area = (MAX_DOT_MM * 72.27 / 25.4) ** 2
    fig, ax = plt.subplots(layout="constrained")
    dots = ax.scatter(
        shown["contrast"].cat.codes,
        shown["set"].cat.codes,
        s=shown["neg_log10_padj"] / peak * max_area,
        c=shown[effect].clip(-p["limit"], p["limit"]),
        cmap=diverging,
        vmin=-p["limit"],
        vmax=p["limit"],
        edgecolors=np.where(shown["significant"], "black", "none"),
        linewidths=0.5,
    )
    ax.set_xticks(range(len(contrasts)), contrasts, rotation=90)
    ax.set_yticks(range(len(top)), top)
    ax.set_xlim(-0.5, len(contrasts) - 0.5)
    ax.set_ylim(len(top) - 0.5, -0.5)
    ax.set_xlabel("contrast")
    fig.colorbar(dots, ax=ax, label=label)
    handles = [
        ax.scatter([], [], s=v / peak * max_area, c="grey", label=f"{v:.1f}")
        for v in (-np.log10(p["padj_max_sets"]), peak)
    ]
    ax.legend(
        handles=handles, title="-log10 padj", loc="lower left", bbox_to_anchor=(1, 0)
    )
    save_panel(fig, shown, figures, "dotplot", 120, height)
    plt.close(fig)

    heat = shown[["set", "contrast", effect, "significant"]].rename(
        columns={effect: label}
    )
    stars = heat.loc[heat["significant"], ["set", "contrast"]].astype(str)
    stars = stars.rename(columns={"set": "row", "contrast": "col"}).assign(star="*")
    plot = plot_heatmap(heat, "set", "contrast", label, limit=p["limit"]) + geom_text(
        aes("col", "row", label="star"), data=stars, inherit_aes=False
    )
    save_panel(plot, heat, figures, "heatmap", 89 + 6 * len(contrasts), height)

    captions = {
        "dotplot": f"Dot plot of {label} for the {len(top)} gene sets with the lowest "
        "adjusted p in any contrast; dot area is -log10 adjusted p, colour saturates "
        f"at ±{p['limit']:g}, black outline marks {rule} (BH within each contrast). "
        f"{order}",
        "heatmap": f"Heatmap of {label} for the same {len(top)} sets by contrast; "
        f"colour saturates at ±{p['limit']:g}, * marks {rule}, blank means the set "
        f"was outside --min-size/--max-size in that contrast. {order}",
    }
    if ranks is None:
        return captions

    tops = table.groupby("contrast", sort=False).head(1).set_index("contrast")["set"]
    curves = pd.concat(
        [
            running_es(
                ranks[ranks["contrast"] == c], sets.loc[sets["set"] == tops[c], "gene"]
            ).assign(set=tops[c])
            for c in contrasts
        ]
    )
    fig, axes = plt.subplots(len(contrasts), 1, layout="constrained", squeeze=False)
    for ax, contrast in zip(axes[:, 0], contrasts, strict=True):
        curve = curves[curves["contrast"] == contrast]
        ax.plot(curve["position"], curve["running_es"], lw=0.8, c=CURVE_COLOUR)
        ax.axhline(0, lw=0.5, c="black")
        # Hit ticks along the bottom 8% of the panel, as in fgsea::plotEnrichment.
        ax.vlines(
            curve.loc[curve["hit"], "position"],
            0,
            0.08,
            transform=ax.get_xaxis_transform(),
            lw=0.3,
            colors="black",
        )
        ax.set_ylabel(f"running ES\n{contrast}")
    axes[-1, 0].set_xlabel(f"gene rank by {p['rank_col']} (high to low)")
    columns = ["contrast", "set", "position", "gene", "rank", "hit", "running_es"]
    save_panel(
        fig, curves[columns], figures, "running_es", 89, 15 + 35 * len(contrasts)
    )
    plt.close(fig)
    pairs = ", ".join(f"{c}: {tops[c]}" for c in contrasts)
    captions["running_es"] = (
        "Running enrichment score along the gene ranking for the top set of each "
        f"contrast ({pairs}); ticks mark the set's genes, the peak is the ES that "
        "fgsea normalises to NES."
    )
    return captions


def analyse(argv, argv_as_run):
    p = parse_args(argv)
    contrasts = check_values(p)
    skillkit.seed_everything(p["seed"])
    cols = [p["rank_col"]] if p["method"] == "fgsea" else ["padj", "log2_fc"]
    read = {name: read_de(path, cols) for name, path in contrasts.items()}
    tables = {name: de for name, (de, _) in read.items()}
    dropped = {name: n for name, (_, n) in read.items()}
    sets = read_gene_sets(Path(p["gene_sets"]))
    set_size = sets.groupby("set").size()
    n_in_resource = sets["gene"].nunique()

    stats_by = {}
    if p["method"] == "fgsea":
        for name, de in tables.items():
            stats_by[name] = {
                "n_genes": len(de),
                "n_sets_tested": sets_in(sets, de["gene"], p, name)["set"].nunique(),
            }
        ranks = pd.concat(tables, names=["contrast", None]).reset_index(level=0)
        ranks = ranks.rename(columns={p["rank_col"]: "rank"})
        table, r_versions = run_fgsea(ranks, sets, p)
        effect, label = "nes", "NES"
    else:
        parts = []
        for name, de in tables.items():
            part, counts = run_ora(de, sets, p, name)
            parts.append(part)
            stats_by[name] = {
                "n_genes": len(de),
                "n_sets_tested": len(part),
                **counts,
            }
        if not any(c["n_foreground"] for c in stats_by.values()):
            fail(
                "no_foreground",
                direction=p["direction"],
                padj_max=p["padj_max"],
                lfc_min=p["lfc_min"],
            )
        table, r_versions, ranks = pd.concat(parts), [], None
        with np.errstate(divide="ignore"):
            table["log2_odds_ratio"] = np.log2(table["odds_ratio"])
        effect, label = "log2_odds_ratio", "log2 odds ratio"

    table["padj"] = table.groupby("contrast")["pvalue"].transform(bh)
    table["significant"] = table["padj"] < p["padj_max_sets"]
    names = list(contrasts)
    table["contrast"] = pd.Categorical(table["contrast"], categories=names)
    table = table.sort_values(["contrast", "padj", "set"], na_position="last")
    table["contrast"] = table["contrast"].astype(str)
    table["set_size"] = table["set"].map(set_size)
    columns = ["contrast", "set", "set_size", "size"]
    columns += ["nes"] if p["method"] == "fgsea" else ["overlap", "odds_ratio"]
    columns += ["pvalue", "padj", "significant", "leading_edge", "status"]
    full = with_status(table, sets, tables, names)

    outdir = skillkit.make_outdir(p["outdir"])
    full[columns].to_csv(
        outdir / "tables" / "enrichment.tsv",
        sep="\t",
        index=False,
        float_format="%.6g",
    )
    plt.switch_backend("Agg")
    captions = draw_figures(
        outdir / "figures", table, p, names, effect, label, ranks, sets
    )

    per_contrast = []
    for name in names:
        rows = table[(table["contrast"] == name) & table["significant"]]
        per_contrast.append(
            {
                "contrast": name,
                **stats_by[name],
                # Resource overlap: distinct set genes found in this DE table.
                "n_matched": int(
                    sets.loc[sets["gene"].isin(tables[name]["gene"]), "gene"].nunique()
                ),
                "n_in_resource": int(n_in_resource),
                "n_significant": len(rows),
                "significant_sets": dict(
                    zip(rows["set"], rows[effect].round(3), strict=True)
                ),
            }
        )
    tested = (
        f"ranked by {p['rank_col']}"
        if p["method"] == "fgsea"
        else f"{p['direction']} genes with padj < {p['padj_max']} and |log2_fc| > "
        f"{p['lfc_min']} against universe {p['universe']}"
    )
    summary = [
        f"{p['method']} on {len(names)} contrasts ({', '.join(names)}), {tested}; "
        f"{sets['set'].nunique()} sets in {Path(p['gene_sets']).name}, those with "
        f"{p['min_size']}-{p['max_size']} genes tested, the rest kept as not_tested. "
        "Set genes found in each DE "
        "table (a missing gene was not tested): "
        + "; ".join(
            f"{c['contrast']} {c['n_matched']} of {c['n_in_resource']}, "
            f"{dropped[c['contrast']]} NA rows dropped"
            for c in per_contrast
        )
        + ".",
        f"Significant sets (padj < {p['padj_max_sets']}, BH within each contrast): "
        + "; ".join(
            f"{c['contrast']} {c['n_significant']} of {c['n_sets_tested']}"
            for c in per_contrast
        )
        + ".",
        f"Top ({label}): "
        + "; ".join(
            f"{c['contrast']} "
            + (
                ", ".join(f"{s} ({v:+.2f})" for s, v in c["significant_sets"].items())
                or "none"
            )
            for c in per_contrast
        ),
    ]
    keys = [*REQUIRED, *METHOD_ARGS[p["method"]]]
    title = f"Gene-set enrichment: {p['method']}, {', '.join(names)}"
    skillkit.write_report(outdir, title, p, keys, summary, captions)
    result = {
        "method": p["method"],
        "n_contrasts": len(names),
        "n_sets_in_file": int(sets["set"].nunique()),
        "contrasts": per_contrast,
        "dropped_per_group": dropped,
        "status_counts": full["status"].value_counts().to_dict(),
        "params": p,
        "seed": p["seed"],
    }
    artifacts = [
        "report.md",
        "tables/enrichment.tsv",
        *[f"figures/{n}.png" for n in captions],
    ]
    skillkit.write_result(
        outdir,
        result,
        summary,
        artifacts,
        context=p["context"],
        reference=p["reference"],
    )
    skillkit.write_reproducibility(
        outdir, __file__, p, keys, argv_as_run, PACKAGES, r_versions
    )
    return outdir


if __name__ == "__main__":
    raise SystemExit(skillkit.run("enrichment", analyse))
