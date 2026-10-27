"""Purpose: one test per feature on sample-level values, corrected across features.
Inputs:  --input <long.tsv> with feature, subject, value and a group or covariate column
Outputs: <outdir>/report.md, result.json, tables/tests.tsv, figures/, reproducibility/
Run:     python scripts/group_tests.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skillkit
from plots import FEW_POINTS, plot_box, plot_scatter_fit
from scipy import stats
from theme import save_panel

SKILL_DIR = Path(__file__).resolve().parents[1]
GROUP_TESTS = ("wilcoxon", "paired-wilcoxon")
REQUIRED = ["input", "test", "value_unit", "adjust", "alpha", "n_plot", "seed"]
OPTIONAL = ["group_col", "levels", "pair_col", "covariate"]
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "input": str(SKILL_DIR / "examples" / "demo.tsv"),
    "test": "paired-wilcoxon",
    "group_col": "visit",
    "levels": "pre,post",
    "pair_col": "donor",
    "covariate": None,
    "value_unit": "a.u.",
    "adjust": "bh",
    "alpha": 0.05,
    "n_plot": 4,
    "seed": 0,
}

ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "needs_args": "--test {test} needs {args}",
    "unused_args": "--test {test} does not use {args}",
    "bad_levels": "--levels must be two different values a,b, got {value!r}",
    "missing_file": "input file not found: {path}",
    "missing_column": "input column not found: {col}",
    "not_numeric": "column {col!r} must be numeric",
    "level_absent": "level {level!r} not found in column {col!r}",
    "duplicate_rows": "{n} rows repeat a ({keys}) combination; one value each",
    "no_tests": "no feature could be tested; check the levels, pairs or covariate",
}

# Fixed plot spec, not analysis thresholds: corpus jitter width (position_jitter(0.2)),
# reference blue and level-b red from palettes.pal_condition.
JITTER_WIDTH = 0.2
LEVEL_COLOURS = ("#4477AA", "#EE6677")
PANEL_MM = 60


def fail(key, **fields):
    skillkit.fail(ERRORS, key, **fields)


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input")
    parser.add_argument("--test", choices=[*GROUP_TESTS, "spearman", "pearson"])
    parser.add_argument("--group-col")
    parser.add_argument("--levels")
    parser.add_argument("--pair-col")
    parser.add_argument("--covariate")
    parser.add_argument("--value-unit")
    parser.add_argument("--adjust", choices=["bh", "bonferroni"])
    parser.add_argument("--alpha", type=float)
    parser.add_argument("--n-plot", type=int)
    parser.add_argument("--seed", type=int)
    return skillkit.parse_required(parser, REQUIRED, DEMO, argv, optional=OPTIONAL)


def check_values(p):
    rules = [
        ("--alpha", p["alpha"], 0 < p["alpha"] < 1, "in (0, 1)"),
        ("--n-plot", p["n_plot"], p["n_plot"] >= 1, "an integer >= 1"),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ]
    for arg, value, ok, rule in rules:
        if not ok:
            fail("bad_value", arg=arg, rule=rule, value=value)
    needs = {
        "wilcoxon": ["group_col", "levels"],
        "paired-wilcoxon": ["group_col", "levels", "pair_col"],
    }.get(p["test"], ["covariate"])
    missing = [k for k in needs if p[k] is None]
    if missing:
        fail("needs_args", test=p["test"], args=" ".join(map(skillkit.flag, missing)))
    unused = [k for k in OPTIONAL if k not in needs and p[k] is not None]
    if unused:
        fail("unused_args", test=p["test"], args=" ".join(map(skillkit.flag, unused)))
    if p["levels"] is None:
        return None
    levels = [lv.strip() for lv in p["levels"].split(",")]
    if len(levels) != 2 or not all(levels) or levels[0] == levels[1]:
        fail("bad_levels", value=p["levels"])
    return levels


def read_values(p, levels):
    """Steps 1-2: read, keep the rows the test uses, refuse repeated keys."""
    path = Path(p["input"])
    if not path.is_file():
        fail("missing_file", path=path)
    # Labels stay strings ("01" is not 1); value and covariate are parsed below.
    data = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    x_col = p["covariate"] or p["group_col"]
    cols = [
        c
        for c in dict.fromkeys(["feature", "subject", "value", x_col, p["pair_col"]])
        if c
    ]
    for col in cols:
        if col not in data.columns:
            fail("missing_column", col=col)
    for col in [c for c in ["value", p["covariate"]] if c]:
        parsed = pd.to_numeric(data[col], errors="coerce")
        if (parsed.isna() & data[col].notna()).any():
            fail("not_numeric", col=col)
        data[col] = parsed

    n_rows = len(data)
    # Drops per group: rows with a missing field, counted by their group level.
    missing = data[data[cols].isna().any(axis=1)]
    dropped = (
        {lv: int((missing[x_col] == lv).sum()) for lv in levels}
        if levels
        else {"all": len(missing)}
    )
    data = data.dropna(subset=cols)
    if levels:
        for level in levels:
            if level not in set(data[x_col]):
                fail("level_absent", level=level, col=x_col)
        data = data[data[x_col].isin(levels)]
    for keys in (["feature", "subject"], ["feature", x_col, p["pair_col"]]):
        if None in keys:
            continue
        repeated = int(data.duplicated(keys, keep=False).sum())
        if repeated:
            fail("duplicate_rows", n=repeated, keys=", ".join(keys))
    return data[cols], n_rows - len(data), dropped


def hodges_lehmann(b, a):
    """Median of all b - a differences: the shift the Mann-Whitney test targets."""
    return float(np.median(np.subtract.outer(b, a)))


def walsh_median(d):
    """Median of the Walsh averages of paired differences (signed-rank target)."""
    i, j = np.triu_indices(len(d))
    return float(np.median((d[i] + d[j]) / 2))


def test_feature(rows, p, levels):
    """Step 3: one row of tests.tsv; pvalue NA when the test cannot run."""
    if p["test"] in ("spearman", "pearson"):
        x, y = rows[p["covariate"]].to_numpy(), rows["value"].to_numpy()
        name = "rho" if p["test"] == "spearman" else "r"
        out = {"n": len(rows), name: np.nan, "pvalue": np.nan}
        # A constant column has no correlation; scipy would warn and return NaN.
        if len(rows) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
            return out
        fit = (stats.spearmanr if p["test"] == "spearman" else stats.pearsonr)(x, y)
        return out | {name: float(fit.statistic), "pvalue": float(fit.pvalue)}

    a_lv, b_lv = levels
    if p["test"] == "paired-wilcoxon":
        # Complete pairs only, as in the corpus (drop_na before wilcox.test).
        wide = rows.pivot(index=p["pair_col"], columns=p["group_col"], values="value")
        wide = wide.reindex(columns=levels).dropna()
        a, b = wide[a_lv].to_numpy(), wide[b_lv].to_numpy()
        out = {"n_pairs": len(wide)}
    else:
        a = rows.loc[rows[p["group_col"]] == a_lv, "value"].to_numpy()
        b = rows.loc[rows[p["group_col"]] == b_lv, "value"].to_numpy()
        out = {f"n_{a_lv}": len(a), f"n_{b_lv}": len(b)}
    out |= {
        f"median_{a_lv}": np.median(a) if len(a) else np.nan,
        f"median_{b_lv}": np.median(b) if len(b) else np.nan,
        "hodges_lehmann": np.nan,
        "statistic": np.nan,
        "pvalue": np.nan,
    }
    if p["test"] == "paired-wilcoxon":
        d = b - a
        # zero_method="wilcox" drops zero differences, as R's wilcox.test does.
        if not np.any(d != 0):
            return out
        res = stats.wilcoxon(d, zero_method="wilcox", method="auto")
        effect = walsh_median(d)
    else:
        if not len(a) or not len(b):
            return out
        res = stats.mannwhitneyu(b, a, alternative="two-sided", method="auto")
        effect = hodges_lehmann(b, a)
    return out | {
        "hodges_lehmann": effect,
        "statistic": float(res.statistic),
        "pvalue": float(res.pvalue),
    }


def adjust(pvalues, method):
    """Step 4: correct across the tested features only; NA stays NA."""
    padj = pd.Series(np.nan, index=pvalues.index)
    tested = pvalues.dropna()
    if method == "bh":
        padj[tested.index] = stats.false_discovery_control(tested, method="bh")
    else:
        padj[tested.index] = np.minimum(tested * len(tested), 1)
    return padj


def draw_figures(figures, data, table, p, levels):
    """Step 6: one panel per top feature; returns {name: caption}."""
    top = table.dropna(subset=["padj"]).head(p["n_plot"])
    family = f"{p['adjust'].upper()} across {table['pvalue'].notna().sum()} features"
    captions = {}
    for rank, row in enumerate(top.itertuples(index=False), start=1):
        name = f"top{rank:02d}"
        rows = data[data["feature"] == row.feature].rename(
            columns={"value": row.feature}
        )
        if levels:
            pal = dict(zip(levels, LEVEL_COLOURS, strict=True))
            pair = p["pair_col"]
            if pair:
                # Complete pairs only, matching the test.
                complete = rows.groupby(pair)[p["group_col"]].transform("nunique") == 2
                rows = rows[complete]
            fig = plot_box(
                rows,
                p["group_col"],
                row.feature,
                pal,
                jitter_width=JITTER_WIDTH,
                seed=p["seed"],
                pair=pair,
                value_label=f"{row.feature} ({p['value_unit']})",
            )
            shown = (
                "lines join the two samples of one pair"
                if pair
                else "points are jittered sideways"
            )
            few = rows[p["group_col"]].value_counts().min() < FEW_POINTS
            form = "Dot plot with medians" if few else "Box plot"
            caption = (
                f"{form} of {row.feature} per {p['group_col']} "
                f"({levels[1]} vs {levels[0]}); each point is one subject, {shown}. "
                f"Hodges-Lehmann shift {row.hodges_lehmann:.3g}, "
                f"padj {row.padj:.2g} ({family})."
            )
        else:
            rows = rows.sort_values([p["covariate"], "subject"])
            fig = plot_scatter_fit(rows, p["covariate"], row.feature)
            fig.axes[0].set_ylabel(f"{row.feature} ({p['value_unit']})")
            estimate = "rho" if p["test"] == "spearman" else "r"
            caption = (
                f"Scatter plot of {row.feature} against {p['covariate']}; each point "
                f"is one subject, the line is a least-squares fit with its 95% "
                f"confidence band. {p['test'].title()} {estimate} "
                f"{getattr(row, estimate):.2f}, padj {row.padj:.2g} ({family})."
            )
        fig.set_layout_engine("constrained")
        save_panel(fig, rows, figures, name, PANEL_MM, PANEL_MM)
        plt.close(fig)
        captions[name] = caption
    return captions


def analyse(argv, argv_as_run):
    p = parse_args(argv)
    levels = check_values(p)
    skillkit.seed_everything(p["seed"])
    data, n_dropped, dropped = read_values(p, levels)

    table = pd.DataFrame(
        [
            {"feature": feature, **test_feature(rows, p, levels)}
            for feature, rows in data.groupby("feature", sort=True)
        ]
    )
    n_tested = int(table["pvalue"].notna().sum())
    if not n_tested:
        fail("no_tests")
    table["padj"] = adjust(table["pvalue"], p["adjust"])
    table["significant"] = table["padj"] < p["alpha"]
    table = table.sort_values(["padj", "feature"], na_position="last")

    outdir = skillkit.make_outdir(p["outdir"])
    table.to_csv(
        outdir / "tables" / "tests.tsv", sep="\t", index=False, float_format="%.6g"
    )
    plt.switch_backend("Agg")
    captions = draw_figures(outdir / "figures", data, table, p, levels)

    n_sig = int(table["significant"].sum())
    contrast = f"{levels[1]} vs {levels[0]}" if levels else f"vs {p['covariate']}"
    sig = table.loc[table["significant"], "feature"].tolist()
    summary = [
        f"{n_tested} of {len(table)} features tested ({p['test']}, {contrast}); "
        f"{n_dropped} input rows dropped (missing value or other level); "
        "missing values per group: "
        + ", ".join(f"{k} {v}" for k, v in dropped.items())
        + ".",
        f"{n_sig} features with padj < {p['alpha']} "
        f"({p['adjust']} across the {n_tested} tested features).",
        f"Significant: {', '.join(sig[:10]) or 'none'}"
        + (f" and {len(sig) - 10} more." if len(sig) > 10 else "."),
    ]
    keys = [*REQUIRED, *OPTIONAL]
    title = f"Group tests: {p['test']}, {contrast}"
    skillkit.write_report(outdir, title, p, keys, summary, captions)
    result = {
        "test": p["test"],
        "n_features": len(table),
        "n_tested": n_tested,
        "n_significant": n_sig,
        "n_rows_dropped": n_dropped,
        "n_missing_dropped_by_group": dropped,
        "significant_features": sig,
        "params": p,
        "seed": p["seed"],
    }
    artifacts = [
        "report.md",
        "tables/tests.tsv",
        *[f"figures/{n}.png" for n in captions],
    ]
    skillkit.write_result(outdir, result, summary, artifacts)
    packages = ["numpy", "pandas", "scipy", "matplotlib"]
    skillkit.write_reproducibility(outdir, __file__, p, keys, argv_as_run, packages)
    return outdir


if __name__ == "__main__":
    raise SystemExit(skillkit.run("group-tests", analyse))
