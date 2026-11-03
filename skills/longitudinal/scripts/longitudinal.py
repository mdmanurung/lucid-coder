"""Purpose: one linear mixed model per feature over visits (lmerTest via Rscript).
Inputs:  --input <long.tsv> with feature, subject, time, value and covariate columns
Outputs: <outdir>/report.md, result.json, tables/model_results.tsv, figures/,
         reproducibility/
Run:     python scripts/longitudinal.py --demo --outdir demo_out
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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skillkit
from scipy import stats
from theme import save_panel

SKILL_DIR = Path(__file__).resolve().parents[1]
ENGINE = Path(__file__).resolve().parent / "lmer.R"
REQUIRED = [
    "context",
    "input",
    "formula",
    "term",
    "time_type",
    "reference_time",
    "value_unit",
    "time_unit",
    "adjust",
    "padj_max",
    "n_plot",
    "seed",
]
# Asked only when they apply: --knots for spline time, --missing for a covariate NA.
OPTIONAL = ["knots", "missing"]
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "context": "synthetic demo: 24 subjects, 3 visits, 30 features, no real tissue",
    "input": str(SKILL_DIR / "examples" / "demo.tsv"),
    "formula": "value ~ time + age + (1|subject)",
    "term": "time3",
    "time_type": "factor",
    "reference_time": "1",
    "value_unit": "a.u.",
    "time_unit": "visit",
    "adjust": "bh",
    "padj_max": 0.05,
    "n_plot": 4,
    "seed": 0,
    "knots": None,
    "missing": None,
}

ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "bad_formula": "--formula must be 'value ~ time + covariates + (1|subject)' with "
    "column names, ':' or '*' interactions and no other random term, got {value!r}",
    "no_subject_term": "--formula needs the random intercept (1|subject) once, "
    "got {value!r}",
    "no_time_term": "--formula must have time as a fixed term, got {value!r}",
    "missing_file": "input file not found: {path}",
    "missing_column": "input column not found: {col}",
    "not_numeric": "column {col!r} must be numeric",
    "reference_absent": "reference time {value!r} not found in column 'time'",
    "duplicate_rows": "{n} rows repeat a (feature, subject, time) combination; "
    "one value each",
    "no_rscript": "Rscript not found; set RSCRIPT or put Rscript on PATH",
    "missing_r_package": "{detail}; ask for it to be installed, do not install it here",
    "term_absent": "--term {term!r} is not a fixed-effect coefficient; "
    "choose one of: {choices}",
    "r_failed": "lmer.R failed with exit {code}: {tail}",
    "no_tests": "no feature could be fitted; first lmer message: {message}",
    "knots_needed": "--time-type spline needs --knots",
    "knots_unused": "--knots is used only with --time-type spline",
    "missing_needed": "formula covariates have missing values ({detail}); "
    "give --missing complete-case or impute",
    "missing_unused": "--missing is used only when a formula covariate has a "
    "missing value",
    "impute_varies": "--missing impute needs covariates fixed within subject; "
    "varies: {cols}",
    "impute_joint": "--missing impute pools one coefficient; the spline joint test "
    "of --term time needs --missing complete-case",
}

NAME = r"[A-Za-z_][\w.]*"
FIXED = rf"{NAME}(?:\s*[:*]\s*{NAME})*"
SUBJECT = r"\(\s*1\s*\|\s*subject\s*\)"
TERM = rf"(?:{FIXED}|{SUBJECT})"
FORMULA = re.compile(rf"\s*value\s*~\s*{TERM}(?:\s*\+\s*{TERM})*\s*")

# Fixed plot spec, not analysis thresholds. save_panel saves without a tight
# bbox, so every figure sets a layout engine.
N_FOREST = 20
MEAN_COLOUR, SUBJECT_COLOUR, OTHER = "#B2182B", "0.6", "0.4"
PANEL_MM = 60


def fail(key, **fields):
    skillkit.fail(ERRORS, key, **fields)


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--context")
    parser.add_argument("--input")
    parser.add_argument("--formula")
    parser.add_argument("--term")
    parser.add_argument("--time-type", choices=["factor", "numeric", "spline"])
    parser.add_argument("--reference-time")
    parser.add_argument("--value-unit")
    parser.add_argument("--time-unit")
    parser.add_argument("--adjust", choices=["bh", "bonferroni"])
    parser.add_argument("--padj-max", type=float)
    parser.add_argument("--n-plot", type=int)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--knots", type=int)
    parser.add_argument("--missing", choices=["complete-case", "impute"])
    return skillkit.parse_required(parser, REQUIRED, DEMO, argv, optional=OPTIONAL)


def check_values(p):
    """Step 1: argument values and the formula; returns the formula's columns."""
    numeric_ref = p["time_type"] == "factor" or is_number(p["reference_time"])
    rules = [
        ("--padj-max", p["padj_max"], 0 < p["padj_max"] <= 1, "in (0, 1]"),
        ("--n-plot", p["n_plot"], p["n_plot"] >= 1, "an integer >= 1"),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
        (
            "--reference-time",
            p["reference_time"],
            numeric_ref,
            f"a number for --time-type {p['time_type']}",
        ),
    ]
    for arg, value, ok, rule in rules:
        if not ok:
            fail("bad_value", arg=arg, rule=rule, value=value)
    spline = p["time_type"] == "spline"
    if spline and p["knots"] is None:
        fail("knots_needed")
    if not spline and p["knots"] is not None:
        fail("knots_unused")
    if spline and p["term"] == "time" and p["missing"] == "impute":
        fail("impute_joint")
    formula = p["formula"]
    if len(re.findall(SUBJECT, formula)) != 1:
        fail("no_subject_term", value=formula)
    if not FORMULA.fullmatch(formula):
        fail("bad_formula", value=formula)
    names = list(dict.fromkeys(re.findall(NAME, formula)))
    if "time" not in names:
        fail("no_time_term", value=formula)
    return names


def is_number(text):
    return pd.notna(pd.to_numeric(text, errors="coerce"))


def read_values(p, names):
    """Step 2: read, keep the formula's columns, drop missing rows, refuse repeats.

    A covariate NA is dropped with its row (complete-case) or kept for mice (impute).
    """
    path = Path(p["input"])
    if not path.is_file():
        fail("missing_file", path=path)
    # Labels stay strings ("01" is not 1); value (and numeric time) parsed below.
    data = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    cols = list(dict.fromkeys(["feature", "subject", "time", "value", *names]))
    for col in cols:
        if col not in data.columns:
            fail("missing_column", col=col)
    for col in ["value"] if p["time_type"] == "factor" else ["value", "time"]:
        parsed = pd.to_numeric(data[col], errors="coerce")
        if (parsed.isna() & data[col].notna()).any():
            fail("not_numeric", col=col)
        data[col] = parsed
    n_rows = len(data)
    covariates = [c for c in cols if c not in ("feature", "subject", "time", "value")]
    has_na = data[covariates].isna()
    if has_na.any().any() and p["missing"] is None:
        detail = ", ".join(f"{c}: {n} rows" for c, n in has_na.sum().items() if n)
        fail("missing_needed", detail=detail)
    if not has_na.any().any() and p["missing"] is not None:
        fail("missing_unused")
    imputed = has_na.columns[has_na.any()].tolist()
    if p["missing"] == "impute":
        per_subject = data.groupby("subject")[imputed].nunique(dropna=False)
        varies = per_subject.columns[(per_subject > 1).any()].tolist()
        if varies:
            fail("impute_varies", cols=", ".join(varies))
    na_rows = data[imputed].isna().any(axis=1)
    imputation = {
        "missing": p["missing"],
        "covariates": imputed,
        "n_rows": int(na_rows.sum()),
        "n_subjects": int(data.loc[na_rows, "subject"].nunique()),
    }
    keep = cols if p["missing"] != "impute" else ["feature", "subject", "time", "value"]
    # Drops per group: rows with a missing field, counted by their visit.
    missing = data.loc[data[keep].isna().any(axis=1), "time"].dropna()
    dropped = {
        str(t): int((missing == t).sum())
        for t in sorted(data["time"].dropna().unique())
    }
    data = data.dropna(subset=keep)[cols]
    if p["time_type"] == "factor" and p["reference_time"] not in set(data["time"]):
        fail("reference_absent", value=p["reference_time"])
    repeated = int(data.duplicated(["feature", "subject", "time"], keep=False).sum())
    if repeated:
        fail("duplicate_rows", n=repeated)
    check_rank(data.dropna(), names, p)
    return data, n_rows - len(data), dropped, imputation


def spline_knots(time, p):
    """Interior knots at time quantiles (as ns() places them) and the time range.

    Fixed once across features, so every feature is fitted on the same basis.
    """
    k = p["knots"]
    interior = np.quantile(time, np.linspace(0, 1, k + 2))[1:-1].tolist()
    boundary = [float(time.min()), float(time.max())]
    n_times = time.nunique()
    ok = (
        1 <= k <= n_times - 2
        and (np.diff([boundary[0], *interior, boundary[1]]) > 0).all()
    )
    if not ok:
        rule = (
            f"an integer from 1 to {n_times - 2} (distinct times - 2), distinct knots"
        )
        fail("bad_value", arg="--knots", rule=rule, value=k)
    return {"interior": interior, "boundary": boundary}


def check_rank(data, names, p):
    """Stop when the fixed terms are aliased (a covariate confounded with time)."""
    fixed = [n for n in names if n not in ("value", "subject")]
    meta = data[["subject", *fixed]].drop_duplicates().copy()
    for col in fixed:
        parsed = pd.to_numeric(meta[col], errors="coerce")
        numeric = parsed.notna().all() and not (
            col == "time" and p["time_type"] == "factor"
        )
        meta[col] = parsed if numeric else meta[col].astype(str)
    skillkit.check_design_rank(meta, fixed)


def run_lmer(data, p, knots):
    """Step 3: lmer.R fits one model per feature; returns its table and versions."""
    rscript = os.environ.get("RSCRIPT") or shutil.which("Rscript")
    if not rscript or not Path(rscript).is_file():
        fail("no_rscript")
    # One BLAS/OpenMP thread: lme4 would otherwise fan out on shared nodes.
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    with tempfile.TemporaryDirectory() as tmp:
        long_tsv, out = Path(tmp) / "long.tsv", Path(tmp) / "fits.tsv"
        data.to_csv(long_tsv, sep="\t", index=False)
        cmd = [
            rscript,
            str(ENGINE),
            str(long_tsv),
            p["formula"],
            p["term"],
            p["time_type"],
            str(p["reference_time"]),
            str(p["seed"]),
            str(out),
            ",".join(map(str, knots["interior"])) if knots else "",
            ",".join(map(str, knots["boundary"])) if knots else "",
            p["missing"] or "complete-case",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False, env=env)
        lines = proc.stderr.strip().splitlines()
        if proc.returncode == 4:
            fail("missing_r_package", detail=lines[-1].strip())
        if proc.returncode == 3:
            choices = next(x for x in lines if x.startswith("coefficients:"))
            fail(
                "term_absent", term=p["term"], choices=choices.split(":", 1)[1].strip()
            )
        if proc.returncode != 0:
            fail("r_failed", code=proc.returncode, tail=" | ".join(lines[-5:]))
        fits = pd.read_csv(out, sep="\t", keep_default_na=False, na_values=[""])
    fits["fit_message"] = fits["fit_message"].fillna("")
    return fits, [f"Rscript {rscript}", *proc.stdout.strip().splitlines()]


def adjust(pvalues, method):
    """Step 4: correct across the fitted features only; NA stays NA."""
    padj = pd.Series(np.nan, index=pvalues.index)
    tested = pvalues.dropna()
    if method == "bh":
        padj[tested.index] = stats.false_discovery_control(tested, method="bh")
    else:
        padj[tested.index] = np.minimum(tested * len(tested), 1)
    return padj


def time_axis(values, p):
    """x positions: visit index for factor time (reference first), else the number."""
    if p["time_type"] != "factor":
        return values.astype(float), None
    ref = p["reference_time"]
    levels = [ref, *sorted(set(values) - {ref})]
    return values.map({lv: i for i, lv in enumerate(levels)}).astype(float), levels


def subject_lines(rows, levels):
    """Subject x visit grid of values; NaN at a missed visit breaks the drawn line.

    The grid is every factor level, or every time seen for numeric time; with
    irregular numeric times a line breaks wherever another subject was seen.
    """
    grid = (
        [float(i) for i in range(len(levels))] if levels else sorted(rows["x"].unique())
    )
    return rows.pivot(index="subject", columns="x", values="value").reindex(
        columns=grid
    )


def draw_trajectories(rows, levels, p, value_label):
    """Per-subject lines that never bridge a missing visit, plus mean and 95% CI."""
    wide = subject_lines(rows, levels)
    means = rows.groupby("x")["value"].agg(["mean", "std", "count"]).reset_index()
    half = stats.t.ppf(0.975, means["count"] - 1) * means["std"]
    half = half / np.sqrt(means["count"])
    fig, ax = plt.subplots(layout="constrained")
    for subject, values in wide.iterrows():
        ax.plot(wide.columns, values, c=SUBJECT_COLOUR, lw=0.4, gid=subject)
    ax.scatter(rows["x"], rows["value"], s=3, c=SUBJECT_COLOUR, zorder=2)
    ax.fill_between(
        means["x"],
        means["mean"] - half,
        means["mean"] + half,
        color=MEAN_COLOUR,
        alpha=0.25,
        lw=0,
        zorder=3,
    )
    ax.plot(means["x"], means["mean"], c=MEAN_COLOUR, lw=1.2, zorder=4)
    if levels:
        ax.set_xticks(range(len(levels)), levels)
    ax.set_xlabel(f"time ({p['time_unit']})")
    ax.set_ylabel(value_label)
    return fig, wide


def draw_figures(figures, data, table, p, family):
    """Step 6: spaghetti per top feature, forest, p-value histogram."""
    captions = {}
    joint = p["time_type"] == "spline" and p["term"] == "time"
    top = table.dropna(subset=["padj"]).head(p["n_plot"])
    _, levels = time_axis(data["time"], p)
    n_subjects = data["subject"].nunique()
    for rank, row in enumerate(top.itertuples(index=False), start=1):
        name = f"top{rank:02d}"
        rows = data.loc[data["feature"] == row.feature, ["subject", "time", "value"]]
        rows = rows.assign(x=time_axis(rows["time"], p)[0]).sort_values(
            ["subject", "x"]
        )
        label = f"{row.feature} ({p['value_unit']})"
        fig, wide = draw_trajectories(rows, levels, p, label)
        # The figure CSV holds the drawn grid: an empty value is a broken line.
        drawn = wide.rename_axis(columns="x").stack(future_stack=True)
        save_panel(
            fig, drawn.rename("value").reset_index(), figures, name, PANEL_MM, PANEL_MM
        )
        plt.close(fig)
        n_gap = int(wide.isna().any(axis=1).sum())
        n_absent = n_subjects - len(wide)
        test = (
            f"joint F {row.f:.3g} on {row.num_df:g} and {row.df:.3g} df"
            if joint
            else f"{p['term']} estimate {row.estimate:.3g}"
        )
        captions[name] = (
            f"Line plot of {row.feature} over time; grey lines join one subject's "
            f"visits and break at a missed visit (subjects with a gap: {n_gap}; "
            f"dropped, no values: {n_absent}). Red: mean with 95% CI (t). Mixed-model "
            f"{test}, padj {row.padj:.2g} ({family})."
        )

    # A joint F test has no single estimate to draw.
    if joint:
        return captions | pvalue_histogram(table, figures, p)

    forest = table.dropna(subset=["padj"]).head(N_FOREST).iloc[::-1]
    half = stats.t.ppf(0.975, forest["df"]) * forest["se"]
    forest = forest.assign(
        ci_low=forest["estimate"] - half, ci_high=forest["estimate"] + half
    )
    flagged = forest["singular"].eq(True) | forest["converged"].eq(False)
    y = np.arange(len(forest))
    fig, ax = plt.subplots(layout="constrained")
    ax.axvline(0, ls="--", lw=0.5, c="black")
    colours = np.where(forest["significant"], MEAN_COLOUR, OTHER)
    ax.hlines(y, forest["ci_low"], forest["ci_high"], colors=colours, lw=0.8)
    ax.scatter(
        forest["estimate"],
        y,
        s=12,
        c=np.where(flagged, "white", colours),
        edgecolors=colours,
        zorder=3,
    )
    ax.set_yticks(y, forest["feature"])
    ax.set_xlabel(f"{p['term']} estimate ({p['value_unit']}, 95% CI)")
    save_panel(fig, forest, figures, "forest", 89, 20 + 3.5 * len(forest))
    plt.close(fig)
    captions["forest"] = (
        f"Forest plot of the {p['term']} estimate with its 95% CI (Satterthwaite df) "
        f"for the {len(forest)} features with the lowest padj; red rows have "
        f"padj < {p['padj_max']}, hollow points are singular or unconverged fits."
    )

    return captions | pvalue_histogram(table, figures, p)


def pvalue_histogram(table, figures, p):
    """Raw p-values of the fitted features in 20 bins; returns its caption."""
    tested = table.dropna(subset=["pvalue"])
    fig, ax = plt.subplots(layout="constrained")
    ax.hist(
        tested["pvalue"], bins=np.linspace(0, 1, 21), color=OTHER, edgecolor="white"
    )
    ax.set_xlabel("raw p-value")
    ax.set_ylabel("features")
    save_panel(fig, tested[["feature", "pvalue"]], figures, "pvalue_hist", 89, 70)
    plt.close(fig)
    return {
        "pvalue_hist": f"Histogram of raw {p['term']} p-values for {len(tested)} "
        "fitted features in 20 bins of width 0.05; a flat floor with a peak near 0 "
        "is expected."
    }


def analyse(argv, argv_as_run):
    p = parse_args(argv)
    names = check_values(p)
    skillkit.seed_everything(p["seed"])
    data, n_dropped, dropped, imputation = read_values(p, names)
    knots = spline_knots(data["time"], p) if p["time_type"] == "spline" else None
    fits, r_versions = run_lmer(data, p, knots)

    n_fitted = int(fits["pvalue"].notna().sum())
    if not n_fitted:
        fail("no_tests", message=fits["fit_message"].iloc[0] or "none")
    fits["padj"] = adjust(fits["pvalue"], p["adjust"])
    fits["significant"] = fits["padj"] < p["padj_max"]
    table = fits.sort_values(["padj", "feature"], na_position="last")
    table = table[
        [
            "feature",
            "estimate",
            "se",
            "df",
            "t",
            "num_df",
            "f",
            "pvalue",
            "padj",
            "significant",
            "singular",
            "converged",
            "n_obs",
            "n_subjects",
            "fit_message",
        ]
    ]

    outdir = skillkit.make_outdir(p["outdir"])
    table.to_csv(
        outdir / "tables" / "model_results.tsv",
        sep="\t",
        index=False,
        float_format="%.6g",
    )
    plt.switch_backend("Agg")
    family = f"{p['adjust'].upper()} across {n_fitted} features"
    captions = draw_figures(outdir / "figures", data, table, p, family)

    sig = table.loc[table["significant"], "feature"].tolist()
    singular = table.loc[table["singular"].eq(True), "feature"].tolist()
    not_converged = table.loc[table["converged"].eq(False), "feature"].tolist()
    failed = table.loc[table["pvalue"].isna(), "feature"].tolist()
    # The missing-covariate choice and the spline knots ride on the first line.
    choice = {
        None: "",
        "complete-case": "complete-case, those rows dropped",
        "impute": "multiple imputation (mice, m = 5, pooled with Rubin's rules)",
    }[p["missing"]]
    notes = []
    if knots:
        notes.append(
            f"natural-spline time, interior knots {knots['interior']}, boundary "
            f"{knots['boundary']}"
        )
    if choice:
        notes.append(
            f"{', '.join(imputation['covariates'])} missing in {imputation['n_rows']} "
            f"rows of {imputation['n_subjects']} subjects: {choice}"
        )
    summary = [
        f"{n_fitted} of {len(table)} features fitted ({p['formula']}, term "
        f"{p['term']}); {n_dropped} input rows dropped (missing value), per time: "
        + ", ".join(f"{k} {v}" for k, v in dropped.items())
        + "".join(f"; {n}" for n in notes)
        + ".",
        f"{len(singular)} singular and {len(not_converged)} unconverged fits, kept "
        f"and flagged; {len(failed)} features without a p-value.",
        f"{len(sig)} features with padj < {p['padj_max']} ({family}): "
        f"{', '.join(sig[:10]) or 'none'}"
        + (f" and {len(sig) - 10} more." if len(sig) > 10 else "."),
    ]
    title = f"Longitudinal mixed models: {p['term']}"
    keys = [*REQUIRED, *OPTIONAL]
    skillkit.write_report(outdir, title, p, keys, summary, captions)
    result = {
        "term": p["term"],
        "n_features": len(table),
        "n_fitted": n_fitted,
        "n_significant": len(sig),
        "n_singular": len(singular),
        "n_not_converged": len(not_converged),
        "n_rows_dropped": n_dropped,
        "n_rows_dropped_by_time": dropped,
        "time_knots": knots,
        "missing_covariates": imputation,
        "significant_features": sig,
        "singular_features": singular,
        "not_converged_features": not_converged,
        "failed_features": failed,
        "params": p,
        "seed": p["seed"],
    }
    artifacts = [
        "report.md",
        "tables/model_results.tsv",
        *["figures/forest.png"] * ("forest" in captions),
        *[f"figures/{n}.png" for n in captions if n.startswith("top")],
    ]
    skillkit.write_result(outdir, result, summary, artifacts, context=p["context"])
    packages = ["numpy", "pandas", "scipy", "matplotlib", "plotnine"]
    skillkit.write_reproducibility(
        outdir, __file__, p, keys, argv_as_run, packages, extra=r_versions
    )
    return outdir


if __name__ == "__main__":
    raise SystemExit(skillkit.run("longitudinal", analyse))
