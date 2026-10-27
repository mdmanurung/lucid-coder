"""Purpose: multi-omics factor analysis with MOFA2 (mofapy2) on per-view sample tables.
Inputs:  repeated --view name=<samples x features TSV>; --metadata <TSV>, group column
Outputs: <outdir>/report.md, result.json, tables/, figures/, data/, reproducibility/
Run:     python scripts/mofa.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""

import os

# Pin BLAS threads before numpy loads: one thread gives bit-identical reruns on CPU.
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
for var in THREAD_VARS:
    os.environ[var] = "1"

import argparse  # noqa: E402
import contextlib  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import h5py  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import skillkit  # noqa: E402
from plots import grey_to_red, plot_box, plot_embedding  # noqa: E402
from scipy import stats  # noqa: E402
from theme import save_panel  # noqa: E402

SKILL_DIR = Path(__file__).resolve().parents[1]
LIKELIHOODS = ("gaussian", "bernoulli", "poisson")
REQUIRED = [
    "view",
    "metadata",
    "group_col",
    "likelihoods",
    "factors",
    "scale_views",
    "convergence",
    "min_var_explained",
    "n_top_features",
    "seed",
    "reference",
]
EXAMPLES = SKILL_DIR / "examples"
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "view": [f"{v}={EXAMPLES / f'demo_{v}.tsv'}" for v in ("rna", "adt", "tf")],
    "metadata": str(EXAMPLES / "demo_meta.tsv"),
    "group_col": "group",
    "likelihoods": "rna=gaussian,adt=gaussian,tf=gaussian",
    "factors": 6,
    "scale_views": "yes",
    "convergence": "fast",
    "min_var_explained": 5.0,
    "n_top_features": 10,
    "seed": 0,
    "reference": "synthetic demo, no genome build",
}

ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "bad_view": "--view must be name=path, got {value!r}",
    "duplicate_view": "view {view!r} is given more than once",
    "bad_likelihoods": "--likelihoods must give one of {choices} for each view "
    "(name=likelihood,...); {problem}",
    "missing_file": "input file not found: {path}",
    "missing_column": "column {col!r} not found in {path}",
    "not_numeric": "view {view!r} has non-numeric values in {n} cells",
    "shared_feature": "{n} feature names appear in more than one view, e.g. "
    "{example!r}; prefix them with the view name",
    "bad_bernoulli": "view {view!r} is bernoulli but holds values other than 0 and 1",
    "bad_poisson": "view {view!r} is poisson but holds negative or non-integer values",
    "unknown_samples": "{n} samples in view {view!r} are not in the metadata, "
    "e.g. {example!r}",
    "missing_group": "{n} modelled samples have no value in {col!r}",
    "too_few_groups": "{col!r} needs at least 2 levels, got {levels}",
    "too_many_labels": "{what} has {n} levels; the fixed palette holds {max}",
    "too_many_factors": "--factors ({factors}) must be below the number of "
    "samples ({n})",
}

# Fixed plot spec, not analysis thresholds: Tol bright colours for groups and views,
# corpus jitter width (ggboxplot add="jitter"), panel sizes in mm.
COLOURS = ("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB")
JITTER_WIDTH = 0.2
PANEL_MM = 60
PACKAGES = ["mofapy2", "h5py", "numpy", "pandas", "scipy", "matplotlib"]


def fail(key, **fields):
    skillkit.fail(ERRORS, key, **fields)


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--view", action="append")
    parser.add_argument("--metadata")
    parser.add_argument("--group-col")
    parser.add_argument("--likelihoods")
    parser.add_argument("--factors", type=int)
    parser.add_argument("--scale-views", choices=["yes", "no"])
    parser.add_argument("--convergence", choices=["fast", "medium", "slow"])
    parser.add_argument("--min-var-explained", type=float)
    parser.add_argument("--n-top-features", type=int)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--reference")
    return skillkit.parse_required(parser, REQUIRED, DEMO, argv)


def check_values(p):
    """Step 1: value ranges, views as {name: path}, likelihoods as {name: lik}."""
    rules = [
        ("--factors", p["factors"], p["factors"] >= 2, "an integer >= 2"),
        (
            "--min-var-explained",
            p["min_var_explained"],
            0 <= p["min_var_explained"] < 100,
            "a percentage in [0, 100)",
        ),
        ("--n-top-features", p["n_top_features"], p["n_top_features"] >= 1, ">= 1"),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ]
    for arg, value, ok, rule in rules:
        if not ok:
            fail("bad_value", arg=arg, rule=rule, value=value)
    views = {}
    for item in p["view"]:
        name, sep, path = item.partition("=")
        if not (sep and name and path):
            fail("bad_view", value=item)
        if name in views:
            fail("duplicate_view", view=name)
        views[name] = Path(path)
    if len(views) > len(COLOURS):
        fail("too_many_labels", what="--view", n=len(views), max=len(COLOURS))
    liks = dict(item.partition("=")[::2] for item in p["likelihoods"].split(","))
    problems = [
        f"missing for {sorted(views.keys() - liks.keys())}",
        f"unknown view {sorted(liks.keys() - views.keys())}",
        f"unknown likelihood {sorted(set(liks.values()) - set(LIKELIHOODS))}",
    ]
    bad = [
        msg
        for msg, ok in zip(
            problems,
            [
                views.keys() <= liks.keys(),
                liks.keys() <= views.keys(),
                set(liks.values()) <= set(LIKELIHOODS),
            ],
            strict=True,
        )
        if not ok
    ]
    if bad:
        fail("bad_likelihoods", choices="/".join(LIKELIHOODS), problem="; ".join(bad))
    return dict(sorted(views.items())), liks


def read_tsv(path):
    if not path.is_file():
        fail("missing_file", path=path)
    # Labels stay strings ("01" is not 1); values are parsed by the caller.
    return pd.read_csv(
        path, sep="\t", dtype=str, keep_default_na=False, na_values=["", "NA", "NaN"]
    )


def read_views(views, liks):
    """Step 2: read each view, check values against its likelihood."""
    tables = {}
    for name, path in views.items():
        raw = read_tsv(path)
        # read_csv renames a repeated column ("x.1"), so check the raw header.
        with path.open() as fh:
            header = fh.readline().rstrip("\n").split("\t")[1:]
        skillkit.check_unique(header, f"feature names in view {name!r}")
        raw = raw.set_index(raw.columns[0])
        skillkit.check_unique(raw.index, f"sample IDs in view {name!r}")
        table = raw.apply(pd.to_numeric, errors="coerce")
        n_bad = int((table.isna() & raw.notna()).sum().sum())
        if n_bad:
            fail("not_numeric", view=name, n=n_bad)
        values = table.to_numpy()[~np.isnan(table.to_numpy())]
        if liks[name] == "bernoulli" and not np.isin(values, (0, 1)).all():
            fail("bad_bernoulli", view=name)
        if liks[name] == "poisson" and ((values < 0) | (values % 1 != 0)).any():
            fail("bad_poisson", view=name)
        tables[name] = table.sort_index(axis=1)
    features = pd.Series([f for t in tables.values() for f in t.columns])
    shared = features[features.duplicated()].unique()
    if len(shared):
        fail("shared_feature", n=len(shared), example=shared[0])
    return tables


def read_metadata(p, tables):
    """Step 3: metadata rows for the modelled samples (sorted union over views)."""
    path = Path(p["metadata"])
    meta = read_tsv(path)
    for col in ("sample", p["group_col"]):
        if col not in meta.columns:
            fail("missing_column", col=col, path=path)
    skillkit.check_unique(meta["sample"], f"sample IDs in {path}")
    meta = meta.set_index("sample")
    for name, table in tables.items():
        unknown = table.index.difference(meta.index)
        if len(unknown):
            fail("unknown_samples", n=len(unknown), view=name, example=unknown[0])
    samples = sorted(set().union(*(t.index for t in tables.values())))
    # Metadata rows in no view are left out; count them per group for the report.
    left_out = meta.loc[~meta.index.isin(samples), p["group_col"]]
    dropped = {str(k): int(v) for k, v in left_out.value_counts().sort_index().items()}
    meta = meta.loc[samples]
    groups = meta[p["group_col"]]
    if groups.isna().any():
        fail("missing_group", n=int(groups.isna().sum()), col=p["group_col"])
    levels = sorted(groups.unique())
    if len(levels) < 2:
        fail("too_few_groups", col=p["group_col"], levels=levels)
    if len(levels) > len(COLOURS):
        fail("too_many_labels", what=p["group_col"], n=len(levels), max=len(COLOURS))
    if p["factors"] >= len(samples):
        fail("too_many_factors", factors=p["factors"], n=len(samples))
    return meta, levels, dropped


def train(p, tables, liks, samples, model_path):
    """Step 4: train MOFA on CPU; mofapy2 seeds numpy inside build() from --seed.

    Views, samples and features are sorted first, because the seeded random
    initialisation maps onto input order.
    """
    from mofapy2.run.entry_point import entry_point

    # mofapy2 prints its log to stdout; keep stdout for the skill's one-line result.
    with contextlib.redirect_stdout(sys.stderr):
        ent = entry_point()
        ent.set_data_options(scale_views=p["scale_views"] == "yes")
        ent.set_data_matrix(
            [[t.reindex(samples).to_numpy(float)] for t in tables.values()],
            likelihoods=[liks[v] for v in tables],
            views_names=list(tables),
            samples_names=[samples],
            features_names=[list(t.columns) for t in tables.values()],
        )
        ent.set_model_options(factors=p["factors"])
        # dropR2=None: no factor is dropped during training; low ones are flagged.
        ent.set_train_options(
            convergence_mode=p["convergence"],
            seed=p["seed"],
            dropR2=None,
            gpu_mode=False,
            quiet=True,
        )
        ent.build()
        ent.run()
        # save() sorts factors by variance explained, so Factor1 explains most.
        ent.save(str(model_path))


def read_model(model_path, views):
    """Step 5: factors, weights and R2 (percent) from the saved, sorted model."""
    with h5py.File(model_path, "r") as h5:
        group = h5["groups/groups"][0].decode()
        z = h5[f"expectations/Z/{group}"][:]
        names = [f"Factor{k + 1}" for k in range(z.shape[0])]
        samples = [s.decode() for s in h5[f"samples/{group}"][:]]
        factors = pd.DataFrame(
            z.T, index=pd.Index(samples, name="sample"), columns=names
        )
        weights = pd.concat(
            [
                pd.DataFrame(
                    h5[f"expectations/W/{v}"][:].T,
                    index=[f.decode() for f in h5[f"features/{v}"][:]],
                    columns=names,
                )
                .rename_axis("feature")
                .reset_index()
                .melt(id_vars="feature", var_name="factor", value_name="weight")
                .assign(view=v)
                for v in views
            ]
        )
        view_order = [v.decode() for v in h5["views/views"][:]]
        r2 = pd.DataFrame(
            h5[f"variance_explained/r2_per_factor/{group}"][:],
            index=view_order,
            columns=names,
        )
        r2_total = dict(
            zip(view_order, h5[f"variance_explained/r2_total/{group}"][:], strict=True)
        )
    return factors, weights, r2, r2_total


def associate(factors, meta):
    """Step 7: each factor against each metadata column; BH across all rows."""
    rows = []
    for col in sorted(meta.columns):
        values = meta[col]
        numeric = pd.to_numeric(values, errors="coerce")
        is_numeric = numeric.notna().sum() == values.notna().sum()
        for factor in factors.columns:
            keep = values.notna().to_numpy()
            x, z = values[keep], factors.loc[keep, factor].to_numpy()
            out = {"factor": factor, "covariate": col, "n": int(keep.sum())}
            if is_numeric:
                out |= {"test": "spearman", "statistic": np.nan, "pvalue": np.nan}
                if numeric[keep].nunique() > 1:
                    fit = stats.spearmanr(numeric[keep], z)
                    out |= {"statistic": fit.statistic, "pvalue": fit.pvalue}
            else:
                out |= {"test": "kruskal-wallis", "statistic": np.nan, "pvalue": np.nan}
                parts = [z[(x == lv).to_numpy()] for lv in sorted(x.unique())]
                if len(parts) > 1:
                    fit = stats.kruskal(*parts)
                    out |= {"statistic": fit.statistic, "pvalue": fit.pvalue}
            rows.append(out)
    table = pd.DataFrame(rows)
    tested = table["pvalue"].notna()
    table["padj"] = np.nan
    if tested.any():
        table.loc[tested, "padj"] = stats.false_discovery_control(
            table.loc[tested, "pvalue"], method="bh"
        )
    table["order"] = table["factor"].str[6:].astype(int)
    return table.sort_values(["order", "covariate"]).drop(columns="order")


def plot_weights(top, pal):
    """Cleveland dot chart of weights; plots.py has no such form."""
    fig, ax = plt.subplots()
    y = np.arange(len(top))
    ax.axvline(0, ls="--", lw=0.5, c="black")
    ax.scatter(top["weight"], y, s=8, c=[pal[v] for v in top["view"]], zorder=3)
    ax.set_yticks(y, top["feature"])
    ax.invert_yaxis()
    handles = [
        plt.Line2D([], [], marker="o", ls="", c=pal[v], label=v)
        for v in top["view"].unique()
    ]
    ax.legend(handles=handles, title="view", loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_xlabel("weight (a.u.)")
    return fig


def draw_figures(figures, p, factors, weights, r2_long, meta, levels, active, assoc):
    """Step 8: heatmap, factor scatter, weights and box per factor above threshold."""
    captions = {}
    thr = p["min_var_explained"]

    def save(fig, data, name, width, height, caption):
        fig.set_layout_engine("tight")
        save_panel(fig, data, figures, name, width, height)
        plt.close(fig)
        captions[name] = caption

    # Sequential, not plots.plot_heatmap's diverging scale: R2 is never negative.
    vmax = float(r2_long["r2"].max())
    matrix = r2_long.pivot(index="factor", columns="view", values="r2").loc[
        factors.columns
    ]
    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(matrix, cmap=grey_to_red, vmin=0, vmax=vmax)
    ax.set_xticks(np.arange(matrix.shape[1]) + 0.5, matrix.columns, rotation=90)
    ax.set_yticks(np.arange(matrix.shape[0]) + 0.5, matrix.index)
    ax.invert_yaxis()
    fig.colorbar(mesh, ax=ax, label="variance explained (%)")
    save(
        fig,
        r2_long,
        "variance_explained",
        PANEL_MM,
        PANEL_MM,
        f"Heatmap of variance explained (R2, %) per factor and view; rows fixed in "
        f"factor order (total R2, as mofapy2 saves them), views sorted by name, not "
        f"clustered; grey-to-red scale from 0 to {vmax:.1f}%. Factors above "
        f"{thr}% in at least one view: {', '.join(active) or 'none'}.",
    )

    group_pal = dict(zip(levels, COLOURS, strict=False))
    scores = factors.join(meta[p["group_col"]]).reset_index()
    fig = plot_embedding(
        scores, "Factor1", "Factor2", p["group_col"], group_pal, axis_name="Factor"
    )
    ax = fig.axes[0]
    ax.collections[0].set_sizes([6])
    ax.set_aspect("auto")
    ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    ax.xaxis.set_major_locator(plt.AutoLocator())
    ax.yaxis.set_major_locator(plt.AutoLocator())
    save(
        fig,
        scores,
        "factor1_vs_factor2",
        PANEL_MM + 20,
        PANEL_MM,
        f"Scatter plot of Factor1 against Factor2 scores; each point is one sample, "
        f"coloured by {p['group_col']}. MOFA with {p['factors']} factors, views "
        f"{'scaled to unit variance' if p['scale_views'] == 'yes' else 'unscaled'}, "
        f"seed {p['seed']}; no downsampling.",
    )

    hits = assoc[assoc["covariate"] == p["group_col"]]
    padj = dict(zip(hits["factor"], hits["padj"], strict=True))
    n_tests = int(assoc["pvalue"].notna().sum())
    view_pal = dict(zip(weights["view"].unique(), COLOURS, strict=False))
    for factor in active:
        top = (
            weights[weights["factor"] == factor]
            .assign(absw=lambda d: d["weight"].abs())
            .sort_values(["view", "absw", "feature"], ascending=[True, False, True])
            .groupby("view", sort=True)
            .head(p["n_top_features"])
            .drop(columns="absw")
        )
        save(
            plot_weights(top, view_pal),
            top,
            f"weights_{factor}",
            PANEL_MM + 20,
            max(PANEL_MM, 3 * len(top) + 15),
            f"Dot chart of the {p['n_top_features']} largest absolute {factor} "
            f"weights per view; unscaled weights, on the "
            f"{'unit-variance view' if p['scale_views'] == 'yes' else 'input'} scale.",
        )
        box = scores[[p["group_col"], factor]]
        save(
            plot_box(
                box,
                p["group_col"],
                factor,
                group_pal,
                jitter_width=JITTER_WIDTH,
                seed=p["seed"],
                value_label=f"{factor} score (a.u.)",
            ),
            box,
            f"box_{factor}",
            PANEL_MM,
            PANEL_MM,
            f"Box plot of {factor} scores per {p['group_col']}; each point is one "
            f"sample, jittered sideways (seed {p['seed']}). Kruskal-Wallis padj "
            f"{padj[factor]:.2g} "
            f"(BH over {n_tests} factor-covariate tests).",
        )
    return captions


def analyse(argv, argv_as_run):
    p = parse_args(argv)
    views, liks = check_values(p)
    skillkit.seed_everything(p["seed"])
    tables = read_views(views, liks)
    meta, levels, dropped = read_metadata(p, tables)
    samples = list(meta.index)

    outdir = skillkit.make_outdir(p["outdir"])
    (outdir / "data").mkdir()
    model_path = outdir / "data" / "model.hdf5"
    train(p, tables, liks, samples, model_path)
    factors, weights, r2, r2_total = read_model(model_path, list(tables))

    # Step 6: flag, never drop. A factor explains a view when its R2 there >= threshold.
    thr = p["min_var_explained"]
    r2_long = (
        r2.rename_axis("view")
        .reset_index()
        .melt(id_vars="view", var_name="factor", value_name="r2")
    )
    r2_long["explains"] = r2_long["r2"] >= thr
    r2_long["active"] = r2_long.groupby("factor")["explains"].transform("any")
    r2_long["order"] = r2_long["factor"].str[6:].astype(int)
    r2_long = r2_long.sort_values(["order", "view"]).drop(columns="order")
    r2_long = r2_long[["factor", "view", "r2", "explains", "active"]]
    active = [
        f
        for f in factors.columns
        if r2_long.loc[r2_long["factor"] == f, "active"].any()
    ]
    explains = {
        f: r2_long.loc[(r2_long["factor"] == f) & r2_long["explains"], "view"].tolist()
        for f in active
    }
    weights["order"] = weights["factor"].str[6:].astype(int)
    weights = weights.sort_values(["view", "feature", "order"])[
        ["view", "feature", "factor", "weight"]
    ]
    assoc = associate(factors, meta)

    tab = outdir / "tables"
    fmt = {"sep": "\t", "index": False, "float_format": "%.6g"}
    r2_long.to_csv(tab / "variance_explained.tsv", **fmt)
    factors.reset_index().to_csv(tab / "factors.tsv", **fmt)
    weights.to_csv(tab / "weights.tsv", **fmt)
    assoc.to_csv(tab / "factor_vs_metadata.tsv", **fmt)

    plt.switch_backend("Agg")
    captions = draw_figures(
        outdir / "figures", p, factors, weights, r2_long, meta, levels, active, assoc
    )

    by_group = assoc[assoc["covariate"] == p["group_col"]].sort_values("padj")
    summary = [
        f"MOFA on {len(samples)} samples and {len(tables)} views "
        f"({', '.join(f'{v} {t.shape[1]}' for v, t in tables.items())} features); "
        f"{p['factors']} factors trained, none dropped; metadata samples in no view "
        f"(left out): {', '.join(f'{g} {n}' for g, n in dropped.items()) or 'none'}.",
        f"{len(active)} factors explain >= {thr}% of a view: "
        + ("; ".join(f"{f} ({', '.join(v)})" for f, v in explains.items()) or "none")
        + ".",
        f"Closest factor to {p['group_col']}: {by_group['factor'].iloc[0]}, "
        f"Kruskal-Wallis padj {by_group['padj'].iloc[0]:.2g} (BH over "
        f"{int(assoc['pvalue'].notna().sum())} factor-covariate tests).",
    ]
    shown = {**p, "view": " ".join(p["view"])}
    skillkit.write_report(
        outdir, "MOFA multi-omics factors", shown, REQUIRED, summary, captions
    )
    threads = [f"{v}={os.environ[v]}" for v in THREAD_VARS]
    result = {
        "n_samples": len(samples),
        "views": {v: t.shape[1] for v, t in tables.items()},
        "n_factors": len(factors.columns),
        "n_active": len(active),
        "active_factors": explains,
        "r2_total": {v: round(float(x), 4) for v, x in r2_total.items()},
        "device": "cpu",
        "threads": threads,
        "dropped_per_group": dropped,
        "params": p,
        "seed": p["seed"],
    }
    artifacts = [
        "report.md",
        "tables/variance_explained.tsv",
        "tables/factor_vs_metadata.tsv",
        *[f"figures/{n}.png" for n in captions],
    ]
    skillkit.write_result(outdir, result, summary, artifacts, reference=p["reference"])
    skillkit.write_reproducibility(
        outdir, __file__, p, REQUIRED, argv_as_run, PACKAGES, ["device cpu", *threads]
    )
    return outdir


if __name__ == "__main__":
    raise SystemExit(skillkit.run("mofa", analyse))
