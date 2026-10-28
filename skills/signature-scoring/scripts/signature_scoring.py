"""Purpose: score gene signatures per cell, summarise per sample, compare to random.
Inputs:  --adata <file.h5ad> with log-normalised X and obs --sample-col, --group-col
         (e.g. sc-embed's data/embedded.h5ad); --gene-sets TSV with columns set, gene
Outputs: <outdir>/report.md, result.json, data/scored.h5ad, tables/, figures/,
         reproducibility/
Run:     python scripts/signature_scoring.py --demo --outdir demo_out
Notes:   Every analysis argument is required; --demo alone pins the demo values.
         Steps run in a fixed order; see SKILL.md. Bad input: "ERROR: ..." and exit 2.
"""

import os

# AUCell runs on numba; one thread keeps runs identical.
for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(var, "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import argparse  # noqa: E402
import re  # noqa: E402
from pathlib import Path  # noqa: E402

import anndata as ad  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
from plots import FEW_POINTS, grey_to_red, plot_box  # noqa: E402
from scipy import sparse  # noqa: E402
from skillkit import (  # noqa: E402
    check_unique,
    fail,
    flag,
    make_outdir,
    parse_required,
    run,
    seed_everything,
    with_status,
    write_report,
    write_reproducibility,
    write_result,
)
from theme import save_panel  # noqa: E402

SKILL_DIR = Path(__file__).resolve().parents[1]

REQUIRED = [
    "adata",
    "gene_sets",
    "sample_col",
    "group_col",
    "method",
    "min_genes",
    "n_random",
    "reference",
    "seed",
    "context",
]
ENGINE_ARGS = {"score_genes": ["ctrl_size", "n_bins"], "aucell": ["n_up"]}
OPTIONAL = ["ctrl_size", "n_bins", "n_up"]
KEYS = [*REQUIRED[:5], *OPTIONAL, *REQUIRED[5:]]
# Demo-only values, not defaults: SKILL.md lists them under the same label.
DEMO = {
    "adata": str(SKILL_DIR / "examples" / "demo.h5ad"),
    "gene_sets": str(SKILL_DIR / "examples" / "demo_gene_sets.tsv"),
    "sample_col": "sample",
    "group_col": "group",
    "method": "score_genes",
    "ctrl_size": 50,
    "n_bins": 25,
    "n_up": None,
    "min_genes": 5,
    "n_random": 50,
    "reference": "synthetic demo genes, no genome build",
    "seed": 0,
    "context": "synthetic demo: 8 samples, healthy vs disease, no real tissue",
}

ERRORS = {
    "bad_value": "{arg} must be {rule}, got {value}",
    "needs_args": "--method {method} needs {args}",
    "unused_args": "--method {method} does not use {args}",
    "missing_file": "input file not found: {path}",
    "missing_column": "obs column not found: {col}",
    "bad_gene_sets": "--gene-sets needs columns set, gene; found {cols}",
    "bad_set_name": "set names may use only letters, digits, _ . -: {names}",
    "too_few_genes": "sets with fewer than {min} genes in the data: {sets}",
    "not_lognorm": "X does not look log-normalised: {why}",
    "split_sample": "samples spanning several values of {col}: {samples}",
    "few_groups": "--group-col {col} needs at least 2 groups, found {n}",
    "many_groups": "--group-col {col} has {n} groups; box plot colours: {max}",
    "same_columns": "--sample-col and --group-col must differ, both are {col}",
}

# Fixed plot spec, not analysis values: Paul Tol bright, in group order; box-plot
# jitter in box widths (group-tests' position_jitter(0.2)); y-axis unit per engine.
GROUP_COLOURS = ("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377")
JITTER_WIDTH = 0.2
PANEL_MM = 60
SCORE_UNIT = {"score_genes": "a.u.", "aucell": "AUC"}


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--adata")
    parser.add_argument("--gene-sets")
    parser.add_argument("--sample-col")
    parser.add_argument("--group-col")
    parser.add_argument("--method", choices=list(ENGINE_ARGS))
    parser.add_argument("--ctrl-size", type=int)
    parser.add_argument("--n-bins", type=int)
    parser.add_argument("--n-up", type=int)
    parser.add_argument("--min-genes", type=int)
    parser.add_argument("--n-random", type=int)
    parser.add_argument("--reference")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--context")
    p = parse_required(parser, REQUIRED, DEMO, argv, optional=OPTIONAL)
    used = ENGINE_ARGS[p["method"]]
    missing = [flag(k) for k in OPTIONAL if k in used and p[k] is None]
    if missing:
        fail(ERRORS, "needs_args", method=p["method"], args=" ".join(missing))
    extra = [flag(k) for k in OPTIONAL if k not in used and p[k] is not None]
    if extra:
        fail(ERRORS, "unused_args", method=p["method"], args=" ".join(extra))
    rules = [
        ("--min-genes", p["min_genes"], p["min_genes"] >= 1, "an integer >= 1"),
        ("--n-random", p["n_random"], p["n_random"] >= 1, "an integer >= 1"),
        ("--seed", p["seed"], p["seed"] >= 0, "an integer >= 0"),
    ] + [(flag(k), p[k], p[k] >= 1, "an integer >= 1") for k in used]
    for arg, value, ok, rule in rules:
        if not ok:
            fail(ERRORS, "bad_value", arg=arg, rule=rule, value=value)
    return p


def load(p):
    """Steps 1-2: read the h5ad and gene sets; check matrix, samples and groups."""
    for key in ("adata", "gene_sets"):
        if not Path(p[key]).is_file():
            fail(ERRORS, "missing_file", path=p[key])
    adata = ad.read_h5ad(p["adata"])
    check_unique(adata.var_names, "gene names")
    check_unique(adata.obs_names, "cell ids")
    if p["sample_col"] == p["group_col"]:
        fail(ERRORS, "same_columns", col=p["sample_col"])
    for col in (p["sample_col"], p["group_col"]):
        if col not in adata.obs.columns:
            fail(ERRORS, "missing_column", col=col)
    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    if (values < 0).any():
        fail(ERRORS, "not_lognorm", why="negative values (scaled?)")
    if np.array_equal(values, np.round(values)):
        fail(ERRORS, "not_lognorm", why="all integers (raw counts?)")
    if p["method"] == "aucell" and p["n_up"] > adata.n_vars:
        fail(
            ERRORS,
            "bad_value",
            arg="--n-up",
            rule=f"<= {adata.n_vars}",
            value=p["n_up"],
        )

    # Categorical order is the group order of every table and figure.
    sample, group = p["sample_col"], p["group_col"]
    groups = adata.obs[group]
    if not isinstance(groups.dtype, pd.CategoricalDtype):
        groups = pd.Categorical(
            groups.astype(str), categories=sorted(set(groups.astype(str)))
        )
    adata.obs[group] = pd.Categorical(groups).remove_unused_categories()
    adata.obs[sample] = adata.obs[sample].astype(str)
    per_sample = adata.obs.groupby(sample, observed=True)[group].nunique()
    if (per_sample > 1).any():
        split = ", ".join(per_sample.index[per_sample > 1])
        fail(ERRORS, "split_sample", col=group, samples=split)
    n_groups = len(adata.obs[group].cat.categories)
    if n_groups < 2:
        fail(ERRORS, "few_groups", col=group, n=n_groups)
    if n_groups > len(GROUP_COLOURS):
        fail(ERRORS, "many_groups", col=group, n=n_groups, max=len(GROUP_COLOURS))

    listed = pd.read_csv(p["gene_sets"], sep="\t", dtype=str)
    if not {"set", "gene"} <= set(listed.columns):
        fail(ERRORS, "bad_gene_sets", cols=", ".join(listed.columns))
    listed = listed[["set", "gene"]].dropna().drop_duplicates()
    listed = listed.sort_values(["set", "gene"])
    # Set names become file names; anything else could escape figures/.
    bad = sorted(n for n in listed["set"].unique() if not re.fullmatch(r"[\w.-]+", n))
    if bad:
        fail(ERRORS, "bad_set_name", names=", ".join(bad))
    listed["found"] = listed["gene"].isin(adata.var_names)
    missing = listed[~listed["found"]].groupby("set")["gene"].agg(",".join)
    found = listed.groupby("set").agg(
        n_listed=("gene", "size"), n_found=("found", "sum")
    )
    found = found.join(missing.rename("missing")).fillna("").reset_index()
    # Sets under --min-genes are kept as not_tested rows; stop only if none is left.
    small = found[found["n_found"] < p["min_genes"]]
    if len(small) == len(found):
        sets = ", ".join(
            f"{s} ({n})" for s, n in zip(small["set"], small["n_found"], strict=True)
        )
        fail(ERRORS, "too_few_genes", min=p["min_genes"], sets=sets)
    sets = {
        name: sorted(d.loc[d["found"], "gene"])
        for name, d in listed.groupby("set", sort=True)
        if name not in set(small["set"])
    }
    return adata, sets, found, set(listed["gene"])


def score(adata, sets, p):
    """Cells x sets score frame from the chosen engine, same arguments every call."""
    if p["method"] == "score_genes":
        scores = {}
        # ponytail: one score_genes call per set re-bins every gene; n_random x
        # n_sets calls on 100k cells is slow. Batch the bins if that bites.
        for name, genes in sets.items():
            sc.tl.score_genes(
                adata,
                genes,
                ctrl_size=p["ctrl_size"],
                n_bins=p["n_bins"],
                random_state=p["seed"],
                # use_raw=False: an h5ad carrying .raw would otherwise be scored.
                use_raw=False,
                score_name="_score",
            )
            scores[name] = adata.obs.pop("_score").to_numpy()
        return pd.DataFrame(scores, index=adata.obs_names)
    import decoupler as dc

    net = pd.DataFrame(
        [(name, gene) for name, genes in sets.items() for gene in genes],
        columns=["source", "target"],
    )
    # empty=False keeps all-zero cells so rows stay aligned with obs.
    dc.mt.aucell(
        adata, net, tmin=p["min_genes"], n_up=p["n_up"], empty=False, verbose=False
    )
    return adata.obsm.pop("score_aucell").loc[adata.obs_names, list(sets)]


def sample_means(adata, scores, p):
    """Mean score per sample x group, long: set, sample, group, value, n_cells."""
    sample, group = p["sample_col"], p["group_col"]
    cells = pd.concat([adata.obs[[sample, group]], scores], axis=1)
    long = cells.melt(id_vars=[sample, group], var_name="set", value_name="value")
    return (
        long.groupby(["set", sample, group], observed=True)["value"]
        .agg(value="mean", n_cells="size")
        .reset_index()
    )


def group_spread(means, group):
    """Per set: max minus min of the group means of the sample means."""
    by_set = means.groupby(["set", group], observed=True)["value"].mean().groupby("set")
    return by_set.max() - by_set.min()


def random_control(adata, sets, means, p):
    """Step 5: n_random size-matched random sets per signature, same engine."""
    rng = np.random.default_rng(p["seed"])
    random_sets = {}
    for name, genes in sets.items():
        pool = np.array(sorted(set(adata.var_names) - set(genes)))
        for i in range(p["n_random"]):
            draw = rng.choice(pool, size=len(genes), replace=False)
            random_sets[f"{name}|{i}"] = sorted(draw)
    random_means = sample_means(adata, score(adata, random_sets, p), p)
    random_spread = group_spread(random_means, p["group_col"]).rename("spread")
    draws = random_spread.reset_index()
    draws[["set", "draw"]] = draws["set"].str.split("|", expand=True)
    draws["draw"] = draws["draw"].astype(int)
    observed = group_spread(means, p["group_col"])
    draws["exceeds"] = draws["spread"] >= draws["set"].map(observed)
    table = (
        draws.groupby("set")
        .agg(
            random_mean=("spread", "mean"),
            random_max=("spread", "max"),
            n_exceed=("exceeds", "sum"),
        )
        .reset_index()
    )
    table = table.assign(
        n_genes=table["set"].map(lambda s: len(sets[s])),
        n_random=p["n_random"],
        observed_spread=table["set"].map(observed),
        # Empirical p with the +1 correction, so it is never 0.
        p_empirical=(1 + table["n_exceed"]) / (p["n_random"] + 1),
    )
    table["beats_random"] = table["observed_spread"] > table["random_max"]
    cols = ["set", "n_genes", "n_random", "observed_spread", "random_mean"]
    table = table[[*cols, "random_max", "p_empirical", "beats_random"]]
    draws = pd.concat(
        [
            draws.drop(columns="exceeds").assign(kind="random"),
            observed.rename("spread").reset_index().assign(kind="observed", draw=-1),
        ]
    )
    return table.sort_values("set"), draws.sort_values(["set", "kind", "draw"])


def draw_figures(figures, adata, scores, means, control, draws, p):
    group = p["group_col"]
    levels = [str(g) for g in adata.obs[group].cat.categories]
    pal = dict(zip(levels, GROUP_COLOURS, strict=False))
    for name in scores.columns:
        rows = means[means["set"] == name].copy()
        rows[group] = rows[group].astype(str)
        plot = plot_box(
            rows,
            group,
            "value",
            pal,
            jitter_width=JITTER_WIDTH,
            seed=p["seed"],
            value_label=f"{name} score, sample mean ({SCORE_UNIT[p['method']]})",
        )
        # 20 mm per box keeps "(n = 12)" tick labels apart.
        width = min(183, 40 + 20 * len(levels))
        save_panel(plot, rows, figures, f"score_by_group_{name}", width, PANEL_MM)

    if "X_umap" in adata.obsm:
        for name in scores.columns:
            umap = pd.DataFrame(
                adata.obsm["X_umap"][:, :2], columns=["UMAP1", "UMAP2"]
            ).assign(cell=adata.obs_names, score=scores[name].to_numpy())
            unit = SCORE_UNIT[p["method"]]
            # High scores on top so a small positive population stays visible.
            umap = umap.sort_values(["score", "cell"])
            fig, ax = plt.subplots()
            dots = ax.scatter(
                umap["UMAP1"],
                umap["UMAP2"],
                c=umap["score"],
                s=1,
                cmap=grey_to_red,
                rasterized=True,
            )
            fig.colorbar(dots, ax=ax, label=f"{name} score ({unit})")
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlabel("UMAP 1")
            ax.set_ylabel("UMAP 2")
            fig.set_layout_engine("tight")
            save_panel(fig, umap, figures, f"score_umap_{name}", 75, PANEL_MM)
            plt.close(fig)

    names = list(control["set"])
    x = draws["set"].map({n: i for i, n in enumerate(names)}).to_numpy(float)
    is_random = (draws["kind"] == "random").to_numpy()
    jitter = np.random.default_rng(p["seed"]).uniform(
        -JITTER_WIDTH, JITTER_WIDTH, len(x)
    )
    fig, ax = plt.subplots()
    ax.scatter(
        x[is_random] + jitter[is_random],
        draws["spread"][is_random],
        s=4,
        c="0.6",
        label="random sets",
    )
    ax.scatter(
        x[~is_random],
        draws["spread"][~is_random],
        s=20,
        c="#B2182B",
        marker="D",
        label="signature",
    )
    ax.set_xticks(range(len(names)), names)
    ax.set_xlabel("gene set")
    ax.set_ylabel(f"spread of {group} means (max - min)")
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
    fig.set_layout_engine("tight")
    width = min(183, max(89, 50 + 15 * len(names)))
    save_panel(fig, draws, figures, "random_control", width, PANEL_MM)
    plt.close(fig)


def analyse(argv, full_argv):
    p = parse_args(argv)
    seed_everything(p["seed"])
    adata, sets, found, listed_genes = load(p)
    scores = score(adata, sets, p)
    means = sample_means(adata, scores, p)
    control, draws = random_control(adata, sets, means, p)

    outdir = make_outdir(p["outdir"])
    tables, figures, data = outdir / "tables", outdir / "figures", outdir / "data"
    data.mkdir()
    sample, group = p["sample_col"], p["group_col"]
    for name in scores.columns:
        adata.obs[f"{p['method']}_{name}"] = scores[name]
    adata.write_h5ad(data / "scored.h5ad")

    cells = pd.concat([adata.obs[[sample, group]], scores], axis=1)
    cells.rename_axis("cell").reset_index().sort_values("cell").to_csv(
        tables / "cell_scores.tsv", sep="\t", index=False
    )
    # group-tests' input names: feature, subject, value; the group keeps its name.
    sample_scores = means.rename(columns={"set": "feature", sample: "subject"})
    sample_scores.sort_values(["feature", "subject"]).to_csv(
        tables / "sample_scores.tsv", sep="\t", index=False
    )
    found.to_csv(tables / "genes_found.tsv", sep="\t", index=False)
    untested = dict.fromkeys(sorted(set(found["set"]) - set(sets)), "not_tested")
    control_all = with_status(control, found["set"], "set", untested).sort_values("set")
    control_all.to_csv(tables / "random_control.tsv", sep="\t", index=False)
    status_counts = control_all["status"].value_counts().sort_index().to_dict()

    plt.switch_backend("Agg")
    draw_figures(figures, adata, scores, means, control, draws, p)

    engine = (
        f"scanpy score_genes, ctrl_size {p['ctrl_size']}, n_bins {p['n_bins']}"
        if p["method"] == "score_genes"
        else f"decoupler AUCell, n_up {p['n_up']}"
    )
    n_samples = means[sample].nunique()
    # plot_box drops the box under FEW_POINTS samples in any group; say which form.
    few = means.groupby(group, observed=True)[sample].nunique().min() < FEW_POINTS
    form = "Dot plot" if few else "Box plot"
    median = " with a bar at each group median" if few else ""
    captions = {}
    for name in scores.columns:
        captions[f"score_by_group_{name}"] = (
            f"{form} of the per-sample mean {name} score by {group}{median}, one "
            f"point per sample ({n_samples} samples); {engine}."
        )
        if "X_umap" in adata.obsm:
            captions[f"score_umap_{name}"] = (
                f"UMAP of {adata.n_obs} cells coloured by the per-cell {name} score, "
                f"highest drawn on top; {engine}; the UMAP is the input's X_umap, "
                "made upstream, not recomputed."
            )
    captions["random_control"] = (
        f"Spread (max minus min) of {group} means of sample-mean scores for each "
        f"signature (red) and {p['n_random']} size-matched random gene sets (grey, "
        f"seed {p['seed']}); a signature beats random when it exceeds every random set."
    )
    # Resource overlap over distinct genes across all sets.
    n_in_resource, n_matched = (
        len(listed_genes),
        len(listed_genes & set(adata.var_names)),
    )
    beats = control.set_index("set")["beats_random"].astype(bool).to_dict()
    missing = {s: m for s, m in found.set_index("set")["missing"].items() if m}
    summary = [
        f"{adata.n_obs} cells in {n_samples} samples and "
        f"{len(adata.obs[group].cat.categories)} groups of {group}; "
        f"{len(sets)} gene sets scored ({engine}).",
        f"Sets not tested (fewer than {p['min_genes']} genes in the data): "
        + (", ".join(untested) or "none")
        + ".",
        "Genes not in the data (not scored; not evidence of absence): "
        + ("; ".join(f"{s}: {m}" for s, m in missing.items()) or "none")
        + ".",
        "Beats size-matched random sets (group spread above all "
        f"{p['n_random']}): "
        + ", ".join(f"{s} {'yes' if b else 'no'}" for s, b in beats.items())
        + ".",
        "Group tests on sample means: hand tables/sample_scores.tsv to group-tests.",
    ]
    keys = [k for k in KEYS if p[k] is not None]
    write_report(outdir, "Gene-signature scores", p, keys, summary, captions)
    write_result(
        outdir,
        {
            "n_cells": adata.n_obs,
            "n_samples": n_samples,
            "n_sets": len(sets),
            "n_matched": n_matched,
            "n_in_resource": n_in_resource,
            "n_found": found.set_index("set")["n_found"].astype(int).to_dict(),
            "missing_genes": {s: m.split(",") for s, m in missing.items()},
            "beats_random": beats,
            "status_counts": status_counts,
            "p_empirical": control.set_index("set")["p_empirical"].to_dict(),
            "method": p["method"],
            "seed": p["seed"],
            "params": p,
        },
        summary,
        ["report.md", "tables/sample_scores.tsv", "figures/random_control.png"],
        context=p["context"],
        reference=p["reference"],
    )
    write_reproducibility(
        outdir,
        __file__,
        p,
        keys,
        full_argv,
        ["anndata", "scanpy", "decoupler", "numpy", "pandas", "scipy", "matplotlib"],
    )
    return outdir


if __name__ == "__main__":
    raise SystemExit(run("signature-scoring", analyse))
