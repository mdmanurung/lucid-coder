"""The seven common plot forms, drawn one fixed way; plots.R is the R twin.

Each function takes a long DataFrame and column names, returns a Figure and
never saves it: the caller imports theme and passes the Figure to save_panel().
No title, subtitle or caption is drawn; the caption goes in the report (fig.text).
Colours and line styles are fixed house style. Every tunable value is a
required keyword argument; a grill asks the user for each one per project:
  padj_max:     which adjusted-p cutoff marks a gene significant?
  lfc_min:      which absolute log2 fold change must a significant gene exceed?
  limit:        at which absolute value should the heatmap colour scale saturate?
  seed:         which seed fixes the box-plot jitter?
  jitter_width: how far (in box widths, each side) may box-plot points jitter?
  max_size:     what dot size (ggplot2 size units, mm) does frac = 1 get?
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import PercentFormatter
from scipy import stats

from palettes import palette_for

diverging = LinearSegmentedColormap.from_list(
    "diverging", ["#2166AC", "white", "#B2182B"]
)
grey_to_red = LinearSegmentedColormap.from_list("grey_to_red", ["#E5E5E5", "#B2182B"])


def check_cols(data, cols, fn):
    missing = [col for col in cols if col not in data.columns]
    if missing:
        raise KeyError(f"{fn}: missing columns: {', '.join(missing)}")


def levels_of(values, pal=None):
    """Palette key order if pal is given, else Categorical order, else first seen."""
    if pal is not None:
        palette_for(values.unique(), pal)
        return [key for key in pal if key in set(values.astype(str))]
    if isinstance(values.dtype, pd.CategoricalDtype):
        return list(values.cat.categories)
    return list(values.unique())


def plot_volcano(data, lfc, padj, *, padj_max, lfc_min, label=None):
    check_cols(data, [lfc, padj] + ([label] if label else []), "plot_volcano")
    significant = (data[padj] < padj_max) & (data[lfc].abs() > lfc_min)
    neg_log_p = -np.log10(data[padj])
    fig, ax = plt.subplots()
    ax.scatter(data[lfc], neg_log_p, s=4, c=np.where(significant, "darkred", "0.6"))
    ax.axhline(-np.log10(padj_max), ls="--", lw=0.5, c="black")
    for x in (-lfc_min, lfc_min):
        ax.axvline(x, ls="--", lw=0.5, c="black")
    if label:
        for text, x, y in zip(
            data.loc[significant, label],
            data.loc[significant, lfc],
            neg_log_p[significant],
            strict=True,
        ):
            ax.annotate(text, (x, y), xytext=(2, 2), textcoords="offset points")
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p")
    return fig


def plot_box(data, group, value, pal, *, jitter_width, seed, pair=None):
    check_cols(data, [group, value] + ([pair] if pair else []), "plot_box")
    levels = levels_of(data[group], pal)
    x = pd.Categorical(data[group].astype(str), categories=levels).codes
    fig, ax = plt.subplots()
    boxes = ax.boxplot(
        [data.loc[x == i, value] for i in range(len(levels))],
        positions=range(len(levels)),
        showfliers=False,
        patch_artist=True,
    )
    for box, colour in zip(boxes["boxes"], palette_for(levels, pal), strict=True):
        box.set_facecolor(colour)
    if pair:
        # Paired points stay on the box centre so the joining lines are readable.
        wide = data.pivot(index=pair, columns=group, values=value)[levels]
        ax.plot(range(len(levels)), wide.T, c="grey", lw=0.5)
        ax.scatter(x, data[value], s=4, c="black", zorder=3)
    else:
        rng = np.random.default_rng(seed)
        jitter = rng.uniform(-jitter_width, jitter_width, len(data))
        ax.scatter(x + jitter, data[value], s=4, c="black", zorder=3)
    counts = data[group].astype(str).value_counts()
    ax.set_xticks(range(len(levels)), [f"{lv} (n = {counts[lv]})" for lv in levels])
    ax.set_xlabel(group)
    ax.set_ylabel(value)
    return fig


def plot_embedding(data, x, y, colour, pal):
    check_cols(data, [x, y, colour], "plot_embedding")
    fig, ax = plt.subplots()
    ax.scatter(data[x], data[y], s=1, c=palette_for(data[colour], pal), rasterized=True)
    handles = [
        plt.Line2D([], [], marker="o", ls="", c=pal[lv], label=lv)
        for lv in levels_of(data[colour], pal)
    ]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    return fig


def plot_heatmap(data, row, col, value, *, limit):
    check_cols(data, [row, col, value], "plot_heatmap")
    rows, cols = levels_of(data[row]), levels_of(data[col])
    matrix = data.pivot(index=row, columns=col, values=value).reindex(
        index=rows, columns=cols
    )
    fig, ax = plt.subplots()
    # Clip to +-limit so outliers saturate instead of stretching the scale (squish).
    mesh = ax.pcolormesh(
        matrix.clip(-limit, limit), cmap=diverging, vmin=-limit, vmax=limit
    )
    ax.set_xticks(np.arange(len(cols)) + 0.5, cols, rotation=90)
    ax.set_yticks(np.arange(len(rows)) + 0.5, rows)
    ax.invert_yaxis()
    fig.colorbar(mesh, ax=ax, label=value)
    return fig


def plot_dotplot(data, feature, group, mean, frac, *, max_size):
    check_cols(data, [feature, group, mean, frac], "plot_dotplot")
    features, groups = levels_of(data[feature]), levels_of(data[group])
    # ggplot2 size is a diameter in mm; scatter s is area in pt^2 (72.27 pt/25.4 mm).
    max_area = (max_size * 72.27 / 25.4) ** 2
    fig, ax = plt.subplots()
    dots = ax.scatter(
        pd.Categorical(data[feature], categories=features).codes,
        pd.Categorical(data[group], categories=groups).codes,
        s=data[frac] * max_area,
        c=data[mean],
        cmap=grey_to_red,
    )
    ax.set_xticks(range(len(features)), features, rotation=90)
    ax.set_yticks(range(len(groups)), groups)
    fig.colorbar(dots, ax=ax, label=mean)
    handles = [
        ax.scatter([], [], s=f * max_area, c="grey", label=f"{f:.0%}")
        for f in (0.25, 0.5, 1)
    ]
    ax.legend(handles=handles, title=frac, loc="upper left", bbox_to_anchor=(1.3, 1))
    return fig


def plot_composition(data, sample, label, count, pal):
    check_cols(data, [sample, label, count], "plot_composition")
    labels = levels_of(data[label], pal)
    props = pd.crosstab(
        data[sample], data[label], values=data[count], aggfunc="sum", normalize="index"
    ).reindex(index=levels_of(data[sample]), columns=labels, fill_value=0)
    fig, ax = plt.subplots()
    bottom = np.zeros(len(props))
    for lv, colour in zip(labels, palette_for(labels, pal), strict=True):
        ax.bar(
            props.index.astype(str), props[lv], bottom=bottom, color=colour, label=lv
        )
        bottom += props[lv].to_numpy()
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_xlabel(sample)
    ax.set_ylabel(f"% of {count}")
    return fig


def plot_scatter_fit(data, x, y, label=None):
    check_cols(data, [x, y] + ([label] if label else []), "plot_scatter_fit")
    xs, ys = data[x].to_numpy(float), data[y].to_numpy(float)
    slope, intercept = np.polyfit(xs, ys, 1)
    # 95% band for the fitted mean, as geom_smooth(method = "lm", se = TRUE), n = 80.
    n = len(xs)
    sigma = np.sqrt(np.sum((ys - intercept - slope * xs) ** 2) / (n - 2))
    grid = np.linspace(xs.min(), xs.max(), 80)
    se = sigma * np.sqrt(
        1 / n + (grid - xs.mean()) ** 2 / np.sum((xs - xs.mean()) ** 2)
    )
    half = stats.t.ppf(0.975, n - 2) * se
    fit = intercept + slope * grid
    fig, ax = plt.subplots()
    ax.fill_between(grid, fit - half, fit + half, color="0.6", alpha=0.4, lw=0)
    ax.plot(grid, fit, c="#3366FF")
    ax.scatter(xs, ys, s=4, c="black")
    if label:
        for text, px, py in zip(data[label], xs, ys, strict=True):
            ax.annotate(text, (px, py), xytext=(2, 2), textcoords="offset points")
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    return fig


if __name__ == "__main__":
    import tempfile
    from pathlib import Path

    from palettes import pal_cell_type, pal_condition
    from theme import save_panel

    plt.switch_backend("Agg")
    rng = np.random.default_rng(1)
    genes = pd.DataFrame(
        {
            "gene": [f"G{i}" for i in range(20)],
            "log2_fc": rng.normal(0, 2, 20),
            "padj": rng.uniform(0.0001, 1, 20),
        }
    )
    donors = pd.DataFrame(
        {
            "donor": [f"D{i}" for i in range(6)] * 2,
            "condition": ["healthy"] * 6 + ["disease"] * 6,
            "score": rng.normal(0, 1, 12),
        }
    )
    cells = pd.DataFrame(
        {
            "umap_1": rng.normal(0, 1, 60),
            "umap_2": rng.normal(0, 1, 60),
            "cell_type": ["t_cell", "b_cell", "monocyte"] * 20,
        }
    )
    grid = pd.DataFrame(
        {
            "gene": ["CD3E", "MS4A1", "CD14"] * 3,
            "cell_type": np.repeat(["t_cell", "b_cell", "monocyte"], 3),
            "z": rng.normal(0, 3, 9),
            "frac": rng.uniform(0, 1, 9),
        }
    )
    counts = pd.DataFrame(
        {
            "sample": np.repeat(["S1", "S2"], 3),
            "cell_type": ["t_cell", "b_cell", "monocyte"] * 2,
            "n_cells": [50, 30, 20, 10, 60, 30],
        }
    )
    figures = {
        "volcano": (
            plot_volcano(
                genes, "log2_fc", "padj", padj_max=0.1, lfc_min=0.5, label="gene"
            ),
            genes,
        ),
        "box": (
            plot_box(
                donors, "condition", "score", pal_condition, seed=1, jitter_width=0.2
            ),
            donors,
        ),
        "box-paired": (
            plot_box(
                donors,
                "condition",
                "score",
                pal_condition,
                pair="donor",
                seed=1,
                jitter_width=0.2,
            ),
            donors,
        ),
        "embedding": (
            plot_embedding(cells, "umap_1", "umap_2", "cell_type", pal_cell_type),
            cells,
        ),
        "heatmap": (plot_heatmap(grid, "cell_type", "gene", "z", limit=2), grid),
        "dotplot": (
            plot_dotplot(grid, "gene", "cell_type", "z", "frac", max_size=4),
            grid,
        ),
        "composition": (
            plot_composition(counts, "sample", "cell_type", "n_cells", pal_cell_type),
            counts,
        ),
        "scatter-fit": (plot_scatter_fit(genes, "log2_fc", "padj"), genes),
    }
    with tempfile.TemporaryDirectory() as run_dir:
        for name, (fig, data) in figures.items():
            save_panel(fig, data, run_dir, name, 89, 60)
            plt.close(fig)
            assert (Path(run_dir) / f"{name}.pdf").exists(), name
    for bad_call in (
        lambda: plot_volcano(genes, "lfc", "padj", padj_max=0.1, lfc_min=0.5),
        lambda: plot_box(
            donors, "condition", "score", pal_cell_type, seed=1, jitter_width=0.2
        ),
    ):
        try:
            plt.close(bad_call())
        except KeyError:
            pass
        else:
            raise SystemExit("a missing column or unknown label did not raise")
