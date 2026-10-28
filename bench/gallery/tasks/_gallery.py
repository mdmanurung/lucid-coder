"""House-style reference figures for the gallery tasks; each reproducer calls one function.

Every function draws in the cwd: reference.pdf, reference.png and reference.csv (source
data), plus reference.R for the R tasks, then writes answer.md as `done`. Python figures
go through templates/plots.py and theme.py; R figures source theme.R, palettes.R and
plots.R. forest, slope-free forms (density, qc violins, feature plot, PCA) and the
ComplexHeatmap panel have no template, so they are drawn here on the same theme.
"""

import subprocess
import sys
from pathlib import Path

TEMPLATES = Path(__file__).resolve().parents[3] / "templates"
sys.path.insert(0, str(TEMPLATES))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
# Constrained layout from creation, so outside legends and colourbars fit the panel.
matplotlib.rcParams["figure.constrained_layout.use"] = True
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from palettes import palette_for  # noqa: E402
from plots import (  # noqa: E402
    grey_to_red,
    plot_composition,
    plot_dotplot,
    plot_embedding,
    plot_volcano,
)
from scipy.stats import gaussian_kde  # noqa: E402
from theme import save_panel  # noqa: E402

# Palettes keyed by label (fig.palette); key order is the plotting order.
PAL_CONDITION = {"ctrl": "#4477AA", "stim": "#EE6677"}
PAL_CELL_TYPE = {
    "CD4 T cells": "#332288",
    "CD8 T cells": "#88CCEE",
    "NK cells": "#44AA99",
    "B cells": "#117733",
    "CD14+ Monocytes": "#CC6677",
    "FCGR3A+ Monocytes": "#882255",
    "Dendritic cells": "#DDCC77",
    "Megakaryocytes": "#999933",
}
PAL_SAMPLE = {"pbmc_1k_v3": "#4477AA", "pbmc3k": "#228833", "pbmc_10k_v3": "#AA3377"}
LIMIT = 2  # heatmap z-score saturation (fig.heatmap)
DIVERGING = ["#2166AC", "white", "#B2182B"]


def done():
    Path("answer.md").write_text("done\n")


def save(fig, data, width_mm, height_mm):
    save_panel(fig, data, ".", "reference", width_mm, height_mm)
    done()


def r_vec(name, pal):
    pairs = ", ".join(f'`{k}` = "{v}"' for k, v in pal.items())
    return f"{name} <- c({pairs})\n"


def run_r(body):
    """Write reference.R (templates sourced, Kang palettes set) and run it."""
    head = "".join(f'source("{TEMPLATES / f}")\n' for f in ("theme.R", "palettes.R", "plots.R"))
    code = head + r_vec("pal_condition", PAL_CONDITION) + r_vec("pal_cell_type", PAL_CELL_TYPE)
    Path("reference.R").write_text(code + body)
    # Plain `Rscript` is the target venv's shim: acumen check puts the venv bin first.
    subprocess.run(["Rscript", "reference.R"], check=True)
    done()


# R tasks ---------------------------------------------------------------------


def box_few(gene, cell_type):
    run_r(f"""
d <- read.delim("box_few.tsv")
d <- d[d$gene == "{gene}" & d$cell_type == "{cell_type}", ]
# n < FEW_POINTS per group, so plot_box draws points and a median, no box.
p <- plot_box(d, "condition", "mean_expr", pal_condition, jitter_width = 0.1, seed = 1,
  value_label = "{gene}, mean log1p(CP10k)")
save_panel(p, d, ".", "reference", 60, 60)
""")


def slope(gene, cell_type):
    run_r(f"""
d <- read.delim("sample_means.tsv")
d <- d[d$gene == "{gene}" & d$cell_type == "{cell_type}", ]
# Paired donors: a line per donor, points, and the median per condition.
p <- plot_box(d, "label", "mean_expr", pal_condition, jitter_width = 0, seed = 1,
  value_label = "{gene}, mean log1p(CP10k)", pair = "replicate")
save_panel(p, d, ".", "reference", 60, 60)
""")


def forest(cell_type):
    run_r(f"""
d <- read.delim("de.tsv")
d <- d[d$cell_type == "{cell_type}", ]
d <- head(d[order(d$padj), ], 15)
d$lower <- d$log2fc - 1.96 * d$se
d$upper <- d$log2fc + 1.96 * d$se
# Genes ordered by estimate (fig.wilke: unordered categories by value).
d$gene <- reorder(factor(d$gene), d$log2fc)
p <- ggplot(d, aes(log2fc, gene)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = lower, xmax = upper), size = 0.2) +
  labs(x = "log2 fold change, stim vs ctrl (95% CI)", y = NULL)
save_panel(p, d, ".", "reference", 89, 80)
""")


def heatmap(cell_type):
    run_r(f"""
d <- read.delim("heatmap_expr.tsv")
d <- d[d$cell_type == "{cell_type}", ]
d$sample <- paste(d$replicate, d$label, sep = "_")
# Row z-score: each gene centred and scaled across this cell type's samples.
d$`z-score` <- ave(d$log2cpm, d$gene, FUN = function(x) (x - mean(x)) / sd(x))
# Order fixed by label: genes by mean stim - ctrl, samples ctrl block then stim.
shift <- tapply(d$log2cpm * ifelse(d$label == "stim", 1, -1), d$gene, mean)
d$gene <- factor(d$gene, levels = names(sort(shift)))
d$sample <- factor(d$sample, levels = unique(d$sample[order(d$label, d$replicate)]))
p <- plot_heatmap(d, "gene", "sample", "z-score", limit = {LIMIT})
save_panel(p, d, ".", "reference", 120, 110)
""")


def pseudobulk_pca(slug):
    run_r(f"""
m <- as.matrix(read.delim("pb_logcpm_{slug}.tsv", row.names = 1, check.names = FALSE))
top <- head(order(apply(m, 1, var), decreasing = TRUE), 500)
pca <- prcomp(t(m[top, ]))
explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
d <- data.frame(
  sample = colnames(m),
  replicate = sub("_[^_]+$", "", colnames(m)),
  label = sub(".*_", "", colnames(m)),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)
# Paired donors joined by a line; colour carries the condition.
p <- ggplot(d, aes(PC1, PC2)) +
  geom_line(aes(group = replicate), colour = "grey60") +
  geom_point(aes(colour = label)) +
  scale_colour_manual(values = pal_condition) +
  labs(
    x = sprintf("PC1 (%.0f%% of variance)", explained[1]),
    y = sprintf("PC2 (%.0f%% of variance)", explained[2]),
    colour = "condition"
  )
save_panel(p, d, ".", "reference", 89, 70)
""")


def complexheatmap(cell_types):
    cts = ", ".join(f'"{ct}"' for ct in cell_types)
    run_r(f"""
suppressPackageStartupMessages({{
  library(ComplexHeatmap)
  library(circlize)
}})
d <- read.delim("heatmap_expr.tsv")
d <- d[d$cell_type %in% c({cts}), ]
# Columns fixed by label: cell type (palette order), then condition, then donor.
d$cell_type <- factor(d$cell_type, levels = intersect(names(pal_cell_type), d$cell_type))
d <- d[order(d$cell_type, d$label, d$replicate), ]
d$column <- factor(paste(d$cell_type, d$replicate, d$label), levels = unique(paste(d$cell_type, d$replicate, d$label)))
m <- tapply(d$log2cpm, list(d$gene, d$column), identity)
z <- t(scale(t(m)))
cols <- unique(d[c("column", "cell_type", "label")])
ha <- HeatmapAnnotation(
  condition = cols$label,
  `cell type` = as.character(cols$cell_type),
  col = list(condition = pal_condition, `cell type` = pal_cell_type[levels(cols$cell_type)])
)
# Rows clustered: Euclidean distance, complete linkage; colours squish beyond +-{LIMIT}.
ht <- Heatmap(
  z,
  name = "z-score",
  col = colorRamp2(c(-{LIMIT}, 0, {LIMIT}), c("{DIVERGING[0]}", "white", "{DIVERGING[2]}")),
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 6),
  top_annotation = ha
)
cairo_pdf("reference.pdf", width = 183 / 25.4, height = 120 / 25.4)
draw(ht)
invisible(dev.off())
png("reference.png", width = 183, height = 120, units = "mm", res = 300)
draw(ht)
invisible(dev.off())
write.csv(d, "reference.csv", row.names = FALSE)
""")


# Python tasks ----------------------------------------------------------------


def composition(labels):
    df = pd.read_csv("celltype_counts.tsv", sep="\t")
    df = df[df["label"].isin(labels)].rename(columns={"n_cells": "cells"})
    order = df.sort_values(["label", "replicate"])["sample"].unique()
    df["sample"] = pd.Categorical(df["sample"], categories=order)
    fig = plot_composition(df, "sample", "cell_type", "cells", PAL_CELL_TYPE)
    fig.axes[0].tick_params(axis="x", labelrotation=90)
    save(fig, df, 183, 90)


def density(cell_types, value):
    xlabel = {
        "total_counts": "UMI counts per cell (log10)",
        "n_genes": "genes detected per cell (log10)",
    }[value]
    df = pd.read_csv("cell_qc.tsv", sep="\t")
    df = df[df["cell_type"].isin(cell_types)].copy()
    df["log10_value"] = np.log10(df[value])
    grid = np.linspace(df["log10_value"].min(), df["log10_value"].max(), 300)
    fig, ax = plt.subplots()
    levels = [ct for ct in PAL_CELL_TYPE if ct in cell_types]
    # Overlapping densities: alpha on the fills (fig.wilke).
    for ct, colour in zip(levels, palette_for(levels, PAL_CELL_TYPE), strict=True):
        values = df.loc[df["cell_type"] == ct, "log10_value"]
        y = gaussian_kde(values)(grid)
        ax.fill_between(grid, y, color=colour, alpha=0.3, lw=0)
        ax.plot(grid, y, c=colour, lw=1, label=f"{ct} (n = {len(values)})")
    ax.set_ylim(bottom=0)
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("density")
    save(fig, df, 120, 60)


def volcano(cell_type):
    df = pd.read_csv("de.tsv", sep="\t")
    df = df[df["cell_type"] == cell_type]
    fig = plot_volcano(df, "log2fc", "padj", padj_max=0.05, lfc_min=1)
    save(fig, df, 89, 75)


def dotplot(genes):
    df = pd.read_csv("marker_summary.tsv", sep="\t")
    df = df[df["gene"].isin(genes)].rename(
        columns={"mean_expr": "mean log1p(CP10k)", "frac_expressing": "cells expressing"}
    )
    df["gene"] = pd.Categorical(df["gene"], categories=genes)
    df["cell_type"] = pd.Categorical(df["cell_type"], categories=list(PAL_CELL_TYPE)[::-1])
    fig = plot_dotplot(
        df, "gene", "cell_type", "mean log1p(CP10k)", "cells expressing", max_size=4
    )
    save(fig, df, 140, 80)


def umap(cell_types):
    df = pd.read_csv("pbmc3k_umap.tsv", sep="\t")
    df = df[df["cell_type"].isin(cell_types)]
    fig = plot_embedding(df, "umap_1", "umap_2", "cell_type", PAL_CELL_TYPE, axis_name="UMAP")
    save(fig, df, 140, 89)


def marker_feature(gene):
    df = pd.read_csv("pbmc3k_markers.tsv", sep="\t")[["cell", "umap_1", "umap_2", gene]]
    df = df.sort_values(gene)  # expressing cells drawn last, on top
    fig, ax = plt.subplots()
    dots = ax.scatter(
        df["umap_1"], df["umap_2"], s=1, c=df[gene], cmap=grey_to_red, rasterized=True
    )
    fig.colorbar(dots, ax=ax, label=f"{gene}, log1p(CP10k)")
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    save(fig, df, 110, 89)


def qc_violin(thresholds):
    """thresholds: {metric: [cut-offs]}; one panel per metric, samples share the x axis."""
    ylabel = {
        "n_genes": "genes detected per cell",
        "total_counts": "UMI counts per cell",
        "pct_mito": "mitochondrial counts (%)",
    }
    df = pd.read_csv("qc_metrics.tsv", sep="\t")
    levels = [s for s in PAL_SAMPLE if s in set(df["sample"])]
    n = df["sample"].value_counts()
    fig, axes = plt.subplots(len(thresholds), 1, sharex=True)
    for ax, (metric, cuts) in zip(np.atleast_1d(axes), thresholds.items(), strict=True):
        parts = ax.violinplot(
            [df.loc[df["sample"] == s, metric] for s in levels],
            positions=range(len(levels)),
            showmedians=True,
            showextrema=False,
        )
        for body, colour in zip(parts["bodies"], palette_for(levels, PAL_SAMPLE), strict=True):
            body.set_facecolor(colour)
            body.set_alpha(0.5)
        parts["cmedians"].set_color("black")
        for cut in cuts:
            ax.axhline(cut, ls="--", lw=0.5, c="black")
        ax.set_ylabel(ylabel[metric])
    ax.set_xticks(range(len(levels)), [f"{s}\n(n = {n[s]})" for s in levels])
    save(fig, df, 89, 45 * len(thresholds))
