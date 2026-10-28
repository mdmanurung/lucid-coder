---
name: figure-style
description: >-
  House figure style for scientific plots in R (ggplot2, ComplexHeatmap) and
  Python (plotnine, plotnine-extra): use it for ANY figure drawn from
  experimental, omics or single-cell data, including plot types not listed
  here (box, violin, slopegraph, forest, volcano, heatmap, dot plot, UMAP,
  PCA, composition bars, QC, survival curves, enrichment, spatial maps).
  Takes precedence over generic charting or dashboard skills for these
  figures. Do NOT fire for web dashboards, slide graphics or UI charts.
---

# figure-style: one look for every scientific figure

You draw figures the way this lab does, also for a plot type you have never
been shown. Rules come first as principles, so they carry over to new forms;
the form table after them shows how the principles land on the forms the user
has already judged. Exact values and their sources are the `fig.*` rows in
`${CLAUDE_PLUGIN_ROOT}/skills/lucid-coder/references/house-style.md`; the
drawing code is `${CLAUDE_PLUGIN_ROOT}/templates/plots.R` and `plots.py`
(source `theme.*` and `palettes.*` first; copy all four into the project).

## The reader

A colleague reads the figure in a paper or a lab meeting, without the code.
They must see what was measured, on how many units, with which test, from the
figure and its caption alone. Every choice below serves that reader.

## Libraries

- **R:** ggplot2, with `ggpubr::stat_pwc` for brackets, patchwork for layout,
  and ComplexHeatmap with circlize for annotated heatmaps.
- **Python:** plotnine with plotnine-extra (`from plotnine_extra import *`
  also gives all of plotnine). Python figures read like their ggplot2 twin,
  layer for layer. plotnine-extra supplies the following:
  - layout: `|`, `/`, `plot_layout()`, `plot_annotation()`;
  - facets: `facet_wrap2()`, `facet_grid2()`;
  - heatmap dendrograms: `guide_dendro()`;
  - spatial scale bars: `guide_axis_scalebar()`;
  - overplotting: `geom_pointdensity()`.
  Brackets come from `plotnine_extra.stat_pwc`, the same interface as
  `ggpubr::stat_pwc`:
  - `method`: Wilcoxon by default, or `"t.test"`;
  - `p_adjust_method`: for example `"BH"` or `"holm"`;
  - `label="p.adj.format"`;
  - for paired data, `paired=True` with `wid="<subject column>"`;
  - add headroom with `scale_y_continuous(expand=(0.05, 0, 0.15, 0))`.
  Do not use matplotlib, seaborn or statannotations for a new figure.
- The `plotnine` and `ggplot2` skills hold the API references. Never invent an
  argument.

## Before you plot: stop and check

Answer these from the data and the analysis before mapping anything to an
aesthetic (HadleyViz, biology-core). If one stays unknown, say so in the
caption instead of filling it with a common default.

1. **The question.** Distribution, group comparison, association, time course,
   composition, or high-dimensional structure.
2. **The units.** The experimental unit (donor, sample) versus the observation
   unit (cell, read). Note any nesting, pairing or repeated measures.
3. **The identities.** Exact identifiers as the source gives them: gene symbol
   plus Ensembl ID where both exist, and variants in HGVS on a stated transcript
   (MANE Select where possible). Never truncate or remap an identifier silently.
4. **The scale.** Note each transform: log, CPM, z-score, CLR. A log axis needs
   a stated treatment of zeros (a pseudocount, or "zeros not shown, n = k").
5. **The context.** Species, tissue, cell type, condition, time point, assay,
   and genome build for genomic panels. These go in the caption.
6. **What was measured versus derived.** Observed values, model scores,
   predictions and transferred labels are labelled as what they are.

## Eleven principles

1. **Name it, then explain it below.** One short descriptive title: the gene,
   the population, or the comparison ("ISG15 in CD14+ monocytes",
   "CD14+ monocyte pseudobulks (top 500 HVGs)"). No subtitle and no caption in
   the plot. The caption sits in the report under the figure, at most 60 words:
   what is shown, then how (n, test, adjustment, reduction, seed, thresholds).
2. **Show the unit of replication.** Each donor or sample is a visible point
   over any summary. Under 10 per group: points plus a median line, no box.
   Write n under each group label on a second line: `ctrl\n(n = 4)`.
3. **Make inference visible and honest.** A comparison shows its p-value
   bracket with the adjusted p (`stat_pwc` in both languages: ggpubr in R,
   plotnine-extra in Python; paired data use the paired test and the subject
   column). The test and the adjustment come from
   the analysis, never a default, and the caption names both. No stars without
   a named test. Error bars or bands name what they are (SD, SE or 95% CI)
   and appear only when n supports them. A figure claims no more than the
   evidence. An association is never drawn or captioned as a cause, and a
   score is never labelled as a probability. A fitted line does not extend
   past the data.
4. **Colour has one job per figure.**
   - Categories: one palette keyed by label, the same colour for the same label
     in every figure; an unknown label is an error, never a fallback colour.
   - Continuous values: viridis, high values drawn on top.
   - Signed values: diverging blue-white-red, centred at 0, symmetric limits
     held in a named constant, values beyond squished.
   - Direction: blue down, red up, grey for not significant or background.
5. **Encode groups twice.** Colour plus shape (or colour plus position), so the
   figure survives greyscale print and colour-vision deficiency.
6. **Axes carry meaning.** Every axis and legend title names the quantity and
   its unit; never a raw column name. Math is typeset: log₂ via
   `expression(log[2]~"fold change")` or `r"$\log_2$ fold change"`. Bars start
   at 0. Crowded category labels rotate 90°. Embeddings have equal aspect, no
   frame and no ticks, with two short corner arrows ("UMAP 1", "UMAP 2"). PCA
   axes state the variance explained.
7. **Order is information.** Natural orders stay (time, dose, ctrl before
   stim); unordered categories sort by value. A heatmap states its order basis
   in the caption: clustered (metric and linkage) or fixed by label.
8. **Mark the design, not decoration.** Thresholds are dashed lines. Lines
   join paired samples only when the design is paired; never show paired data
   as two independent boxes. Heatmaps carry an annotation bar for each design
   variable (condition, cell type), and clustering is used only when it has a
   biological reason that the caption states. Proportions share one stated
   denominator ("% of cells per sample"). Labels go on the top N items only
   when asked, N per direction.
9. **Distributions match their size.** n < 10: points and median. n ≥ 10:
   violin with a median line and no extrema whiskers, or a box with
   `outlier.shape = NA` under jittered points. Overlapping densities get alpha
   fills. Jitter takes a seed. Too many points to read: rasterise, then use
   alpha, then `geom_pointdensity()`, hexbin or a 2D density.
10. **Missing is not zero.** A gene that was filtered, undetected, below the
    detection limit or not measured gets its own mark or colour, labelled "not
    tested" or "not detected", and never a 0 or the colour of "no change".
11. **One build, three files.** `theme_bw(base_size = 8)` (Python: the plotnine
    theme in `theme.py`); size in mm (89 single column, 183 double);
    `cairo_pdf` plus a PNG for review plus the source-data CSV of the plotted
    frame, all through `save_panel()`. Never a JPEG, never 3D.

## Forms the user has judged

| Form | How the principles land | Template |
|---|---|---|
| Box, few donors | Points + median, n under labels, bracket with adjusted p, gene as title | `plot_box(..., test=, p_adjust=)` |
| Paired slopegraph | Lines join each donor, no jitter, bracket with paired test | `plot_box(..., pair=)` |
| Forest | Point + named interval, dashed line at 0 (or 1 on ratio scales), sorted by estimate, log₂ typeset | same theme, no template |
| Heatmap | Diverging centred at 0, symmetric limits, rotated labels, condition annotation bar, order basis in caption | `plot_heatmap(..., col_annotation=)`, ComplexHeatmap with `colorRamp2(c(-L, 0, L))` |
| PCA | Colour and shape by condition, variance on axes, pair lines only when paired | same theme |
| Composition | Stacked proportions 0–100% per sample, label-keyed palette, rotated sample labels | `plot_composition` |
| Distributions | Violin + median, no extrema | `plot_violin` |
| Volcano | Blue down, red up, grey rest, dashed cut-offs, top N per direction labelled only when asked | `plot_volcano(..., n_label=)` |
| Dot plot | Area = fraction expressing with a grey-circle size legend at 25/50/75/100%, viridis mean with a unit colourbar | `plot_dotplot(..., mean_label=)` |
| UMAP | No frame, corner arrows, label-keyed palette, rasterised points | `plot_embedding(..., axis_name="UMAP")` |
| Feature on UMAP | viridis, high on top, unit in colourbar | `plot_feature` |

## A form not in the table

Run the stop-and-check above first, then work the form out from the data,
not from the form's usual look:

1. **What is one point?** Donor, sample, cell, gene. Show that unit (2, 9).
2. **What varies?** Categories, a continuous value, a signed effect, a
   direction. Pick the colour job (4) and encode groups twice (5).
3. **Is there a comparison?** Then show the test result and name it (3).
4. **Is there a design structure?** Pairing, batch, time: draw it (8).
5. **What is the order?** State it (7).
6. **Is anything missing or below detection?** Give it its own mark (10).
7. **Write the title, axes with units and the caption** (1, 6), and save the
   three files (11).

Worked examples of that procedure:

- **Kaplan–Meier curve.** Step function per group, label-keyed colours plus
  line type; a number-at-risk table under the x axis; CI band with alpha;
  censor marks; log-rank or Cox p in the plot, the test in the caption; x axis
  in days with its unit.
- **Correlation matrix.** Diverging centred at 0 with limits ±1; order by
  clustering, metric and linkage in the caption; the diagonal left blank.
- **Ridgeline of a marker per cluster.** Clusters in a stated order, alpha
  fills, a shared x axis with the unit, a median tick per ridge.
- **Spatial map.** Like an embedding: equal aspect, no frame, a scale bar in
  µm instead of axis arrows; continuous values in viridis.
- **GSEA running score.** Enrichment curve, a barcode rug of set members, NES
  and adjusted p in the plot, the gene set and ranking metric in the caption.

## The model will want to... Do not

- **…accept the library default title, legend title or axis label.** Do not:
  default labels are column names. Write each one.
- **…add a subtitle with n or the test.** Do not: n goes under the group
  labels and the test in the caption.
- **…use a rainbow or the default categorical cycle.** Do not: use the label-keyed
  palette or viridis.
- **…draw a box for 4 donors, or a bar with an error bar.** Do not: points
  plus a median.
- **…pick the test for the brackets itself.** Do not: the test and adjustment
  come from the analysis plan; ask if they are not stated.
- **…label every significant gene.** Do not: label the top N per direction
  only when asked.
- **…put a log axis on data with zeros, or fill missing cells with 0.** Do
  not: state the pseudocount or the dropped count, and draw "not tested" in its
  own colour.
- **…shorten an identifier to fit the axis.** Do not: keep the full symbol,
  accession or HGVS name, and rotate or wrap the labels instead.
- **…reach for matplotlib or seaborn in Python.** Do not: use plotnine with
  plotnine-extra, so the figure matches its R twin.
- **…caption a UMAP distance or a cluster gap as a biological distance.** Do
  not: embedding axes have no units, and the caption says so.

## Where the rules came from

The principles generalise the `fig.*` rows of the house style. The installed
`hadley-viz` skill gives the stop-and-check, the identities, log axes with
zeros, overplotting, denominators and the clustering rationale. The installed
`biology-core` skill gives measured versus derived, non-detection versus
absence, scores that are not probabilities, and context in the caption. The
other sources are the immgenT figure site for minimal in-plot text, Wilke's *Fundamentals of Data
Visualization* and data-to-viz for chart choice, the biology-skills review for
n and order, and the user's picks and notes on the figure gallery
(2026-10-28; all 13 forms picked the house-style reference, with the
refinements above). The bench result that motivated a separate skill: with
only the lucid-coder persona installed, agents loaded Claude Code's generic
dataviz skill for every figure task and never applied the house style.
