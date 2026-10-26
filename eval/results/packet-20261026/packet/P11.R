# Purpose: Plot cell-type proportion by condition, one point per donor, in lab colours.
# Inputs:  data/proportions.tsv (donor_id, condition, cell_type, proportion) [gap: producer unknown]
#          data/palette.tsv (label, hex) [gap: producer unknown]
# Outputs: <out_dir>/proportions_by-condition.pdf
#          <out_dir>/proportions_by-condition.png
#          <out_dir>/proportions_by-condition.csv (source data for the figure)
#          <out_dir>/session_info.txt
# Run:     Rscript scripts/03_plot-proportions.R results/run1
# Notes:   Proportion is the share of each donor's cells (denominator: all cells of that donor).
#          Points are coloured by condition; the bar marks the condition median.

options(warn = 1)

# Setup ------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)

# healthy first so it is the reference on the left of each panel
CONDITION_LEVELS <- c("healthy", "disease")
# lineage order: lymphoid, then myeloid
CELL_TYPE_LEVELS <- c("t_cell", "b_cell", "nk_cell", "monocyte", "dendritic_cell")
CELL_TYPE_LABELS <- c(
  t_cell = "T cell",
  b_cell = "B cell",
  nk_cell = "NK cell",
  monocyte = "Monocyte",
  dendritic_cell = "Dendritic cell"
)
# fixed seed so the horizontal jitter of donor points is identical on every rerun
JITTER_SEED <- 1
FIG_WIDTH_MM <- 183
FIG_HEIGHT_MM <- 70
PNG_DPI <- 300

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("args: expected one run directory; pass e.g. results/run1.", call. = FALSE)
}
out_dir <- here::here(args[[1]])
if (dir.exists(out_dir)) {
  stop("out_dir: ", out_dir, " already exists; choose a new run directory.", call. = FALSE)
}
dir.create(out_dir, recursive = TRUE)

theme_set(theme_bw(base_size = 8))

# Load inputs ------------------------------------------------------------------

proportions_raw <- read_tsv(here::here("data", "proportions.tsv"), show_col_types = FALSE)
palette <- read_tsv(here::here("data", "palette.tsv"), show_col_types = FALSE)

unknown_labels <- setdiff(
  c(proportions_raw$condition, proportions_raw$cell_type),
  c(CONDITION_LEVELS, CELL_TYPE_LEVELS)
)
if (length(unknown_labels) > 0) {
  stop("proportions: unknown labels ", toString(unknown_labels), "; add to levels.", call. = FALSE)
}

# Colour lookup ----------------------------------------------------------------

condition_colours <- setNames(palette$hex, palette$label)[CONDITION_LEVELS]
if (anyNA(condition_colours)) {
  missing <- CONDITION_LEVELS[is.na(condition_colours)]
  stop("palette: no colour for ", toString(missing), "; add it to data/palette.tsv.", call. = FALSE)
}
names(condition_colours) <- CONDITION_LEVELS

# Build plot data --------------------------------------------------------------

plot_data <- proportions_raw |>
  mutate(
    condition = factor(condition, levels = CONDITION_LEVELS),
    cell_type = factor(cell_type, levels = CELL_TYPE_LEVELS),
    hex = unname(condition_colours[as.character(condition)])
  ) |>
  arrange(cell_type, condition, donor_id)

n_donors <- plot_data |>
  distinct(donor_id, condition) |>
  count(condition)
caption <- paste0(
  "Proportion = share of each donor's cells. One point per donor (",
  paste0(n_donors$condition, " n = ", n_donors$n, collapse = ", "),
  "); bar = median."
)

# Plot -------------------------------------------------------------------------

p <- ggplot(plot_data, aes(x = condition, y = proportion, colour = condition)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.3, colour = "grey40") +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = JITTER_SEED), size = 1.5) +
  facet_wrap(~cell_type, nrow = 1, scales = "free_y", labeller = as_labeller(CELL_TYPE_LABELS)) +
  scale_colour_manual(values = condition_colours, guide = "none") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Proportion of donor's cells", caption = caption)

# Write outputs ----------------------------------------------------------------

ggsave(
  file.path(out_dir, "proportions_by-condition.pdf"),
  plot = p, width = FIG_WIDTH_MM, height = FIG_HEIGHT_MM, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(out_dir, "proportions_by-condition.png"),
  plot = p, width = FIG_WIDTH_MM, height = FIG_HEIGHT_MM, units = "mm", dpi = PNG_DPI
)
write_csv(
  select(plot_data, donor_id, condition, cell_type, proportion, hex),
  file.path(out_dir, "proportions_by-condition.csv")
)

# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
