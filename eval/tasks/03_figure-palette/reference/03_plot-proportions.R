# Purpose: plot cell-type proportion by condition, one point per donor, in the lab palette.
# Inputs:  data/proportions.tsv   [gap: provided with the task]
#          data/palette.tsv       [gap: lab palette, keyed by label]
# Outputs: <run_dir>/proportions_by-condition.{pdf,png,csv}, <run_dir>/session_info.txt
# Run:     Rscript scripts/03_plot-proportions.R results/run1
# Notes:   [Stochastic] none. Unit: donor. Denominator: all of a donor's cells.

# ---- Setup ----
library(ggplot2)

# Healthy first, so it sits left of disease in every cell-type group.
CONDITION_LEVELS <- c("healthy", "disease")
PANEL <- "proportions_by-condition"
# Single-column journal width.
PANEL_WIDTH_MM <- 89
PANEL_HEIGHT_MM <- 60
REQUIRED_COLS <- c("donor_id", "condition", "cell_type", "proportion")

theme_set(ggpubr::theme_pubr(base_size = 8))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: Rscript scripts/03_plot-proportions.R <run_dir>", call. = FALSE)
}
run_dir <- args[[1]]
if (dir.exists(run_dir)) {
  stop(
    "run_dir ",
    run_dir,
    " already exists; pick a new name so nothing is overwritten.",
    call. = FALSE
  )
}
dir.create(run_dir, recursive = TRUE)

# ---- Load inputs ----
proportions <- read.delim(here::here("data", "proportions.tsv"))
missing_cols <- setdiff(REQUIRED_COLS, names(proportions))
if (length(missing_cols) > 0) {
  stop(
    "proportions.tsv lacks columns ",
    toString(missing_cols),
    ".",
    call. = FALSE
  )
}
palette_table <- read.delim(here::here("data", "palette.tsv"))
palette <- setNames(palette_table$hex, palette_table$label)
# Colours by label, never by position: an unknown condition must stop the script.
unknown_labels <- setdiff(proportions$condition, names(palette))
if (length(unknown_labels) > 0) {
  stop(
    "palette.tsv has no colour for ",
    toString(unknown_labels),
    ".",
    call. = FALSE
  )
}
proportions$condition <- factor(
  proportions$condition,
  levels = CONDITION_LEVELS
)
if (anyNA(proportions$condition)) {
  stop(
    "proportions.tsv has conditions outside CONDITION_LEVELS; recode them at load.",
    call. = FALSE
  )
}

# ---- Plot ----
p_proportions <- ggplot(
  proportions,
  aes(cell_type, proportion, colour = condition)
) +
  geom_point(position = position_dodge(width = 0.6), size = 1) +
  scale_colour_manual(values = palette[CONDITION_LEVELS]) +
  labs(
    x = NULL,
    y = "Proportion of the donor's cells",
    colour = NULL,
    caption = "One point per donor; denominator: all cells of that donor."
  )

# ---- Write outputs ----
stem <- file.path(run_dir, PANEL)
ggsave(
  paste0(stem, ".pdf"),
  p_proportions,
  width = PANEL_WIDTH_MM,
  height = PANEL_HEIGHT_MM,
  units = "mm",
  device = cairo_pdf
)
ggsave(
  paste0(stem, ".png"),
  p_proportions,
  width = PANEL_WIDTH_MM,
  height = PANEL_HEIGHT_MM,
  units = "mm",
  dpi = 300
)
# Source data: the same frame passed to ggplot().
write.csv(proportions, paste0(stem, ".csv"), row.names = FALSE)

# ---- Session info ----
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "session_info.txt")
)
