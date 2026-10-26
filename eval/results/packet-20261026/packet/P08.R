# Purpose: Plot cell-type proportion by condition, one point per donor, in lab colours.
# Inputs:  data/proportions.tsv (donor_id, condition, cell_type, proportion) [gap: producer unknown]
#          data/palette.tsv (label, hex; one row per cell type and per condition) [gap]
# Outputs: <run_dir>/proportions_by-condition.pdf
#          <run_dir>/proportions_by-condition.png
#          <run_dir>/proportions_by-condition.csv (source data for the panel)
#          <run_dir>/session_info.txt
# Run:     Rscript scripts/03_plot-proportions.R results/run1
# Notes:   Points are coloured by condition; the crossbar is the median per condition.
#          Each point is one donor, so n = 4 donors per condition per cell type.
#          janitor is not in the allowed package set, so column names are used as-is.

# Setup ------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)

# Healthy first so it is the reference group on the left of each pair.
CONDITION_LEVELS <- c("healthy", "disease")
# Fixed jitter seed so donor point positions are identical on every rerun.
JITTER_SEED <- 1
# Small horizontal jitter keeps donors visible without crossing into the other condition.
JITTER_WIDTH <- 0.08
DODGE_WIDTH <- 0.6

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("args: expected one run directory; run `Rscript scripts/03_plot-proportions.R results/run1`.",
    call. = FALSE
  )
}
run_dir <- here::here(args[[1]])
if (dir.exists(run_dir)) {
  stop("run_dir: ", run_dir, " already exists; choose a new run directory.", call. = FALSE)
}

theme_set(theme_bw(base_size = 8))

# Load inputs ------------------------------------------------------------------

proportions <- read_tsv(
  here::here("data", "proportions.tsv"),
  col_types = cols(
    donor_id = col_character(),
    condition = col_character(),
    cell_type = col_character(),
    proportion = col_double()
  )
)
missing_cols <- setdiff(c("donor_id", "condition", "cell_type", "proportion"), names(proportions))
if (length(missing_cols) > 0) {
  stop("proportions: missing columns ", toString(missing_cols), "; check the input file.",
    call. = FALSE
  )
}

palette <- read_tsv(
  here::here("data", "palette.tsv"),
  col_types = cols(label = col_character(), hex = col_character())
)
if (anyDuplicated(palette$label) > 0) {
  stop("palette: duplicated labels; keep one hex per label.", call. = FALSE)
}
palette_hex <- setNames(palette$hex, palette$label)

# Prepare plot data ------------------------------------------------------------

# An unlisted condition would otherwise become NA and silently vanish from the plot.
unknown_conditions <- setdiff(unique(proportions$condition), CONDITION_LEVELS)
if (length(unknown_conditions) > 0) {
  stop("proportions: unknown condition ", toString(unknown_conditions),
    "; add it to CONDITION_LEVELS.",
    call. = FALSE
  )
}
# Every plotted condition needs a lab colour; ggplot would otherwise draw it grey.
unpaletted_conditions <- setdiff(CONDITION_LEVELS, names(palette_hex))
if (length(unpaletted_conditions) > 0) {
  stop("palette: no colour for ", toString(unpaletted_conditions), "; add it to palette.tsv.",
    call. = FALSE
  )
}

# Cell types ordered by overall median proportion, most abundant first.
cell_type_order <- proportions |>
  group_by(cell_type) |>
  summarise(median_proportion = median(proportion), .groups = "drop") |>
  arrange(desc(median_proportion)) |>
  pull(cell_type)

plot_data <- proportions |>
  mutate(
    condition = factor(condition, levels = CONDITION_LEVELS),
    cell_type = factor(cell_type, levels = cell_type_order)
  ) |>
  arrange(cell_type, condition, donor_id)

# Plot -------------------------------------------------------------------------

n_donors <- plot_data |>
  distinct(donor_id, condition) |>
  count(condition)
condition_labels <- setNames(
  paste0(n_donors$condition, " (n = ", n_donors$n, " donors)"),
  n_donors$condition
)

p <- ggplot(plot_data, aes(x = cell_type, y = proportion, colour = condition)) +
  stat_summary(
    fun = median, fun.min = median, fun.max = median,
    geom = "crossbar", width = 0.4, linewidth = 0.3,
    position = position_dodge(width = DODGE_WIDTH),
    show.legend = FALSE
  ) +
  geom_point(
    size = 1.2,
    position = position_jitterdodge(
      jitter.width = JITTER_WIDTH, jitter.height = 0,
      dodge.width = DODGE_WIDTH, seed = JITTER_SEED
    )
  ) +
  scale_colour_manual(
    values = palette_hex[CONDITION_LEVELS],
    breaks = CONDITION_LEVELS,
    labels = condition_labels[CONDITION_LEVELS]
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL,
    y = "Proportion of donor's cells",
    colour = NULL,
    caption = "Points: donors. Bars: median per condition."
  ) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust = 1))

# Write outputs ----------------------------------------------------------------

dir.create(run_dir, recursive = TRUE)

write_csv(plot_data, file.path(run_dir, "proportions_by-condition.csv"))
ggsave(
  file.path(run_dir, "proportions_by-condition.pdf"),
  plot = p, width = 89, height = 70, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(run_dir, "proportions_by-condition.png"),
  plot = p, width = 89, height = 70, units = "mm", dpi = 300
)

# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
