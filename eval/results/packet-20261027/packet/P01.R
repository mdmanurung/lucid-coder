#$ ===== Reviewer: type after each colon, then add #$ comments anywhere below =====
#$ q1 one read (1-5; can you say what this script does after one read?): 
#$ q2 find decisions (1-5; can you find each threshold, filter, model choice and path?): 
#$ q3 change checks (yes/no; would you add, remove or change a check?): 
#$ q3 which (one sentence): 
#$ ---- your #$ comments: on their own line above the code, or after the code on its line ----
# Plot cell-type proportions by condition, one point per donor, in lab colours.
#
# Usage: Rscript scripts/03_plot-proportions.R <run_dir>
#
# Inputs:  data/proportions.tsv, data/palette.tsv
# Outputs: <run_dir>/proportions_by-condition.{pdf,png,csv}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(forcats)
  library(ggplot2)
  library(here)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript scripts/03_plot-proportions.R <run_dir>", call. = FALSE)
}
run_dir <- args[[1]]
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Read inputs ------------------------------------------------------------

proportions <- read_tsv(
  here("data", "proportions.tsv"),
  col_types = cols(
    donor_id = col_character(),
    condition = col_character(),
    cell_type = col_character(),
    proportion = col_double()
  )
)

palette <- read_tsv(
  here("data", "palette.tsv"),
  col_types = cols(label = col_character(), hex = col_character())
)

# ---- Checks -----------------------------------------------------------------

conditions <- c("healthy", "disease")
stopifnot(
  !anyNA(proportions),
  all(proportions$condition %in% conditions),
  !anyDuplicated(proportions[c("donor_id", "cell_type")]),
  # each donor belongs to exactly one condition
  all(count(distinct(proportions, donor_id, condition), donor_id)$n == 1)
)

# Colours are looked up by label, not by row order in palette.tsv
missing_labels <- setdiff(conditions, palette$label)
if (length(missing_labels) > 0) {
  stop("palette.tsv has no colour for: ", paste(missing_labels, collapse = ", "))
}
condition_colours <- setNames(palette$hex, palette$label)[conditions]

donor_sums <- proportions |>
  group_by(donor_id) |>
  summarise(total = sum(proportion), .groups = "drop")
if (any(abs(donor_sums$total - 1) > 0.01)) {
  warning("Some donors' proportions do not sum to 1 (tolerance 0.01).")
}

# ---- Prepare plot data ------------------------------------------------------

plot_data <- proportions |>
  mutate(
    # healthy is the reference level, shown first
    condition = factor(condition, levels = conditions),
    # panels ordered from most to least abundant cell type
    cell_type = fct_reorder(cell_type, proportion, .fun = mean, .desc = TRUE),
    cell_type_label = fct_relabel(cell_type, \(x) gsub("_", " ", x))
  ) |>
  arrange(cell_type, condition, donor_id)

condition_medians <- plot_data |>
  group_by(cell_type_label, condition) |>
  summarise(median_proportion = median(proportion), .groups = "drop")

# ---- Plot -------------------------------------------------------------------

p <- ggplot(plot_data, aes(x = condition, y = proportion, colour = condition)) +
  geom_crossbar(
    data = condition_medians,
    aes(y = median_proportion, ymin = median_proportion, ymax = median_proportion),
    width = 0.5,
    linewidth = 0.4,
    colour = "grey40"
  ) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 1),
    size = 2.5,
    alpha = 0.9
  ) +
  facet_wrap(~cell_type_label, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = condition_colours, name = "Condition") +
  scale_y_continuous(
    labels = \(x) paste0(round(x * 100), "%"),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = NULL,
    y = "Proportion of donor's cells",
    caption = "One point per donor (n = 4 per condition); horizontal bar = median."
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    plot.caption = element_text(colour = "grey40", hjust = 0)
  )

# ---- Save -------------------------------------------------------------------

fig_width <- 9
fig_height <- 3.6

ggsave(
  file.path(run_dir, "proportions_by-condition.pdf"),
  p,
  width = fig_width,
  height = fig_height,
  device = cairo_pdf
)
ggsave(
  file.path(run_dir, "proportions_by-condition.png"),
  p,
  width = fig_width,
  height = fig_height,
  dpi = 300
)

plot_data |>
  transmute(
    donor_id,
    condition = as.character(condition),
    cell_type = as.character(cell_type),
    proportion,
    colour = unname(condition_colours[condition])
  ) |>
  write_csv(file.path(run_dir, "proportions_by-condition.csv"))

message("Wrote figure and data to ", normalizePath(run_dir))
