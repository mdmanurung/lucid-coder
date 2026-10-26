# Plot cell-type proportions by condition, one point per donor, in lab colours.
#
# Usage: Rscript scripts/03_plot-proportions.R <run_dir>
# Writes proportions_by-condition.{pdf,png,csv} into <run_dir>.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(magrittr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript scripts/03_plot-proportions.R <run_dir>", call. = FALSE)
}
run_dir <- args[[1]]
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

proportions <- read_tsv(here::here("data", "proportions.tsv"), show_col_types = FALSE)
palette <- read_tsv(here::here("data", "palette.tsv"), show_col_types = FALSE)

# The palette file is not in a fixed order, so look colours up by label.
conditions <- c("healthy", "disease")
missing <- setdiff(c(conditions, unique(proportions$cell_type)), palette$label)
if (length(missing) > 0) {
  stop("No palette colour for: ", paste(missing, collapse = ", "), call. = FALSE)
}
unknown <- setdiff(unique(proportions$condition), conditions)
if (length(unknown) > 0) {
  stop("Unexpected condition(s): ", paste(unknown, collapse = ", "), call. = FALSE)
}
colours <- setNames(palette$hex, palette$label)

# Order cell types by overall mean proportion, most abundant first.
cell_type_order <- proportions %>%
  group_by(cell_type) %>%
  summarise(mean_proportion = mean(proportion), .groups = "drop") %>%
  arrange(desc(mean_proportion)) %>%
  pull(cell_type)

plot_data <- proportions %>%
  mutate(
    condition = factor(condition, levels = conditions),
    cell_type = factor(cell_type, levels = cell_type_order),
    colour = unname(colours[as.character(condition)])
  ) %>%
  group_by(cell_type, condition) %>%
  mutate(median_proportion = median(proportion)) %>%
  ungroup() %>%
  arrange(cell_type, condition, donor_id)

write_csv(plot_data, file.path(run_dir, "proportions_by-condition.csv"))

medians <- distinct(plot_data, cell_type, condition, median_proportion)

p <- ggplot(plot_data, aes(x = condition, y = proportion, colour = condition)) +
  geom_crossbar(
    data = medians,
    aes(y = median_proportion, ymin = median_proportion, ymax = median_proportion),
    width = 0.5, linewidth = 0.4, colour = "grey40"
  ) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 1),
    size = 2.5
  ) +
  facet_wrap(~cell_type, nrow = 1, scales = "free_y",
             labeller = as_labeller(function(x) gsub("_", " ", x))) +
  scale_colour_manual(values = colours[conditions], breaks = conditions) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = NULL,
    y = "Proportion of donor's cells",
    colour = "Condition",
    caption = "Points: donors (n = 4 per condition); bars: median."
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  )

ggsave(file.path(run_dir, "proportions_by-condition.pdf"), p, width = 9, height = 3.5)
ggsave(file.path(run_dir, "proportions_by-condition.png"), p, width = 9, height = 3.5, dpi = 300)
