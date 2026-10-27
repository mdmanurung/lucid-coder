#$ ===== Reviewer: type after each colon, then add #$ comments anywhere below =====
#$ q1 one read (1-5; can you say what this script does after one read?): 4
#$ q2 find decisions (1-5; can you find each threshold, filter, model choice and path?): 4
#$ q3 change checks (yes/no; would you add, remove or change a check?): no
#$ q3 which (one sentence): 
#$ ---- your #$ comments: on their own line above the code, or after the code on its line ----
# Purpose: Plot cell-type proportion by condition, one point per donor.
# Inputs: data/proportions.tsv [gap], data/palette.tsv [lab palette]
# Outputs: <run_dir>/proportions_by-condition.{pdf,png,csv}
# Run: Rscript scripts/03_plot-proportions.R results/run1

library(dplyr)
library(readr)
library(forcats)
library(ggplot2)

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
condition_levels <- c("healthy", "disease")
jitter_seed <- 1

theme_set(theme_bw(base_size = 8))

out_stem <- file.path(run_dir, "proportions_by-condition")
stopifnot(
  "run dir already holds this figure" = !file.exists(paste0(out_stem, ".pdf"))
)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

# Load data ----
proportions <- read_tsv(here::here("data/proportions.tsv"), show_col_types = FALSE) |>
  mutate(
    condition = fct(condition, levels = condition_levels),
    # Order panels by overall abundance so the largest populations read first
    cell_type = fct_reorder(cell_type, proportion, .fun = median, .desc = TRUE)
  )

palette <- read_tsv(here::here("data/palette.tsv"), show_col_types = FALSE)
condition_colours <- setNames(palette$hex, palette$label)[condition_levels]
stopifnot("every condition has a palette colour" = !anyNA(condition_colours))

donors_per_condition <- proportions |>
  distinct(donor_id, condition) |>
  count(condition)

# Plot ----
p <- proportions |>
  ggplot(aes(x = condition, y = proportion, colour = condition)) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.5,
    linewidth = 0.3,
    colour = "grey40"
  ) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = jitter_seed),
    size = 1.5
  ) +
  facet_wrap(vars(cell_type), nrow = 1) +
  # Breaks = condition_levels keeps legend order equal to x-axis order
  scale_colour_manual(values = condition_colours, breaks = condition_levels) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    x = NULL,
    y = "Proportion of donor's cells (fraction)",
    colour = "Condition",
    caption = sprintf(
      "Points: donors (%s). Bar: median per condition.",
      paste0(donors_per_condition$condition, " n = ", donors_per_condition$n, collapse = ", ")
    )
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save ----
ggsave(
  paste0(out_stem, ".pdf"),
  plot = p,
  width = 183,
  height = 70,
  units = "mm",
  device = cairo_pdf
)
ggsave(paste0(out_stem, ".png"), plot = p, width = 183, height = 70, units = "mm", dpi = 300)
write_csv(p$data, paste0(out_stem, ".csv"))

writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
