# Plot cell-type proportions by condition, one point per donor, in lab colours.
# Usage: Rscript scripts/03_plot-proportions.R <run_dir>

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(magrittr)
  library(here)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript scripts/03_plot-proportions.R <run_dir>")
run_dir <- args[[1]]
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

props <- read_tsv(here("data", "proportions.tsv"), show_col_types = FALSE)
palette <- read_tsv(here("data", "palette.tsv"), show_col_types = FALSE)
pal <- set_names(palette$hex, palette$label)

conditions <- c("healthy", "disease")
missing <- setdiff(c(conditions, unique(props$cell_type)), names(pal))
if (length(missing) > 0) stop("No palette colour for: ", paste(missing, collapse = ", "))
if (!all(props$condition %in% conditions)) stop("Unexpected condition values")

# Order cell types by overall mean proportion, most abundant first.
cell_order <- props %>%
  group_by(cell_type) %>%
  summarise(m = mean(proportion), .groups = "drop") %>%
  arrange(desc(m)) %>%
  pull(cell_type)

plot_data <- props %>%
  mutate(
    condition = factor(condition, levels = conditions),
    cell_type = factor(cell_type, levels = cell_order)
  ) %>%
  arrange(cell_type, condition, donor_id)

summary_data <- plot_data %>%
  group_by(cell_type, condition) %>%
  summarise(mean_proportion = mean(proportion), .groups = "drop")

write_csv(plot_data, file.path(run_dir, "proportions_by-condition.csv"))

p <- ggplot(plot_data, aes(condition, proportion, colour = condition)) +
  geom_crossbar(
    data = summary_data,
    aes(y = mean_proportion, ymin = mean_proportion, ymax = mean_proportion),
    width = 0.5, linewidth = 0.4
  ) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), size = 2.2) +
  facet_wrap(~cell_type, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = pal[conditions], breaks = conditions) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, NA)) +
  labs(x = NULL, y = "Proportion of donor's cells", colour = "Condition",
       caption = "Points: donors (n = 4 per condition); bar: mean") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", colour = "black"),
    legend.position = "bottom"
  )

# Colour each facet strip with that cell type's lab colour.
g <- ggplotGrob(p)
strips <- grep("^strip", g$layout$name)
strip_order <- order(g$layout$l[strips])
for (i in seq_along(strips)) {
  idx <- strips[strip_order[i]]
  ct <- cell_order[i]
  rect <- which(vapply(g$grobs[[idx]]$grobs[[1]]$children, function(x) inherits(x, "rect"), logical(1)))
  g$grobs[[idx]]$grobs[[1]]$children[[rect[1]]]$gp$fill <- pal[[ct]]
}

w <- 10; h <- 4
ggsave(file.path(run_dir, "proportions_by-condition.pdf"), g, width = w, height = h)
ggsave(file.path(run_dir, "proportions_by-condition.png"), g, width = w, height = h, dpi = 300)
