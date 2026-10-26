# Naive answer: makes the figure but forgets the plotted-data CSV.
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

df <- read.delim("data/proportions.tsv")
pal <- read.delim("data/palette.tsv")
cols <- setNames(pal$hex, pal$label)

p <- ggplot(df, aes(cell_type, proportion, colour = condition)) +
  geom_point(position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = cols) +
  theme_bw()
ggsave(file.path(out, "proportions_by-condition.pdf"), p, width = 6, height = 4)
ggsave(file.path(out, "proportions_by-condition.png"), p, width = 6, height = 4, dpi = 150)
