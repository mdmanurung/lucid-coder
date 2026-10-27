# Purpose: fit the planned disease-vs-healthy contrast per gene on donor-level expression.
# Inputs:  data/donor_expression.tsv   [gap: toy data built inline until the producer exists]
# Outputs: <run_dir>/contrasts.tsv, <run_dir>/top-gene_by-condition.{pdf,png,csv},
#          <run_dir>/session_info.txt
# Run:     Rscript scripts/01_fit-contrasts.R results/2026-10-26_contrasts
# Notes:   [Stochastic] none. Inference unit: donor. One row per donor x gene.

# Setup ----
library(dplyr)
library(forcats)
library(ggplot2)
library(purrr)
library(tidyr)
# align copies palettes.R and theme.R to the project root.
source(here::here("palettes.R"))
source(here::here("theme.R"))

# Reference = healthy, so log2_fc reads as disease minus healthy.
condition_levels <- c("healthy", "disease")
# Multiple-testing family: every gene tested in this script, adjusted together.
p_adjust_method <- "BH"

args <- commandArgs(trailingOnly = TRUE)
stopifnot("usage: Rscript scripts/01_fit-contrasts.R <run_dir>" = length(args) == 1)
run_dir <- args[[1]]
stopifnot("run_dir already exists; pick a new name" = !dir.exists(run_dir))
dir.create(run_dir, recursive = TRUE)

# Load inputs ----
# Toy data stands in for data/donor_expression.tsv; replace with the read call.
# Per gene: four healthy donors, then four disease donors.
expression <- expand_grid(gene = c("IFI27", "ACTB", "CD3E"), donor_id = sprintf("D%02d", 1:8)) |>
  relocate(donor_id) |>
  mutate(
    condition = fct(rep(rep(condition_levels, each = 4), times = 3), levels = condition_levels),
    log2_expr = c(
      c(5.1, 4.8, 5.4, 5.0, 8.2, 7.9, 8.5, 8.1), # IFI27
      c(10.1, 9.9, 10.2, 10.0, 10.0, 10.3, 9.8, 10.1), # ACTB
      c(6.2, 6.0, 6.4, 6.1, 5.3, 5.6, 5.1, 5.4) # CD3E
    )
  )

# Fit contrasts ----
de_results <- expression |>
  split(~gene) |>
  map(\(gene_expression) {
    fit <- lm(log2_expr ~ condition, data = gene_expression)
    disease_term <- coef(summary(fit))["conditiondisease", ]
    tibble(log2_fc = disease_term[["Estimate"]], p_value = disease_term[["Pr(>|t|)"]])
  }) |>
  list_rbind(names_to = "gene") |>
  mutate(q_value = p.adjust(p_value, method = p_adjust_method)) |>
  arrange(q_value)

# Plot the top gene ----
# Donor points, because the inference unit is the donor.
top_gene <- de_results$gene[[1]]
top_gene_expression <- filter(expression, gene == top_gene)
p_top_gene <- ggplot(top_gene_expression, aes(condition, log2_expr, colour = condition)) +
  geom_point() +
  scale_colour_manual(values = palette_for(condition_levels, pal_condition)) +
  labs(title = top_gene, y = "log2 expression (one point per donor)") +
  theme(legend.position = "none")

# Write outputs ----
write.table(
  de_results,
  file.path(run_dir, "contrasts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
save_panel(p_top_gene, top_gene_expression, run_dir, "top-gene_by-condition", 45, 45)

# Session info ----
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
