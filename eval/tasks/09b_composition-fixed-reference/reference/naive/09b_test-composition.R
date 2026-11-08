# Naive answer: task 09's reference, crumblr + limma with no reference type, so the Enterocyte
# expansion pulls every other type's centred log-ratio down (closure).
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(limma)
  library(crumblr)
})

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
dir.create(run_dir, recursive = TRUE)

cells <- read_tsv(here::here("data", "cells.tsv"), show_col_types = FALSE)
counts <- unclass(table(cells$sample_id, cells$cell_type))
samples <- cells |>
  distinct(sample_id, condition) |>
  mutate(
    condition = forcats::fct(condition, levels = c("Control", "Salmonella"))
  ) |>
  arrange(match(sample_id, rownames(counts)))
stopifnot(
  "samples match counts" = identical(samples$sample_id, rownames(counts))
)

fit <- eBayes(lmFit(crumblr(counts), model.matrix(~condition, samples)))
composition <- topTable(fit, coef = "conditionSalmonella", number = Inf) |>
  as_tibble(rownames = "cell_type") |>
  select(cell_type, log_fc = logFC, p_value = P.Value, q_value = adj.P.Val)
write_tsv(composition, file.path(run_dir, "composition.tsv"))
