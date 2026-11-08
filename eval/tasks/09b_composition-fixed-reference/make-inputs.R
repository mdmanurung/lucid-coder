# Purpose: freeze the task 09b inputs and truth: one simulated draw of per-mouse cell counts in
#          which one type expands 10-fold and a stable reference type is known.
# Run:     module load statistical/R/4.5.2; Rscript eval/tasks/09b_composition-fixed-reference/make-inputs.R
# Inputs:  none (simulated; the simulator is sim_counts from the 09b design simulation)
# Outputs: eval/tasks/09b_composition-fixed-reference/inputs/data/cells.tsv  one row per cell
#          eval/tasks/09b_composition-fixed-reference/truth/cell_types.tsv    role per cell type
# Notes:   [Stochastic] counts are negative binomial, one draw. The seed is the smallest one for
#          which (i) edgeR QL with offset log(Goblet count) and the t-test on log(count / Goblet
#          count) call no fixed type down at q < 0.05 and call Enterocyte up at q < 0.05, and
#          (ii) the proportion t-test, the quasi-binomial GLM and crumblr + limma each call at
#          least one fixed type down at q < 0.05 (BH over the types each method tests).
#          Search result: seed 1 (the first seed tried). Fixed types down at q < 0.05, and
#          Enterocyte up at q < 0.05, per method on this draw:
#            prop t-test        2 (Stem q 7.2e-05, Endocrine q 0.015)   Enterocyte up
#            quasi-binomial     3 (Stem, TA, Endocrine; q 3.6e-05-0.037) Enterocyte up
#            crumblr + limma    1 (Stem q 0.032)                         Enterocyte up
#            edgeR QL offset    0 (smallest fixed q 0.35)                Enterocyte q 4.5e-05
#            log-ratio t-test   0 (smallest fixed q 0.48)                Enterocyte q 0.026

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(edgeR)
  library(limma)
  library(crumblr)
})

# Parameters --------------------------------------------------------------------

# Mean cells per mouse in Control; Goblet is the stable reference, listed first.
MEAN_CONTROL <- c(
  Goblet = 300,
  Enterocyte = 100,
  Stem = 150,
  TA = 150,
  Endocrine = 150,
  Tuft = 150
)
# Enterocyte expands 10-fold under Salmonella; every other type keeps its mean.
FOLD_UP <- 10
# Negative binomial dispersion phi; size = 1 / phi.
PHI <- 0.1
# Mice per group.
N_MICE <- 4
# Significance level for the seed search, as in the gate.
Q_CUTOFF <- 0.05
# Seeds searched, smallest first.
SEEDS <- 1:499

REFERENCE <- names(MEAN_CONTROL)[[1]]
UP <- "Enterocyte"
FIXED <- setdiff(names(MEAN_CONTROL), c(REFERENCE, UP))
task_dir <- here::here("eval", "tasks", "09b_composition-fixed-reference")

# Simulate one draw -------------------------------------------------------------

# sim_counts from the design simulation: mice x types, group B is Salmonella.
sim_counts <- function(seed) {
  set.seed(seed)
  condition <- rep(c("Control", "Salmonella"), each = N_MICE)
  mu <- matrix(MEAN_CONTROL, 2 * N_MICE, length(MEAN_CONTROL), byrow = TRUE)
  mu[condition == "Salmonella", names(MEAN_CONTROL) == UP] <- MEAN_CONTROL[[
    UP
  ]] *
    FOLD_UP
  counts <- matrix(rnbinom(length(mu), mu = mu, size = 1 / PHI), nrow(mu))
  dimnames(counts) <- list(paste0("B", seq_len(nrow(mu))), names(MEAN_CONTROL))
  list(
    counts = counts,
    condition = factor(condition, levels = c("Control", "Salmonella"))
  )
}

# Call each type with five methods ----------------------------------------------

# One row per method and type: estimate (Salmonella vs Control) and BH q over that method's types.
call_types <- function(draw) {
  y <- draw$counts
  condition <- draw$condition
  design <- model.matrix(~condition)
  share <- y / rowSums(y)
  log_ratio <- log(y[, colnames(y) != REFERENCE] / y[, REFERENCE])
  t_by_type <- \(m) {
    map(colnames(m), \(k) {
      t <- t.test(m[, k] ~ condition)
      tibble(cell_type = k, est = diff(t$estimate), p = t$p.value)
    }) |>
      list_rbind()
  }
  # Proportion t-test, closure-unaware.
  prop_t <- t_by_type(share)
  # Quasi-binomial GLM on per-mouse counts, closure-unaware.
  quasi <- map(colnames(y), \(k) {
    fit <- glm(
      cbind(y[, k], rowSums(y) - y[, k]) ~ condition,
      family = quasibinomial
    )
    coefs <- coef(summary(fit))["conditionSalmonella", ]
    tibble(cell_type = k, est = coefs[["Estimate"]], p = coefs[["Pr(>|t|)"]])
  }) |>
    list_rbind()
  # crumblr + limma with no reference, closure-unaware.
  crumblr_fit <- eBayes(lmFit(crumblr(y), design))
  crumblr_tab <- topTable(
    crumblr_fit,
    coef = 2,
    number = Inf,
    sort.by = "none"
  ) |>
    as_tibble(rownames = "cell_type") |>
    select(cell_type, est = logFC, p = P.Value)
  # edgeR QL with offset log(reference count), no library-size normalisation.
  dge <- DGEList(t(y[, colnames(y) != REFERENCE]))
  dge$offset <- matrix(log(y[, REFERENCE]), nrow(dge), ncol(dge), byrow = TRUE)
  dge <- estimateDisp(dge, design)
  edger_tab <- topTags(
    glmQLFTest(glmQLFit(dge, design), coef = 2),
    n = Inf,
    sort.by = "none"
  )$table |>
    as_tibble(rownames = "cell_type") |>
    select(cell_type, est = logFC, p = PValue)
  # t-test on log(count / reference count).
  ratio_t <- t_by_type(log_ratio)
  list(
    prop_t = prop_t,
    quasi_binomial = quasi,
    crumblr = crumblr_tab,
    edger_offset = edger_tab,
    log_ratio_t = ratio_t
  ) |>
    list_rbind(names_to = "method") |>
    mutate(q = p.adjust(p, "BH"), .by = method)
}

# Search for the seed -----------------------------------------------------------

# Per method: fixed types called down, and whether Enterocyte is called up.
summarise_calls <- function(calls) {
  calls |>
    summarise(
      fixed_down = sum(cell_type %in% FIXED & est < 0 & q < Q_CUTOFF),
      up_called = any(cell_type == UP & est > 0 & q < Q_CUTOFF),
      .by = method
    )
}
seed <- detect(SEEDS, \(s) {
  calls <- summarise_calls(call_types(sim_counts(s)))
  aware <- filter(calls, method %in% c("edger_offset", "log_ratio_t"))
  unaware <- filter(calls, method %in% c("prop_t", "quasi_binomial", "crumblr"))
  all(aware$fixed_down == 0 & aware$up_called) && all(unaware$fixed_down > 0)
})
if (is.null(seed)) {
  stop("no seed below 500 meets both conditions", call. = FALSE)
}
draw <- sim_counts(seed)
message("seed ", seed)
print(as.data.frame(summarise_calls(call_types(draw))))
print(as.data.frame(call_types(draw)))
print(draw$counts)

# Expand counts to cells --------------------------------------------------------

mice <- tibble(
  sample_id = rownames(draw$counts),
  condition = as.character(draw$condition)
)
cells <- as_tibble(draw$counts, rownames = "sample_id") |>
  pivot_longer(-sample_id, names_to = "cell_type", values_to = "n") |>
  uncount(n) |>
  left_join(
    mice,
    join_by(sample_id),
    relationship = "many-to-one",
    unmatched = "error"
  ) |>
  mutate(
    cell_id = sprintf("%s_%05d", sample_id, row_number()),
    .by = sample_id
  ) |>
  select(cell_id, sample_id, condition, cell_type)
roles <- tibble(cell_type = names(MEAN_CONTROL)) |>
  mutate(
    role = case_when(
      cell_type == REFERENCE ~ "reference",
      cell_type == UP ~ "up",
      .default = "unchanged"
    )
  )

# Write inputs and truth --------------------------------------------------------

dir.create(
  file.path(task_dir, "inputs", "data"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(file.path(task_dir, "truth"), showWarnings = FALSE)
write_tsv(cells, file.path(task_dir, "inputs", "data", "cells.tsv"))
write_tsv(roles, file.path(task_dir, "truth", "cell_types.tsv"))
