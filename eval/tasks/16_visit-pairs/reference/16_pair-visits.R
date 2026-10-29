# Purpose: pair each donor's baseline and day 7 IL-6 from the latest runs and test the change.
# Run:     Rscript scripts/16_pair-visits.R results/run1
# Inputs:  data/runs.tsv  [gap] one row per assay run: run_id, donor_id, visit, run_date, il6
# Outputs: <run_dir>/pairs.tsv         per donor with both visits: il6_d0, il6_d7, change
#          <run_dir>/paired_test.tsv   paired t-test of change: n_donors, mean_change, p_value
#          <run_dir>/session_info.txt  R session

# Setup ------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})
run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
stopifnot("run dir is new" = !dir.exists(run_dir))
dir.create(run_dir, recursive = TRUE)

# Load data --------------------------------------------------------------------
runs <- readr::read_tsv(here::here("data", "runs.tsv"), show_col_types = FALSE)

# Keep the latest run per sample -----------------------------------------------
# A rerun replaces earlier runs of the same donor and visit; a tie on run_date stops below.
latest_runs <- runs |>
  filter(visit %in% c("D0", "D7")) |>
  slice_max(run_date, n = 1, by = c(donor_id, visit))
stopifnot("one run per donor and visit" = !anyDuplicated(latest_runs[c("donor_id", "visit")]))

# Pair baseline and day 7 ------------------------------------------------------
visit_wide <- latest_runs |>
  select(donor_id, visit, il6) |>
  pivot_wider(names_from = visit, values_from = il6, names_glue = "il6_{tolower(visit)}")
unpaired <- visit_wide |>
  filter(is.na(il6_d0) | is.na(il6_d7))
if (nrow(unpaired) > 0) {
  message("Dropped, no D0 or D7 run: ", toString(unpaired$donor_id))
}
pairs <- visit_wide |>
  drop_na(il6_d0, il6_d7) |>
  mutate(change = il6_d7 - il6_d0) |>
  arrange(donor_id)

# Test the change --------------------------------------------------------------
change_test <- t.test(pairs$change)
paired_test <- tibble(
  n_donors = nrow(pairs),
  mean_change = unname(change_test$estimate),
  p_value = change_test$p.value
)

# Write outputs ----------------------------------------------------------------
readr::write_tsv(pairs, file.path(run_dir, "pairs.tsv"))
readr::write_tsv(paired_test, file.path(run_dir, "paired_test.tsv"))
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
