# Purpose: test each TF's chromVAR score for a trend over time since symptom onset, with the
#          person (pid) as the unit of replication.
# Run:     Rscript scripts/07c_model-tf-time.R results/run1
# Inputs:  data/samples.tsv    [gap: one row per sample; pid, centred time, age, sex]
#          data/tf_scores.tsv  [gap: chromVAR deviation z-scores, samples x TFs]
# Outputs: <run_dir>/tf_time.tsv          one row per TF: time slope, p, BH q
#          <run_dir>/session_info.txt     R session
# Notes:   [Stochastic] none. The 69 samples come from 18 pids with 3 to 5 samples each, and age
#          and sex are constant within a pid, so samples are repeated visits of one person and
#          not independent: a random intercept per pid. Fit by ML, as the MOCHA manuscript did;
#          Satterthwaite df from lmerTest.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(lmerTest)
})

# Parameters --------------------------------------------------------------------

# Age and sex are between-person; time and its square are within-person.
MODEL <- score ~ age + sex + time + time_sqrd + (1 | pid)
# Multiple-testing family: the time slope across every TF in tf_scores.tsv.
P_ADJUST_METHOD <- "BH"

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
if (dir.exists(run_dir)) {
  stop("run_dir ", run_dir, " exists; pick a new name.", call. = FALSE)
}
dir.create(run_dir, recursive = TRUE)

# Load data ---------------------------------------------------------------------

samples <- read_tsv(
  here::here("data", "samples.tsv"),
  col_types = cols(pid = "c")
)
tf_scores <- read_tsv(
  here::here("data", "tf_scores.tsv"),
  show_col_types = FALSE
)
stopifnot(
  "samples match tf_scores" = identical(samples$sample_id, tf_scores$sample_id)
)

# Fit one mixed model per TF ----------------------------------------------------

tf_long <- tf_scores |>
  pivot_longer(-sample_id, names_to = "tf", values_to = "score") |>
  left_join(samples, join_by(sample_id), relationship = "many-to-one")
tf_time <- tf_long |>
  split(~tf) |>
  map(\(df) {
    time_row <- coef(summary(lmer(MODEL, data = df, REML = FALSE)))["time", ]
    tibble(estimate = time_row[["Estimate"]], p_value = time_row[["Pr(>|t|)"]])
  }) |>
  list_rbind(names_to = "tf") |>
  mutate(q_value = p.adjust(p_value, method = P_ADJUST_METHOD)) |>
  arrange(p_value)

# Write outputs -----------------------------------------------------------------

write_tsv(tf_time, file.path(run_dir, "tf_time.tsv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "session_info.txt")
)
