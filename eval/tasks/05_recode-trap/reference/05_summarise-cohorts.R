# Purpose: label each donor's age group from its cohort and summarise the score per group.
# Inputs:  data/samples.tsv   [gap: recruitment sample sheet]
#          data/scores.tsv    [gap: plasma score per donor]
# Outputs: <run_dir>/samples_labelled.tsv, <run_dir>/cohort_summary.tsv, <run_dir>/session_info.txt
# Run:     Rscript scripts/05_summarise-cohorts.R results/run1
# Notes:   [Stochastic] none. One cohort maps to one age group; any other cohort label stops
#          the script instead of joining a group.

# ---- Setup ----
library(dplyr)

# Cohort to age group. Drives both the recode and the factor levels.
AGE_GROUP <- c(BR1 = "young", BR2 = "older")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: Rscript scripts/05_summarise-cohorts.R <run_dir>", call. = FALSE)
}
run_dir <- args[[1]]
if (dir.exists(run_dir)) {
  stop(
    "run_dir ",
    run_dir,
    " already exists; pick a new name so nothing is overwritten.",
    call. = FALSE
  )
}
dir.create(run_dir, recursive = TRUE)

# ---- Load inputs ----
samples <- read.delim(here::here("data", "samples.tsv"))
scores <- read.delim(here::here("data", "scores.tsv"))
if (
  !all(c("sample_id", "cohort") %in% names(samples)) ||
    !all(c("sample_id", "score") %in% names(scores)) ||
    anyDuplicated(samples$sample_id)
) {
  stop(
    "samples.tsv/scores.tsv: missing columns or duplicated sample_id; fix the sheet.",
    call. = FALSE
  )
}

# ---- Label age groups ----
samples$cohort <- trimws(samples$cohort)
samples_labelled <- samples |>
  mutate(age_group = factor(AGE_GROUP[cohort], levels = AGE_GROUP)) |>
  select(sample_id, cohort, age_group)
unmapped <- samples_labelled |> filter(is.na(age_group))
if (nrow(unmapped) > 0) {
  stop(
    "samples$cohort: unmapped labels ",
    toString(sprintf("%s (%s)", unmapped$cohort, unmapped$sample_id)),
    "; fix the sheet or extend AGE_GROUP.",
    call. = FALSE
  )
}

# ---- Summarise per age group ----
cohort_summary <- samples_labelled |>
  inner_join(scores, by = "sample_id", relationship = "one-to-one") |>
  group_by(age_group) |>
  summarise(n = n(), mean_score = mean(score), .groups = "drop")

# ---- Write outputs ----
readr::write_tsv(samples_labelled, file.path(run_dir, "samples_labelled.tsv"))
readr::write_tsv(cohort_summary, file.path(run_dir, "cohort_summary.tsv"))

# ---- Session info ----
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
