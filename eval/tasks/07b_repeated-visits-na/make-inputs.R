# Purpose: freeze the task 07b inputs and truth: task 07's data with planted missing TF scores.
# Run:     Rscript eval/tasks/07b_repeated-visits-na/make-inputs.R
# Inputs:  eval/tasks/07_repeated-visits/inputs/data/samples.tsv    [task 07 make-inputs.R]
#          eval/tasks/07_repeated-visits/inputs/data/tf_scores.tsv  [same, 69 samples x 30 TFs]
#          eval/tasks/07_repeated-visits/truth/tf_time.tsv          [same, trap/filler role]
# Outputs: eval/tasks/07b_repeated-visits-na/inputs/data/samples.tsv    copied unchanged
#          eval/tasks/07b_repeated-visits-na/inputs/data/tf_scores.tsv  with NA scores planted
#          eval/tasks/07b_repeated-visits-na/truth/tf_time.tsv          per-TF lmer on its own rows
# Notes:   One filler TF misses 17 samples (a motif absent from some chromVAR runs); two more
#          fillers miss 1 and 2 other samples, so 20 samples are incomplete. The truth is task
#          07's model, exp ~ age + sex + time + time_sqrd + (1 | donor), REML = FALSE, lmerTest
#          Satterthwaite p, fitted on each TF's available rows. On the full data this refit
#          matches the MOCHA repo's p-values within 3e-8 log10 units. Data are GPL-3.0 (NOTICE.md).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(lmerTest)
})

# Parameters --------------------------------------------------------------------

# Missing samples per TF: one heavy gap, two stray ones; the three sample sets are disjoint.
N_MISSING <- c(HOXB2 = 17, TET1 = 1, RUNX1 = 2)

src_dir <- here::here("eval", "tasks", "07_repeated-visits")
task_dir <- here::here("eval", "tasks", "07b_repeated-visits-na")

# Load task 07 ------------------------------------------------------------------

samples <- read_tsv(
  file.path(src_dir, "inputs", "data", "samples.tsv"),
  col_types = cols(.default = "c")
)
tf_text <- read_tsv(
  file.path(src_dir, "inputs", "data", "tf_scores.tsv"),
  col_types = cols(.default = "c")
)
role <- read_tsv(
  file.path(src_dir, "truth", "tf_time.tsv"),
  col_types = cols(.default = "c")
) |>
  select(tf, role)
stopifnot("samples match tf_scores" = identical(samples$sample_id, tf_text$sample_id))

# Plant missing scores ----------------------------------------------------------

set.seed(71)
missing_rows <- sample(nrow(tf_text), sum(N_MISSING)) |>
  split(rep(names(N_MISSING), N_MISSING))
for (tf in names(N_MISSING)) {
  tf_text[[tf]][missing_rows[[tf]]] <- NA
}

# Fit the truth -----------------------------------------------------------------

tf_long <- type.convert(tf_text, as.is = TRUE) |>
  pivot_longer(-sample_id, names_to = "tf", values_to = "score") |>
  filter(!is.na(score)) |>
  left_join(
    type.convert(samples, as.is = TRUE) |>
      mutate(donor_id = as.character(donor_id)),
    join_by(sample_id),
    relationship = "many-to-one"
  )
truth <- tf_long |>
  split(~tf) |>
  map(\(df) {
    fit <- lmer(score ~ age + sex + time + time_sqrd + (1 | donor_id), df, REML = FALSE)
    time_row <- coef(summary(fit))["time", ]
    tibble(
      estimate = time_row[["Estimate"]],
      p_value = time_row[["Pr(>|t|)"]],
      n_samples = nobs(fit)
    )
  }) |>
  list_rbind(names_to = "tf") |>
  left_join(role, join_by(tf), relationship = "one-to-one") |>
  relocate(role, .after = tf)

# Write inputs and truth --------------------------------------------------------

dir.create(file.path(task_dir, "inputs", "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(task_dir, "truth"), showWarnings = FALSE)
write_tsv(samples, file.path(task_dir, "inputs", "data", "samples.tsv"))
write_tsv(tf_text, file.path(task_dir, "inputs", "data", "tf_scores.tsv"), na = "NA")
write_tsv(truth, file.path(task_dir, "truth", "tf_time.tsv"))
