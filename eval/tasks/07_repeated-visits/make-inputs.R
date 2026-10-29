# Purpose: freeze the task 07 inputs and truth from the MOCHA manuscript repository, once.
# Run:     Rscript eval/tasks/07_repeated-visits/make-inputs.R <SuppFig17_ZeroInflatedGLMMs dir>
# Inputs:  <dir>/metadata.csv                [aifimmunology/MOCHA_Manuscript, 69 samples, 18 donors]
#          <dir>/tf_matrix.csv               [same, 69 samples x 870 chromVAR TF z-scores]
#          <dir>/tf_usage_model_results.csv  [same, the repo's per-TF lmer fit, the answer]
# Outputs: eval/tasks/07_repeated-visits/inputs/data/samples.tsv    one row per sample (visit)
#          eval/tasks/07_repeated-visits/inputs/data/tf_scores.tsv  samples x 30 TFs
#          eval/tasks/07_repeated-visits/truth/tf_time.tsv          the repo's time slope and p
# Notes:   Source: https://github.com/aifimmunology/MOCHA_Manuscript at
#          48c75cd78075bc59a022d68f4090b0e831cd6a11, GPL-3.0; these files are a subset of it.
#          The repo fits exp ~ Age + sex + time + time_sqrd + (1 | PTID), REML = FALSE, per TF
#          (Fig5/Fig5_chromVarModeling.R). Trap TFs: the 8 distinct TFs whose time p-value moves
#          most when the donor term is dropped (lm), among those with repo p < 1e-3; 22 seeded
#          fillers. Values are copied as text, never re-rounded.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

N_TRAP <- 8
N_FILLER <- 22

src_dir <- commandArgs(trailingOnly = TRUE)[[1]]
task_dir <- here::here("eval", "tasks", "07_repeated-visits")

# Load source -------------------------------------------------------------------

meta <- read_csv(
  file.path(src_dir, "metadata.csv"),
  col_types = cols(.default = "c")
)
tf_text <- read_csv(
  file.path(src_dir, "tf_matrix.csv"),
  col_types = cols(.default = "c")
)
repo <- read_csv(
  file.path(src_dir, "tf_usage_model_results.csv"),
  show_col_types = FALSE
)
stopifnot("tf rows match metadata" = identical(tf_text[[1]], meta$Sample))
tf <- as.matrix(mutate(tf_text[, -1], across(everything(), as.numeric)))

# Pick TFs ----------------------------------------------------------------------

# Clones (identical motifs) and Ensembl-only names are dropped so each TF is distinct.
keep <- !duplicated(t(round(tf, 8))) & !grepl("^ENSG", colnames(tf))
meta_num <- type.convert(meta, as.is = TRUE)
lm_p <- apply(tf, 2, \(y) {
  coef(summary(lm(y ~ Age + Sex + time + time_sqrd, meta_num)))["time", 4]
})
shift <- tibble(tf = colnames(tf), keep, lm_p) |>
  left_join(
    select(repo, tf = TF, repo_p = time_pval),
    join_by(tf),
    relationship = "one-to-one"
  ) |>
  filter(keep) |>
  mutate(shift = log10(lm_p) - log10(repo_p))
trap_tfs <- shift |>
  filter(repo_p < 1e-3) |>
  slice_max(shift, n = N_TRAP) |>
  pull(tf)
set.seed(7)
filler_tfs <- sample(setdiff(shift$tf, trap_tfs), N_FILLER)
picked <- c(trap_tfs, filler_tfs)

# Write inputs and truth --------------------------------------------------------

samples <- meta |>
  select(
    sample_id = Sample,
    donor_id = PTID,
    time,
    time_sqrd,
    age = Age,
    sex = Sex
  )
tf_scores <- tf_text[, c(names(tf_text)[[1]], picked)] |>
  rename(sample_id = 1)
truth <- repo |>
  filter(TF %in% picked) |>
  transmute(
    tf = TF,
    role = if_else(tf %in% trap_tfs, "trap", "filler"),
    estimate = time,
    p_value = time_pval
  )

dir.create(
  file.path(task_dir, "inputs", "data"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(file.path(task_dir, "truth"), showWarnings = FALSE)
write_tsv(samples, file.path(task_dir, "inputs", "data", "samples.tsv"))
write_tsv(tf_scores, file.path(task_dir, "inputs", "data", "tf_scores.tsv"))
write_tsv(truth, file.path(task_dir, "truth", "tf_time.tsv"))
