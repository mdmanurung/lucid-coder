# Purpose: write the seeded inputs and gate truth for task 11 (contrast direction).
# Run:     Rscript eval/tasks/11_contrast-direction/make-data.R
# Inputs:  none
# Outputs: eval/tasks/11_contrast-direction/inputs/data/samples.tsv     sample_id, group
#          eval/tasks/11_contrast-direction/inputs/data/expression.tsv  log2 levels, proteins x samples
#          eval/tasks/11_contrast-direction/truth/proteins.tsv          planted disease - healthy log2fc
# Notes:   12 samples (6 disease, 6 healthy) x 40 proteins. "disease" sorts before "healthy", so
#          a default factor makes disease the reference and flips every sign. 10 proteins are
#          planted up in disease, 10 down, 20 null. Reruns are byte-identical.

# Setup ------------------------------------------------------------------------
set.seed(1111)
task_dir <- here::here("eval", "tasks", "11_contrast-direction")
write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Plant the effects ------------------------------------------------------------
proteins <- sprintf("PROT%02d", 1:40)
samples <- data.frame(
  sample_id = sprintf("S%02d", 1:12),
  group = sample(rep(c("disease", "healthy"), each = 6))
)
planted_log2fc <- setNames(rep(0, length(proteins)), proteins)
planted <- sample(proteins, 20)
planted_log2fc[planted[1:10]] <- 1.5
planted_log2fc[planted[11:20]] <- -1.5
is_disease <- samples$group == "disease"
expression <- rnorm(length(proteins), mean = 6, sd = 1.5) +
  outer(planted_log2fc, is_disease) +
  matrix(
    rnorm(length(proteins) * nrow(samples), sd = 0.3),
    nrow = length(proteins)
  )
expression <- round(expression, 3)
dimnames(expression) <- list(proteins, samples$sample_id)

# The trap must hold: a default factor flips the sign of every planted protein.
default_fit <- limma::eBayes(limma::lmFit(
  expression,
  model.matrix(~ factor(samples$group))
))
default_log2fc <- default_fit$coefficients[, 2]
stopifnot(all(sign(default_log2fc[planted]) == -sign(planted_log2fc[planted])))

# Write inputs and truth -------------------------------------------------------
write_tsv(samples, file.path(task_dir, "inputs", "data", "samples.tsv"))
write_tsv(
  data.frame(protein = proteins, expression, check.names = FALSE),
  file.path(task_dir, "inputs", "data", "expression.tsv")
)
write_tsv(
  data.frame(
    protein = proteins,
    planted_log2fc = planted_log2fc,
    direction = c("down", "none", "up")[sign(planted_log2fc) + 2]
  ),
  file.path(task_dir, "truth", "proteins.tsv")
)
