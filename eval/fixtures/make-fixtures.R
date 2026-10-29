# Purpose: write the seeded input data and gate truth for every eval task.
# Inputs:  eval/fixtures/00_legacy.R   [hand-written legacy script for task 04]
# Outputs: eval/tasks/<id>/inputs/data/*, eval/tasks/<id>/truth/*,
#          eval/tasks/04_refactor-legacy/inputs/scripts/00_legacy.R
# Run:     Rscript eval/fixtures/make-fixtures.R
# Notes:   [Stochastic] one seed per task, so editing one task leaves the others unchanged.
#          Reruns are byte-identical; eval/tests/gates.test.js checks this.
#          Tasks 07-16 are not built here: each has its own make-inputs.*/make-data.R in its folder.

# ---- Setup ----
TASKS_DIR <- here::here("eval", "tasks")
FIXTURES_DIR <- here::here("eval", "fixtures")

write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---- 01_pseudobulk-de ----
# Six donors x 200 cells x 50 genes. Two genes carry the test:
#   PSEUDOREP_GENE: donor effects differ a lot but the condition means barely differ, so a test
#     that treats cells as replicates calls it and a donor-level test does not.
#   TRUE_GENE: a large, consistent donor-level effect that any donor-level test calls.
# Neutral gene names, so the data does not hint which gene is which; truth/genes.tsv maps them.
set.seed(101)
task_dir <- file.path(TASKS_DIR, "01_pseudobulk-de")
CELLS_PER_DONOR <- 200
genes <- sprintf("GENE%02d", 1:50)
PSEUDOREP_GENE <- "GENE17"
TRUE_GENE <- "GENE33"
donors <- data.frame(
  donor_id = sprintf("donor%02d", 1:6),
  condition = rep(c("healthy", "disease"), each = 3)
)

# Share of each cell's counts per gene, before donor effects.
gene_share <- setNames(rlnorm(length(genes), meanlog = 0, sdlog = 0.8), genes)
gene_share <- gene_share / sum(gene_share)
# High expression for PSEUDOREP_GENE keeps Poisson noise small, so the cell-level test is sharp.
gene_share[[PSEUDOREP_GENE]] <- 0.05
gene_share[[TRUE_GENE]] <- 0.01

# log2 donor effects, genes x donors. Null genes vary a little between donors.
donor_effect <- matrix(
  rnorm(length(genes) * nrow(donors), sd = 0.15),
  nrow = length(genes),
  dimnames = list(genes, donors$donor_id)
)
# Healthy mean 0, disease mean 0.3, donor SD about 0.3: donor-level p about 0.3.
donor_effect[PSEUDOREP_GENE, ] <- c(-0.35, 0.10, 0.25, 0.35, -0.05, 0.60)
# Disease 2 log2 units up, equal spread in both groups so a Welch test keeps about 4 df.
donor_effect[TRUE_GENE, ] <- c(-0.10, 0.00, 0.10, 1.90, 2.00, 2.10)

cells <- data.frame(
  cell_id = sprintf("cell%04d", seq_len(nrow(donors) * CELLS_PER_DONOR)),
  donor_id = rep(donors$donor_id, each = CELLS_PER_DONOR)
)
cells$condition <- donors$condition[match(cells$donor_id, donors$donor_id)]
library_size <- round(rlnorm(nrow(cells), meanlog = log(3000), sdlog = 0.25))
donor_share <- gene_share * 2^donor_effect
donor_share <- sweep(donor_share, 2, colSums(donor_share), "/")
expected <- donor_share[, cells$donor_id] * rep(library_size, each = length(genes))
cell_counts <- matrix(
  rnbinom(length(expected), mu = expected, size = 20),
  nrow = length(genes),
  dimnames = list(genes, cells$cell_id)
)

write_tsv(
  data.frame(gene = genes, cell_counts, check.names = FALSE),
  file.path(task_dir, "inputs", "data", "counts.tsv")
)
write_tsv(cells, file.path(task_dir, "inputs", "data", "cells.tsv"))
write_tsv(
  data.frame(role = c("pseudoreplicated", "true_effect"), gene = c(PSEUDOREP_GENE, TRUE_GENE)),
  file.path(task_dir, "truth", "genes.tsv")
)

# The trap must hold: a cell-level test calls PSEUDOREP_GENE, donor-level tests do not, and
# donor-level tests call TRUE_GENE. BH over the 50 genes, as the gate reads it.
bh_q <- function(p_values, gene) p.adjust(p_values, method = "BH")[[gene]]
is_disease_cell <- cells$condition == "disease"
cell_log_norm <- log1p(sweep(cell_counts, 2, colSums(cell_counts), "/") * 1e4)
cell_t_p <- apply(cell_log_norm, 1, function(e) {
  t.test(e[is_disease_cell], e[!is_disease_cell])$p.value
})
cell_wilcox_p <- apply(cell_log_norm, 1, function(e) {
  wilcox.test(e[is_disease_cell], e[!is_disease_cell], exact = FALSE)$p.value
})

donor_sums <- sapply(donors$donor_id, function(d) rowSums(cell_counts[, cells$donor_id == d]))
donor_log_cpm <- log2(sweep(donor_sums, 2, colSums(donor_sums), "/") * 1e6 + 1)
donor_means <- sapply(donors$donor_id, function(d) rowMeans(cell_log_norm[, cells$donor_id == d]))
is_disease_donor <- donors$condition == "disease"
welch_p <- function(m) {
  apply(m, 1, function(e) t.test(e[is_disease_donor], e[!is_disease_donor])$p.value)
}
pooled_p <- function(m) {
  apply(m, 1, function(e) {
    t.test(e[is_disease_donor], e[!is_disease_donor], var.equal = TRUE)$p.value
  })
}
p_by_method <- list(
  "cell t-test" = cell_t_p,
  "cell Wilcoxon" = cell_wilcox_p,
  "donor sum logCPM, Welch" = welch_p(donor_log_cpm),
  "donor sum logCPM, pooled t" = pooled_p(donor_log_cpm),
  "donor mean log-norm, Welch" = welch_p(donor_means),
  "donor mean log-norm, pooled t" = pooled_p(donor_means)
)
trap <- data.frame(
  method = names(p_by_method),
  q_pseudorep = sapply(p_by_method, bh_q, gene = PSEUDOREP_GENE),
  q_true = sapply(p_by_method, bh_q, gene = TRUE_GENE),
  row.names = NULL
)
is_cell_level <- startsWith(trap$method, "cell")
message("task 01 trap (BH q over 50 genes):")
message(paste(capture.output(print(trap, digits = 3)), collapse = "\n"))
# Margins well clear of the gate's 0.05, so any reasonable donor-level method lands the same side.
stopifnot(
  all(trap$q_pseudorep[is_cell_level] < 1e-6),
  all(trap$q_pseudorep[!is_cell_level] > 0.2),
  all(trap$q_true[!is_cell_level] < 0.01)
)

# ---- 02_metadata-join ----
# Counts columns and metadata rows are in different random orders, so a positional join
# mislabels most samples.
set.seed(202)
task_dir <- file.path(TASKS_DIR, "02_metadata-join")
samples <- data.frame(
  sample_id = sprintf("S%02d", 1:12),
  condition = rep(c("healthy", "disease"), each = 6),
  batch = rep(c("batch1", "batch2", "batch3"), times = 4)
)
genes <- sprintf("GENE%03d", 1:200)
sample_depth <- round(runif(nrow(samples), min = 0.5, max = 1.5), 2)
sample_counts <- matrix(
  rnbinom(
    length(genes) * nrow(samples),
    mu = outer(rlnorm(length(genes), 4, 1), sample_depth),
    size = 5
  ),
  nrow = length(genes),
  dimnames = list(genes, samples$sample_id)
)
samples$library_size <- colSums(sample_counts)

count_order <- sample(samples$sample_id)
metadata_order <- sample(nrow(samples))
metadata <- samples[metadata_order, c("sample_id", "condition", "batch")]
positional <- samples[match(count_order, samples$sample_id), ]
n_mislabelled <- sum(
  positional$condition != metadata$condition | positional$batch != metadata$batch
)
message("task 02: a positional join mislabels ", n_mislabelled, " of 12 samples")
stopifnot(n_mislabelled >= 6)

write_tsv(
  data.frame(gene = genes, sample_counts[, count_order], check.names = FALSE),
  file.path(task_dir, "inputs", "data", "counts.tsv")
)
write_tsv(metadata, file.path(task_dir, "inputs", "data", "metadata.tsv"))
write_tsv(samples, file.path(task_dir, "truth", "samples.tsv"))

# ---- 03_figure-palette ----
set.seed(303)
task_dir <- file.path(TASKS_DIR, "03_figure-palette")
cell_types <- c("t_cell", "b_cell", "nk_cell", "monocyte", "dendritic_cell")
donors <- data.frame(
  donor_id = sprintf("donor%02d", 1:8),
  condition = rep(c("healthy", "disease"), each = 4)
)
# Dirichlet draws via normalised gammas; disease shifts share from t cells to monocytes.
type_weight <- c(t_cell = 40, b_cell = 15, nk_cell = 12, monocyte = 25, dendritic_cell = 3)
proportions <- do.call(
  rbind,
  lapply(seq_len(nrow(donors)), function(i) {
    weight <- type_weight
    if (donors$condition[[i]] == "disease") {
      weight[c("t_cell", "monocyte")] <- c(28, 37)
    }
    draw <- rgamma(length(weight), shape = weight)
    data.frame(
      donor_id = donors$donor_id[[i]],
      condition = donors$condition[[i]],
      cell_type = cell_types,
      proportion = round(draw / sum(draw), 4)
    )
  })
)
# Palette rows in a different order from the data, keyed by label.
palette <- data.frame(
  label = c("monocyte", "healthy", "nk_cell", "t_cell", "disease", "dendritic_cell", "b_cell"),
  hex = c("#AA3377", "#4477AA", "#66CCEE", "#228833", "#EE6677", "#BBBBBB", "#CCBB44")
)

write_tsv(proportions, file.path(task_dir, "inputs", "data", "proportions.tsv"))
write_tsv(palette, file.path(task_dir, "inputs", "data", "palette.tsv"))
write_tsv(
  proportions[c("donor_id", "cell_type")],
  file.path(task_dir, "truth", "donor_cell_types.tsv")
)

# ---- 04_refactor-legacy ----
# Ten samples in three groups; one control sample is shallow, so the legacy library-size
# filter (colSums > 20000) removes it and changes every output.
set.seed(404)
task_dir <- file.path(TASKS_DIR, "04_refactor-legacy")
samples <- data.frame(
  sample_id = sprintf("S%02d", 1:10),
  group = c(rep("ctrl", 4), rep("treatA", 3), rep("treatB", 3))
)
genes <- sprintf("GENE%03d", 1:300)
gene_mean <- rlnorm(length(genes), meanlog = 4.5, sdlog = 1.5)
group_effect <- cbind(
  ctrl = 1,
  treatA = 2^rnorm(length(genes), sd = 0.8),
  treatB = 2^rnorm(length(genes), sd = 0.8)
)
sample_depth <- c(1.0, 1.1, 0.9, 0.04, 1.2, 0.8, 1.0, 1.1, 0.9, 1.0)
legacy_counts <- matrix(
  rnbinom(
    length(genes) * nrow(samples),
    mu = gene_mean * group_effect[, samples$group] * rep(sample_depth, each = length(genes)),
    size = 10
  ),
  nrow = length(genes),
  dimnames = list(genes, samples$sample_id)
)
library_sizes <- colSums(legacy_counts)
message("task 04 library sizes: ", toString(library_sizes))
stopifnot(sum(library_sizes <= 20000) == 1, min(library_sizes[library_sizes > 20000]) > 40000)

write_tsv(
  data.frame(gene = genes, legacy_counts, check.names = FALSE),
  file.path(task_dir, "inputs", "data", "counts.tsv")
)
write_tsv(samples, file.path(task_dir, "inputs", "data", "samples.tsv"))
dir.create(file.path(task_dir, "inputs", "scripts"), recursive = TRUE, showWarnings = FALSE)
invisible(file.copy(
  file.path(FIXTURES_DIR, "00_legacy.R"),
  file.path(task_dir, "inputs", "scripts", "00_legacy.R"),
  overwrite = TRUE
))

# ---- 05_recode-trap ----
# Two cohorts of 12 donors. Two sheet rows carry a malformed cohort label ("BR 1", "br2"), so a
# catch-all recode arm (TRUE ~ "older") silently files one young donor under "older".
set.seed(505)
task_dir <- file.path(TASKS_DIR, "05_recode-trap")
AGE_GROUP <- c(BR1 = "young", BR2 = "older")
cohort_samples <- data.frame(
  sample_id = sprintf("S%02d", 1:24),
  cohort = rep(names(AGE_GROUP), each = 12)
)
cohort_samples$age_group <- AGE_GROUP[cohort_samples$cohort]
cohort_samples$age <- ifelse(
  cohort_samples$age_group == "young",
  sample(22:35, nrow(cohort_samples), replace = TRUE),
  sample(62:79, nrow(cohort_samples), replace = TRUE)
)
cohort_samples$malformed <- cohort_samples$sample_id %in% c("S07", "S19")
sheet <- cohort_samples
sheet$cohort[sheet$sample_id == "S07"] <- "BR 1"
sheet$cohort[sheet$sample_id == "S19"] <- "br2"
stopifnot(sum(!sheet$cohort %in% names(AGE_GROUP)) == 2)
scores <- data.frame(
  sample_id = sample(cohort_samples$sample_id),
  score = round(rnorm(nrow(cohort_samples), mean = 50, sd = 8), 2)
)

write_tsv(
  sheet[c("sample_id", "cohort", "age")],
  file.path(task_dir, "inputs", "data", "samples.tsv")
)
write_tsv(scores, file.path(task_dir, "inputs", "data", "scores.tsv"))
write_tsv(
  cohort_samples[c("sample_id", "cohort", "age_group", "malformed")],
  file.path(task_dir, "truth", "samples.tsv")
)

# ---- 06_grep-trap ----
# 300 genes x 10 samples: 20 ERCC spike-ins and no MT- genes, so `x[-grep("^MT-", x$gene), ]`
# drops every row (integer(0) index), while `!grepl()` keeps them all.
set.seed(606)
task_dir <- file.path(TASKS_DIR, "06_grep-trap")
genes <- c(sprintf("GENE%03d", 1:280), sprintf("ERCC-%05d", 1:20))
genes <- sample(genes)
control_samples <- sprintf("S%02d", 1:10)
control_counts <- matrix(
  rnbinom(length(genes) * length(control_samples), mu = rlnorm(length(genes), 4, 1), size = 5),
  nrow = length(genes),
  dimnames = list(genes, control_samples)
)
control_table <- data.frame(gene = genes, control_counts, check.names = FALSE)
stopifnot(
  sum(grepl("^MT-", genes)) == 0,
  sum(grepl("^ERCC-", genes)) == 20,
  nrow(control_table[-grep("^MT-", control_table$gene), ]) == 0
)

write_tsv(control_table, file.path(task_dir, "inputs", "data", "counts.tsv"))
write_tsv(
  data.frame(gene = genes, gene_set = ifelse(grepl("^ERCC-", genes), "ERCC", "none")),
  file.path(task_dir, "truth", "genes.tsv")
)
