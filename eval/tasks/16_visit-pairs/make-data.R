# Purpose: write the seeded inputs and gate truth for task 16 (visit and pair selection).
# Run:     Rscript eval/tasks/16_visit-pairs/make-data.R
# Inputs:  none
# Outputs: eval/tasks/16_visit-pairs/inputs/data/runs.tsv  one row per assay run
#          eval/tasks/16_visit-pairs/truth/pairs.tsv       per donor: paired or not, the right runs
# Notes:   12 donors x visits D0, D7, D90. P09 has no D0 run, so it has no pair. Three samples
#          (P03 D0, P06 D7, P11 D0) were run twice; the first run reads 2.5 high and the later
#          rerun is the one to keep. Rows are shuffled. Reruns are byte-identical.

# Setup ------------------------------------------------------------------------
set.seed(1616)
task_dir <- here::here("eval", "tasks", "16_visit-pairs")
write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Plant visits, a missing baseline and reruns ---------------------------------
donors <- sprintf("P%02d", 1:12)
visit_day <- c(D0 = 0, D7 = 7, D90 = 90)
enrolled <- setNames(
  as.Date("2024-01-08") + sample(0:120, length(donors)),
  donors
)
baseline_il6 <- setNames(rnorm(length(donors), mean = 3, sd = 0.6), donors)
runs <- expand.grid(
  visit = names(visit_day),
  donor_id = donors,
  stringsAsFactors = FALSE
)
runs$il6 <- baseline_il6[runs$donor_id] +
  c(D0 = 0, D7 = 1, D90 = 0.2)[runs$visit] +
  rnorm(nrow(runs), sd = 0.3) * (runs$visit != "D0")
runs$run_date <- enrolled[runs$donor_id] + visit_day[runs$visit] + 3
runs <- runs[!(runs$donor_id == "P09" & runs$visit == "D0"), ]

rerun_keys <- c("P03 D0", "P06 D7", "P11 D0")
first_runs <- runs[paste(runs$donor_id, runs$visit) %in% rerun_keys, ]
first_runs$il6 <- first_runs$il6 + 2.5
runs$run_date[
  paste(runs$donor_id, runs$visit) %in% rerun_keys
] <- first_runs$run_date + 37
runs <- rbind(runs, first_runs)
runs$il6 <- round(runs$il6, 2)
runs <- runs[sample(nrow(runs)), ]
runs$run_id <- sprintf("R%03d", seq_len(nrow(runs)))
runs$run_date <- format(runs$run_date)

# Truth: the latest run at D0 and D7 for every donor that has both.
latest <- runs[order(runs$run_date, decreasing = TRUE), ]
latest <- latest[!duplicated(paste(latest$donor_id, latest$visit)), ]
d0 <- latest[latest$visit == "D0", ]
d7 <- latest[latest$visit == "D7", ]
truth <- data.frame(donor_id = donors)
truth$run_d0 <- d0$run_id[match(donors, d0$donor_id)]
truth$run_d7 <- d7$run_id[match(donors, d7$donor_id)]
truth$il6_d0 <- d0$il6[match(donors, d0$donor_id)]
truth$il6_d7 <- d7$il6[match(donors, d7$donor_id)]
truth$paired <- !is.na(truth$run_d0) & !is.na(truth$run_d7)
stopifnot(
  sum(!truth$paired) == 1,
  sum(duplicated(paste(runs$donor_id, runs$visit))) == 3,
  all(first_runs$il6 > 0)
)

# Write inputs and truth -------------------------------------------------------
write_tsv(
  runs[c("run_id", "donor_id", "visit", "run_date", "il6")],
  file.path(task_dir, "inputs", "data", "runs.tsv")
)
write_tsv(truth, file.path(task_dir, "truth", "pairs.tsv"))
