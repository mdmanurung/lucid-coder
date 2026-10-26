# Naive answer: a catch-all recode arm files every non-BR1 label under "older".
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim("data/samples.tsv")
scores <- read.delim("data/scores.tsv")

labelled <- samples %>%
  mutate(age_group = case_when(cohort == "BR1" ~ "young", TRUE ~ "older")) %>%
  select(sample_id, cohort, age_group)
res <- labelled %>%
  left_join(scores, by = "sample_id") %>%
  group_by(age_group) %>%
  summarise(n = n(), mean_score = mean(score))

write.table(labelled, file.path(out, "samples_labelled.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res, file.path(out, "cohort_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
