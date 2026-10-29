# Naive answer: the first visit per donor is taken as baseline and duplicate runs are kept, so
# P09 (no D0) pairs D7 with itself and rerun samples give two rows per donor.
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

runs <- read.delim("data/runs.tsv")

base <- runs %>%
  group_by(donor_id) %>%
  filter(visit == sort(visit)[1]) %>%
  ungroup() %>%
  select(donor_id, il6_d0 = il6)
d7 <- runs %>%
  filter(visit == "D7") %>%
  select(donor_id, il6_d7 = il6)
pairs <- inner_join(base, d7, by = "donor_id") %>%
  mutate(change = il6_d7 - il6_d0)
tt <- t.test(pairs$il6_d7, pairs$il6_d0, paired = TRUE)
res <- data.frame(n_donors = nrow(pairs), mean_change = unname(tt$estimate), p_value = tt$p.value)

write.table(pairs, file.path(out, "pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res, file.path(out, "paired_test.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
