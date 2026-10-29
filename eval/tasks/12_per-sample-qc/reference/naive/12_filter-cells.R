# Naive answer: pools the libraries and applies one UMI and gene floor to every cell.
library(scuttle)
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

samples <- c("SRR28679756", "SRR28679757", "SRR28679759")
sce <- DropletUtils::read10xCounts(file.path("data", samples), sample.names = samples)
qc <- perCellQCMetrics(sce)
keep <- qc$sum >= 100 & qc$detected >= 50

res <- data.frame(sample = sce$Sample, barcode = sce$Barcode, keep = keep)
summ <- aggregate(keep ~ sample, res, function(k) c(n_cells = length(k), n_kept = sum(k)))
summ <- data.frame(sample = summ$sample, summ$keep)
write.table(res, file.path(out, "cells_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summ, file.path(out, "qc_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
