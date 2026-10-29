# Naive answer: pools cells per condition and tests each type's share (cells as the unit).
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

cells <- read.delim("data/cells.tsv")
tab <- table(cells$condition, cells$cell_type)
total <- rowSums(tab)

res <- do.call(
  rbind,
  lapply(colnames(tab), function(ct) {
    pt <- prop.test(
      c(tab["Salmonella", ct], tab["Control", ct]),
      c(total["Salmonella"], total["Control"])
    )
    share <- tab[, ct] / total
    data.frame(
      cell_type = ct,
      log_fc = log2(share["Salmonella"] / share["Control"]),
      p_value = pt$p.value
    )
  })
)
res$q_value <- p.adjust(res$p_value, method = "BH")
write.table(res, file.path(out, "composition.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
