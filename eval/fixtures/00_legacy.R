args = commandArgs(T)
out = args[1]
dir.create(out, showWarnings = F, recursive = T)

x = read.delim("data/counts.tsv", check.names = F, stringsAsFactors = F)
s = read.delim("data/samples.tsv", stringsAsFactors = F)
rownames(x) = x$gene
x = x[, -1]
#x = x[, 1:8]
s = s[match(colnames(x), s$sample_id), ]

# remove bad ones
keep = colSums(x) > 20000
x = x[, keep]
s = s[keep, ]

x = x[rowSums(x >= 10) >= 3, ]
#x = x[rowSums(x) > 50, ]

cpm = t(t(x) / colSums(x)) * 1e6
l = log2(cpm + 1)
#l = log(cpm + 1)

# ctrl
a = l[, s$group == "ctrl"]
m1 = rowMeans(a)
sd1 = apply(a, 1, sd)
n1 = ncol(a)

# treatA
b = l[, s$group == "treatA"]
m2 = rowMeans(b)
sd2 = apply(b, 1, sd)
n2 = ncol(b)

# treatB
c2 = l[, s$group == "treatB"]
m3 = rowMeans(c2)
sd3 = apply(c2, 1, sd)
n3 = ncol(c2)

res = data.frame(gene = rownames(l), ctrl_mean = m1, ctrl_sd = sd1, treatA_mean = m2, treatA_sd = sd2, treatB_mean = m3, treatB_sd = sd3, stringsAsFactors = F)
res$lfcA = res$treatA_mean - res$ctrl_mean
res$lfcB = res$treatB_mean - res$ctrl_mean
res$dirA = ifelse(res$lfcA > 1, "up", ifelse(res$lfcA < -1, "down", "ns"))
res$dirB = ifelse(res$lfcB > 1, "up", ifelse(res$lfcB < -1, "down", "ns"))
res = res[order(-abs(res$lfcA)), ]
#print(head(res))
write.table(res, paste0(out, "/group_means.tsv"), sep = "\t", quote = F, row.names = F)

tmp = res[res$dirA != "ns" | res$dirB != "ns", ]
tmp = tmp[, c("gene", "lfcA", "lfcB", "dirA", "dirB")]
write.table(tmp, paste0(out, "/top_genes.tsv"), sep = "\t", quote = F, row.names = F)

df = data.frame(sample_id = colnames(x), group = s$group, lib = colSums(x), ngenes = colSums(x > 0), stringsAsFactors = F)
#df$lib2 = df$lib / 1e6
df2 = data.frame(group = c("ctrl", "treatA", "treatB"), n = c(n1, n2, n3), stringsAsFactors = F)
write.table(df, paste0(out, "/sample_qc.tsv"), sep = "\t", quote = F, row.names = F)
write.table(df2, paste0(out, "/group_sizes.tsv"), sep = "\t", quote = F, row.names = F)
