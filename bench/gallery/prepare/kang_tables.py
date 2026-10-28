"""Write the Kang 2018 tables the gallery tasks start from.

Reads the plots bench's kang2018.h5ad (raw counts), else the local pertpy file that
common/prepare/kang2018.py copies. Outputs, all tab-separated, in the cache dir:
  sample_means.tsv      replicate, label, cell_type, gene, mean_expr (slope)
  box_few.tsv           sample, condition, cell_type, gene, mean_expr (box_few):
                        ctrl from the first 4 donors, stim from the other 4, so unpaired
  de.tsv                cell_type, gene, log2fc, se, pvalue, padj (forest, volcano)
  heatmap_expr.tsv      cell_type, replicate, label, gene, log2cpm (heatmap, complexheatmap)
  pb_logcpm_<ct>.tsv    gene x sample log2 CPM, one per DE cell type (pseudobulk_pca)
  celltype_counts.tsv   sample, replicate, label, cell_type, n_cells (composition)
  cell_qc.tsv           cell_type, label, total_counts, n_genes (density)
  marker_summary.tsv    cell_type, gene, mean_expr, frac_expressing (dotplot)
"""

import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

out = Path(sys.argv[1])
local = Path(__file__).resolve().parents[2] / "plots" / "cache" / "kang2018.h5ad"
source = os.environ.get(
    "LUCID_KANG_H5AD",
    "/exports/para-lipg-hpc/mdmanurung/.cache/cellimo-docs-data/kang_2018.h5ad",
)
adata = ad.read_h5ad(local if local.is_file() else source)
obs = adata.obs[["label", "replicate", "cell_type"]].astype(str)
obs["sample"] = obs["replicate"] + "_" + obs["label"]
counts = adata.X.tocsr().astype("float64")

DE_TYPES = ["B cells", "CD4 T cells", "CD8 T cells", "CD14+ Monocytes", "NK cells"]
ISGS = ["ISG15", "IFI6", "MX1", "IFIT3", "ISG20", "IFITM3", "LY6E"]
MARKERS = ["CD3E", "IL7R", "CD8A", "NKG7", "GNLY", "MS4A1", "CD79A", "CD14", "LYZ",
           "FCGR3A", "MS4A7", "FCER1A", "CST3", "PPBP", "PF4", "HBB"]  # fmt: skip
slug = {ct: ct.lower().replace("+", "").replace(" ", "_") for ct in DE_TYPES}


def write(df, name):
    df.to_csv(out / name, sep="\t", index=False)
    print(name, df.shape)


# per-cell log-normalised expression (log1p CP10k)
lognorm = adata.copy()
lognorm.X = counts.copy()
sc.pp.normalize_total(lognorm, target_sum=1e4)
sc.pp.log1p(lognorm)

# per-sample mean of ISGs per cell type
isg = pd.DataFrame(lognorm[:, ISGS].X.toarray(), columns=ISGS, index=obs.index)
isg = pd.concat([obs[["replicate", "label", "cell_type"]], isg], axis=1)
means = (
    isg.groupby(["replicate", "label", "cell_type"])[ISGS]
    .mean()
    .reset_index()
    .melt(["replicate", "label", "cell_type"], var_name="gene", value_name="mean_expr")
)
write(means, "sample_means.tsv")

# unpaired: ctrl from the first four donors, stim from the last four
donors = sorted(obs["replicate"].unique())
keep = ((means["label"] == "ctrl") & means["replicate"].isin(donors[:4])) | (
    (means["label"] == "stim") & means["replicate"].isin(donors[4:])
)
box = means[keep].rename(columns={"replicate": "sample", "label": "condition"})
write(box, "box_few.tsv")

# pseudobulk log2 CPM and a paired t-test per gene (stim - ctrl within donor)
de_rows, heat_rows = [], []
for ct in DE_TYPES:
    cells = (obs["cell_type"] == ct).to_numpy()
    pb = pd.DataFrame(
        counts[cells].toarray(), columns=adata.var_names, index=obs.index[cells]
    ).groupby(obs.loc[cells, "sample"].to_numpy()).sum()
    pb = pb.loc[:, pb.sum() >= 10]
    logcpm = np.log2(pb.div(pb.sum(axis=1), axis=0) * 1e6 + 1)
    logcpm.T.rename_axis("gene").reset_index().to_csv(
        out / f"pb_logcpm_{slug[ct]}.tsv", sep="\t", index=False, float_format="%.4f"
    )
    diff = (
        logcpm.loc[[f"{d}_stim" for d in donors]].to_numpy()
        - logcpm.loc[[f"{d}_ctrl" for d in donors]].to_numpy()
    )
    se = diff.std(axis=0, ddof=1) / np.sqrt(len(donors))
    ok = se > 0
    lfc = diff.mean(axis=0)[ok]
    p = 2 * stats.t.sf(np.abs(lfc / se[ok]), len(donors) - 1)
    de_rows.append(
        pd.DataFrame(
            {
                "cell_type": ct,
                "gene": logcpm.columns[ok],
                "log2fc": lfc,
                "se": se[ok],
                "pvalue": p,
                "padj": stats.false_discovery_control(p),
            }
        )
    )
    heat_rows.append((ct, logcpm))
de = pd.concat(de_rows)
write(de, "de.tsv")

# 30 heatmap genes kept in every cell type: the top 6 by padj in each, then the next best
common = set.intersection(*(set(logcpm.columns) for _, logcpm in heat_rows))
ranked = de[de["gene"].isin(common)].sort_values("padj")
top = ranked.groupby("cell_type").head(6)["gene"].unique().tolist()
rest = [g for g in ranked["gene"].unique() if g not in top]
genes = (top + rest)[:30]
heat = []
for ct, logcpm in heat_rows:
    long = logcpm[genes].rename_axis("sample").reset_index()
    long = long.melt("sample", var_name="gene", value_name="log2cpm")
    long[["replicate", "label"]] = long["sample"].str.rsplit("_", n=1, expand=True)
    long["cell_type"] = ct
    heat.append(long[["cell_type", "replicate", "label", "gene", "log2cpm"]])
write(pd.concat(heat), "heatmap_expr.tsv")

# cells per sample and cell type
comp = obs.groupby(["sample", "replicate", "label", "cell_type"]).size()
write(comp.rename("n_cells").reset_index(), "celltype_counts.tsv")

# per-cell library size and detected genes
qc = obs[["cell_type", "label"]].copy()
qc["total_counts"] = np.asarray(counts.sum(axis=1)).ravel().astype(int)
qc["n_genes"] = np.asarray((counts > 0).sum(axis=1)).ravel()
write(qc, "cell_qc.tsv")

# marker mean expression and fraction expressing per cell type
mk = lognorm[:, MARKERS].X.toarray()
grp = obs["cell_type"].to_numpy()
mean = pd.DataFrame(mk, columns=MARKERS).groupby(grp).mean()
frac = pd.DataFrame(mk > 0, columns=MARKERS).groupby(grp).mean()
dot = mean.stack().rename("mean_expr").to_frame().join(
    frac.stack().rename("frac_expressing")
)
write(dot.rename_axis(["cell_type", "gene"]).reset_index(), "marker_summary.tsv")
