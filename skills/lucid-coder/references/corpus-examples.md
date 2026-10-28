# Corpus examples for the lucid-coder habits

Every fragment is verbatim from `mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings/`; paths are relative to that folder. `cell N` is the raw index in `nb.cells`, markdown cells included. Blocks marked **Rewrite** are not corpus code. Read 2026-10-28.

## Habit 1: Open cold

### Example (R)
immgenT-GP-analysis `repos/immgenT-GP-analysis/script/Figure3.R:36-40` (inputs with producer or `[gap]`)
```r
#   igt1_96_..._ADTonly.Rds                  [primary input Seurat object]
#   L_pm_filtered.rds, F_pm_filtered.rds     [code/pipeline/01b_filter_cells.R]
#   umap_result.rds                          [gap, no producer script here]
#   level_1_AUC_list_figure_no_thymocytes_healthy.rds [code/pipeline/02_compute_auc.R]
#   mean_shifted_log_expr.rds                [gap, no producer script here]
```
immgenT-GP-analysis `repos/immgenT-GP-analysis/script/Figure3.R:53-60` (source lines name the functions used; divider names the step)
```r
source("code/R/plot_utils.R") # lineage_colors()
source("code/R/lineage_plots.R") # plot_gp_swarm(), plot_loadings_on_mde()
source("code/R/volcano_helpers.R") # plot_gp_signature_volcano(), normalize_maxabs()
source("code/R/cross_gp_helpers.R") # plot_cross_gp_heatmap()

# ============================================================
# Load data
# ============================================================
```

### Example (Python)
multiple-myeloma `repos/multiple-myeloma/analyses-rna/01-preprocessing/02-doublets_and_label_transfer/01-predict_doublets.py:16-23`
```python
# %% [markdown]
# # 2025-07-10: Doublet detection
# ### By [Aishwarya Chander](aishwarya.chander@alleninstitute.org), High Resolution Translational Immunology, Allen Institute for Immunology
#
# **Main aim**: Run scrublet for doublet detection across all samples.

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## 1 Imports and functions
```
Outputs appear as a per-panel list in the same header (`Figure3.R:8-16`, `# Panels produced:` with one line per panel).

### Rejected
cns-haematopoiesis `repos_new/cns-haematopoiesis/analysis/scRNA_seq/02_Integration/13_standard_reintegation_script.py:1-6`
```python
#!/usr/bin/env python3
import os, sys, argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
```
Nothing says what the script makes or what it reads, so a cold reader must run it to find out.

## Habit 2: Knobs at the top

### Example (R)
immgenT-GP-analysis `repos/immgenT-GP-analysis/script/Figure1.R:114-121` (the only UPPER_SNAKE block in the R corpus)
```r
N_TOP   <- 5                                      # top up-genes taken per GP
SIG_THR <- 0.1                                    # min per-GP-normalized loading
GP_HIGHLIGHTS <- c(                               # GP -> highlight color
  GP68  = "pink2",  GP58  = "orange2", GP35  = "purple", GP171 = "blue",
  GP1   = "cyan2",  GP56  = "red2",    GP161 = "brown",  GP6   = "green2",
  GP7   = "green3", GP196 = "yellow3")
GP_COL   <- "darkgrey"                            # non-highlighted GP nodes
GENE_COL <- "#C7A76C"                             # gene nodes (tan)
```
MOCHA_Manuscript `repos/MOCHA_Manuscript/Fig5/Fig5_chromVarModeling.R:49-50` (formula as an object)
```r
linFormula <- exp ~ Age + sex_at_birth + time + (1|PTID)
quadFormula <- exp ~ Age + sex_at_birth + time + time_sqrd + (1|PTID)
```

### Example (Python)
sound-life-scrna-analysis `repos/sound-life-scrna-analysis/02-reference_labeling/05-Python_subset_qc_filtering.ipynb` cell 4
```python
max_mito = 10
min_genes = 200
max_genes = 5000
```

### Rejected
immgenT-GP-analysis `repos/immgenT-GP-analysis/script/Figure1.R:66-70`
```r
  inner_join(seurat_meta_filtered %>% select(cellID, annotation_level1), by = "cellID") %>%
  filter(annotation_level1 != "thymocyte")
max_total <- 1000000
min_per_group <- 300
cap_per_group <- 20000
```
Three subsampling knobs sit mid-panel after a join, with no reason beside them.

## Habit 3: Say it once, inline until it repeats

### Example (R)
dreamlet_analysis `repos_new/dreamlet_analysis/Kfoury_CancerCell_2021/Kfoury_CancerCell_2021.Rmd:181-186` (one named vector drives one map; the key becomes a column)
```r
names(contrasts) %>%
  map_df(function(x) 
    topTable(res.dl, coef = x, number=Inf) %>%
    as_tibble %>%
    mutate(coef = x)) %>%
    write.table(file, quote=FALSE, sep='\t', row.names=FALSE)
```
immgenT-GP-analysis `repos/immgenT-GP-analysis/script/Figure3.R:348-350`
```r
gp_letter <- c("GP3" = "3J", "GP29" = "3K", "GP22" = "3L")
for (gp_name in names(gp_letter)) {
  gp_num <- as.numeric(sub("^GP", "", gp_name))
```
aifi-healthy-pbmc-reference `repos/aifi-healthy-pbmc-reference/08-Differential_expression/36-R_pseudobulk_deseq2.ipynb` cell 20
```r
cohort_res <- map2_dfr(
    deseq_res, names(deseq_res),
    format_deseq_results,
    group_by = "cohort.cohortGuid",
    fg = "BR1", bg = "BR2"
)
```

### Example (Python)
sound-life-scrna-analysis `repos/sound-life-scrna-analysis/04-file-sets/Python_assemble_deg_results.ipynb` cell 12, then cell 16
```python
contrast_sets = {
    'cohort.cohortGuid': { 'fg': 'Older Adult', 'bg': 'Younger Adult' },
    'CMV': { 'fg': 'Positive', 'bg': 'Negative' },
    'subject.biologicalSex': {'fg': 'Male', 'bg': 'Female' },
    'sample.visitName': {'fg': 'Flu Year 1 Day 7', 'bg': 'Flu Year 1 Day 0'}
}
```
```python
for contrast in contrasts:
    df = deg.loc[deg['contrast'] == contrast].copy()
    contrast_name = contrast_to_name[contrast]

    df['fg'] = contrast_sets[contrast]['fg']
    df['bg'] = contrast_sets[contrast]['bg']
```

### Rejected
MOCHA_Manuscript `repos/MOCHA_Manuscript/Fig3/Fig3E_DiscriminativePower.R:40-44` (the first of eight pasted `mclapply` blocks, differing only in `allFeatures[[n]]`)
```r
subClustListsMocha <- mclapply(1:1000, function(x){
                subMat <- t(tileMat[rownames(tileMat) %in% sample(allFeatures[[1]],50),])
                 kClust <- kmeans(subMat, centers=2, iter.max = 50, nstart = 20 )$cluster
                findConcordance(meta, kClust)
    }, mc.cores = 60)
```
sound-life-scrna-analysis `repos/sound-life-scrna-analysis/02-reference_labeling/05-Python_subset_qc_filtering.ipynb` cell 18, then the pasted copy at cell 27
```python
search_df = search_df[search_df['name'].str.contains('.h5ad')]
```
```python
search_df = ps_df[ps_df['name'].str.contains('.h5ad')]
```
The copy filters `ps_df` instead of `search_df`, so the `search_id` filter applied in cell 25 is lost.

## Habit 4: Trust the package

### Example (R)
dreamlet_analysis `repos_new/dreamlet_analysis/Nathan_NatImm_2021/Nathan_NatImm_2021.Rmd:153-159` (the package filters, weights and loops over cell types)
```r
form = ~ TB_status + (1|donor) + (1|batch) + (1|season) + (1|sex) + scale(age) + I(scale(age^2)) + scale(prop_NAT)
res.proc = processAssays( pb, form,
              min.samples = 10,
              min.cells=5,
              min.count = 1,
              assays = assayNames(pb), 
              BPPARAM = SnowParam(6))
```
immgenT-GP-analysis `repos/immgenT-GP-analysis/code/R/centered_mean_heatmap.R:67-72` (a count matrix stays a matrix)
```r
  labels <- droplevels(factor(as.character(labels)))
  group_sums <- rowsum(L_mat, group = labels, reorder = TRUE)
  group_counts <- as.integer(table(labels)[rownames(group_sums)])

  list(
    matrix = t(sweep(group_sums, 1L, group_counts, "/")),
```

### Example (Python)
multiple-myeloma `repos/multiple-myeloma/analyses-rna/01-preprocessing/02-doublets_and_label_transfer/01-predict_doublets.py:60-63`
```python
def scrublet_processing(file_name):
    result = read_h5_anndata(file_name)
    sc.external.pp.scrublet(result)
    return result.obs[['barcodes','predicted_doublet','doublet_score']]
```

### Rejected
skin_fibroblast_atlas `repos_new/skin_fibroblast_atlas/figures/supp_fig1_DEGs.ipynb` cell 8
```python
for group in categories:
    genes = adata.uns["de_leiden"]["names"][group]
    logFCs = adata.uns["de_leiden"]["logfoldchanges"][group]
    pvals = adata.uns["de_leiden"]["pvals"][group]
    
    gene_list.extend(genes)
    logFC_list.extend(logFCs)
    pval_list.extend(pvals)
```
This rebuilds `sc.get.rank_genes_groups_df()` by hand and drops the adjusted p-value the package computed.

## Habit 5: Check where data enters, then trust it

### Example (R)
immgenT-GP-analysis `repos/immgenT-GP-analysis/script/ExtendedDataTable1_GP_summary.R:75-81`
```r
# ---- Manual comments (curated input; blank where a GP has no comment) ----
# Keyed by the GP*n* output labels. Insisting on an exact GP1..GP200 match means
# a renamed or dropped row fails loudly here instead of silently blanking a
# comment in the published table.
manual_comments <- read.csv(paste0(curation_path, "GP_manual_annotations.csv"),
                            stringsAsFactors = FALSE, colClasses = "character")
stopifnot(identical(manual_comments$GP, gp_labels))
```
Metadata meets the matrix by name at `Figure3.R:66` above (`seurat_meta[rownames(L_pm_filtered), ]`).

### Example (Python)
immgenT-GP-analysis `repos/immgenT-GP-analysis/script/Figure8b_rematch.py:85-87`
```python
    fr.close()
    if not np.array_equal(obs, r_obs):
        raise ValueError("EBMF and RQVI cell orders differ")
```

### Rejected
immgenT-GP-analysis `repos/immgenT-GP-analysis/code/R/centered_mean_heatmap.R:55-57`
```r
  if (ncol(L) != 200L || nrow(L) != nrow(meta) || anyNA(L)) {
    stop("The healthy non-thymocyte loading matrix has unexpected dimensions or missing values.")
  }
```
An expected-count audit; a wrong `ncol` would already fail at the next `[`, so this line adds reading and no safety. Same for the column asserts at cns `01_processing/05_basic_qc_filtering_add_metadata.py:26-28`: a missing column already raises `KeyError`.

No corpus join states its cardinality; the only argument seen is MOCHA's `multiple='first'` (see What goes wrong). Closest fragment is Figure1.R:66 above (`inner_join(..., by = "cellID")`).
**Rewrite** in the persona's style:
```r
plot_df <- df_mde_a |>
  inner_join(select(seurat_meta_filtered, cellID, annotation_level1),
             by = "cellID", relationship = "one-to-one")
```

No corpus run writes to a fresh run dir; outputs carry the date in the file name instead (sound-life `05` cell 36: `'output/diha_PBMC_qc_filter_counts_{d}.csv'.format(d = date.today())`).
**Rewrite**:
```r
run_dir <- file.path("runs", format(Sys.time(), "%Y%m%d-%H%M%S"))
stopifnot(dir.create(run_dir, recursive = TRUE))  # FALSE, so a stop, if it exists
```

## Habit 6: Scripts compute, notebooks show

### Example (R)
GTestimate-Paper `repos_new/GTestimate-Paper/notebooks/run_analysis.R:9-19` (numbered run order in one driver; `pub_figures.Rmd` only reads)
```r
library(rmarkdown)
render('read_counts_overview.Rmd')
render('select_cells.Rmd')
render('umi_correction.Rmd')
render('downsampling_results.Rmd')
render('cta_seq_results.Rmd')
render('downstream_analysis_results.Rmd')
render('spatial_transcriptomics_results.Rmd')
render('toy_example.Rmd')

render('pub_figures.Rmd')
```
IFN_manuscript `repos/IFN_manuscript/Figures/Main/IFN_Fig1.ipynb` cells 13-14 (the figure notebook reads a saved object and writes the source data from the plotted object)
```r
umap_df <- p1[[1]]$data
umap_df$cell_id <- rownames(umap_df)
```
```r
umap_df %>% fwrite(file.path(wd, "figures_source_data", "Fig1b.csv"))
```

### Example (Python)
multiple-myeloma `repos/multiple-myeloma/analyses-rna/01-preprocessing/02-doublets_and_label_transfer/01-predict_doublets.py:67-74` (a numbered step writes one table named for what it holds)
```python
# max_workers=30; adjust to match available cores (check with `nproc`)
results = []
with ThreadPoolExecutor(max_workers=30) as executor:  
    for result in executor.map(scrublet_processing, metadata['manual.file_paths']):
        results.append(result)

final_result = pd.concat(results, ignore_index=True)
final_result.to_parquet(file_path_csvs+'all_doublet_scores.parquet')
```
Batch steps run under a scheduler: cns `01_processing/bash_scripts/01_submit_generate_h5ad_bone_atlas.sh` wraps one `.py` in `bsub`.

### Rejected
sound-life-scrna-analysis `repos/sound-life-scrna-analysis/02-reference_labeling/10a-Python_cluster_L3_b_cell_data.ipynb` cell 27
```python
    out_file = 'output/diha_clustered_celltypist_L3_{ct}_{d}.h5ad'.format(ct = out_type, d = date.today())
    if os.path.isfile(out_file):
        print('Previously processed {ct}; Skipping.'.format(ct = out_type))
        out_files.append(out_file)
    else:
        adata = process_adata(adata, resolution = 2)
        adata.write_h5ad(out_file)
```
A notebook clusters and writes the object, and `10a`-`10g` are seven copies of it per cell class; a parameterised script would run once.

## Habit 7: Every number reruns

### Example (R)
GTestimate-Paper `repos_new/GTestimate-Paper/notebooks/downstream_analysis_results.Rmd:22-24` (seed in the setup chunk; the repo also commits `notebooks/renv.lock`, the only lockfile in the corpus)
```r
umap_dot_size <- 0.15

set.seed(1)
```
dreamlet_analysis `repos_new/dreamlet_analysis/Nathan_NatImm_2021/Nathan_NatImm_2021.Rmd:84-85` (explicit factor levels; the file ends with `sessionInfo()` at line 473)
```r
sce$cluster_name = factor(sce$cluster_name, unique(tab$cluster_name))
sce$cluster_ids = factor(sce$cluster_ids, unique(tab$cluster_ids))
```
The family is named where it is adjusted: `mutate(fdr = p.adjust(p.value, "fdr"))` after `group_by(assay, name)` (dreamlet `PsychAD_r0/intersection_QTL.Rmd:29-31`).

### Example (Python)
multiple-myeloma `repos/multiple-myeloma/analyses-rna/00-utilities/functions/python/process_scrna_data.py:92-93`
```python
    # UMAP min_dist=0.45 preserves local structure; random_state=0 for reproducibility
    sc.tl.umap(adata, min_dist=0.45, random_state=0, n_components=2)
```
sound-life-scrna-analysis `repos/sound-life-scrna-analysis/02-reference_labeling/05-Python_subset_qc_filtering.ipynb` cells 44-45 (the run ends by recording its inputs, then the session)
```python
hisepy.upload.upload_files(
    files = out_list,
    study_space_id = study_space_uuid,
    title = title,
    input_file_ids = in_files,
    destination = search_id
)
```
```python
import session_info
session_info.show()
```

### Rejected
dreamlet_analysis `repos_new/dreamlet_analysis/Nathan_NatImm_2021/Nathan_NatImm_2021.Rmd:59-61`
```r
path = "/sc/arion/projects/CommonMind/hoffman/scRNAseq_data/Nathan_NatImm_2021/"
file = paste0(path, "/Nathan_NatImm_2021.h5ad")
sce = readH5AD(file, use_hdf5=TRUE)   
```
An absolute home path runs on one machine only; IHA `Figure3/02_FigB_C.ipynb` cell 1 reads from `/home/workspace/` the same way.

`suppressPackageStartupMessages` appears in 12 of 17 repos (for example immgenT `script/Figure1.R:111-113`; dreamlet `PsychAD_r0/PsychAD_r0_analysis.Rmd:28-31`); the native pipe in 2 (immgenT `script/Figure6.R:477`; multiple-myeloma `analyses-rna/04-pseudobulk/l3-pseudobulk/12-plot-deseq2-bars.R:157`), and `%>%` is the norm.

## Habit 8: Figures show the data

### Example (R)
immgenT-GP-analysis `repos/immgenT-GP-analysis/script/Figure4.R:522-530` (bars from 0 via `expansion(mult = c(0, 0.02))` at line 501; source data from the plotted object, named by panel)
```r
ggsave(
  paste0(figure_path, "4c.pdf"), p_4c,
  width = bar_figure_width, height = bar_figure_height, dpi = 300
)

# Every bar segment of all seven rows, raw and as a share of its cluster's
# total; drawn_in_4c marks the row this panel shows.
bar_record <- do.call(rbind, lapply(bar_data, `[[`, "record"))
write.csv(bar_record, file.path(record_path, "4c_mean_loading_by_cluster.csv"), row.names = FALSE)
```
ALTRA-manuscript `repos/ALTRA-manuscript/Figures/Figure3_IC/IC_flow_bcell_data_analysis_final.R:106-109` (each subject a point over the summary; units on the axis)
```r
  ggplot(perTNFastat, aes(status, PTNFa, fill=status))  + 
    geom_boxplot(width=0.5, color="black", position="dodge") + 
    geom_jitter(aes(status, PTNFa), size=8, height=0, width=0.03) +
    ylab("TNFa+ Naive B cells (%)") +
```
The interval is named where drawn: `geom_errorbar(aes(ymin= rho - 1.96*rho.se, ymax= rho + 1.96*rho.se), width=0)` (dreamlet `Velmeshev_2019/Velmeshev_2019.Rmd:295`).

### Example (Python)
IHA-Figure `repos/IHA-Figure/Figure3/02_FigB_C.ipynb` cell 9 (panels share one axis grid; size in mm)
```python
fig, axes = plt.subplots(
    nrows, ncols,
    figsize=(60 / 25.4, 18 / 25.4),
    gridspec_kw={'wspace': 0, 'hspace': 0}
)

plot_idx = 0
for visit in visit_list:
    br1_sub = subject_adata[subject_adata.obs["sample.visitName"] == visit].copy()
```

### Rejected
dreamlet_analysis `repos_new/dreamlet_analysis/Nathan_NatImm_2021/Nathan_NatImm_2021.Rmd:428-429`
```r
  scale_fill_manual(values=c("blue", "red")) +
  scale_color_manual(values=c("blue", "red"))+
```
An unnamed palette assigns colours by factor order, so a relevelled factor swaps CASE and CONTROL silently.

No corpus palette fails on an unknown label; `scale_*_manual` with a named vector silently drops to NA grey. Closest is immgenT `code/R/plot_utils.R:49-55` (`lineage_colors()` as a named vector).
**Rewrite**:
```r
lineage_colors <- c(CD4 = "blue", CD8 = "darkorange2", Treg = "purple")
stopifnot(all(plot_df$Lineage %in% names(lineage_colors)))
p <- p + scale_color_manual(values = lineage_colors)
```

## What goes wrong: evidence

### Donors are the replicates; cells are the measurements
dreamlet_analysis `repos_new/dreamlet_analysis/Nathan_NatImm_2021/Nathan_NatImm_2021.Rmd:87-91`
```r
# create pseudobulk 
pb <- aggregateToPseudoBulk(sce,
    assay = "counts",     
    cluster_id = "cluster_name",  
    sample_id = "donor_batch")
```
Same shape at variancePartition_analysis `repos_new/variancePartition_analysis/yazar_2022.Rmd:67-68`: `# Compute pseudobulk by donor` then `aggregateToPseudoBulk(...)`.

### Position is a permutation waiting for a reorder
immgenT-GP-analysis `repos/immgenT-GP-analysis/script/README.md:618-622`
```text
lives in one named vector per script (`fig6_gating` in `Figure6.R`,
`figs5_gating` in `FigureS5.R`) and the loops iterate over its *names*, so a GP
cannot be drawn under another GP's letter. Earlier versions kept the GP list and
the letters in two separate vectors and assigned the letters positionally, which
silently permuted three published panels -- do not reintroduce that shape.
```
Python form: `df_lat.reindex(adata_full.obs_names.astype(str))` before a positional copy of the latent space (cns `02_Integration/13_standard_reintegation_script.py:97-100`).

### A silent filter is a result nobody chose
sound-life-scrna-analysis `repos/sound-life-scrna-analysis/02-reference_labeling/05-Python_subset_qc_filtering.ipynb` cell 9 (each drop is a data row), then cells 35-36
```python
    # Filter doublets
    counts['n_start'] = [adata.shape[0]]
    counts['n_doublets'] = [sum(adata.obs['predicted_doublet'] == True)]
    adata = adata[adata.obs['predicted_doublet'] == False]
    counts['n_singlets'] = [adata.shape[0]]
```
```python
all_filter_counts = pd.concat(filter_counts)
all_filter_counts
```
```python
counts_file = 'output/diha_PBMC_qc_filter_counts_{d}.csv'.format(d = date.today())
all_filter_counts.to_csv(counts_file)
```
Counter: the lost `search_id` filter at cell 27 (Habit 3 Rejected) removed nothing visible and reported nothing.

### A plausible wrong number
Many-to-many join, silenced: MOCHA_Manuscript `repos/MOCHA_Manuscript/SuppFig5_Mouse_Data_Analysis/evaluateResults.R:260-262`
```r
tile_counts = dplyr::left_join(tile_counts, cells_per_sample, multiple='first')
tss_counts = dplyr::left_join(tss_counts, cells_per_sample, multiple='first')
ctcf_counts = dplyr::left_join(ctcf_counts, cells_per_sample, multiple='first')
```
`multiple='first'` keeps one arbitrary match per row instead of failing.

Label that fell into "Other": multiple-myeloma `repos/multiple-myeloma/analyses-rna/02-object-building/02-pbmc/11-subsample-500k-build_umap.py:55-59`
```python
    'Progenitor cell': '#A31001',
    'Platelet': '#5B3930',
    'Erythrocyte': '#1F1515',
    'Other': '#1F1515'
}
```
"Other" shares a hex with Erythrocyte, so a label that fell through is invisible on the plot.

Mutated view: sound-life-scrna-analysis `repos/sound-life-scrna-analysis/02-reference_labeling/05-Python_subset_qc_filtering.ipynb` cell 9
```python
    adata = adata[adata.obs['predicted_doublet'] == False]
    counts['n_singlets'] = [adata.shape[0]]
    
    # Compute fraction mitochondrial
    adata.var["mito"] = adata.var_names.str.startswith("MT-")
```
`adata` is a view after the `[`; the assignment to `.var` works only because AnnData copies on write with a warning.

Muted warning: skin_fibroblast_atlas `repos_new/skin_fibroblast_atlas/integration/run_nichecompass.py:152` and multiple-myeloma `repos/multiple-myeloma/manuscript-figures/utilities/functions/gex-violins-colors.R:14`
```python
warnings.filterwarnings("ignore")
```
```r
options(warn = -1)
```
sound-life silences by category instead (`05` cell 2: `warnings.simplefilter(action='ignore', category=FutureWarning)`).

### A number without a script is a rumour
immgenT-GP-analysis `repos/immgenT-GP-analysis/script/Figure3.R:38` marks the input it cannot rebuild: `#   umap_result.rds                          [gap, no producer script here]`.

MOCHA_Manuscript `repos/MOCHA_Manuscript/Fig3/Fig3D_MOCHA_Pathways.R:19-22`
```r
DAP_Reactome$TopLevel = PathAnnotDAP_Reactome$SecondLevel = c(rep('Innate Immune System', times = 10), 
		rep('Cytokine Signaling in Immune System', n = 1), 
		rep('Innate Immune System', times = 6), 
		rep('Adaptive Immune System', times = 2), 
```
Pathway labels were typed by hand against a row order that lives nowhere in the code, so the table cannot be regenerated. immgenT `code/R/doc_code.R:4-7` records the same lesson for line-range citations, which "silently went stale every time a script was edited".
