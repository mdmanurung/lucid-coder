# Corpus patterns: concise code and project shape

Source: three Opus agents read the code itself, not the summaries, on 2026-10-27. They covered about 140 scripts and notebooks in 17 repos under `mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings/repos/` and `repos_new/`. The user named immgenT-GP-analysis as the R exemplar, and sound-life-scrna-analysis and aifi-healthy-pbmc-reference as the Python exemplars. Paths below are relative to those two folders. Corpus code is evidence, not a template. The last section lists where the plugin overrides it.

## Why scripts got long

A 2026-10-27 demo showed lucid-coder doubling script length against baseline. The added lines were column checks at load, a `message()` block per drop, long headers and constants for one-use regexes. Nothing in the corpus supports those. The corpus does support three things: an alignment check before a positional copy (`stopifnot(identical(manual_comments$GP, gp_labels))`, immgenT `script/ExtendedDataTable1_GP_summary.R:81`), seeds, and inspection output in notebooks.

Expert files leave out:

- **Load checks.** IHA has `stop(` in 7 of about 287 files, all inside helpers. Hoffman reads an h5ad and goes straight to analysis (`dreamlet_analysis/Nathan_NatImm_2021/Nathan_NatImm_2021.Rmd:59-61`).
- **Drop logging in prose.** IHA and IFN have no `message()` calls. sound-life logs removals as data rows, one table per notebook (`02-.../05:9`).
- **Long headers.** MOCHA uses one line, and Rmds use YAML plus a sentence. immgenT's long headers in `script/` record a 2026 refactor, so they are not a template.
- **Named one-use values.** Cutoffs and regexes stay inline (`FDR < 0.2`, MOCHA `Fig3D_MOCHA_Pathways.R:6`). Names appear for values that are reused or vary per run (`contrasts`, `Kfoury_CancerCell_2021.Rmd:148`).
- **CLI scaffolding.** `argparse`, `__main__` and `assert` are absent from aifi, sound-life, reflex and single_cell_aging. IHA, ALTRA and MOCHA have 3 `argparse` or `commandArgs` files between them.

## Brevity idioms with two or more witnesses

### R

- **Split, map, combine once.** `type_meta <- split(meta, meta$AIFI_L3)`, then `map(type_meta, function(meta) mat[keep_genes, meta$barcodes])` (aifi `07-Pseudobulk/35-R_pseudobulk_seurat.ipynb:10`). Also `lapply(names(structure_plot_panels), bar_means)`, then one `do.call(rbind, ...)` (immgenT `script/Figure4.R:494-529`). All five R repos do this.
- **Stack results per contrast with the key as a column.** `names(contrasts) |> map_df(\(x) topTable(res, coef = x, number = Inf) |> mutate(coef = x))` (`Kfoury_CancerCell_2021.Rmd:181-187`; aifi `36-R:20` `map2_dfr`).
- **Let the package loop.** Write one formula and call `processAssays(pb, form)` per cell type (dreamlet, `Nathan_NatImm_2021.Rmd:153-159`; variancePartition).
- **Count in `summarise()`.** `nDE = sum(adj.P.Val < 0.05)` (`Nathan_NatImm_2021.Rmd:202-208`; single_cell_aging `05/02:32`; aifi `36-R:8`).
- **Fit a model per group inside the pipe.** `group_modify(~ tidy(cor.test(.x$a, .x$b)))` (`PsychAD_r0/intersection_QTL.Rmd:29-31`); `group_modify(~ freq_glm(.x, formula = form))` (ALTRA `RA_06...Rmd:201-204`); `reframe()` (MM `frequency-volcanos.R:28-38`).
- **Keep groups above a size in two lines.** `tab <- table(sce$cell_type) > 1000`, then `keep <- names(tab)[tab]` (`variancePartition_analysis/yazar_2022.Rmd:64-65`; `write_pb_for_qtl.Rmd:50-52`).
- **Use a named vector as both lookup and loop set.** `gp_letter <- c(GP3 = "3J", ...)`, then `for (gp in names(gp_letter))` (immgenT `Figure3.R:348-349`; ALTRA, MM).
- **Group means in one expression.** `rowsum(L, labels)`, then `sweep(..., "/")` (immgenT `code/R/centered_mean_heatmap.R:68-72`).
- **Select columns by pattern.** `starts_with("NMF")` (IHA) and `select(starts_with("subject"), ...) |> distinct()` (sound-life `15-R:12`; single_cell_aging).
- **Align by names, not joins.** `meta[rownames(L), ]` (immgenT `02_compute_auc.R:52`). `keep <- intersect(colnames(counts), meta$barcode)`, then subset both (`Kfoury_..._get_data.Rmd:115-117`).
- **Grids.** `expand.grid()` with `with()` (dreamlet `theory_sims.Rmd:218-221`); `expand.grid()` with `anti_join()` to complete pairs (IFN `IFN_Fig1.ipynb` cell 25).

### Python

- **One boolean mask, one subset.** `mask = adata.var.mt | adata.var.ribo | adata.var_names.isin(extra)`, then `adata[:, ~mask]` (`skin_fibroblast_atlas/integration/run_scvi_fbs.py:60-70`; spatial_skin_atlas; cns-haematopoiesis).
- **Per-group thresholds without a loop.** `obs.groupby("sample")[m].transform(f)` (`cns-haematopoiesis/.../05_basic_qc_filtering_add_metadata.py:92-98`). Use `np.select()` for multi-way labels (cns; spatial `3_integrate_by_lineage.py:115-130`).
- **Composition in one chain.** `groupby("cluster")["tissue"].value_counts(normalize=True).unstack(fill_value=0)` or `pd.crosstab(a, b, normalize="index")` (skin, cns, spatial).
- **Named aggregation, with a boolean mean as a proportion.** `df.assign(is_x=df.niche == "X").groupby("id").agg(prop_x=("is_x", "mean"))` (spatial `1a_numbers_and_colors.ipynb` cell 135).
- **Map, then concat once.** `pd.concat(executor.map(f, items), ignore_index=True)` (aifi `02/04:11-12`; IHA; MM `01-predict_doublets.py:70-73`). For AnnData, `sc.concat(adatas, label="sample", keys=names)` (cns `helper_functions.py:122-128`).
- **Use the package getter.** `sc.pl.dotplot(..., return_fig=True).dot_size_df` (aifi `07`; sound-life `12a`). `sc.get.rank_genes_groups_df` replaces 25 hand-written lines in `skin_fibroblast_atlas/figures/supp_fig1_DEGs.ipynb` cell 8.
- **Config as a dict keyed by the loop variable.** `contrast_sets = {"CMV": {"fg": "Positive", "bg": "Negative"}}`, then a `.items()` loop (sound-life `04-file-sets/...deg_results.ipynb:12`).

## Larger patterns

1. **Three layers.** Numbered compute steps write `.rds`, `.h5ad` or parquet. Each figure or table gets one script or notebook that only reads them. Helpers live in a file only when two or more steps use them. immgenT has `code/pipeline/NN_*.R`, `script/FigureN.R` and `code/R/<topic>.R`, with the helper rule stated at `plot_utils.R:1-3`. MM and IHA have the same shape.
2. **A script reads top to bottom.** First paths, then `source()` with a comment naming the functions used. Then load, then one block per panel: compute, plot, save (immgenT `Figure3.R:53-56`, `Figure1.R:63-88`). In notebooks, markdown headings replace banner comments (aifi `07/33`).
3. **Code runs inline until it repeats.** A block used twice becomes a local function above its first use (immgenT `Figure4.R:140-162`; sound-life `10a:11` `process_adata()` once it loops). Functions take data and return one value, or a named list for several (`centered_mean_heatmap.R:71-74`). A writer returns `invisible(file)`.
4. **Repeated analyses are one loop over names, never a pasted block.** The corpus's worst lengths are copies: seven `mclapply` blocks (MOCHA `Fig3E_DiscriminativePower.R:40-99`), and six (aifi `36-R`). The short form sits in the same corpus.
5. **Parameters are plain variables in the first cell or at the top.** Formulas are objects passed in (`linFormula <- exp ~ Age + sex + (1|PTID)`, MOCHA `Fig5_chromVarModeling.R:49`). No repo uses YAML, papermill or Rmd `params:`.
6. **Source data comes from the plotted object.** `umap_df <- p1[[1]]$data`, then `fwrite(..., "Fig1b.csv")` (IFN `IFN_Fig1.ipynb` cells 13-14; immgenT `Figure4.R:530`; IHA).
7. **Use caching for heavy compute.** knitr `cache = TRUE` (`Nathan_NatImm_2021.Rmd:22`). Another form ends a chunk with `saveRDS` and starts the next with `readRDS` (GTestimate `umi_correction.Rmd:152-181`). sound-life and aifi skip finished items with `os.path.exists()`.

## Notebooks against scripts

- **Repo origin sets the split, not the kind of step.** Allen Institute repos are nearly all `.ipynb`: aifi 82, sound-life 82, IHA 286, ALTRA 163. MOCHA is 108 `.R`. Hoffman's repos are Rmd. immgenT computes in `.R` and uses Rmd only to display.
- **Where a split exists, it follows runtime.** Long or batch jobs are scripts run by `.sh`/SLURM wrappers (single_cell_aging, cns, MM `02-object-building`). Judgement steps are numbered notebooks: annotation, QC review and figures. MM keeps 17 `.py` files in jupytext percent format, so they open as notebooks.
- **Notebooks never call notebooks.** No repo uses `%run`, papermill or Rmd `params`. Notebooks reach shared code through `source()`, or `sys.path` plus import. GTestimate renders its Rmds in order from one driver, `notebooks/run_analysis.R:9-19`.
- **Hand-off is through files named after the producing step.** Examples are CSV (IHA), `.rds` (immgenT, reflex), h5ad or parquet (MM, cns `.obs` parquet), and a data platform by UUID (aifi, sound-life).
- **Numbering gives run order.** Names follow `NN-Language_task` when R and Python alternate (aifi `35-R`, then `36-R`). Variants are copies (`10a`-`10g`, ALTRA `05a/b/c`) where a parameter would do. The plugin does not copy that.
- **Bare `dim()`, `value_counts()` or `crosstab` cells are the audit trail in a notebook.** They do not carry over to scripts, where nobody reads the output.

## Where lucid-coder overrides the corpus

| Corpus habit | Plugin rule | Why |
|---|---|---|
| `%>%` in 16 of 17 repos | `\|>` | base R 4.1+, no dependency; user rule |
| Positional coefficients and columns (`coefficients[2:4, ]`, `dat[, 2]`) | select by name | IHA, MOCHA and MM all relabel by position |
| `$` partial matching (`resList[[CT]]$d` reads `dispObj`) | `[[ ]]` with the full name | works only by accident |
| `x[-grep()]` | `!grepl()` | empty result when nothing matches |
| Absolute paths, `setwd()`, sourcing a stale `_V1` path | root-relative paths, `here::here()` | every repo has them |
| `warnings.filterwarnings("ignore")`, `options(warn = -1)` | silence by category only | hides real problems |
| Helpers pasted into each notebook | helper file once 2+ use it | the copies drift (cns redefines `save_obs`) |
| Copied notebooks per subset (`05a/b/c`) | one parameter, one loop | copies drift (sound-life `05` cell 27 lost a filter) |
| Each plot pasted twice for `pdf()` | `p <- ...; ggsave()` | halves the file |
| `[not x for x in mask]` | `~mask` | vectorised |
| Bare `try/except: print()` | no bare except | errors swallowed |
| Printed `identical()` in scripts | one named `stopifnot()` | a print in a script is never read |

There is one open conflict. immgenT keeps `script/verify_*.R` files that re-derive published thresholds, while core bans `check_*.R` files. The plugin keeps the ban for now and lists this case here.
