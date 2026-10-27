# Analysis and plot catalogue: cc_learnings corpus vs ClawBio determinism

Written 2026-10-27. This file indexes 213 analysis rows and 174 plot rows. Each row has its own ID. The rows come from 17 repos, 16 counted as distinct because reflex v1 and v2 merge.
Every row lives in a part file under `parts/`. The part file carries the exact `path:line` citation.

## How to read a citation

Citations are relative to `$B=/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`.
Notebook line numbers count raw `.ipynb` JSON lines, not cell lines.

| Part | Repos | Abbreviations in its Source column |
|---|---|---|
| `parts/A_skin.md` | spatial_skin_atlas, skin_fibroblast_atlas | SA, FB = `repos_new/<repo>` |
| `parts/B_pbmc.md` | aifi-healthy-pbmc-reference, sound-life-scrna-analysis, IHA-Figure, MOCHA_Manuscript | AIFI, SL, IHA, MO = `repos/<repo>` |
| `parts/C_disease.md` | ALTRA-manuscript, IFN_manuscript, multiple-myeloma, single_cell_aging | full paths |
| `parts/D_reflex.md` | reflex_manuscript (+v2), immgenT-GP-analysis, cns-haematopoiesis | R1, R2, I, C, CF (legend on line 3) |
| `parts/E_methods.md` | dreamlet_analysis, variancePartition_analysis, GTestimate-Paper | full paths |
| `parts/F_clawbio_patterns.md` | ClawBio 0.3.0, Mycelium 0.7.2, mycelium-extra 0.16.1 | paths relative to each plugin root |

**Unverified rows.** These came from the summary or from grep only, so recheck them before you build on them:
- A-an-07, A-an-11, A-an-24
- B-an-04, B-an-10, B-an-38, B-an-43, B-an-55, B-pl-02, B-pl-10
- D-pl-22, plus the Fig2c-heatmap, ridge and lollipop rows in D

## Analysis families

"Repos" counts distinct repos out of 16. "Rows" lists part-table IDs. A merged row in a part file can cite up to 4 files.

| Family | Repos | Rows | Usual calls |
|---|---|---|---|
| QC and doublets | 15 | A-an-01..03; B-an-01,02,21,33,38; C-an-01..03,45; D-an-02,03,07,12,30,34..36; E-an-02,26..28 | `sc.pp.calculate_qc_metrics`, Scrublet, Seurat `subset()`, MAD outliers |
| Embed, cluster, integrate | 12 | A-an-04..06,08..10; B-an-06,07; C-an-05,06,09,10; D-an-04,05,37..40; E-an-24 | HVG → PCA → Harmony/scVI/CCA → neighbors → UMAP → Leiden |
| Enrichment | 11 | A-an-26,28; B-an-22,23,50; C-an-25..28,38,39; D-an-15,17,27; E-an-08 | fgsea, GSVA, Enrichr/ORA, decoupler, zenith |
| Annotation and reference mapping | 10 | A-an-07,33,34; B-an-03..05,08,10; C-an-07,08; D-an-01,06,41..43 | CellTypist (AIFI L1/L2/L3), scANVI/scArches, FAISS kNN |
| Composition / differential abundance | 10 | A-an-16,17; B-an-15,16,18,27; C-an-20,22..24; D-an-09; E-an-10 | CLR + Wilcoxon/lm/lmer, scCODA, milopy, crumblr |
| Pseudobulk DE | 9 | A-an-13; B-an-11..14; C-an-12..14,16; E-an-01,04,05,15 | DESeq2 / PyDESeq2, dreamlet, lme4 |
| Signatures and gene programs | 8 | A-an-14,15; B-an-24,28; C-an-29,30; D-an-16,18,20..25,31,32 | `score_genes`, UCell, NMF/starCAT, flashier EBMF |
| Marker DE | 7 | A-an-12; B-an-09; C-an-19; D-an-14 | `rank_genes_groups`, `FindAllMarkers` |
| Clinical / association tests | 6 | A-an-18,36; B-an-17,19,25,34..36; C-an-37,42..44,46..50 | Wilcoxon, Spearman, Fisher combination, lm |
| Single-cell and proteomic DE | 5 | B-an-20; C-an-17,18,21; D-an-13,26 | MAST, Seurat DE, Olink NPX lm |
| Trajectory / pseudotime | 5 | A-an-11,29..32; B-an-26,29,52; C-an-35 | scVelo, CellRank, PAGA, Monocle3, Slingshot |
| Variance partitioning | 4 | C-an-41; D-an-28,29; E-an-06,07,16 | `fitVarPart`, `fitExtractVarPartModel`, ANOVA η² |
| Benchmark / robustness | 4 | B-an-43..47; E-an-11..14,17,18,20..23,25,29 | simulation, downsampling, ARI, leave-one-out |
| TF activity / GRN | 4 | A-an-27; B-an-30,51; D-an-33 | CollecTRI + ULM, chromVAR, motif enrichment |
| Normalisation | 4 | C-an-04,34; E-an-03,19 | SCTransform, DSB, voom, GTestimate |
| Repertoire, flow, other assays | 4 | B-an-32,37; C-an-11,32,36; D-an-08,10,11,19 | scRepertoire, FlowSOM/scyan, MOFA2, infercnvpy |
| Cell-cell communication | 3 | A-an-25; B-an-31; C-an-31 | LIANA, CellPhoneDB, LARIS |
| Spatial | 2 | A-an-19..24,35 | squidpy nhood/co-occurrence, NicheCompass |
| Chromatin (scATAC) | 2 | B-an-39..42,48,49,53..55; C-an-33 | MOCHA, ArchR, Signac, MACS2 |
| Cross-cell-type sharing | 2 | C-an-15,40; E-an-09 | mashr, Spearman of logFC |

These analyses are absent from the corpus: survival/KM/Cox, SCENIC, forest-plot meta-analysis and ROC/FDR-calibration benchmarks. `parts/C_disease.md` and `parts/E_methods.md` both report this in their "Recurring defaults" sections.

## Plot families

| Family | Repos | Rows | Usual calls |
|---|---|---|---|
| Box / violin / swarm (+ stats) | 15 | A-pl-10,12,13,16,25; B-pl-04,09,12,13,28,32; C-pl-05..07; D-pl-06,14,22,25,30; E-pl-15,16,20,35 | `geom_boxplot(outlier.shape=NA)` + jitter, `stat_compare_means`, seaborn |
| Embedding by label | 14 | A-pl-01; B-pl-01,10,40; C-pl-01; D-pl-01,04,18,27; E-pl-12,33,34 | `sc.pl.umap`, `DimPlot`, scattermore |
| Heatmap | 12 | A-pl-03,06,08,23; B-pl-05,14,15,24,25,39; C-pl-03,14; D-pl-09,10,28,34..36,40; E-pl-08,09,13,17,31,36,40 | ComplexHeatmap, pheatmap, `geom_tile`, seaborn |
| Scatter + fit / labels | 12 | A-pl-15; B-pl-22,35,38,44; C-pl-09,16,18; D-pl-13,26,31,32; E-pl-18,19,32 | `geom_smooth(method="lm")`, ggrepel |
| Line / time course / benchmark curve | 10 | A-pl-07,24; B-pl-19,23,34; C-pl-10,23; D-pl-39; E-pl-21..26,28..30,38 | `geom_line` + `geom_ribbon` |
| Marker dotplot | 10 | A-pl-02,11; B-pl-03; C-pl-02; D-pl-07,19 | `sc.pl.dotplot`, Seurat `DotPlot` |
| Composition bar | 9 | A-pl-05; B-pl-08; C-pl-08; D-pl-08,20,29,38; E-pl-14 | `geom_bar(position="fill")`, pandas stacked bar |
| Volcano | 7 | B-pl-11,37; C-pl-04; D-pl-23; E-pl-05,06 | `plotVolcano`, ggplot + ggrepel |
| Enrichment dot / bar | 7 | A-pl-14,21; B-pl-36; C-pl-12,20; D-pl-37 | NES dotplot, `geom_col` |
| Ridge / density / histogram | 6 | B-pl-20,41; C-pl-17; D-pl-15,17,24 | ggridges, `geom_histogram` |
| Embedding feature / density | 5 | A-pl-17; B-pl-02; D-pl-02,03,05,41 | `FeaturePlot`, `sc.pl.embedding_density` |
| Graph / trajectory | 5 | A-pl-09,18..20; B-pl-21,27; C-pl-21; D-pl-33 | PAGA, milopy nhood graph, velocity stream |
| Alluvial / Sankey | 4 | B-pl-18; C-pl-11; D-pl-11; E-pl-37 | ggalluvial, `clonalCompare` |
| Set overlap | 3 | B-pl-31; C-pl-13; D-pl-16 | UpSetR, eulerr |
| Variance fractions | 3 | C-pl-15; E-pl-01..03,27 | `plotVarPart`, `plotPercentBars` |
| Forest / effect | 2 | C-pl-19; E-pl-07,10,11 | `plotForest`, scCODA effect plot |
| Spatial map | 2 | A-pl-04,22 | `sq.pl.spatial_scatter` |
| Other forms | 7 | B-pl-06,07,16,17,26,30,33,43; D-pl-12,21; E-pl-04,39 | sunburst, treemap, circos, ternary, pie, genome track, lollipop |
| Layout / export | 4 | B-pl-29,42; C-pl-22; E-pl-41 | patchwork, cowplot, mm-sized export |

## Where repos disagree

The same analysis gets different defaults in different repos. A deterministic skill has to fix each of these choices once.

- **DE call thresholds differ:**
  - IHA uses `padj < 0.05` and `abs(log2FC) > 0.1` (B-an-13).
  - multiple-myeloma uses `Qvalue < 0.1` and `|log2FC| > 0.5`. Its source is `repos/multiple-myeloma/analyses-rna/04-pseudobulk/l3-pseudobulk/13-plots-deseq2-volcano_fixed.ipynb:79-80`.
  - dreamlet uses FDR < 0.05 (E-an-04).
  - The skin atlases set no cutoff (A-an-12, A-an-13).
- **Pseudobulk gene filters differ:**
  - B-an-14 keeps genes detected in ≥ 10% of cells.
  - A-an-13 uses PyDESeq2 with `min_cells=20, min_counts=1000`.
  - E-an-03 uses dreamlet `min.cells`/`min.counts`, from 5/10 down to 2.
- **Seeds:** 3030 (AIFI, ALTRA, IHA), 42, 1, 123, 2022 and 6173 all appear. Most training and embedding steps are unseeded:
  - scVI and NicheCompass are unseeded (`parts/A_skin.md` Recurring defaults).
  - The reflex notebooks set no seed (`parts/D_reflex.md`).
  - The dreamlet scripts set no seed (`parts/E_methods.md`).
- **Graph and embedding parameters differ:**
  - `n_neighbors` takes 20, 25 and 50.
  - UMAP `min_dist` takes 0.05, 0.1 and 0.45.
  - Leiden resolution runs from 0.1 to 3.
- **Figure size and export differ:**
  - `theme_classic(base_size=6)` with mm sizes (IHA).
  - 7×7 `pdf()` defaults (MOCHA).
  - `ggsave(dpi=300)` (immgenT).
  - `dpi_save=300, pdf.fonttype=42` (skin).
  - A hard-coded `~/www/test.pdf` (dreamlet).
- **Paths:** every repo hard-codes absolute lab paths. Only cns-haematopoiesis uses argparse.

## Which ClawBio mechanisms fix which family

`parts/F_clawbio_patterns.md` ranks 8 gaps between ClawBio and Mycelium. The root gap is F-01 + F-28. ClawBio has one tested, parameterised script per analysis, and the agent may not re-implement it. Mycelium fixes the process (grill, gate, verify) but leaves code, thresholds, seeds and plot specs to the agent.

The top families above are the ones to build first. They appear in ≥ 9 of 16 repos, and they use standard tools.

| Rank | Candidate skill or convention | Covers families | ClawBio mechanisms to port |
|---|---|---|---|
| 1 | Shared plot-spec module (R `theme.R` + Python `theme.py`) for the 7 common forms: box/violin+stats, embedding, heatmap, scatter+fit, dotplot, composition bar, volcano | plot families with ≥ 7 repos | F-10, F-11, F-13, F-14 (fixed size, dpi, palette, drawn thresholds) |
| 2 | Pseudobulk DE script: aggregate → filter → DESeq2 or dreamlet → `de_results.csv` + volcano + heatmap | Pseudobulk DE, Volcano, Heatmap | F-01, F-05..07 (input contract), F-23 (named thresholds), F-08/09 (output tree) |
| 3 | scRNA QC + embed recipe script with seeded Scrublet, neighbors, UMAP and Leiden | QC and doublets, Embed/cluster, Embedding plots | F-01, F-16..20 (seed, commands.sh, env, checksums) |
| 4 | Composition / DA script: CLR + mixed model or crumblr → composition bar + effect plot | Composition / DA, Composition bar, Forest | F-01, F-23, F-03/04 (demo data + golden test) |
| 5 | Enrichment script: preranked fgsea or zenith on a DE table → NES dotplot | Enrichment, Enrichment dot | F-01, F-05..07, F-28 |

Ranks 2–5 take the DE table, counts and metadata layout of `clawbio:rnaseq-de` as their template. That skill is the most complete ClawBio example the F agent read (`parts/F_clawbio_patterns.md`, section 1).

Skipped: skills for families in ≤ 4 repos (spatial, scATAC, CCC, repertoire). Add one when a new project needs it twice.
