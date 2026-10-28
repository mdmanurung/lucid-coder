# Pipelines that can become skills

Written 2026-10-27. Each pipeline below is an ordered chain of families from `CATALOGUE.md`.

- **Repos:** the number of the 16 corpus repos that run at least one of its steps.
- **Rows:** the catalogue IDs. The `path:line` sources are in `parts/`.
- **Plot sources:** a plot marked with row IDs exists in the corpus. A plot marked **std** is standard practice but missing from the corpus.
  - Single-cell plots: Heumos et al. 2023, *Nat Rev Genet* 24:550 (sc-best-practices.org).
  - Bulk and pseudobulk DE plots: the DESeq2 and limma vignettes.
- **Plots module:** forms marked `plots.*` already exist in `lucid-coder/templates/plots.{R,py}`.
- **Parameters:** every cutoff a skill needs is a required argument, settled by a grill question (`grill/references/analysis-decisions.md`, "Plot parameters" row). No skill borrows the corpus values.

## Build first (≥ 9 repos)

### 1. scRNA preprocessing: QC → doublets → normalise → embed → cluster

- **Repos:** 15
- **Rows:** A-an-01..05,10; B-an-01,02,06,07; C-an-01,02,05,06,10; D-an-02,03,34..38,40
- **Steps:**
  - QC metrics, then thresholds (a grill question);
  - Scrublet with a seed;
  - filter;
  - `normalize_total` + `log1p`;
  - HVG;
  - PCA;
  - optional Harmony or scVI by batch;
  - neighbours, UMAP and Leiden, each with a seed.
- **ClawBio template:** `scrna-orchestrator` and `scrna-embedding`.

| Plot | Source |
|---|---|
| QC violins: genes, counts, % mito | B-pl-04, D-pl-06, D-pl-22, A-pl-13 (`plots.*` box) |
| Counts vs genes scatter, coloured by % mito, with threshold lines | **std** |
| Doublet-score histogram with the cut-off | **std** |
| Cells kept per sample, before and after QC | E-pl-14 (`plots.*` composition) |
| PCA variance (elbow / scree) | C-pl-16 |
| UMAP by batch and by cluster | A-pl-01, B-pl-01, D-pl-18 (`plots.*` embedding) |
| Integration metrics (scIB) | D-an-39 |

### 2. Cell-type annotation: reference mapping → marker check → label hierarchy

- **Repos:** 10
- **Rows:** A-an-07,34; B-an-03..05,08,10; C-an-07,08; D-an-01,06,41..43
- **Steps:**
  - CellTypist, scANVI or kNN transfer;
  - majority voting;
  - a marker dotplot to confirm labels;
  - L1/L2/L3 cleanup.

| Plot | Source |
|---|---|
| UMAP by label | B-pl-01, C-pl-01 (`plots.*` embedding) |
| Marker dotplot | A-pl-02, B-pl-03, C-pl-02, D-pl-07, D-pl-19 (`plots.*` dotplot) |
| Label vs cluster flow | B-pl-18, D-pl-11 (alluvial) |
| Prediction-confidence violin per label | **std** |
| Label × cluster confusion heatmap | **std** (`plots.*` heatmap) |

### 3. Pseudobulk differential expression

- **Repos:** 9
- **Rows:** A-an-13; B-an-11..14; C-an-12..14,16; E-an-01..05,15
- **Steps:**
  - sum counts per sample × cell type;
  - drop small groups and low-detection genes (a grill question);
  - DESeq2 or dreamlet with the design (grill: covariates and the replication unit);
  - write `de_results.tsv`.
- **ClawBio template:** `rnaseq-de`, which is the most complete ClawBio example.

| Plot | Source |
|---|---|
| Sample PCA before testing, by condition and batch | **std** |
| Raw p-value histogram | **std** |
| MA plot | **std** |
| Mean-variance trend | E-pl-04 (voom) |
| Volcano with drawn cut-offs | B-pl-11, C-pl-04, E-pl-05, E-pl-06 (`plots.*` volcano) |
| Top-gene heatmap, gene × cell type | E-pl-08, E-pl-17 (`plots.*` heatmap) |
| DE genes per cell type | E-pl-10 |
| One gene across cell types (forest) | E-pl-07 |
| Expression of top genes per donor | B-pl-12 (`plots.*` box) |

### 4. Composition / differential abundance

- **Repos:** 10
- **Rows:** A-an-16,17; B-an-15,16,18,27; C-an-20,22..24; D-an-09; E-an-10
- **Steps:**
  - count cells per sample × label;
  - CLR or crumblr;
  - a mixed model or Wilcoxon on the donor unit;
  - BH correction.
- **Optional route:** Milo (milopy / miloR) at the neighbourhood level.

| Plot | Source |
|---|---|
| Stacked composition bar per sample | A-pl-05, C-pl-08, D-pl-08 (`plots.*` composition) |
| Per-donor proportion by group | E-pl-15 (`plots.*` box, `pair` for repeated visits) |
| Effect size with interval per cell type | C-pl-19 (forest) |
| Milo neighbourhood graph and beeswarm | A-pl-09, B-pl-27 |

### 5. Enrichment on a DE result

- **Repos:** 11
- **Rows:** A-an-26,28; B-an-22,23,50; C-an-25..28; D-an-15,17,27; E-an-08
- **Steps:**
  - rank the genes or take a gene set (grill question: rank statistic and gene-set source and version);
  - fgsea, zenith, ORA, or decoupler ULM/MLM.
- **ClawBio template:** `pathway-enricher`.

| Plot | Source |
|---|---|
| NES dotplot | C-pl-12, D-pl-37 |
| Pathway × cell type NES heatmap with stars | D-pl-10, E-pl-09 (`plots.*` heatmap) |
| Top-pathway bar | B-pl-36 |
| Running enrichment-score curve for a named pathway | **std** |

## Second wave (4–8 repos)

| Pipeline | Repos | Rows | Plots in corpus | Missing standard plot |
|---|---|---|---|---|
| Signature scoring and gene programs (`score_genes`, UCell, NMF, flashier) | 8 | A-an-14; B-an-24,28; C-an-29,30; D-an-16,20..25 | score violin (A-pl-10), program heatmap (D-pl-28), score on UMAP (B-pl-02) | score vs random-gene-set control |
| Group-comparison and association tests (scores, frequencies, Olink) | 6 | B-an-17,19,25; C-an-46,47,49 | box + stats (B-pl-12, C-pl-06), scatter + fit (C-pl-09) | the paired-line form (`plots.*` box `pair`) |
| Trajectory (PAGA, scVelo, CellRank, Monocle3, Slingshot) | 5 | A-an-11,29..32; B-an-26,29; C-an-35 | PAGA graph (A-pl-18, C-pl-21), velocity stream (A-pl-19), pseudotime ridge (B-pl-20) | gene trend along pseudotime with a band |
| Variance partitioning | 4 | C-an-41; D-an-29; E-an-06,07,16 | fraction violin (E-pl-01, C-pl-15), percent bars (E-pl-02) | none |
| TF / pathway activity (CollecTRI, PROGENy, chromVAR) | 4 | A-an-27,28; B-an-30,51 | TF violin (A-pl-25) | activity heatmap (`plots.*` heatmap) |

## Leave as reference skills (≤ 3 repos)

These pipelines appear in too few repos to fix their defaults from the corpus. They are CCC (A-an-25, B-an-31, C-an-31), spatial (A-an-19..24), scATAC (B-an-39..55, C-an-33), repertoire (D-an-08..11) and benchmarks (E-an-11..23). Installed reference skills already cover the APIs: `liana-ccc`, `cellrank`, `milo-r`, `crumblr-compositional`, `dreamlet` and `decoupler`.

Those skills document how to call each tool. They do not fix the run, so they are not deterministic in the ClawBio sense. Add a pipeline skill when a second project needs one.
