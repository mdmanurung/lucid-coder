---
name: sc-pipelines
description: >-
  Router to tested, deterministic scripts for single-cell and bulk-omics analyses whose
  exact numeric answer depends on parameter and formula choices an ad-hoc script is likely
  to get slightly wrong: QC and doublet filtering, normalisation, embedding, clustering,
  cell-type annotation, marker genes (rank_genes_groups), pseudobulk differential
  expression, cell-type composition and differential abundance, group-comparison or
  correlation tests (Wilcoxon, Spearman, Pearson), longitudinal mixed models, gene-set
  enrichment, gene-program discovery (consensus NMF, cNMF), gene-signature scoring, TF or
  pathway activity, cell-cell communication, trajectory and pseudotime, RNA velocity, MOFA
  multi-omics factor analysis, mashr effect sharing, and variance partitioning. Fire on any
  of those analyses, even when the task already spells out its own step-by-step method --
  check whether a fixed script covers it before writing one by hand.
  (acumen bench 2026-10-27: this wording raised script use from 30/48 to 39/48 runs.)
---

# sc-pipelines

Each analysis below has a tested script that fixes every step. Do not write your own
analysis code when a pipeline covers the task.

1. Pick the pipeline from the table.
2. Read its `SKILL.md`, at `.claude/skills/sc-pipelines/<pipeline>/SKILL.md`.
3. Run its script as `python .claude/skills/sc-pipelines/<pipeline>/scripts/<script>`
   with every required argument. Take each value from the task; never invent one.
4. Read the answer from the script's `result.json` and `tables/`. Do not recompute it
   by hand, even for a question about one intermediate step.

| Pipeline | Script | Use it for |
|---|---|---|
| `ccc` | `ccc.py` | Ligand-receptor cell-cell communication between cell types from a log-normalised single-cell h5ad with LIANA+ rank_aggregate (CellPhoneDB, Connectome, log2FC, NATMI, SingleCellSignalR, or a subset), a bundled or exported ligand-receptor resource TSV, a fixed pass rule on magnitude and specificity ranks, a source x target count heatmap, an LR dot plot, a report and a reproducibility bundle. |
| `composition` | `composition.py` | Cell-type composition and differential abundance with crumblr: counts cells per sample and label (from a cell-level h5ad or a counts TSV), drops rare labels, runs crumblr's precision-weighted CLR and variancePartition dream (fixed and random effects, e.g. (1|donor)), and writes tables, a stacked composition bar, per-label proportion boxes, an effect forest, a report and a reproducibility bundle. |
| `enrichment` | `enrichment.py` | Gene-set enrichment on one or more DE tables (contrasts or cell types): preranked fgsea in R, or over-representation (one-sided Fisher / hypergeometric) of the significant genes against a stated universe, BH within each contrast, with an NES dot plot, an NES heatmap, running enrichment-score curves, a report and a reproducibility bundle. |
| `gene-programs` | `gene_programs.py` | Gene-program discovery by consensus NMF (cNMF-style) on raw single-cell counts: HVG selection, repeated seeded sklearn NMF, a local-density filter and k-means consensus of the replicate components, per-cell program usage, gene scores, top genes, mean usage, fraction of active cells and one-vs-rest AUC per group, a stability table, four fixed figures, a report and a reproducibility bundle. |
| `group-tests` | `group_tests.py` | Group-comparison and association tests on sample-level values (cell-type frequencies or CLR values, signature scores, Olink NPX): one Wilcoxon rank-sum, paired Wilcoxon signed-rank, Spearman or Pearson test per feature, BH or Bonferroni across features, an effect size for every test, a table, box or scatter figures for the top features, a report and a reproducibility bundle. |
| `longitudinal` | `longitudinal.py` | Longitudinal mixed models per feature on repeated sample-level values (CLR cell-type frequencies, pseudobulk vst or log CPM, signature scores, Olink NPX): one lmerTest `lmer(value ~ time + covariates + (1|subject))` per feature, the named coefficient tested with Satterthwaite df, BH or Bonferroni across features, singular and unconverged fits flagged per feature, subject trajectory, forest and p-value figures, a report and a reproducibility bundle. |
| `mashr-sharing` | `mashr_sharing.py` | Cross-cell-type effect sharing with mashr: takes one DE table per cell type (gene, logFC, SE; a pseudobulk-de de_results.tsv works as is), aligns genes, fits mashr with canonical or canonical plus data-driven (ED) covariances, and writes posterior means, posterior SDs and lfsr per gene and cell type, counts of genes significant in 1..K cell types, mashr pairwise sharing, raw-logFC Spearman per pair, three figures, a report and a reproducibility bundle. |
| `milo-da` | `milo_da.py` | Milo differential abundance on a kNN graph with pertpy (pt.tl.Milo, pydeseq2 solver): seeded index-cell sampling, counts per sample per neighbourhood, NB GLM per neighbourhood, SpatialFDR, annotation with a Mixed label. |
| `mofa` | `mofa.py` | Multi-omics factor analysis with MOFA2 (Python mofapy2) on per-view sample tables (pseudobulk RNA, ADT, Olink, chromVAR TF scores): a seeded CPU fit with a fixed number of factors, variance explained per factor and view, factor scores, feature weights, factor-metadata tests, a model file, figures, a report and a reproducibility bundle. |
| `pseudobulk-de` | `pseudobulk_de.py` | Pseudobulk differential expression for one cell type from a single-cell h5ad: sums raw counts per sample, filters samples and genes, runs DESeq2 (PyDESeq2 or R) or dreamlet (repeated measures, random effects) with a design formula and contrast, and writes tables, five fixed figures, a report and a reproducibility bundle. |
| `sc-annotate` | `sc_annotate.py` | Cell-type annotation of a single-cell query h5ad by reference label transfer: CellTypist (a given model, or one trained on a labelled reference h5ad) or scanpy ingest kNN transfer, then a marker dotplot to check the labels. |
| `sc-de` | `sc_de.py` | Per-feature differential expression in two modes: a MAST hurdle model on the log-normalised cells of one cell type (seeded per-donor downsampling, optional cellular detection rate covariate, optional (1|donor) random effect), or a per-feature linear model with covariates on sample-level values such as Olink NPX. |
| `sc-embed` | `sc_embed.py` | Embed and cluster a QC'd single-cell h5ad: normalise, pick highly variable genes, PCA, optional Harmony integration by batch, kNN graph, UMAP and igraph Leiden, every step seeded. |
| `sc-markers` | `sc_markers.py` | Marker genes per cluster or cell-type label from a log-normalised single-cell h5ad: scanpy rank_genes_groups (Wilcoxon or t-test, each group against the rest, Benjamini-Hochberg), a fixed pass rule on adjusted p, log fold change and fraction expressing, a top-marker dot plot and heatmap, a report and a reproducibility bundle. |
| `sc-normalize` | `sc_normalize.py` | Normalise a single-cell count h5ad one of three fixed ways: log1p of library-size-scaled counts, analytic Pearson residuals (the model SCTransform v2 approximates), or dsb for ADT counts using empty droplets and isotype controls. |
| `sc-qc` | `sc_qc.py` | scRNA-seq cell QC and doublet removal from a raw-count h5ad with a sample column: computes per-cell QC metrics (genes, counts, % mito, % ribo), scores doublets with Scrublet per sample, applies fixed or per-sample MAD thresholds, drops failed cells and predicted doublets, and writes a filtered h5ad, tables, four fixed figures, a report and a reproducibility bundle. |
| `signature-scoring` | `signature_scoring.py` | Gene-signature scores per cell from a log-normalised single-cell h5ad (scanpy score_genes or decoupler AUCell), the mean per sample and group in a long table ready for group-tests, a size-matched random-gene-set control per signature, box, UMAP and control figures, a report and a reproducibility bundle. |
| `tf-activity` | `tf_activity.py` | TF or pathway activity from a prior-knowledge network (CollecTRI, PROGENy or any source-target-weight TSV) with decoupler 2.x ULM or MLM. |
| `trajectory` | `trajectory.py` | Trajectory and pseudotime from an embedded single-cell h5ad: PAGA group connectivity and diffusion pseudotime (scanpy DPT) from a fixed root cell, with tables, four fixed figures, a report and a reproducibility bundle. |
| `variance-partition` | `variance_partition.py` | Variance partitioning of pseudobulk gene expression across sample covariates with R variancePartition: voomWithDreamWeights, then fitExtractVarPartModel with categorical covariates as random effects (1|x) and continuous ones as fixed, after a canCorPairs collinearity gate. |
| `velocity` | `velocity.py` | RNA velocity and fate probabilities from a single-cell h5ad with spliced and unspliced layers: scVelo velocity (stochastic or dynamical), velocity pseudotime, and CellRank 2 VelocityKernel + GPCCA fate probabilities towards named or inferred terminal states, with tables, four fixed figures, a report and a reproducibility bundle. |
