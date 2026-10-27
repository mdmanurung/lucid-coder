---
name: sc-pipelines
description: >-
  Router to tested, deterministic scripts for single-cell and bulk-omics analyses
  whose exact numeric answer depends on parameter and formula choices an ad-hoc
  script is likely to get slightly wrong: QC and doublet filtering, normalisation,
  embedding, clustering, cell-type annotation, marker genes (rank_genes_groups),
  pseudobulk differential expression, cell-type composition and differential
  abundance, group-comparison or correlation tests (Wilcoxon, Spearman, Pearson),
  longitudinal mixed models, gene-set enrichment, gene-program discovery (consensus
  NMF, cNMF), gene-signature scoring, TF or pathway activity, cell-cell
  communication, trajectory and pseudotime, RNA velocity, MOFA multi-omics factor
  analysis, mashr effect sharing, and variance partitioning. Fire on any of those
  analyses, even when the task already spells out its own step-by-step method --
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
| `composition` | `composition.py` | Cell-type composition and differential abundance with crumblr. |
| `enrichment` | `enrichment.py` | Gene-set enrichment on one or more DE tables (contrasts or cell types). |
| `gene-programs` | `gene_programs.py` | Gene-program discovery by consensus NMF (cNMF-style) on raw single-cell counts. |
| `group-tests` | `group_tests.py` | Group-comparison and association tests on sample-level values (cell-type frequencies or CLR values, signature scores, Olink NPX). |
| `longitudinal` | `longitudinal.py` | Longitudinal mixed models per feature on repeated sample-level values (CLR cell-type frequencies, pseudobulk vst or log CPM, signature scores, Olink NPX). |
| `mashr-sharing` | `mashr_sharing.py` | Cross-cell-type effect sharing with mashr. |
| `milo-da` | `milo_da.py` | Milo differential abundance on a kNN graph with pertpy (pt.tl.Milo, pydeseq2 solver). |
| `mofa` | `mofa.py` | Multi-omics factor analysis with MOFA2 (Python mofapy2) on per-view sample tables (pseudobulk RNA, ADT, Olink, chromVAR TF scores). |
| `pseudobulk-de` | `pseudobulk_de.py` | Pseudobulk differential expression for one cell type from a single-cell h5ad. |
| `sc-annotate` | `sc_annotate.py` | Cell-type annotation of a single-cell query h5ad by reference label transfer. |
| `sc-de` | `sc_de.py` | Per-feature differential expression in two modes. |
| `sc-embed` | `sc_embed.py` | Embed and cluster a QC'd single-cell h5ad. |
| `sc-markers` | `sc_markers.py` | Marker genes per cluster or cell-type label from a log-normalised single-cell h5ad. |
| `sc-normalize` | `sc_normalize.py` | Normalise a single-cell count h5ad one of three fixed ways. |
| `sc-qc` | `sc_qc.py` | scRNA-seq cell QC and doublet removal from a raw-count h5ad with a sample column. |
| `signature-scoring` | `signature_scoring.py` | Gene-signature scores per cell from a log-normalised single-cell h5ad (scanpy score_genes or decoupler AUCell), the mean per sample and group in a long table ready for group-tests, a size-matched random-gene-set control per signature, box, UMAP and control figures, a report and a reproducibility bundle. |
| `tf-activity` | `tf_activity.py` | TF or pathway activity from a prior-knowledge network (CollecTRI, PROGENy or any source-target-weight TSV) with decoupler 2.x ULM or MLM. |
| `trajectory` | `trajectory.py` | Trajectory and pseudotime from an embedded single-cell h5ad. |
| `variance-partition` | `variance_partition.py` | Variance partitioning of pseudobulk gene expression across sample covariates with R variancePartition. |
| `velocity` | `velocity.py` | RNA velocity and fate probabilities from a single-cell h5ad with spliced and unspliced layers. |
