---
name: sc-annotate
description: >-
  Cell-type annotation of a single-cell query h5ad by reference label transfer:
  CellTypist (a given model, or one trained on a labelled reference h5ad) or
  scanpy ingest kNN transfer, then a marker dotplot to check the labels. Writes
  per-cell labels with a confidence, a label-by-cluster table, four fixed
  figures, an annotated h5ad, a report and a reproducibility bundle. Fire on
  "annotate cell types", "label transfer from a reference", "run CellTypist",
  "train a CellTypist model on my reference", "scanpy ingest labels",
  "map query cells to reference labels". Do NOT fire for marker-only manual
  annotation without a reference or model, scANVI/scArches/Seurat reference
  mapping, clustering or UMAP of the query, or label hierarchy cleanup.
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - single-cell
    - annotation
    - label transfer
    - CellTypist
    - scanpy ingest
  inputs:
    - name: query
      type: file
      format:
        - h5ad
      description: cells to label; raw integer counts in X, obsm X_umap, a cluster column
      required: true
    - name: markers
      type: file
      format:
        - tsv
      description: columns label and gene; drives the marker dotplot
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: parameters, summary and captioned figures
    - name: result
      type: file
      format:
        - json
      description: counts per label, median confidence, params, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - anndata
      - numpy
      - pandas
      - scipy
      - matplotlib
      - scanpy
      - celltypist
      - scikit-learn<1.8
  demo_data:
    - path: examples/demo_query.h5ad
      description: synthetic, 400 cells of 4 known types, 8 clusters, 158 genes
    - path: examples/demo_reference.h5ad
      description: synthetic, 300 labelled cells of the same 4 types
    - path: examples/demo_markers.tsv
      description: 2 marker genes per type
  endpoints:
    cli: python scripts/sc_annotate.py --query {query} --markers {markers} ... --outdir {outdir}
  trigger_keywords:
    - cell-type annotation
    - label transfer
    - CellTypist
    - reference mapping
    - scanpy ingest
---

# Cell-type annotation

You are the dispatcher for **sc-annotate**. The script does the analysis; you
collect its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/sc_annotate.py`. Never re-implement the transfer, change a
setting, or pick a value the user has not given. If a required argument is
unknown, ask the user (these are grill questions), or run `--demo` only when the
user asks for a demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--query` | `.h5ad` | raw counts (non-negative integers) in `X`; `obsm["X_umap"]`; the `--over-clustering` obs column | `examples/demo_query.h5ad` |
| `--reference-adata` | `.h5ad` | raw counts in `X`; the `--reference-label` obs column with at least 2 labels | `examples/demo_reference.h5ad` |
| `--model` | `.pkl` | a CellTypist model saved by `model.write` or downloaded beforehand | `data/celltypist_model.pkl` of a demo run |
| `--markers` | `.tsv` | exactly the columns `label` and `gene`, one row per marker | `examples/demo_markers.tsv` |

Gene names must match between query and reference or model (symbols with
symbols). The UMAP is the query's own; the script does not compute or project one.

## Parameters (all required, no defaults)

Each row is a grill question. Ask the user; never fill one in.

| Argument | Grill question | Demo-only value |
|---|---|---|
| `--query` | Which h5ad holds the cells to label? | `examples/demo_query.h5ad` |
| `--markers` | Which marker genes check each label (TSV of label, gene)? | `examples/demo_markers.tsv` |
| `--over-clustering` | Which obs column holds the query clusters (used for voting and the label-by-cluster table)? | `cluster` |
| `--engine` | CellTypist (`celltypist`) or scanpy ingest kNN (`ingest`)? | `celltypist` |
| `--seed` | Which seed? | `0` |
| `--reference` | Which genome build and annotation release named the query genes (e.g. `GRCh38, GENCODE 44`)? Recorded in `result.json`. | `synthetic demo genes, no genome build` |
| `--model` | celltypist only: which saved CellTypist model? Give this or `--reference-adata`. | absent |
| `--reference-adata` | Which labelled reference h5ad? (celltypist: trains a model; ingest: required) | `examples/demo_reference.h5ad` |
| `--reference-label` | Which reference obs column holds the labels? | `cell_type` |
| `--majority-voting` | celltypist only: refine labels by majority vote within each cluster (`yes`/`no`)? | `yes` |
| `--n-pcs` | ingest only: how many reference PCs? | `10` |
| `--n-neighbors` | ingest only: how many reference neighbours vote? | `15` |
| `--outdir` | Which new directory gets the outputs? It must not exist. | `sc_annotate_demo` |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them.
With `--demo`, only `--engine` and `--outdir` may also be given; `--demo
--engine ingest` uses the ingest demo values. An argument of the other engine is
an error, so a run never silently ignores one.

## Workflow (fixed order)

1. **Validate** arguments and inputs; bad input gives `ERROR: <message>`, exit 2.
   Gene names and cell ids must be unique in the query and the reference
   (`skillkit.check_unique`).
2. **Normalise** copies of query and reference: `normalize_total(1e4)`, then
   `log1p`. `X` in the saved h5ad stays raw counts.
3. **Transfer** labels with `--engine`:
   - `celltypist`: load `--model`, or train on `--reference-adata` by
     `--reference-label` (`celltypist.train`, CellTypist's default solver,
     `n_jobs=1`, `random_state=--seed`). Then `celltypist.annotate` in
     best-match mode, with majority voting over `--over-clustering` when
     `--majority-voting yes`.
   - `ingest`: keep shared genes, then PCA (`--n-pcs`) and a kNN graph
     (`--n-neighbors`) on the reference, both with `random_state=--seed`. Then
     scanpy `Ingest` projects the query into the reference PCA and takes the
     kNN mode label. Ties go to the label first in sorted order.
4. **Confidence** per cell:
   - celltypist: the model probability of the assigned label (after voting, a
     cell outvoted by its cluster can score near 0);
   - ingest: the fraction of the k reference neighbours carrying the label.
5. **Check markers**: every `--markers` label must be a model or reference
   label, and every gene must be in the query.
6. **Tabulate** labels per `--over-clustering` cluster.
7. **Draw** the four figures and write the report.

## Output structure

`tables/labels.tsv` is sorted by cell; `tables/label_by_cluster.tsv` by cluster
(obs category order), then label, and lists every cluster × label pair,
zeros included. `n_matched` of `n_in_resource` is the gene overlap: query genes
found among the model's features or the reference's genes (`n_shared_genes`
repeats `n_matched`).

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # n_cells, n_labels, n_shared_genes, n_matched, n_in_resource, label_counts, median_confidence, params, engine, seed, reference, chat_summary_lines, preferred_artifacts
├── data/
│   ├── annotated.h5ad         # the query, raw X, plus obs transferred_label and transferred_confidence
│   └── celltypist_model.pkl   # only when celltypist trained on --reference-adata
├── tables/
│   ├── labels.tsv             # cell, label, confidence
│   └── label_by_cluster.tsv   # cluster, label, n_cells, frac_of_cluster
├── figures/
│   ├── umap_label.pdf
│   ├── umap_label.png
│   ├── umap_label.csv
│   ├── marker_dotplot.pdf
│   ├── marker_dotplot.png
│   ├── marker_dotplot.csv
│   ├── confidence_box.pdf
│   ├── confidence_box.png
│   ├── confidence_box.csv
│   ├── label_cluster_heatmap.pdf
│   ├── label_cluster_heatmap.png
│   └── label_cluster_heatmap.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python and package versions
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

## Gotchas

- **The model will want to download a CellTypist model** (`models.download_models`,
  `Immune_All_Low.pkl`) when none is given. Do not. Ask the user for a model
  path or a reference; the script never touches the network.
- **The model will want to turn majority voting on** because the corpus does.
  Do not choose for the user. Voting overwrites each cell's own prediction with
  its cluster's majority, so it depends on the clustering; ask.
- **The model will want to call confidence a probability that the label is
  right.** Do not. A score is not a probability: CellTypist's value is a model
  output from one classifier, and ingest's is a neighbour fraction. Neither is
  calibrated against truth; report it as "confidence".
- **The model will want to read a marker the query lacks in a cell as absence.**
  Do not. Non-detection is not absence: a zero in a sparse cell is often
  dropout, so the dot plot's fraction expressing is a detection rate.
- **The model will want to guess `--reference`** from the gene names. Do not.
  Ask for the build and annotation release; it is recorded, never inferred.
  `--reference-adata` is the labelled h5ad, a different argument.
- **The model will want to re-label low-confidence cells** as "Unknown" with a
  cutoff it picked. Do not. The script reports confidence; any cutoff is the
  user's, applied downstream.
- **The model will want to drop a marker gene or label** that fails the marker
  check. Do not. The error names the absent genes or labels; ask the user to fix
  the TSV or the gene names.
- **The model will want to feed log-normalised data** to fix a `must hold raw
  counts` error. Do not. Point the user to the raw counts; the script
  normalises them itself.
- **The model will want to compute a UMAP or clusters** for a query that lacks
  them. Do not, inside this skill. That is an upstream step with its own
  parameters.

## CLI reference

```bash
python scripts/sc_annotate.py \
  --query data/query.h5ad --markers config/pbmc_markers.tsv \
  --over-clustering leiden --engine celltypist \
  --model models/ref_pbmc_L2.pkl --majority-voting yes \
  --seed 0 --reference 'GRCh38, GENCODE 44' --outdir results/2026-10-27_annotate

python scripts/sc_annotate.py \
  --query data/query.h5ad --markers config/pbmc_markers.tsv \
  --over-clustering leiden --engine ingest \
  --reference-adata data/reference.h5ad --reference-label cell_type \
  --n-pcs 30 --n-neighbors 15 --seed 0 --reference 'GRCh38, GENCODE 44' --outdir results/2026-10-27_ingest

# Demo (only when the user asks for one)
python scripts/sc_annotate.py --demo --outdir demo_out
```

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
