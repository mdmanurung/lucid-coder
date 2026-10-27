---
name: sc-embed
description: >-
  Embed and cluster a QC'd single-cell h5ad: normalise, pick highly variable
  genes, PCA, optional Harmony integration by batch, kNN graph, UMAP and
  igraph Leiden, every step seeded. Writes data/embedded.h5ad, cluster and
  batch-by-cluster tables, four fixed figures, a report and a
  reproducibility bundle. Fire on "embed and cluster", "normalise, HVG, PCA,
  UMAP, Leiden", "run harmony on batch", "integrate batches and cluster",
  "UMAP coloured by batch and cluster", "Leiden clustering of my h5ad".
  Do NOT fire for QC or doublet filtering (sc-qc), scVI/scANVI integration,
  cell-type annotation, marker genes, or differential expression.
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - single-cell
    - normalisation
    - highly variable genes
    - PCA
    - Harmony
    - UMAP
    - Leiden
  inputs:
    - name: adata
      type: file
      format:
        - h5ad
      description: QC'd cell-level AnnData with raw integer counts in X
      required: true
  outputs:
    - name: embedded
      type: file
      format:
        - h5ad
      description: all genes, log-normalised X, counts layer, PCA/Harmony/UMAP, leiden
    - name: report
      type: file
      format:
        - md
      description: params, summary and captioned figures
    - name: result
      type: file
      format:
        - json
      description: cluster sizes, params, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - scanpy
      - anndata
      - harmonypy
      - igraph
      - scikit-misc
      - numpy
      - pandas
      - scipy
      - matplotlib
  demo_data:
    - path: examples/demo.h5ad
      description: synthetic, 2 batches x 3 groups x 100 cells, 300 genes
  endpoints:
    cli: python scripts/sc_embed.py --adata {adata} ... --outdir {outdir}
  trigger_keywords:
    - embed and cluster
    - normalise HVG PCA
    - harmony integration
    - UMAP Leiden
    - batch correction
---

# sc-embed

You are the dispatcher for **sc-embed**. The script does the analysis; you
collect its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/sc_embed.py`. Never re-implement the analysis, change a
parameter, or pick a value the user has not given. If a required argument is
unknown, ask the user (these are grill questions), or run `--demo` only when the
user asks for a demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--adata` | `.h5ad` | QC'd cells, raw counts (non-negative integers) in `X` | sc-qc `data/qc_filtered.h5ad`, `examples/demo.h5ad` |
| batch column | obs column | `--batch-key`, no missing values; 2+ levels for Harmony | `batch` with `b1`, `b2` |

## Parameters (all required, no defaults)

Each row is a grill question. Ask the user; never fill one in.

| Argument | Grill question | Demo-only value |
|---|---|---|
| `--adata` | Which QC'd h5ad holds the cells? | `examples/demo.h5ad` |
| `--target-sum` | To what total count is each cell scaled before `log1p`? | `10000` |
| `--n-top-genes` | How many highly variable genes? | `100` |
| `--flavor` | Which HVG flavor: `seurat`, `cell_ranger` or `seurat_v3`? | `seurat` |
| `--batch-key` | Which obs column holds the batch (sample, donor, run), or `none`? | `batch` |
| `--n-comps` | How many principal components to compute? | `20` |
| `--integration` | Integrate by batch with Harmony (`harmony`) or not (`none`)? | `harmony` |
| `--n-neighbors` | How many neighbours in the kNN graph? | `15` |
| `--n-pcs` | How many PCs (from the elbow) feed the graph? | `10` |
| `--min-dist` | Which UMAP `min_dist`? | `0.3` |
| `--resolution` | Which Leiden resolution? | `0.3` |
| `--seed` | Which seed for PCA, Harmony, neighbours, UMAP and Leiden? | `0` |
| `--reference` | Which genome build and annotation release named the genes (e.g. `GRCh38, GENCODE 44`)? Recorded in `result.json`. | `synthetic demo genes, no genome build` |
| `--outdir` | Which new directory gets the outputs? It must not exist. | `sc_embed_demo` |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them.
With `--demo`, only `--integration` and `--outdir` may also be given; any other
analysis argument is an error.

## Workflow (fixed order)

1. **Validate** arguments and input; bad input gives `ERROR: <message>`, exit 2.
   X must hold raw counts (non-negative integers); gene names and cell ids must
   be unique (`skillkit.check_unique`). `--batch-key none` adds a constant obs
   column `none`.
2. **Keep raw counts** in `layers["counts"]`.
3. **HVG**: `--flavor seurat_v3` ranks genes on the raw counts before step 4;
   `seurat` and `cell_ranger` rank on the log-normalised data after step 4.
   `--n-top-genes` genes, `batch_key=--batch-key` unless `none`.
4. **Normalise**: `normalize_total(target_sum=--target-sum)`, then `log1p`.
   X holds the log-normalised values for all genes.
5. **PCA** on a copy subset to the HVGs, scaled to zero mean and unit variance
   with no clipping: `--n-comps` components, `svd_solver="arpack"`, seeded.
6. **Integrate**: `--integration harmony` runs `harmonypy.run_harmony` on the
   PCs by `--batch-key` (20 rounds at most, seeded) into `X_pca_harmony`;
   `none` keeps `X_pca`.
7. **Neighbours** on the first `--n-pcs` columns of that representation,
   `--n-neighbors` neighbours, seeded.
8. **UMAP** with `--min-dist`, seeded.
9. **Leiden**: `flavor="igraph"`, `n_iterations=2`, undirected, `--resolution`,
   seeded, into `obs["leiden"]`, levels in numeric order.

BLAS and numba run on one thread (`threadpoolctl`, `sc.settings.n_jobs = 1`).

## Output structure

`cluster_sizes.tsv` is sorted by cluster number; `batch_by_cluster.tsv` by
cluster number, then batch.

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # n_cells, n_genes, n_hvg, n_clusters, cluster_sizes, use_rep, integration, seed, params, reference, chat_summary_lines, preferred_artifacts
├── data/
│   └── embedded.h5ad          # all genes; X log-normalised, layers["counts"], var highly_variable, obsm X_pca (+ X_pca_harmony), X_umap, obs leiden
├── tables/
│   ├── cluster_sizes.tsv      # cluster, n_cells, frac
│   └── batch_by_cluster.tsv   # cluster, batch, n_cells, frac_of_cluster (zeros kept)
├── figures/
│   ├── pca_variance.pdf
│   ├── pca_variance.png
│   ├── pca_variance.csv
│   ├── umap_batch.pdf
│   ├── umap_batch.png
│   ├── umap_batch.csv
│   ├── umap_cluster.pdf
│   ├── umap_cluster.png
│   ├── umap_cluster.csv
│   ├── batch_composition.pdf
│   ├── batch_composition.png
│   └── batch_composition.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python and package versions
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

## Gotchas

- **The model will want to re-cluster in-session** (`sc.tl.leiden` at another
  resolution) to "refine" clusters. Do not. Rerun the script with the user's
  new `--resolution` into a new outdir, so each clustering has provenance.
- **The model will want to switch to Harmony** because the batch UMAP looks
  split. Do not. Integration is the user's call: a split may be biology, such
  as a disease-only cell state. Show both runs and ask.
- **The model will want to pick `--n-pcs` from the elbow itself.** Do not.
  Show `pca_variance.png` and ask the user; the elbow is a judgement.
- **The model will want to drop `--seed` or use the legacy `leidenalg` flavor**
  to match an old notebook. Do not. The igraph flavor with a fixed seed and
  `n_iterations=2` is what makes reruns identical.
- **The model will want to pass log-normalised data** to fix the
  `counts in X must be non-negative integers` error. Do not. Ask where the raw
  counts are; this skill normalises them itself.
- **The model will want to call a non-HVG gene "not expressed"** or absent.
  Do not. Non-detection is not absence: genes outside the HVG set were not used
  for the embedding, and stay in X for later tests.
- **The model will want to read UMAP distances** between clusters as
  similarity. Do not. The caption gives `min_dist`, neighbours and seed because
  the layout depends on them; compare clusters in PC space or by markers.
- **The model will want to guess `--reference`** from the gene names. Do not.
  Ask for the build and annotation release; it is recorded, never inferred.
- **The model will want to name clusters as cell types** from the UMAP. Do not.
  Leiden numbers are not labels; annotation is a separate skill.

## Engine scope

Harmony (`harmonypy`) or no integration. scVI and scANVI (corpus rows A-an-05,
D-an-37) are out of scope: they need a GPU-sized training step and a model
directory. Add them later as `--integration scvi` with a seeded trainer.
CCA (D-an-04) is Seurat-only and out of scope.

## CLI reference

```bash
python scripts/sc_embed.py \
  --adata results/2026-10-27_qc/data/qc_filtered.h5ad \
  --target-sum 10000 --n-top-genes 2000 --flavor seurat_v3 --batch-key sample \
  --n-comps 50 --integration harmony --n-neighbors 20 --n-pcs 30 \
  --min-dist 0.3 --resolution 1 --seed 0 --reference 'GRCh38, GENCODE 44' \
  --outdir results/2026-10-27_embed

# Demo (only when the user asks for one); add --integration none to compare
python scripts/sc_embed.py --demo --outdir demo_out
```

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
