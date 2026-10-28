---
name: sc-normalize
description: >-
  Normalise a single-cell count h5ad one of three fixed ways: log1p of
  library-size-scaled counts, analytic Pearson residuals (the model SCTransform
  v2 approximates), or dsb for ADT counts using empty droplets and isotype
  controls. Writes data/normalized.h5ad with the input counts kept, per-cell
  and per-feature tables, mean-variance, depth-vs-PC1 and value-distribution
  figures, a report and a reproducibility bundle. Fire on "normalise my
  counts", "Pearson residuals", "SCTransform-style normalisation", "dsb
  normalise the ADT", "denoise CITE-seq protein with empty droplets", "does
  depth still drive PC1". Do NOT fire for QC filtering (sc-qc), embedding and
  clustering (sc-embed, which does its own log1p), pseudobulk voom
  (pseudobulk-de), or R SCTransform itself.
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - single-cell
    - normalisation
    - Pearson residuals
    - SCTransform
    - dsb
    - CITE-seq
  inputs:
    - name: context
      type: string
      description: species, tissue or cell type, condition and dataset
      required: true
    - name: adata
      type: file
      format:
        - h5ad
      description: raw integer counts in X (RNA, or ADT of cells for dsb)
      required: true
    - name: empty_droplets
      type: file
      format:
        - h5ad
      description: ADT counts of empty droplets, same features (dsb only)
      required: false
  outputs:
    - name: normalized
      type: file
      format:
        - h5ad
      description: X and layers[--layer-out] normalised, layers["counts"] input
    - name: report
      type: file
      format:
        - md
      description: params, summary and captioned figures
    - name: result
      type: file
      format:
        - json
      description: PC1-depth correlation, sizes, params, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - scanpy
      - anndata
      - muon
      - scikit-learn
      - numpy
      - pandas
      - scipy
      - matplotlib
      - plotnine-extra
      - threadpoolctl
  demo_data:
    - path: examples/demo.h5ad
      description: synthetic RNA, 3 groups x 200 cells, 300 genes, wide depth spread
    - path: examples/demo_adt.h5ad
      description: ADT of the same 600 cells, 4 markers and 3 isotypes
    - path: examples/demo_empty.h5ad
      description: ADT of 2000 empty droplets, same 7 features
  endpoints:
    cli: python scripts/sc_normalize.py --adata {adata} --method {method} ... --outdir {outdir}
  trigger_keywords:
    - normalise counts
    - Pearson residuals
    - SCTransform
    - dsb ADT
    - CITE-seq protein normalisation
---

# sc-normalize

You are the dispatcher for **sc-normalize**. The script does the analysis; you
collect its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/sc_normalize.py`. Never re-implement the analysis, change a
parameter, or pick a value the user has not given. If a required argument is
unknown, ask the user (these are grill questions), or run `--demo` only when the
user asks for a demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

**Fire this skill when the user says any of:**
- "normalise my QC'd counts with Pearson residuals, theta 100"
- "compare log1p and Pearson residuals: does depth still drive PC1?"
- "dsb-normalise the ADT with these empty droplets and isotypes"

**Do NOT fire when:**
- the cells still need QC (use sc-qc first);
- the user wants a UMAP or clusters (sc-embed normalises by log1p itself);
- the data are pseudobulk counts for DE (pseudobulk-de runs voom);
- the user asks for R `SCTransform()` with `vars.to.regress` (out of scope).

### Engine scope

- `log1p`: scanpy `normalize_total` then `log1p`, the sc-embed baseline.
- `pearson`: `scanpy.experimental.pp.normalize_pearson_residuals`, the analytic
  Pearson residuals of Lause, Berens and Kobak (2021) with one fixed overdispersion
  `theta` for all genes. SCTransform v2 (Choudhary and Satija 2022) fits a
  regularised NB GLM per gene and is close to this model; this skill does **not**
  run SCTransform. The corpus "port" (`single_cell_aging/.../sctransform.py`,
  from sctoolkit) is an rpy2 wrapper around R `sctransform`, not a Python
  re-implementation, and also regresses out other covariates; neither is here.
- `dsb`: muon's `mu.prot.pp.dsb` (muon 0.2.0), the method of Mulè, Martins and
  Tsang (2022): log(x + 10), standardise each protein by the mean and SD of the
  empty droplets, then per cell take the background mean from a 2-component GMM,
  combine it with the isotype values into one technical component (first PC)
  and regress it out of every protein. The corpus (ALTRA TEA-seq) calls the
  same function.
- Out of scope: voom/log2 CPM on pseudobulk (E-an-03, in pseudobulk-de) and
  GTestimate (E-an-19, an R Good-Turing estimator with no Python port).

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--adata` | `.h5ad` | raw counts (non-negative integers) in `X`; no cell with zero counts. RNA for `log1p`/`pearson`, ADT of cells for `dsb` | sc-qc `data/qc_filtered.h5ad`, `examples/demo.h5ad` |
| `--empty-droplets` | `.h5ad` | ADT counts of empty droplets in `X`, the same features in the same order as `--adata`, no barcode shared with `--adata` | `examples/demo_empty.h5ad` |
| gene ids | `var[--gene-id-col]` | log1p/pearson: Ensembl gene ids, one per symbol in `var_names`; no id shared by two symbols | `gene_ids` in `examples/demo.h5ad` |
| isotypes | `var_names` | every `--isotype-controls` name is a feature of `--adata` | `IgG1,IgG2a,IgG2b` |

## Parameters (all required, no defaults)

Each row is a grill question. Ask the user; never fill one in. Method rows are
required for that method and refused for the others.

| Argument | Grill question | Demo-only value |
|---|---|---|
| `--adata` | Which h5ad holds the raw counts (ADT of cells for dsb)? | `examples/demo.h5ad`; dsb: `examples/demo_adt.h5ad` |
| `--method` | Normalise by `log1p`, `pearson` residuals or ADT `dsb`? | `pearson` |
| `--layer-out` | Which layer name gets the normalised values (not `counts`)? | `normalized` |
| `--seed` | Which seed for the dsb GMM and the diagnostic PCA? | `0` |
| `--reference` | Which genome build and annotation release named the genes (e.g. `GRCh38, GENCODE 44`)? For dsb, which antibody panel and release? Recorded in `result.json`. | `synthetic demo features, no genome build` |
| `--context` | Species, tissue or cell type, condition and dataset? | `synthetic demo: 3 groups x 200 cells, RNA and ADT, no real tissue` |
| `--target-sum` | log1p: to what total is each cell scaled first? | `10000` |
| `--theta` | pearson: which NB overdispersion theta (100 in Lause et al.)? | `100` |
| `--gene-id-col` | log1p/pearson: which `var` column holds the Ensembl gene ids (`var_names` are the symbols)? `none` if the source has no ids; `gene_id` is then NA. | `gene_ids` |
| `--clip` | pearson: clip residuals to +-what (scanpy suggests sqrt(n cells))? | `25` |
| `--empty-droplets` | dsb: which h5ad holds the ADT counts of empty droplets, and how were they chosen? | `examples/demo_empty.h5ad` |
| `--isotype-controls` | dsb: which features are isotype controls (comma list)? | `IgG1,IgG2a,IgG2b` |
| `--outdir` | Which new directory gets the outputs? It must not exist. | `sc_normalize_demo` |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them.
With `--demo`, only `--method` and `--outdir` may also be given; `--demo --method
dsb` switches to the ADT demo files.

## Workflow (fixed order)

1. **Validate** arguments and inputs; bad input gives `ERROR: <message>`, exit 2.
   Nothing is written before every check passes. Feature names and cell ids
   must be unique in `--adata` and `--empty-droplets` (`skillkit.check_unique`).
   For log1p/pearson, `skillkit.check_gene_ids` stops on a gene id shared by
   two symbols.
2. **Drop all-zero features** (zero counts in every cell); the count goes in
   `result.json` and the summary, and each stays in `gene_mean_variance.tsv`
   with `status` `not_detected` and NA statistics. Pearson residuals are 0/0 for them.
3. **Keep counts** in `layers["counts"]`.
4. **Normalise X** by `--method`:
   - `log1p`: `normalize_total(target_sum=--target-sum)`, `log1p`;
   - `pearson`: `normalize_pearson_residuals(theta=--theta, clip=--clip)`;
   - `dsb`: `mu.prot.pp.dsb(cells, empty, pseudocount=10, denoise_counts=True,
     isotype_controls=..., scale_factor="standardize", random_state=--seed)`.
5. **Copy** the dense normalised X into `layers[--layer-out]`.
6. **Diagnostics**: per-cell library size, normalised total and detected
   features; per-feature mean and variance before and after; PC1 from a
   2-component arpack PCA (seeded, all features, centred, unscaled) and its
   |Spearman rho| with library size.
7. **Figures** and report.

BLAS runs on one thread (`threadpoolctl`). Fixed constants: dsb pseudocount 10,
6 features in the value plot, 20 quantile bins for the trend line.

## Output structure

`per_cell_summary.tsv` is sorted by cell; `gene_mean_variance.tsv` by feature.
`gene_mean_variance.tsv` keeps every input feature: `status` is `tested`, or
`not_detected` for an all-zero feature (NA statistics, left out of the
figures). `gene_id` (log1p/pearson only) is NA under `--gene-id-col none`.

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # method, n_cells, n_features, n_features_dropped, status_counts, gene_id_col, pc1_depth_abs_spearman, layer_out, seed, params, reference, context, chat_summary_lines (first: "Context: ..."), preferred_artifacts
├── data/
│   └── normalized.h5ad        # X normalised (dense), layers["counts"] input, layers[--layer-out] = X
├── tables/
│   ├── per_cell_summary.tsv   # cell, library_size, normalised_total, n_detected, pc1
│   └── gene_mean_variance.tsv # feature, gene_id (not dsb), mean_counts, var_counts, mean_normalised, var_normalised, status
├── figures/
│   ├── mean_variance.pdf
│   ├── mean_variance.png
│   ├── mean_variance.csv
│   ├── depth_pc1.pdf
│   ├── depth_pc1.png
│   ├── depth_pc1.csv
│   ├── top_variable_values.pdf
│   ├── top_variable_values.png
│   ├── top_variable_values.csv
│   ├── isotype_vs_background.pdf   # dsb only
│   ├── isotype_vs_background.png   # dsb only
│   └── isotype_vs_background.csv   # dsb only
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python and package versions
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

Figures come from the vendored `plots.py` (`plot_scatter_fit`) and matplotlib in
the same theme, saved by `save_panel`; captions sit in `report.md`.

- `mean_variance` (E-pl-04, voom-like): log10 variance against log10 mean count
  per feature, counts (left, Poisson line) and normalised (right), binned-median
  trend.
- `depth_pc1`: PC1 against log10 library size, linear fit and 95% CI.
- `top_variable_values`: box plots of normalised values for the 6 features with
  the highest normalised variance.
- `isotype_vs_background` (dsb): per-cell median protein log(count + 10)
  against mean isotype log(count + 10). The median stands in for dsb's GMM
  background mean, which muon does not return.

`scripts/skillkit.py`, `plots.py`, `theme.py` and `palettes.py` are byte copies
of `templates/`; a test fails when they drift.

## Gotchas

- **The model will want to pick `--theta` or `--clip` itself** ("scanpy uses
  100"). Do not. They are grill questions; quote Lause et al. and scanpy's
  sqrt(n) suggestion, then ask.
- **The model will want to call this SCTransform** in a methods section. Do not.
  Write "analytic Pearson residuals (Lause et al. 2021)", and note that
  SCTransform v2 approximates the same model.
- **The model will want to pick empty droplets by a UMI range itself** (the
  corpus used `empty_counts_range=(1.2, 2.8)`). Do not. Which barcodes are empty
  is the user's call; ask for the file and how it was made.
- **The model will want to drop `--isotype-controls`** when the panel has none.
  Do not. Stop and tell the user; dsb without isotypes is a different
  correction and is not offered here.
- **The model will want to feed normalised values back in** to fix
  `counts in X ... must be non-negative integers`. Do not. Ask where the raw
  counts are; this skill normalises them itself.
- **The model will want to copy between `.raw` and `.X`** to "restore" counts
  or stash a normalised copy. Do not. Raw counts live in `layers["counts"]`
  only; a log-normalised `X` cannot give counts back without its size factors.
- **The model will want to guess what an unknown h5ad holds.** Do not. Check a
  few rows first: integers mean counts; `np.expm1(X[:5]).sum(1)` near one total
  means log1p-normalised; negative values mean residuals or scaled data. Then
  ask for the raw counts if `X` is not counts (source: kai,
  davidfischerlab/kai `prompt_manager.py`).
- **The model will want to call a dropped all-zero feature "not expressed".**
  Do not. Non-detection is not absence: the feature was not measured above
  zero in these cells, so downstream it is "not tested".
- **The model will want to guess `--reference`** from the gene names. Do not.
  Ask for the build and annotation release (or the ADT panel); it is recorded,
  never inferred.
- **The model will want to say depth "drives" PC1.** Do not. |rho| is a rank
  correlation, so write "PC1 is associated with library size" (house-style
  `claims.strength`).
- **The model will want to read a high |rho| as a failure and switch method.**
  Do not. Report both runs; a cell type with more RNA can also tie PC1 to depth.

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
