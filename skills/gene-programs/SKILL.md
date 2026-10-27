---
name: gene-programs
description: >-
  Gene-program discovery by consensus NMF (cNMF-style) on raw single-cell counts:
  HVG selection, repeated seeded sklearn NMF, a local-density filter and k-means
  consensus of the replicate components, per-cell program usage, gene scores,
  top genes, mean usage, fraction of active cells and one-vs-rest AUC per group,
  a stability table, four fixed figures, a report and a reproducibility bundle.
  Fire on "gene programs", "cNMF", "consensus NMF", "gene expression programs",
  "GEP usage", "NMF programs per cell type", "which programs mark each cluster",
  "structure plot of program usage".
  Do NOT fire for scoring a known gene set (score_genes, UCell), for projecting
  onto a fixed reference (starCAT, NMFproj), for choosing k by a sweep, for
  marker genes per cluster (use sc-markers), or for log-normalised input.
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - single-cell
    - gene programs
    - NMF
    - cNMF
    - structure plot
  inputs:
    - name: adata
      type: file
      format:
        - h5ad
      description: cell-level AnnData, raw integer counts in X or a layer
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
      description: counts, top genes per program, stability, params, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - anndata
      - scanpy
      - scikit-learn
      - numpy
      - pandas
      - scipy
      - matplotlib
  demo_data:
    - path: examples/demo.h5ad
      description: synthetic counts, 3 groups x 60 cells, 300 genes, 3 planted programs of 10 genes
  endpoints:
    cli: python scripts/gene_programs.py --adata {adata} ... --outdir {outdir}
  trigger_keywords:
    - gene programs
    - cNMF
    - consensus NMF
    - program usage
    - structure plot
---

# Gene programs by consensus NMF

You are the dispatcher for **gene-programs**. The script does the analysis; you
collect its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/gene_programs.py`. Never re-implement the analysis, change a
threshold, or pick a value the user has not given. If a required argument is
unknown, ask the user (these are grill questions), or run `--demo` only when the
user asks for a demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--adata` | `.h5ad` | raw integer counts in `X`, or in `--layer`; QC-filtered cells; no cell with zero counts over the HVGs | sc-qc's filtered h5ad; `examples/demo.h5ad` |
| `--groupby` | obs column | at least 2 groups, for the program-by-group summary | `leiden`, `cell_type`; demo `group` |

The script refuses negative or non-integer values (log-normalised data). A
categorical `--groupby` keeps its category order; any other column is cast to
strings and sorted.

## Parameters (all required, no defaults)

Each row is a grill question. Ask the user; never fill one in.

| Argument | Meaning | Grill question | Demo-only value |
|---|---|---|---|
| `--adata` | input h5ad | Which h5ad holds the raw counts? | `examples/demo.h5ad` |
| `--layer` (optional) | count matrix | Which layer holds the raw counts? Absent means `X`. | absent (`X`) |
| `--groupby` | obs column of groups | Which obs column holds the clusters or labels to summarise programs by? | `group` |
| `--reference` | genome build and annotation | Which genome build and annotation release name the genes (for example `GRCh38, GENCODE 44`)? | `synthetic demo, no genome build` |
| `--k` | number of programs | How many programs (k)? Choose it before the run; this skill does not sweep k. | `3` |
| `--n-hvg` | genes factorised | How many highly variable genes go into NMF? | `200` |
| `--n-replicates` | NMF runs | How many NMF replicates (at least 4)? cNMF uses 100-200 for real data. | `20` |
| `--density-threshold` | outlier filter | Above which mean distance to its nearest replicate components is a component dropped? `2` keeps all. | `0.5` |
| `--active-threshold` | active-cell rule | Above which fraction of the program's maximum usage does a cell count as active? | `0.1` |
| `--n-top-genes` | top genes per program | How many top genes per program go in the table and heatmap? | `10` |
| `--seed` | seed | Which seed? It derives the replicate seeds and fixes k-means and the jitter. | `0` |
| `--outdir` | new output directory | Which new directory gets the outputs? It must not exist. | `gene_programs_demo` |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them.
With `--demo`, only `--outdir` may also be given; any other argument is an error.

```bash
python scripts/gene_programs.py \
  --adata results/2026-10-27_qc/data/filtered.h5ad --layer counts --groupby leiden \
  --reference "GRCh38, GENCODE 44" \
  --k 12 --n-hvg 2000 --n-replicates 100 --density-threshold 0.1 \
  --active-threshold 0.1 --n-top-genes 20 --seed 0 \
  --outdir results/2026-10-27_programs

# Demo (only when the user asks for one)
python scripts/gene_programs.py --demo --outdir demo_out
```

## Workflow (fixed order)

1. **Validate** arguments and input; bad input gives `ERROR: <message>`, exit 2.
   Gene names and cell ids must be unique (`skillkit.check_unique`). The matrix
   must be non-negative integers and `--groupby` must have 2+ groups.
2. **Select HVGs**: the top `--n-hvg` expressed genes by Seurat dispersion on a
   `normalize_total` + `log1p` copy. The NMF input is the HVG counts scaled to
   unit variance per gene, not centred (cNMF `prepare`).
3. **Factorise** `--n-replicates` times with sklearn `NMF` (Frobenius, `cd`
   solver, random init, 1000 iterations); replicate seeds are drawn from
   `numpy.random.default_rng(--seed)`.
4. **Consensus**: L2-normalise the `k × n-replicates` components; drop each
   whose mean Euclidean distance to its `int(0.3 × n-replicates)` nearest
   components is not below `--density-threshold`; k-means (`k`, 10 starts,
   `random_state=--seed`); the median of each cluster, scaled to sum 1, is a
   consensus spectrum. Fewer than `k` kept components stops the run.
5. **Fit usage** with the consensus spectra fixed (`non_negative_factorization`,
   `update_H=False`), then scale each cell's usage to sum 1 (cNMF `load_results`).
6. **Name programs** `GP1`..`GPk` by total usage across cells, largest first.
7. **Score genes** for every expressed gene: OLS of per-gene z-scored TPM on
   usage (cNMF `gene_spectra_score`). Top genes are the highest scores.
8. **Summarise by group**: mean usage; fraction of active cells, where usage
   divided by the program's maximum over cells exceeds `--active-threshold`
   (immgenT Figure2.R); and one-vs-rest AUC of usage (`roc_auc_score`, fixed
   direction: above 0.5 means higher in the group, below 0.5 lower).
9. **Draw** four figures, then write tables, h5ad, report, result and bundle.

## Output structure

Sort keys: `usage.tsv` by program then cell; `gene_loadings.tsv` by program,
score descending, gene; `top_genes.tsv` likewise with a 1-based `rank`;
`program_by_group.tsv` and `stability.tsv` by program then group order.

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # n_cells, n_genes, n_hvg, n_genes_scored, k, n_components_kept, n_components_total, mean_silhouette, top_genes, seed, params, reference, chat_summary_lines, preferred_artifacts
├── tables/
│   ├── usage.tsv              # cell, <groupby>, program, usage (cell sums to 1), usage_max_scaled, active
│   ├── gene_loadings.tsv      # program, gene, score (z-score units), hvg; expressed genes only
│   ├── top_genes.tsv          # program, rank, gene, score
│   ├── program_by_group.tsv   # program, <groupby>, n_cells, mean_usage, frac_active, auc
│   └── stability.tsv          # program, n_components, n_replicates, silhouette (mean over its components)
├── figures/
│   ├── usage_heatmap.pdf
│   ├── usage_heatmap.png
│   ├── usage_heatmap.csv
│   ├── usage_structure.pdf
│   ├── usage_structure.png
│   ├── usage_structure.csv
│   ├── top_gene_heatmap.pdf
│   ├── top_gene_heatmap.png
│   ├── top_gene_heatmap.csv
│   ├── auc_by_group.pdf
│   ├── auc_by_group.png
│   └── auc_by_group.csv
├── data/
│   └── programs.h5ad          # input object + obsm["usage"], varm["gene_scores"], var["highly_variable"]
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python and package versions
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

### Fixed plot spec

Figures come from the vendored `scripts/plots.py` (heatmaps, box) or matplotlib
in the same theme (structure plot), saved by `save_panel` as PDF, PNG and source
CSV. No titles; captions sit in `report.md`.

- `usage_heatmap` (`plot_heatmap`, D-pl-28, D-pl-34): program × group mean
  usage, centred and scaled per program across groups.
- `usage_structure` (D-pl-29): one stacked column per cell (`stackplot`, so the
  PDF stays small), cells ordered by group, then dominant program, then its usage.
- `top_gene_heatmap` (`plot_heatmap`): program × top-gene scores, scaled per
  gene across programs.
- `auc_by_group` (`plot_box`, D-pl-30): one point per program per group, dashed
  line at AUC 0.5. Under 10 programs (`FEW_POINTS`) it draws points and a
  median, no box.

Every caption names k, the HVG count, the replicate count and the seed. Both
heatmaps are fixed by label, not clustered: programs by name (GP1 has the
largest total usage), groups by `--groupby` order, genes by program then rank.

Display constants, not analysis parameters: heatmaps saturate at ±2
(`HEATMAP_LIMIT`, house-style `fig.heatmap`); jitter ±0.2 box widths
(`JITTER_WIDTH`). `scripts/skillkit.py`, `plots.py`, `theme.py` and
`palettes.py` are byte copies of `templates/`; a test fails when they drift.

## Gotchas

- **The model will want to feed log-normalised or scaled data** to get past the
  `does not look like raw counts` error. Do not. NMF here models counts; point
  `--layer` at the raw counts or ask the user where they are.
- **The model will want to change `--k` after seeing the programs**, or run
  several k and keep the nicest. Do not. k is the user's answer; a k sweep with
  a stability-versus-error plot is a separate task the user asks for.
- **The model will want to raise `--density-threshold`** after `too_few_kept`.
  Do not. Report the message; the threshold is the user's choice.
- **The model will want to name programs** ("interferon program") from the top
  genes and write that into the tables. Do not. Labels `GP1`..`GPk` stay; any
  naming is the user's interpretation, stated as such.
- **The model will want to read AUC below 0.5 as no signal.** Do not. The AUC is
  directional: 0.2 means the program is lower in that group.
- **The model will want to call a high mean usage a group difference with a
  p-value.** Do not. Cells are not replicates; per-donor tests go to group-tests
  or pseudobulk-de through `tables/usage.tsv`.
- **The model will want to read usage or AUC as a probability.** Do not. Usage
  is a cell's share of its fitted expression, and AUC ranks cells; neither is
  the chance that a cell is in a program or a group.
- **The model will want to say an unexpressed gene is "not in the program".**
  Do not. A gene with no counts gets no score: it is not tested, so absence
  from `gene_loadings.tsv` says nothing about it.
- **The model will want to guess `--reference`** from the gene names. Do not.
  The build and annotation release are the user's answer.
- **The model will want to add the `cnmf` package or another NMF library.** Do
  not. The engine is fixed: sklearn NMF plus the consensus steps above.

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
