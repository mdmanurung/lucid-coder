---
name: sc-markers
description: >-
  Marker genes per cluster or cell-type label from a log-normalised single-cell
  h5ad: scanpy rank_genes_groups (Wilcoxon or t-test, each group against the
  rest, Benjamini-Hochberg), a fixed pass rule on adjusted p, log fold change and
  fraction expressing, a top-marker dot plot and heatmap, a report and a
  reproducibility bundle. Fire on "marker genes per cluster", "find markers",
  "rank_genes_groups", "FindAllMarkers in scanpy", "top markers per Leiden
  cluster", "marker dotplot", "what genes define each cell type".
  Do NOT fire for DE between conditions or donors (use pseudobulk-de), for
  scoring a known gene set, for annotation from a reference, or for raw counts
  that are not yet normalised.
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - single-cell
    - marker genes
    - rank_genes_groups
    - scanpy
    - dotplot
  inputs:
    - name: adata
      type: file
      format:
        - h5ad
      description: cell-level AnnData, log-normalised values in X or a layer
      required: true
    - name: gene_id_col
      type: string
      description: var column of Ensembl gene ids, or none when the source has none
      required: true
    - name: context
      type: string
      description: species, tissue or cell type, condition and dataset
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
      description: counts, top markers per group, params, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - anndata
      - scanpy
      - numpy
      - pandas
      - scipy
      - matplotlib
      - plotnine-extra
  demo_data:
    - path: examples/demo.h5ad
      description: synthetic, 3 clusters x 50 cells, 200 genes, 5 planted markers each
  endpoints:
    cli: python scripts/sc_markers.py --adata {adata} ... --outdir {outdir}
  trigger_keywords:
    - marker genes
    - cluster markers
    - find markers
    - rank_genes_groups
    - marker dotplot
---

# Marker genes

You are the dispatcher for **sc-markers**. The script does the analysis; you
collect its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/sc_markers.py`. Never re-implement the analysis, change a
threshold, or pick a value the user has not given. If a required argument is
unknown, ask the user (these are grill questions), or run `--demo` only when the
user asks for a demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--adata` | `.h5ad` | log-normalised values (`normalize_total` + `log1p`) in `X`, or in `--layer`; not raw counts, not scaled | sc-embed's `data/embedded.h5ad`; `examples/demo.h5ad` |
| `--groupby` | obs column | at least 2 groups, each with at least 2 cells | `leiden`, `cell_type`; demo `cluster` |
| `--gene-id-col` | var column | Ensembl gene ids (versioned as the source gives them); `var_names` are the symbols | `gene_ids`; `none` when the source has no ids |

The script refuses a matrix with negative values (scaled) or only integers (raw
counts). A categorical `--groupby` keeps its category order; any other column
is cast to strings and sorted. An h5ad `.raw` is ignored (`use_raw=False`).
Cell ids and gene ids must be unique, and each symbol must map to one gene id
(`skillkit.check_gene_ids`, house-style `data.ids`); a duplicate or an ambiguous
map stops the run before any test. `--gene-id-col none` keeps symbols only:
`gene_id` is written as NA, never guessed.

## Parameters (all required, no defaults)

Each row is a grill question. Ask the user; never fill one in.

| Argument | Meaning | Grill question | Demo-only value |
|---|---|---|---|
| `--adata` | input h5ad | Which h5ad holds the log-normalised cells? | `examples/demo.h5ad` |
| `--layer` (optional) | matrix to test | Which layer holds the log-normalised values? Absent means `X`. | absent (`X`) |
| `--groupby` | obs column of groups | Which obs column holds the clusters or labels? | `cluster` |
| `--gene-id-col` | var column of Ensembl ids, or `none` | Which var column holds the Ensembl gene ids (or `none` if the source has none)? | `gene_ids` |
| `--method` | `wilcoxon` or `t-test` | Wilcoxon rank-sum or Welch t-test? | `wilcoxon` |
| `--padj-max` | adjusted-p cutoff | Which BH-adjusted p cutoff marks a marker? | `0.05` |
| `--lfc-min` | log2 fold change floor | Which log2 fold change (group vs rest) must a marker exceed? | `1.0` |
| `--min-pct` | fraction expressing in the group | In what fraction of the group's cells must a marker be detected? | `0.25` |
| `--n-top` | markers per group in the figures | How many top markers per group go in the figures? | `5` |
| `--limit` | heatmap colour saturation | At which absolute z-score should the heatmap colours saturate? | `2.0` |
| `--reference` | genome build and annotation release of the gene names | Which genome build and annotation were the genes counted against (e.g. `GRCh38, GENCODE 44`)? | `synthetic demo genes, no genome build` |
| `--context` | biological context, printed first in the report and summary | Species, tissue or cell type, condition and dataset? | `synthetic demo: 3 clusters x 50 cells, no real tissue` |
| `--seed` | seed | Which seed? (Recorded; neither test draws random numbers.) | `0` |
| `--outdir` | new output directory | Which new directory gets the outputs? It must not exist. | `sc_markers_demo` |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them.
With `--demo`, only `--outdir` may also be given; any other argument is an error.

```bash
python scripts/sc_markers.py \
  --adata results/2026-10-27_embed/data/embedded.h5ad --groupby leiden \
  --gene-id-col gene_ids --context "human PBMC, healthy donors, Leiden clusters" \
  --method wilcoxon --padj-max 0.05 --lfc-min 1 --min-pct 0.25 \
  --n-top 5 --limit 2 --reference "GRCh38, GENCODE 44" --seed 0 \
  --outdir results/2026-10-27_markers

# Demo (only when the user asks for one)
python scripts/sc_markers.py --demo --outdir demo_out
```

## Workflow (fixed order)

1. **Validate** arguments and input; bad input gives `ERROR: <message>`, exit 2.
   Cell ids must be unique; `check_gene_ids` stops on an ambiguous symbol-to-id map.
2. **Check** the matrix is log-normalised and `--groupby` has at least 2 groups
   of at least 2 cells.
3. **Test** each group against the rest with `sc.tl.rank_genes_groups`
   (`--method`, `reference="rest"`, `corr_method="benjamini-hochberg"`,
   `pts=True`, `use_raw=False`), every gene. The BH family is one group's genes.
4. **Call markers**: `passes` when `padj < --padj-max`, `logfc > --lfc-min` and
   `pct_in >= --min-pct`. One-sided: a marker is up in its group. `logfc` is
   scanpy's log2 ratio of `expm1` group means, so it assumes natural `log1p` data.
5. **Pick top markers**: the first `--n-top` passing genes per group by score.
   A group with fewer passing genes gets fewer, never padded. Zero passing genes
   in every group stops the run (`no_markers`).
6. **Draw** the dot plot and heatmap of the top markers; a gene topping two
   groups is drawn once, in its first group's block.
7. **Write** report, result and reproducibility bundle.

## Output structure

`markers.tsv` holds every gene for every group, sorted by group (category
order), then score descending, then gene symbol. `top_markers.tsv` keeps the same
order with a 1-based `rank` per group. `gene_id` sits beside `gene_symbol` (NA with
`--gene-id-col none`). `status` (house-style `data.status`) is `tested` for every
row: each group tests every gene, with no filter; `result.json` counts it in
`status_counts`.

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # n_cells, n_groups, n_genes_tested, status_counts, gene_id_col, n_pass, top_markers, method, seed, params, context, reference, chat_summary_lines, preferred_artifacts
├── tables/
│   ├── markers.tsv            # group, gene_symbol, gene_id, score, logfc, pval, padj, pct_in, pct_out, passes, status
│   └── top_markers.tsv        # group, rank, gene_symbol, gene_id, score, logfc, pval, padj, pct_in, pct_out
├── figures/
│   ├── marker_dotplot.pdf
│   ├── marker_dotplot.png
│   ├── marker_dotplot.csv
│   ├── marker_heatmap.pdf
│   ├── marker_heatmap.png
│   └── marker_heatmap.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python and package versions
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

## Gotchas

- **The model will want to run `rank_genes_groups` on raw counts or on scaled
  data** to get past the `does not look log-normalised` error. Do not. Point
  `--layer` at the log-normalised values or ask the user where they are.
- **The model will want to call these markers "differentially expressed between
  conditions", or read `padj` as donor-level evidence.** Do not. Cells are the
  unit only when the condition was assigned to each cell (a per-cell
  perturbation). Clusters are not assigned per cell, so the donor is the unit and
  these p-values are descriptive: they rank genes within one dataset. Condition
  DE goes to pseudobulk-de.
- **The model will want to call a gene with `pct_in = 0` "absent" from a group.**
  Do not. Non-detection is not absence: the gene was not detected in those cells.
  `score` is a test statistic, not a probability.
- **The model will want to say a marker "defines", "drives" or "regulates" its
  cluster.** Do not. A marker is a gene associated with the group by a test of
  that group against the rest, in this dataset; it is not a mechanism or a cell
  identity (house-style `claims.strength`).
- **The model will want to relax `--padj-max`, `--lfc-min` or `--min-pct`**
  after a group comes back with no markers. Do not. Report the group as having
  none; thresholds are the user's answers, fixed before the run.
- **The model will want to pad a group's top markers** with genes that do not
  pass, so every group shows `--n-top`. Do not. The figures show passing genes only.
- **The model will want to add hand-picked canonical markers** (CD3E, MS4A1) to
  the dot plot. Do not. The figures show this run's markers; a curated panel is a
  separate plot the user asks for.
- **The model will want to fill in a missing argument** with a "standard" value
  (0.05, 0.25, 5). Do not. Ask, or run `--demo` only when the user asked for one.

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
