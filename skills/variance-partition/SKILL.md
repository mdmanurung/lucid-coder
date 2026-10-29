---
name: variance-partition
description: >-
  Variance partitioning of pseudobulk gene expression across sample covariates
  with R variancePartition: voomWithDreamWeights, then fitExtractVarPartModel
  with categorical covariates as random effects (1|x) and continuous ones as
  fixed, after a canCorPairs collinearity gate. Writes per-gene variance
  fractions, a median summary, the canonical-correlation matrix, a violin plot,
  percent bars and a heatmap. Fire on "variance partition", "variancePartition",
  "fitVarPart", "plotVarPart", "how much variance does donor explain",
  "which covariates drive expression", "covariate selection before DE",
  "canCorPairs". Do NOT fire for differential expression (use pseudobulk-de),
  for PCA-only covariate screens, for cell-level variance (no pseudobulk), or
  for cell-type composition.
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - pseudobulk
    - variance partitioning
    - variancePartition
    - mixed models
    - covariates
  inputs:
    - name: counts
      type: file
      format:
        - tsv
      description: genes x samples summed raw counts (pseudobulk-de tables/pseudobulk_counts.tsv)
      required: false
    - name: meta
      type: file
      format:
        - tsv
      description: one row per sample, formula columns (pseudobulk-de tables/sample_meta.tsv)
      required: false
    - name: adata
      type: file
      format:
        - h5ad
      description: cell-level AnnData with raw integer counts, summed per sample
      required: false
    - name: gene_id_col
      type: string
      description: var column of Ensembl gene ids, or none (--adata input)
      required: false
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
      description: top covariate, median fractions, top genes per term, params
  dependencies:
    python: ">=3.11"
    packages:
      - numpy
      - pandas
      - matplotlib
      - plotnine-extra
      - anndata
      - scipy
    r:
      - variancePartition
      - edgeR
      - lme4
  demo_data:
    - path: examples/demo_counts.tsv
      description: synthetic, 12 donors x ctrl/stim, 120 genes, planted donor/condition/noise genes
    - path: examples/demo_meta.tsv
      description: sample, donor, condition, rin
  endpoints:
    cli: python scripts/variance_partition.py --counts {counts} --meta {meta} ... --outdir {outdir}
  trigger_keywords:
    - variance partition
    - variancePartition
    - fitVarPart
    - fitExtractVarPartModel
    - plotVarPart
    - canCorPairs
---

# Variance partition

You are the dispatcher for **variance-partition**. The script does the analysis;
you collect its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/variance_partition.py`. Never re-implement the analysis, change a
threshold, or pick a value the user has not given. If a required argument is
unknown, ask the user (these are grill questions), or run `--demo` only when the
user asks for a demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

## Input contract

Give one of two inputs: `--counts` with `--meta`, or `--adata` with
`--gene-id-col`, `--sample-col`, `--donor-col` and `--min-cells`. Giving both is an error.

| Input | Format | Required content | Example |
|---|---|---|---|
| `--counts` | TSV | columns `gene_id` (Ensembl; empty when the source has none) and `gene_symbol`, then one column per sample; non-negative integer summed counts | pseudobulk-de `tables/pseudobulk_counts.tsv`; `examples/demo_counts.tsv` |
| `--meta` | TSV | first column sample IDs, the same set as the `--counts` columns; every formula term | pseudobulk-de `tables/sample_meta.tsv`; `examples/demo_meta.tsv` |
| `--adata` | `.h5ad` | raw integer counts in `X`, or in `--layer`; one cell per row | one cell type's cells |
| `--gene-id-col` | var column | Ensembl gene ids; `var_names` are the symbols; `none` when the source has no ids | `gene_ids` |
| `--sample-col` | obs column | the sample each cell belongs to; every formula term is constant within a sample | `sample` |
| `--donor-col` | obs column | the donor each sample comes from; constant within a sample | `donor` |

A `(1|x)` term is categorical, even when its values look numeric. A plain term
must be numeric (a continuous covariate); a categorical plain term is an error,
because variancePartition needs categorical covariates as random effects.

Required checks (docs/pipeline-skills.md, adopted 2026-10-27):
- **Unique ids**: `skillkit.check_unique` stops on a repeated sample id in
  `--counts` (header read as written) and `--meta`, and on a repeated cell id in
  `--adata`.
- **Gene ids** (house-style `data.ids`): `skillkit.check_gene_ids` stops on a
  repeated gene id or a symbol with several ids. With no ids (`--gene-id-col
  none`, or an empty `gene_id` column) symbols must be unique and `gene_id` is
  written as NA, never guessed.
- **Design rank**: the plain (fixed) terms go through `skillkit.check_design_rank`;
  aliased numeric covariates stop the run. Random terms are checked by
  `--max-cancor`.
- **Donor versus sample** (`--adata`): `skillkit.check_replication` stops when a
  donor gives several samples and the formula has no donor term. `result.json`
  reports `n_donors` (null for `--counts`, whose donors the script cannot see:
  put `(1|donor)` in the formula yourself) and `n_samples`.
- **Reference record**: `--reference` (genome build and annotation release) is
  required and written to `result.json`.
- **Drops per group**: there is no contrast group here; dropped samples are
  listed by name (`dropped_samples`).

## Parameters (all required, no defaults)

Each row is a grill question. Ask the user; never fill one in.

| Argument | Meaning | Grill question | Demo-only value |
|---|---|---|---|
| `--counts` | summed counts TSV | Which pseudobulk counts table? | `examples/demo_counts.tsv` |
| `--meta` | sample metadata TSV | Which sample table holds the covariates? | `examples/demo_meta.tsv` |
| `--adata` | cell-level h5ad (instead of the two above) | Which h5ad holds the cells, already subset to one cell type? | absent |
| `--layer` (optional, `--adata` only) | raw-count layer | Which layer holds the raw counts? Absent means `X`. | absent |
| `--gene-id-col` (`--adata` only) | var column of Ensembl ids, or `none` | Which var column holds the Ensembl gene ids (or `none` if the source has none)? | absent |
| `--sample-col` (`--adata` only) | obs column of samples | Which obs column names the sample? | absent |
| `--donor-col` (`--adata` only) | obs column of donors | Which obs column names the donor? Several samples per donor need a donor term in `--formula`. | absent |
| `--min-cells` (`--adata` only) | sample floor | How many cells must a sample have to be kept? | absent |
| `--formula` | model | Which covariates enter, and which are categorical (`(1|x)`) or continuous (`x`)? | `~ (1|donor) + (1|condition) + rin` |
| `--min-frac` | gene filter | In what fraction of cells (`--adata`) or samples (`--counts`) must a gene be detected? | `0.1` |
| `--n-top` | genes per term in the percent bars | How many top genes per term go in the percent bars? | `5` |
| `--max-cancor` | collinearity ceiling | Above which canonical correlation between two terms should the run stop? | `0.9` |
| `--reference` | genome build and annotation | Which genome build and annotation release are the genes from (e.g. `GRCh38, GENCODE 44`)? | `synthetic demo genes, no genome build` |
| `--context` | biological context, printed first in the report and summary | Species, tissue or cell type, condition and dataset? | `synthetic demo: 12 donors, ctrl vs stim pseudobulk, no real tissue` |
| `--seed` | seed | Which seed? (Recorded; the lme4 fits draw no random numbers.) | `0` |
| `--outdir` | new output directory | Which new directory gets the outputs? It must not exist. | `variance_partition_demo` |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them.
With `--demo`, only `--outdir` may also be given; any other argument is an error.

```bash
python scripts/variance_partition.py \
  --counts results/2026-10-27_mono-de/tables/pseudobulk_counts.tsv \
  --meta results/2026-10-27_mono-de/tables/sample_meta.tsv \
  --formula "~ (1|donor) + (1|condition) + (1|sex) + age" \
  --min-frac 0.1 --n-top 5 --max-cancor 0.9 --reference "GRCh38, GENCODE 44" \
  --context "human PBMC CD14 monocytes, IFN-beta vs control, 8 donors" \
  --seed 0 --outdir results/2026-10-27_mono-varpart

# From cells
python scripts/variance_partition.py --adata data/mono.h5ad --layer counts \
  --gene-id-col gene_ids --sample-col sample --donor-col donor --min-cells 10 \
  --context "human PBMC CD14 monocytes, IFN-beta vs control, 8 donors" \
  --formula "~ (1|donor) + (1|condition)" --min-frac 0.1 --n-top 5 \
  --max-cancor 0.9 --reference "GRCh38, GENCODE 44" --seed 0 --outdir results/mono-vp

# Demo (only when the user asks for one)
python scripts/variance_partition.py --demo --outdir demo_out
```

## Workflow (fixed order)

1. **Validate** arguments and input; bad input gives `ERROR: <message>`, exit 2.
   The formula accepts only `name` and `(1|name)` terms joined by `+`, at most 23 (the largest colour-blind-safe palette).
   Unique gene and sample ids; fewer than 3 `--counts` samples stops the run.
2. **Pseudobulk** (`--adata` only): sum raw counts per `--sample-col`, then drop
   samples with fewer than `--min-cells` cells, as in pseudobulk-de. Fewer than
   3 samples left stops the run; then the donor-versus-sample check.
3. **Keep genes** detected (count > 0) in at least `--min-frac` of the kept
   samples' cells (`--adata`) or of the samples (`--counts`). The design rank of
   the plain terms is checked here.
4. **Collinearity**: `canCorPairs(formula, meta)`. Any pair above
   `--max-cancor` stops the run before the outdir exists; the error lists every
   pair's value. Every refusal leaves no outdir.
5. **Normalise**: `DGEList`, `calcNormFactors` (TMM), `voomWithDreamWeights`
   with the formula.
6. **Fit** `fitExtractVarPartModel(vobj, formula, meta)`, one mixed model per gene.
7. **Summarise**: median and mean fraction per term; terms ordered as
   `sortCols` does (median descending, Residuals last). Top genes per term: the
   `--n-top` largest fractions, ties by gene name.
8. **Draw** the three figures, then write report, result and reproducibility.

## Output structure

`variance_fractions.tsv` has one row per input gene, sorted by gene symbol:
`gene_symbol`, `gene_id`, the terms in formula order, Residuals, `status`. A
tested row sums to 1. `status` (house-style `data.status`) is `tested`, or
`not_detected` (zero everywhere) or `not_tested` (under `--min-frac`) with NA
fractions; `status_counts` in `result.json` counts each. The summary table and
figures use tested genes only. `summary.tsv` is sorted by median fraction
descending. `cancor.tsv` keeps formula order on both axes.

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # n_samples, n_donors, n_genes, status_counts, gene_id_col (--adata), dropped_samples, context, reference, top_covariate, median_fractions, max_cancor, top_genes, params, seed, chat_summary_lines, preferred_artifacts
├── tables/
│   ├── variance_fractions.tsv # gene_symbol, gene_id, one column per term, Residuals, status
│   ├── summary.tsv            # term, median_fraction, mean_fraction
│   └── cancor.tsv             # term x term canonical correlations
├── figures/
│   ├── variance_violin.pdf
│   ├── variance_violin.png
│   ├── variance_violin.csv
│   ├── percent_bars.pdf
│   ├── percent_bars.png
│   ├── percent_bars.csv
│   ├── cancor_heatmap.pdf
│   ├── cancor_heatmap.png
│   └── cancor_heatmap.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python packages; R, variancePartition, limma, edgeR, lme4
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

## Gotchas

- **The model will want to write a categorical covariate as a fixed effect**
  (`~ donor + condition`) because DE formulas look like that. Do not.
  variancePartition needs `(1|donor)`; the script refuses the fixed form.
- **The model will want to raise `--max-cancor` or drop a term** when the
  collinearity check stops the run. Do not. Report the pairs and values from
  the error message; which term to drop is the user's call.
- **The model will want to read a variance fraction as an effect size or a
  test.** Do not. A fraction is a share of variance, not an effect size, a fold
  change or a probability: a small, consistent shift can explain a large share.
  It has no p-value; DE between conditions goes to pseudobulk-de.
- **The model will want to say a term "drives" a gene's expression** because
  it explains a large fraction. Do not. A fraction is variance attributed to the
  term under this formula, an association, not a cause or a mechanism
  (house-style `claims.strength`).
- **The model will want to guess `--reference`** from gene symbols. Do not. Ask
  the user for the genome build and annotation release.
- **The model will want to compare fractions across runs with different
  formulas.** Do not. Fractions are shares of one model's total; adding a term
  changes every other share.
- **The model will want to pass log-normalised values or CPM** to `--counts`.
  Do not. The script needs raw summed counts; voom does the normalisation.
- **The model will want to fill in a missing argument** with a "standard" value
  (0.1, 10 cells, 0.9). Do not. Ask, or run `--demo` only when asked for one.

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
