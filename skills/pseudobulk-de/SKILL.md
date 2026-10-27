---
name: pseudobulk-de
description: >-
  Pseudobulk differential expression for one cell type from a single-cell h5ad:
  sums raw counts per sample, filters samples and genes, runs DESeq2 (PyDESeq2 or
  R) or dreamlet (repeated measures, random effects) with a design formula and contrast, and writes tables, five fixed figures,
  a report and a reproducibility bundle. Fire on "pseudobulk DE", "pseudo-bulk
  differential expression", "DESeq2 per cell type", "DE between conditions in
  T cells", "sum counts per donor/sample and test", "paired by patient DE",
  "dreamlet", "repeated measures pseudobulk", "random effect for donor".
  Do NOT fire for cell-level DE (Wilcoxon or t-test between cells, marker genes,
  rank_genes_groups), bulk RNA-seq without cells (a genes x samples matrix
  already), or composition questions (cell-type proportions, differential
  abundance).
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - single-cell
    - pseudobulk
    - differential expression
    - DESeq2
    - PyDESeq2
    - dreamlet
  inputs:
    - name: adata
      type: file
      format:
        - h5ad
      description: cell-level AnnData with raw integer counts in X or a layer
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
      description: counts, top genes, params, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - anndata
      - numpy
      - pandas
      - scipy
      - matplotlib
      - pydeseq2
  demo_data:
    - path: examples/demo.h5ad
      description: synthetic, 7 patients x ctrl/stim, 2 cell types, 300 genes
  endpoints:
    cli: python scripts/pseudobulk_de.py --adata {adata} ... --outdir {outdir}
  trigger_keywords:
    - pseudobulk
    - pseudo-bulk
    - pseudobulk differential expression
    - DESeq2 per cell type
    - sum counts per sample
    - donor-level DE
    - dreamlet
    - repeated measures
---

# Pseudobulk DE

You are the dispatcher for **pseudobulk-de**. The script does the analysis; you
collect its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/pseudobulk_de.py`. Never re-implement the analysis, change a
threshold, or pick a value the user has not given. If a required argument is
unknown, ask the user (these are grill questions), or run `--demo` only when the
user asks for a demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--adata` | `.h5ad` | raw counts (non-negative integers) in `X`, or in `--layer` | `examples/demo.h5ad` |
| cell-level `obs` | columns | `--cell-type-col`, `--sample-col`, `--donor-col`, every `--design` term | `cell_type`, `sample`, `patient`, `condition` |
| design terms, donor | obs columns | constant within each sample (sample-level) | `patient`, `condition` |
| ids | `var_names`, `obs_names` | unique | `G000`, `cell0000` |
| contrast factor | obs column | a `--design` term holding both contrast levels | `condition` with `stim`, `ctrl` |

Code a categorical factor as strings in `obs`. A numeric column enters the
design as a continuous covariate in every engine. Only additive designs
(`~ a + b`) are accepted, plus random intercepts `(1|x)` for dreamlet
(`~ condition + (1|patient)`). The contrast factor must be a fixed term.

## Parameters (all required, no defaults)

Each row is a grill question. Ask the user; never fill one in.

| Argument | Grill question | Demo-only value |
|---|---|---|
| `--adata` | Which h5ad holds the cells? | `examples/demo.h5ad` |
| `--layer` (optional) | Which layer holds the raw counts? Absent means `X`. | absent (`X`) |
| `--cell-type-col` | Which obs column holds the cell-type labels? `none` when the input is already one population: all cells are tested, no subset. | `cell_type` |
| `--cell-type` | Which cell type do we test? With `--cell-type-col none`, the label the report gives that population. | `T` |
| `--sample-col` | Which obs column names the sample (one library per donor and condition)? | `sample` |
| `--donor-col` | Which obs column names the donor? A donor with several kept samples needs a donor term in `--design` (fixed or `(1|donor)`). | `patient` |
| `--design` | Which design formula? Does a subject give more than one sample (then `(1|subject)`, dreamlet), and which covariates enter? | `~ patient + condition`; with `--engine dreamlet`, `~ condition + (1|patient)` |
| `--contrast` | Which factor, numerator and denominator (`factor,numerator,denominator`)? | `condition,stim,ctrl` |
| `--min-cells` | How many cells must a sample have in this cell type to be kept? | `10` |
| `--min-frac` | In what fraction of this cell type's cells must a gene be detected? | `0.1` |
| `--padj-max` | Which adjusted-p cutoff marks a gene significant? | `0.05` |
| `--lfc-min` | Which absolute log2 fold change must a significant gene exceed? | `1.0` |
| `--seed` | Which seed? (Recorded; neither DESeq2 path draws random numbers.) | `0` |
| `--engine` | PyDESeq2 (`python`), R DESeq2 (`r`) or dreamlet (`dreamlet`)? Omit only with a `(1|x)` design, which selects dreamlet. | `python` |
| `--reference` | Which genome build and annotation release were the counts made with (e.g. `GRCh38, GENCODE 44`)? Never guess it. | `synthetic demo, no genome` |
| `--outdir` | Which new directory gets the outputs? It must not exist. | `pseudobulk_de_demo` |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them.
With `--demo`, only `--engine` and `--outdir` may also be given; any other
analysis argument is an error. `--demo --engine dreamlet` swaps the design to
`~ condition + (1|patient)`.

## Workflow (fixed order)

1. **Validate** arguments and input; bad input gives `ERROR: <message>`, exit 2.
2. **Subset** to `--cell-type` in `--cell-type-col`; skipped when `--cell-type-col none`.
3. **Sum** raw counts per `--sample-col` value.
4. **Drop** samples with fewer than `--min-cells` cells, recorded per contrast
   level (`dropped_by_level`).
5. **Check the design** on the kept samples: fixed-effect terms must not be
   aliased (`check_design_rank`), and a donor with several samples needs a
   donor term (`check_replication`). Either failure stops the run.
6. **Keep genes** detected (count > 0) in at least `--min-frac` of that cell
   type's cells in the samples kept in step 4. A filtered gene is not tested.
7. **Test** with `--design` and `--contrast`. `--engine python` = PyDESeq2,
   `--engine r` = DESeq2 via Rscript; both use independent filtering with
   `alpha = --padj-max` and no LFC shrinkage. `--engine dreamlet` = dreamlet via
   Rscript (`processAssays`, `dreamlet`, `topTable(coef=, number=Inf)`); it has
   no independent filtering, so `--padj-max` acts only in step 8.
8. **Call significance**: `padj < --padj-max` and `|log2FC| > --lfc-min`.
   A gene with `padj` NA (independent filtering, outlier) is not significant.

## Output structure

`n_genes_tested` is the number of rows in `de_results.tsv`, the genes passing
step 6. A gene the engine could not fit stays as a row with NA statistics. Tables and results are sorted by `padj`, then gene, NA last.

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # n_samples, n_donors, dropped_samples, dropped_by_level, reference, n_genes_tested, n_significant, top_genes, params, engine, seed, chat_summary_lines, preferred_artifacts
├── tables/
│   ├── pseudobulk_counts.tsv  # genes x kept samples, summed raw counts
│   ├── sample_meta.tsv        # one row per kept sample, design terms
│   └── de_results.tsv         # gene, base_mean, log2_fc, lfc_se, stat, pvalue, padj, significant
├── figures/
│   ├── sample_pca.pdf
│   ├── sample_pca.png
│   ├── sample_pca.csv
│   ├── pvalue_hist.pdf
│   ├── pvalue_hist.png
│   ├── pvalue_hist.csv
│   ├── ma.pdf
│   ├── ma.png
│   ├── ma.csv
│   ├── volcano.pdf
│   ├── volcano.png
│   ├── volcano.csv
│   ├── top_genes_heatmap.pdf
│   ├── top_genes_heatmap.png
│   └── top_genes_heatmap.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python and package versions; R and DESeq2 (r) or dreamlet, variancePartition, limma, edgeR, lme4 (dreamlet)
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

## Gotchas

- **The model will want to fake a cell-type column** (pass `--cell-type-col
  condition`, or write a tagged copy of the h5ad) when the input holds one
  population. Do not. Pass `--cell-type-col none` and name the population in
  `--cell-type`.
- **The model will want to re-run DE in-session** (a quick PyDESeq2 or
  `scanpy` call) to "check" or extend a result. Do not. The script's tables are
  the result; a second analysis has none of the fixed filters or provenance.
- **The model will want to drop the patient term** when a paired design gives
  few hits or a convergence warning. Do not. Dropping it changes the estimand
  from a within-patient to a between-group effect. Ask the user.
- **The model will want to relax `--padj-max` or `--lfc-min`** after seeing few
  significant genes. Do not. Thresholds are the user's answers, fixed before the
  run; report the count as it is.
- **The model will want to use a cell-level Wilcoxon test** (`rank_genes_groups`)
  because it is faster or gives more hits. Do not. Cells are not independent
  replicates; the sample is the unit of inference.
- **The model will want to fill in a missing argument** with a "standard" value
  (0.05, 10 cells, 10%). Do not. Every analysis argument is a grill question;
  ask, or run `--demo` only when the user asked for a demo.
- **The model will want to recount by hand** for a question about one
  intermediate step (samples kept, genes passing the filter). Do not: run the
  script and read `n_samples` and `n_genes_tested` from `result.json`. In one
  acumen run the agent counted detection over pseudobulk samples instead of
  cells and answered 11/9299; the expected answer was 11/456.
- **The model will want to fit a repeated-measures design with DESeq2** by
  passing `--engine python` or `r` with `(1|patient)`, or by turning the random
  term into a fixed one to dodge the error. Do not. Repeated measures default to
  dreamlet; a fixed patient term is the user's call, so ask.
- **The model will want to read a DE gene as a change inside each cell.** Do
  not. A pseudobulk change is a bulk change: a shift in subtype mix within the
  cell type gives the same signal. Run `composition` next before calling it
  cell-intrinsic.
- **The model will want to call a gene absent** because it is missing from
  `de_results.tsv`. Do not. Non-detection is not absence: a gene below
  `--min-frac` was filtered and is "not tested".
- **The model will want to drop `--donor-col` or the donor term** to get past
  the `more than one sample` error. Do not. Samples from one donor are not
  independent; ask the user for a fixed or `(1|donor)` term.
- **The model will want to guess `--reference`** from the gene names. Do not;
  ask which genome build and annotation release made the counts.
- **The model will want to feed normalised or log values** to "fix" a
  `counts ... must be non-negative integers` error. Do not. Point `--layer` at the
  raw counts or ask the user where they are.

## Engine rule: repeated measures use dreamlet

When a subject contributes more than one sample (repeated measures), dreamlet
is the default method. The script enforces this through the design:

- A random-effect term `(1|x)` in `--design` selects `--engine dreamlet` when
  `--engine` is omitted.
- `--engine python` or `--engine r` with a random-effect term fails:
  `--design has the random-effect term (1|x): repeated measures need --engine dreamlet, not --engine python`.
- A fixed-effect design (`~ patient + condition`) still needs an explicit `--engine`.

dreamlet fits a precision-weighted linear mixed model per gene (voom weights,
`variancePartition::dream`); see Hoffman et al. (2023) and the vignette
<https://diseaseneurogenomics.github.io/dreamlet/articles/dreamlet.html>. The
local `dreamlet` skill (`~/.claude/skills/dreamlet/SKILL.md`,
`references/api-reference.md`, `references/core-workflow.md`) documents the
same calls and the `~ condition + (1|donor_id)` form.

How `scripts/dreamlet.R` runs it:

- **Input object.** Our counts are already summed per sample. The script
  rebuilds a cell-level SingleCellExperiment that sums back to them (per sample:
  the pseudobulk vector plus `n_cells - 1` zero cells) and calls
  `aggregateToPseudoBulk`. That sets the cell counts per sample
  (`int_colData(pb)$n_cells`, read by `cellCounts`) and the metadata
  `processAssays` reads, through the public API.
- **Filters off.** Steps 4 and 5 stay authoritative, so `processAssays` runs
  with `min.cells = 0` (keep every sample), `min.count = 0`, `min.prop = 0` and
  `min.total.count = 0` (`filterByExpr` keeps every gene), and
  `min.samples = 2` (never binds; step 4 already needs 2 samples per level).
  edgeR still drops an all-zero gene; it stays in the table as an NA row.
  `dreamlet()` gets no `min.cells`: dreamlet 1.8.0 accepts it but does not use it.
- **Contrast.** The contrast factor is relevelled to the denominator, so the
  coefficient is `<factor><numerator>` (e.g. `conditionstim`).
- **Columns.** `base_mean` = `AveExpr` (mean log2 CPM, not a normalised count),
  `log2_fc` = `logFC`, `lfc_se` = `logFC / t`, `stat` = `t`,
  `pvalue` = `P.Value`, `padj` = `adj.P.Val` (BH).

## CLI reference

```bash
python scripts/pseudobulk_de.py \
  --adata data/cells.h5ad --layer counts \
  --cell-type-col cell_type --cell-type "CD14 Mono" --sample-col sample --donor-col patient \
  --design "~ patient + condition" --contrast condition,stim,ctrl \
  --min-cells 10 --min-frac 0.1 --padj-max 0.05 --lfc-min 1 \
  --seed 0 --engine python --reference "GRCh38, GENCODE 44" --outdir results/2026-10-27_mono-de

# Repeated measures: (1|patient) selects dreamlet
python scripts/pseudobulk_de.py ... --design "~ condition + (1|patient)" ...

# Demo (only when the user asks for one)
python scripts/pseudobulk_de.py --demo --outdir demo_out
```

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
