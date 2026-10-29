---
name: sc-de
description: >-
  Per-feature differential expression in two modes: a MAST hurdle model on the
  log-normalised cells of one cell type (seeded per-donor downsampling, optional
  cellular detection rate covariate, optional (1|donor) random effect), or a
  per-feature linear model with covariates on sample-level values such as Olink
  NPX. Writes one sorted DE table, a volcano and a top-feature figure, a report
  and a reproducibility bundle. Fire on "MAST DE", "single-cell DE with MAST",
  "zlm hurdle model", "cell-level DE with CDR", "Olink NPX lm with covariates",
  "linear model per protein adjusted for age and sex". Do NOT fire for
  pseudobulk DE from raw counts (use pseudobulk-de), marker genes per cluster
  (use sc-markers), covariate-free two-group tests or correlations (use
  group-tests), or composition questions.
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - single-cell
    - differential expression
    - MAST
    - linear model
    - Olink
  inputs:
    - name: adata
      type: file
      format:
        - h5ad
      description: cell-level AnnData, log-normalised expression (mast mode)
      required: false
    - name: input
      type: file
      format:
        - tsv
      description: long table, one value per feature and sample (lm mode)
      required: false
    - name: gene_id_col
      type: string
      description: var column of Ensembl gene ids, or none (mast mode)
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
      description: counts, significant features, params, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - anndata
      - numpy
      - pandas
      - scipy
      - matplotlib
      - plotnine-extra
    r:
      - MAST
      - lme4
      - dplyr
      - purrr
  demo_data:
    - path: examples/demo.h5ad
      description: synthetic, 6 donors x ctrl/stim, T and Mono, 200 genes, 15 planted DE in T
    - path: examples/demo_npx.tsv
      description: synthetic NPX, 40 samples x 50 proteins, 8 planted status effects
  endpoints:
    cli: python scripts/sc_de.py --mode {mast|lm} ... --outdir {outdir}
  trigger_keywords:
    - MAST
    - zlm
    - hurdle model
    - cellular detection rate
    - Olink linear model
    - NPX lm
---

# Single-cell and proteomics DE

You are the dispatcher for **sc-de**. The script does the analysis; you collect
its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/sc_de.py`. Never re-implement the model, change a threshold, or
pick a value the user has not given. If a required argument is unknown, ask the
user (these are grill questions), or run `--demo` only when the user asks for a
demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
In `mast` mode, **state the unit of inference**: the p-values treat cells as
replicates, so they overstate the evidence when donors are the replication unit.
Point the user to `pseudobulk-de` for a donor-level test. If the script prints
`ERROR: ...` and exits 2, tell the user the message and the argument it names.
Do not work around it.

**Fire this skill when the user says any of:**
- "run MAST on the T cells, stim vs ctrl, with CDR"
- "cell-level DE with a donor random effect"
- "which Olink proteins differ by status, adjusted for age, sex and BMI"

**Do NOT fire when:** the data are raw counts to sum per sample (`pseudobulk-de`);
the user wants cluster markers (`sc-markers`); the test has no covariates
(`group-tests`); the question is about cell-type proportions.

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--adata` (mast) | `.h5ad` | log-normalised expression (`normalize_total` + `log1p`, or log2) in `X` or `--layer`; non-negative, not all integers | `examples/demo.h5ad` |
| `--gene-id-col` (mast) | var column | Ensembl gene ids (versioned as the source gives them); `var_names` are the symbols | `gene_ids`; `none` when the source has no ids |
| cell `obs` (mast) | columns | `--cell-type-col`, `--donor-col`, every `--formula` column | `cell_type`, `donor`, `condition` |
| `--input` (lm) | TSV, long | columns `sample`, `feature`, `value` (numeric) and every `--formula` column | `examples/demo_npx.tsv` |

Raw counts and scaled (negative) data are refused. A categorical column must
hold non-numeric labels; a numeric column enters as a continuous covariate.
In lm mode a row with a missing value is dropped and counted; one row per
(feature, sample). Mast: each symbol must map to one gene id
(`skillkit.check_gene_ids`, house-style `data.ids`), checked before anything
else; `--gene-id-col none` writes `gene_id` as NA, never a guessed id. Lm
features are proteins or assays named by the input, so lm takes no gene ids.

## Parameters (all required, no defaults)

Ask these through the analysis contract: `${CLAUDE_PLUGIN_ROOT}/skills/lucid-coder/references/contracts.md`.

Each row is a grill question. Ask the user; never fill one in.

| Argument | Grill question | Demo-only value (mast / lm) |
|---|---|---|
| `--mode` | Cells with MAST (`mast`) or samples with a linear model (`lm`)? | `mast` / `lm` |
| `--formula` | Which right-hand side, e.g. `~ condition` or `~ status + age + sex`? A `(1\|donor)` term (mast) fits a glmer random effect: correct for donors but slow (minutes per hundred genes). | `~ condition` / `~ status + age + sex` |
| `--contrast-term` (mast), `--term` (lm) | Which term do we test: `column,level,reference` for a factor, or a numeric `column`? | `condition,stim,ctrl` / `status,disease,healthy` |
| `--padj-max` | Which BH-adjusted p cutoff marks a feature significant? | `0.05` |
| `--lfc-min` | Which absolute effect (log2FC, or NPX for lm) must a significant feature exceed? | `0.5` |
| `--seed` | Which seed fixes downsampling and plot jitter? | `0` |
| `--context` | Species, tissue or cell type, condition and dataset? Printed first in the report and summary. | `synthetic demo: 6 donors, ctrl vs stim T cells, no real tissue` / `synthetic demo: 40 samples, disease vs healthy NPX, no real cohort` |
| `--reference` | Which genome build and annotation release name the features, e.g. `GRCh38, GENCODE 44` (or the Olink panel version)? Recorded in `result.json`; never guessed. | `synthetic demo, no genome` |
| `--adata` (mast) | Which h5ad holds the cells? | `examples/demo.h5ad` |
| `--gene-id-col` (mast) | Which var column holds the Ensembl gene ids (or `none` if the source has none)? | `gene_ids` |
| `--layer` (mast, optional) | Which layer holds the log-normalised values? Absent means `X`. | absent |
| `--log-base` (mast) | Is the expression natural log (`e`, scanpy `log1p`) or `2`? | `e` |
| `--cell-type-col`, `--cell-type` (mast) | Which obs column and which cell type? | `cell_type`, `T` |
| `--donor-col` (mast) | Which obs column names the donor? | `donor` |
| `--n-cells-per-donor` (mast) | Downsample to how many cells per donor (per term level), or `all`? | `50` |
| `--min-frac` (mast) | In what fraction of the kept cells must a gene be detected? | `0.1` |
| `--cdr` (mast) | Add the cellular detection rate (`cngeneson`) as a covariate, `yes` or `no`? | `yes` |
| `--input`, (lm) | Which long TSV holds the values? | `examples/demo_npx.tsv` |
| `--outdir` | Which new directory gets the outputs? It must not exist. | given by the user |

The last column is **demo-only, not defaults**. `--demo` alone pins the mast
values; `--demo --mode lm` pins the lm values. Beside `--demo` only `--mode` and
`--outdir` may be given.

## Workflow (fixed order)

1. **Validate** arguments: `--mode` decides which arguments it needs; an extra
   one is an error. The term column must appear in `--formula`.
2. **Read** the input; check the columns and the term (both levels present, or
   numeric). Mast: subset to `--cell-type` and refuse raw counts or negative
   values; divide by ln 2 when `--log-base e`, so every logFC is log2.
3. **Downsample** (mast): per `--donor-col` value (and per term level when the
   term is a factor), keep `--n-cells-per-donor` cells drawn with
   `numpy.random.default_rng(--seed)`; groups at or below the size keep all.
4. **CDR** (mast, `--cdr yes`): `cngeneson` = genes detected per kept cell over
   all genes, centred and scaled; `+ cngeneson` is appended to the formula.
5. **Filter genes** (mast): keep genes detected (> 0) in at least `--min-frac`
   of the kept cells.
6. **Fit** in R through Rscript. Mast: `zlm` (bayesglm with empirical Bayes; or
   `method = "glmer", ebayes = FALSE` when the formula has `|`), then
   `lrTest(fit, CoefficientHypothesis(term))`: stat and p are the hurdle
   (discrete + continuous) likelihood-ratio test, log2_fc and se come from
   `getLogFC`. Lm: `lm(value ~ formula)` per feature; estimate, se, t and p of
   the term coefficient. A factor term is relevelled so the reference is the
   given one. An unestimable term gives NA. Before the fit, cell ids must be
   unique and each gene symbol must map to one gene id (mast), and `skillkit.check_design_rank` runs on the
   fixed-effect columns of the fitted rows (the formula without `(1|x)`, plus
   `cngeneson`); an aliased design stops the run.
7. **Adjust**: BH across the features with a p-value; NA stays NA.
8. **Call** significance: `padj < --padj-max` and `|effect| > --lfc-min`.
   Sort by padj, then feature, NA last.
9. **Draw** the volcano and the top-padj feature: a box per term level, or a
   scatter with a least-squares line for a numeric term.

Lm mode runs in R `lm`, not Python statsmodels: R builds the design matrix
from the formula (factor contrasts, relevel) exactly as the corpus did, R is
already needed for MAST, and `lm` is a deterministic QR fit with no new Python
dependency.

## Output structure

`de_results.tsv` has one row per input feature: mast `gene_symbol,
gene_id, log2_fc, se, stat, pvalue, padj, significant, status`; lm
`feature, estimate, se, stat, pvalue, padj, significant, status`. Sorted by padj,
then feature, NA last. `status` (house-style `data.status`) is `tested` for a
feature with a p-value; an untested feature keeps NA statistics and says why:
`not_detected` (mast: zero in every kept cell), `not_tested` (mast: under
`--min-frac`; lm: every row dropped or the term not estimable).

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # mode, formula_fitted, n_cells + n_donors (mast) or n_samples (lm), dropped_by_level (lm rows lost to a missing value, per term level), context, reference, n_features, n_tested, status_counts, gene_id_col (mast), n_significant, significant_features, params, seed, chat_summary_lines, preferred_artifacts
├── tables/
│   └── de_results.tsv
├── figures/
│   ├── volcano.pdf
│   ├── volcano.png
│   ├── volcano.csv
│   ├── top_feature.pdf
│   ├── top_feature.png
│   └── top_feature.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python packages, Rscript path, R, MAST, lme4 versions
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

Figures come from the vendored `scripts/plots.py` (`plot_volcano`, `plot_box`,
`plot_scatter_fit`), saved by `save_panel` at 60 x 60 mm. No titles; captions
sit in `report.md`. In the mast box each point is one cell; the reference level
is blue, the other red, jitter ±0.2 box widths with `--seed`. `skillkit.py`,
`plots.py`, `theme.py` and `palettes.py` are byte copies of `templates/`.

## Gotchas

- **The model will want to run mast mode on donor-level conditions.** Do not.
  Cells are the unit only when the condition was assigned to each cell (a
  per-cell perturbation, such as a CRISPR guide). When the condition belongs to
  the donor or sample (disease, stimulation of a whole sample), the donor is the
  unit: run `pseudobulk-de`, or at least add `(1|donor)`.
- **The model will want to read a filtered gene as unexpressed.** Do not. A gene
  under `--min-frac` was not tested; non-detection is not absence.
- **The model will want to report MAST p-values as donor-level evidence.** Do
  not. Cells of one donor are not independent, so cell-level p-values overstate
  the evidence when donors are the replication unit. The demo shows it: G141 is
  null but called significant from donor noise. State the unit, and point to
  `pseudobulk-de` (or a `(1|donor)` formula) for a donor-level answer.
- **The model will want to drop `cngeneson` or the `(1|donor)` term** to get
  more hits or a faster run. Do not. Both are the user's answers; ask.
- **The model will want to pass raw counts or scaled data** to fix the
  `log-normalised` error. Do not. MAST models log expression; point `--layer` at
  the log-normalised values or ask the user where they are.
- **The model will want to call `FindMarkers(test.use = "MAST")` or `zlm`
  in-session** to check or extend a result. Do not. A second analysis has none
  of the fixed downsampling, filters or provenance.
- **The model will want to relax `--padj-max` or `--lfc-min`** after few hits.
  Do not; report the count as it is.
- **The model will want to read an lm coefficient as causal** or omit a
  covariate the user named. Do not. The estimate is adjusted only for the
  formula's covariates.
- **The model will want to say the condition "drives" or "regulates" a
  significant gene.** Do not. A DE call says the gene's expression is associated
  with the term in this design; without a perturbation it is not a mechanism
  (house-style `claims.strength`).

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
