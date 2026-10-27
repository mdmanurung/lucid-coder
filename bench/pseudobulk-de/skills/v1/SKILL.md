---
name: pseudobulk-de
description: >-
  Pseudobulk differential expression for one cell type from a single-cell h5ad:
  sums raw counts per sample, filters samples and genes, runs DESeq2 (PyDESeq2 or
  R) with a design formula and contrast, and writes tables, five fixed figures,
  a report and a reproducibility bundle. Fire on "pseudobulk DE", "pseudo-bulk
  differential expression", "DESeq2 per cell type", "DE between conditions in
  T cells", "sum counts per donor/sample and test", "paired by patient DE".
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

## Trigger

**Fire this skill when the user says any of:**
- "pseudobulk DE for <cell type> between <conditions>"
- "run DESeq2 on pseudobulk counts per donor"
- "which genes change in monocytes after stimulation, paired by patient"
- "sum the counts per sample and test stim vs ctrl"

**Do NOT fire when:**
- the test is between cells (Wilcoxon, t-test, `rank_genes_groups`, markers);
- the data are bulk RNA-seq without cells (use a bulk DE skill);
- the question is about cell-type proportions or differential abundance.

## Why this exists

- **Without it**: the agent writes a fresh pseudobulk script each time. Filters,
  the design and the cutoffs drift from run to run, and cell-level tests creep in.
- **With it**: one tested script with fixed steps, fixed file names and fixed
  figures. Every analysis choice is a named argument in `commands.sh`.

## Workflow (fixed order)

1. **Validate** arguments and input; bad input gives `ERROR: <message>`, exit 2.
2. **Subset** to `--cell-type` in `--cell-type-col`.
3. **Sum** raw counts per `--sample-col` value.
4. **Drop** samples with fewer than `--min-cells` cells.
5. **Keep genes** detected (count > 0) in at least `--min-frac` of that cell
   type's cells in the samples kept in step 4.
6. **DESeq2** with `--design` and `--contrast` (`--engine python` = PyDESeq2,
   `--engine r` = DESeq2 via Rscript). Independent filtering uses
   `alpha = --padj-max` in both engines. No LFC shrinkage.
7. **Call significance**: `padj < --padj-max` and `|log2FC| > --lfc-min`.
   A gene with `padj` NA (independent filtering, outlier) is not significant.

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--adata` | `.h5ad` | raw counts (non-negative integers) in `X`, or in `--layer` | `examples/demo.h5ad` |
| cell-level `obs` | columns | `--cell-type-col`, `--sample-col`, every `--design` term | `cell_type`, `sample`, `patient`, `condition` |
| design terms | obs columns | constant within each sample (sample-level) | `patient`, `condition` |
| contrast factor | obs column | a `--design` term holding both contrast levels | `condition` with `stim`, `ctrl` |

Code a categorical factor as strings in `obs`. A numeric column enters the
design as a continuous covariate in both engines. Only additive designs
(`~ a + b`) are accepted.

## Parameters (all required, no defaults)

Each row is a grill question. Ask the user; never fill one in.

| Argument | Grill question | Demo-only value |
|---|---|---|
| `--adata` | Which h5ad holds the cells? | `examples/demo.h5ad` |
| `--layer` (optional) | Which layer holds the raw counts? Absent means `X`. | absent (`X`) |
| `--cell-type-col` | Which obs column holds the cell-type labels? | `cell_type` |
| `--cell-type` | Which cell type do we test? | `T` |
| `--sample-col` | Which obs column names the sample (one library per donor and condition)? | `sample` |
| `--design` | Which design formula? Is the contrast paired by patient, and which covariates enter? | `~ patient + condition` |
| `--contrast` | Which factor, numerator and denominator (`factor,numerator,denominator`)? | `condition,stim,ctrl` |
| `--min-cells` | How many cells must a sample have in this cell type to be kept? | `10` |
| `--min-frac` | In what fraction of this cell type's cells must a gene be detected? | `0.1` |
| `--padj-max` | Which adjusted-p cutoff marks a gene significant? | `0.05` |
| `--lfc-min` | Which absolute log2 fold change must a significant gene exceed? | `1.0` |
| `--seed` | Which seed? (Recorded; neither DESeq2 path draws random numbers.) | `0` |
| `--engine` | PyDESeq2 (`python`) or R DESeq2 (`r`)? | `python` |
| `--outdir` | Which new directory gets the outputs? It must not exist. | `pseudobulk_de_demo` |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them.
With `--demo`, only `--engine` and `--outdir` may also be given; any other
analysis argument is an error.

## CLI reference

```bash
python scripts/pseudobulk_de.py \
  --adata data/cells.h5ad --layer counts \
  --cell-type-col cell_type --cell-type "CD14 Mono" --sample-col sample \
  --design "~ patient + condition" --contrast condition,stim,ctrl \
  --min-cells 10 --min-frac 0.1 --padj-max 0.05 --lfc-min 1 \
  --seed 0 --engine python --outdir results/2026-10-27_mono-de

# Demo (only when the user asks for one)
python scripts/pseudobulk_de.py --demo --outdir demo_out
```

## Demo

The demo data come from `examples/make_demo.py` (seed 0): 7 patients, `ctrl` and
`stim` samples for P1-P6, one `ctrl` sample for P7 with 5 T cells, 300 genes.
Genes `G000`-`G029` are truly DE in T cells (`var["true_de"]`, log2FC ±1.5); the
last 40 genes are near-silent. With `--engine python` the demo keeps 12 samples
(P7_ctrl dropped), tests 250 genes and calls 28 significant, all truly DE.
`--engine r` gives the same 28.

## Output structure

`n_genes_tested` is the number of rows in `de_results.tsv`, the genes passing
step 5. Tables and results are sorted by `padj`, then gene, NA last.

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # n_samples, n_genes_tested, n_significant, top_genes, params, engine, seed, chat_summary_lines, preferred_artifacts
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
    ├── environment.txt        # Python and package versions; R and DESeq2 for --engine r
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

## Fixed plot spec

Figures come from the vendored `scripts/plots.py` (`plot_volcano`,
`plot_heatmap`) and plain matplotlib in the same theme, saved by `save_panel` as
PDF, PNG and source CSV. No titles; captions sit in `report.md`.

- `sample_pca`: log2 CPM, top 500 variable genes, numerator red, denominator blue.
- `pvalue_hist`: raw p-values, 20 bins of width 0.05.
- `ma`: log2FC against mean normalised count (log x), dashed at ±`--lfc-min`.
- `volcano`: genes with a `padj`, dashed at both cutoffs.
- `top_genes_heatmap`: 20 genes with the lowest `padj`, row z-scored log2 CPM,
  colours saturate at ±2 (`HEATMAP_LIMIT`, a display constant, not a threshold).

`scripts/plots.py`, `theme.py` and `palettes.py` are byte copies of
`templates/`; a test fails when they drift.

## Running the R engine

The script finds Rscript from `RSCRIPT`, else `Rscript` on `PATH`. On this HPC
the minimal working environment is:

```bash
export RSCRIPT=/share/software/statistical/R/R-4.5.2/gcc-8.5.0/R-4.5.2/bin/Rscript
export LD_LIBRARY_PATH=/home/mdmanurung/.local/lib/conda-gcc-runtime  # libstdc++ with GLIBCXX_3.4.29
export R_LIBS_USER=/exports/para-lipg-hpc/mdmanurung/R/4.5           # holds DESeq2 1.50.2
```

`R_HOME` is not needed. `module load statistical/R/4.5.2/gcc.8.5.0` sets the
first two; `~/.Renviron` sets `R_LIBS_USER` only when `HOME` is the user's.

## Gotchas

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
- **The model will want to feed normalised or log values** to "fix" a
  `counts ... must be non-negative integers` error. Do not. Point `--layer` at the
  raw counts or ask the user where they are.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`:

- Aggregation per sample and the 10% detection filter (`min_detection <- 0.1`):
  `repos/aifi-healthy-pbmc-reference/07-Pseudobulk/35-R_pseudobulk_seurat.ipynb:57`.
- Minimum cells and counts per pseudobulk sample (`min_cells=20, min_counts=1000`):
  `repos_new/spatial_skin_atlas/figs/3_Fig3/mintflow_tdc_degs.ipynb:516-523`.
- Additive DESeq2 design with covariates
  (`~ cohort.cohortGuid + subject.biologicalSex + subject.cmv`):
  `repos/aifi-healthy-pbmc-reference/08-Differential_expression/36-R_pseudobulk_deseq2.ipynb:200-213`.

## ClawBio mechanisms used

IDs from `docs/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data | `examples/demo.h5ad`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-11, F-13, F-14 | Fixed plot spec, thresholds drawn, style, Agg | vendored `plots.py`, `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` |
| F-23, F-26 | Thresholds documented, closed option sets | parameter table; `--engine` choices; no defaults |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", "Trigger" |
| F-28 | Golden-value demo test | `tests/test_pseudobulk_de.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Citations

- Love, Huber, Anders (2014). DESeq2. Genome Biology 15:550.
- Muzellec et al. (2023). PyDESeq2. Bioinformatics 39:btad547.

## Where this skill lives

This skill is installed at `.claude/skills/pseudobulk-de/`. Run its script from the working
directory as `python .claude/skills/pseudobulk-de/scripts/pseudobulk_de.py`.
