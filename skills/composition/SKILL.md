---
name: composition
description: >-
  Cell-type composition and differential abundance with crumblr: counts cells per
  sample and label (from a cell-level h5ad or a counts TSV), drops rare labels,
  runs crumblr's precision-weighted CLR and variancePartition dream (fixed and
  random effects, e.g. (1|donor)), and writes tables, a stacked composition bar,
  per-label proportion boxes, an effect forest, a report and a reproducibility
  bundle. Fire on "cell-type composition", "differential abundance of cell
  types", "do cell-type proportions change", "crumblr", "CLR of cell
  frequencies", "which populations expand after stimulation".
  Do NOT fire for neighbourhood-level DA on a KNN graph (Milo), gene expression
  DE (pseudobulk or cell-level), scCODA / Bayesian composition models, or a
  plain Wilcoxon on frequencies (group-tests).
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - single-cell
    - composition
    - differential abundance
    - crumblr
    - dream
  inputs:
    - name: context
      type: string
      description: species, tissue or cell type, condition and dataset
      required: true
    - name: adata
      type: file
      format:
        - h5ad
      description: cell-level AnnData; only obs is read
      required: false
    - name: counts
      type: file
      format:
        - tsv
      description: samples x labels cell counts, with a --meta sample table
      required: false
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
      description: labels tested, dropped, significant with sign, params, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - anndata
      - numpy
      - pandas
      - matplotlib
      - plotnine-extra
    r:
      - crumblr
      - variancePartition
      - limma
  demo_data:
    - path: examples/demo.h5ad
      description: synthetic, 8 donors x ctrl/stim, 7 cell types, obs only
  endpoints:
    cli: python scripts/composition.py --adata {adata} ... --outdir {outdir}
  trigger_keywords:
    - composition
    - differential abundance
    - cell-type proportions
    - crumblr
    - CLR frequencies
---

# Composition

You are the dispatcher for **composition**. The script does the analysis; you
collect its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/composition.py`. Never re-implement the analysis, change a
threshold, or pick a value the user has not given. If a required argument is
unknown, ask the user (these are grill questions), or run `--demo` only when the
user asks for a demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

**Do NOT fire when:**
- the question is DA of KNN neighbourhoods, not of labels (use milo-da);
- the question is gene expression (use pseudobulk-de or sc-de);
- the user asks for scCODA or another Bayesian composition model (later engine);
- the user wants one rank test per frequency without a model (use group-tests).

## Input contract

Give exactly one input mode.

| Input | Format | Required content | Example |
|---|---|---|---|
| `--adata` | `.h5ad` | cell-level `obs` with unique cell ids, `--sample-col`, `--label-col`, `--donor-col`, every `--formula` term and `--group-col`; `X` is not read | `examples/demo.h5ad` |
| `--counts` | `.tsv` | first column sample ID, one column per label, non-negative integer cell counts | `sample  B  Mono ...` |
| `--meta` | `.tsv` | first column the same unique sample IDs as `--counts`, `--donor-col`, every `--formula` term and `--group-col` | `sample  donor  condition` |

Formula terms, `--donor-col` and `--group-col` must be constant within a sample. A numeric
column enters as a continuous covariate; code a factor as strings. Only additive
formulas (`~ a + b`) with random intercepts `(1|x)` are accepted. The reference
level of a factor is its first value in alphabetical order, as in R, so the
coefficient is `<column><level>` (e.g. `conditionstim`).

## Parameters (all required, no defaults)

Ask these through the analysis contract: `${CLAUDE_PLUGIN_ROOT}/skills/lucid-coder/references/contracts.md`.

Each row is a grill question. Ask the user; never fill one in.

| Argument | Grill question | Demo-only value |
|---|---|---|
| `--adata` + `--sample-col` + `--label-col` | Which h5ad holds the cells, which obs column names the sample (one library), which the cell-type label? | `examples/demo.h5ad`, `sample`, `cell_type` |
| `--counts` + `--meta` | Or: which TSV holds samples x labels counts, which TSV the sample metadata? | not used |
| `--context` | Species, tissue or cell type, condition and dataset? | `synthetic demo: 8 donors, ctrl vs stim, 7 PBMC-like labels, no real tissue` |
| `--donor-col` | Which sample column names the donor (person or animal)? A donor with several samples needs a donor term in `--formula`. | `donor` |
| `--formula` | Which model? Does a donor give more than one sample (then `(1\|donor)`), and which covariates enter? | `~ condition + (1\|donor)` |
| `--coef` | Which coefficient is the effect of interest (model.matrix name)? | `conditionstim` |
| `--group-col` | Which sample column orders and colours the figures? | `condition` |
| `--tau` | Which Dirichlet-multinomial overdispersion for crumblr? (crumblr's own default is 1; the weights barely depend on it.) | `1` |
| `--min-count` | How many cells, summed over all samples, must a label have to be tested? | `50` |
| `--padj-max` | Which BH-adjusted p cutoff marks a label significant? | `0.05` |
| `--seed` | Which seed? (Recorded; crumblr and dream draw no random numbers; it seeds the box jitter of unpaired designs.) | `0` |
| `--outdir` | Which new directory gets the outputs? It must not exist. | `composition_demo` |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them;
only `--outdir` may be given beside it.

## Steps (fixed order)

1. **Validate** arguments and the input mode; bad input gives `ERROR: <message>`, exit 2.
   Cell ids (`--adata`) and sample ids (`--counts`, `--meta`) must be unique
   (`skillkit.check_unique`).
2. **Count** cells per sample and label (`pd.crosstab`), or read `--counts`;
   build one metadata row per sample from obs or `--meta`. Then check the design:
   the fixed-effect terms of `--formula` must give a full-rank design
   (`skillkit.check_design_rank`; batch confounded with condition stops the run),
   and a donor with several samples must appear in `--formula`, fixed or
   `(1|donor)` (`skillkit.check_replication`). `n_donors` and `n_samples` go to
   `result.json`.
3. **Drop** labels with fewer than `--min-count` cells over all samples; they are
   listed in `result.json` (`dropped_labels`) and the report, with the cells they
   remove per `--group-col` level (`dropped_cells_by_level`). No sample is
   dropped. Proportions and CLR use the kept labels only.
4. **Test** in `scripts/crumblr.R`: `crumblr(counts, tau = --tau)` (samples as
   rows), `dream(cobj, --formula, meta)`, `eBayes`, then
   `topTable(coef = --coef, number = Inf, confint = TRUE)`. A `--coef` not in the
   fit fails with the list of coefficients.
5. **Call significance**: `padj < --padj-max`, BH across the labels tested.
6. **Draw** the stacked composition bar, one proportion box per label and the
   effect forest.

Calls follow the `crumblr-compositional` reference skill (the reference skill
named for composition in `docs/pipeline-skills.md`)
(`~/.claude/skills/crumblr-compositional/references/getting-started.md:15-28,66-76`;
`integration-pseudobulk.md` for the counts orientation). `se` is `logFC / t`;
`ci_low`/`ci_high` are limma's `CI.L`/`CI.R`, the 95% interval on the moderated t.

**Later engines, not built.**
- Tree test (`--tree`): `buildClusterTreeFromPB(pb)` (dreamlet) needs a
  pseudobulk expression object that neither input mode carries, so
  `treeTest(fit, cobj, hcl, coef=)` (`tree-testing.md:9-21`) waits for an
  expression-aware input.
- scCODA (pertpy): Bayesian HMC, a reference cell type and an FDR on credible
  effects (corpus row C-an-23); not deterministic enough for this contract.

## Output structure

Tables are wide with one row per sample, sorted by sample. `da_results.tsv` is
sorted by padj, then label, and keeps every input label: `status` is `tested`, or
`not_tested` with NA statistics for a label under `--min-count`. Figures use the
tested labels only. Box figures are named `proportion_<label>` with
non-word characters replaced by `_`; the tree below lists the demo's six.

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # n_samples, n_donors, dropped_cells_by_level, n_labels_tested, dropped_labels, status_counts, n_significant, significant (label: up/down), paired_by, engine, params, seed, context, chat_summary_lines (first: "Context: ..."), preferred_artifacts
├── tables/
│   ├── counts.tsv             # sample x kept labels, cells
│   ├── proportions.tsv        # counts / sample total over kept labels
│   ├── clr.tsv                # crumblr's CLR values (cobj$E, pseudocount 0.5)
│   └── da_results.tsv         # label, logFC, se, t, ci_low, ci_high, pvalue, padj, significant, status
├── figures/
│   ├── composition.pdf
│   ├── composition.png
│   ├── composition.csv
│   ├── proportion_B.pdf
│   ├── proportion_B.png
│   ├── proportion_B.csv
│   ├── proportion_CD4_T.pdf
│   ├── proportion_CD4_T.png
│   ├── proportion_CD4_T.csv
│   ├── proportion_CD8_T.pdf
│   ├── proportion_CD8_T.png
│   ├── proportion_CD8_T.csv
│   ├── proportion_DC.pdf
│   ├── proportion_DC.png
│   ├── proportion_DC.csv
│   ├── proportion_Mono.pdf
│   ├── proportion_Mono.png
│   ├── proportion_Mono.csv
│   ├── proportion_NK.pdf
│   ├── proportion_NK.png
│   ├── proportion_NK.csv
│   ├── effect_forest.pdf
│   ├── effect_forest.png
│   └── effect_forest.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python packages; R, crumblr, variancePartition, limma, lme4
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

**Fixed plot spec.** No titles and no stars; captions with padj sit in
`report.md`.
- `composition`: `plots.plot_composition`, one stacked bar per sample, samples
  ordered by `--group-col` then sample, labels coloured by `categorical_palette` (colour-blind safe, up to 23).
- `proportion_<label>`: `plots.plot_box` of each sample's share per
  `--group-col` level (Tol blue, red, ...), y title "<label> (fraction of
  sample's cells)". Under 10 samples per level (the demo has 8) it draws points
  and the group median, no box, and the caption says so. When a random-effect term of
  `--formula` has exactly one sample per group level (paired donors), points sit
  on the box centre and a line joins each donor (`paired_by` in `result.json`);
  otherwise points jitter by ±0.2 with `--seed`.
- `effect_forest`: matplotlib in the house theme, logFC with its 95% CI per
  label, sorted by logFC; red when significant, dashed line at 0.

`scripts/plots.py`, `theme.py`, `palettes.py` and `skillkit.py` are byte copies
of `templates/`; a test fails when they drift.

## Gotchas

- **The model will want to run a Wilcoxon or t-test on proportions** to "check"
  crumblr. Do not. Proportions are compositional and share one denominator; the
  script's model is the result. A rank test per label is group-tests' job, and a
  different question.
- **The model will want to read a label's logFC, or a shift in its proportion,
  as an absolute change in cell numbers.** Do not. Proportions are compositional:
  they share one denominator, so one expanding label lowers every other label's
  share. The logFC is on the CLR scale, relative to the geometric mean of all
  kept labels. Say "relative abundance".
- **The model will want to say the condition recruits, expands or kills a cell
  type**, or that a shifted label drives a phenotype. Do not. crumblr gives a
  relative abundance associated with `--coef`; a mechanism needs a perturbation
  or other data (house-style `claims.strength`).
- **The model will want to treat each sample as independent** when a donor gave
  several (e.g. ctrl and stim). Do not. The script stops unless `--donor-col`
  appears in `--formula`; ask whether it enters as `(1|donor)` or fixed.
- **The model will want to lower `--min-count`** to keep a rare label it cares
  about. Do not. The cutoff is the user's answer; report the dropped labels and
  ask.
- **The model will want to drop `(1|donor)`** when the fit warns of a singular
  random effect or gives few hits. Do not. Repeated samples per donor are not
  independent; ask the user.
- **The model will want to guess `--coef`** from the condition values. Do not.
  Ask, or run once and read the coefficient list from the `--coef` error.
- **The model will want to pass percentages or normalised frequencies as
  `--counts`.** Do not. crumblr's precision weights need raw cell counts; the
  script refuses non-integers.
- **The model will want to fill in a missing argument** with a "standard" value
  (tau 1, 0.05, 10 cells). Do not. Ask, or run `--demo` only when the user asked
  for a demo.

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
