---
name: mofa
description: >-
  Multi-omics factor analysis with MOFA2 (Python mofapy2) on per-view sample
  tables (pseudobulk RNA, ADT, Olink, chromVAR TF scores): a seeded CPU fit
  with a fixed number of factors, variance explained per factor and view,
  factor scores, feature weights, factor-metadata tests, a model file,
  figures, a report and a reproducibility bundle. Fire on "run MOFA",
  "MOFA2 on RNA plus protein", "multi-omics factors", "integrate views with
  factor analysis", "which factor separates disease", "variance explained per
  view". Do NOT fire for single-view PCA or NMF, cell-level integration
  (scVI, totalVI, WNN, Harmony), MEFISTO or other covariate-smoothed models,
  mofaflex, or for regression of factors on several covariates at once.
license: MIT
metadata:
  version: "0.1.0"
  domain: integration
  tags:
    - mofa
    - mofa2
    - mofapy2
    - multi-omics
    - factor analysis
  inputs:
    - name: context
      type: string
      description: species, tissue or cell type, condition and dataset
      required: true
    - name: view
      type: file
      format:
        - tsv
      description: one samples x features table per view, repeated
      required: true
    - name: metadata
      type: file
      format:
        - tsv
      description: one row per sample, a sample column and a group column
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
      description: factors above threshold, views they explain, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - mofapy2
      - h5py
      - numpy
      - pandas
      - scipy
      - matplotlib
      - plotnine-extra
  demo_data:
    - path: examples/demo_rna.tsv
      description: synthetic, 60 samples, 3 views, 2 shared and 1 view-specific factor
  endpoints:
    cli: python scripts/mofa.py --view {name}={path} ... --metadata {metadata} --outdir {outdir}
  trigger_keywords:
    - mofa
    - mofa2
    - multi-omics factor analysis
    - variance explained per view
    - factor weights
---

# MOFA

You are the dispatcher for **mofa**. The script does the analysis; you collect
its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/mofa.py`. Never re-implement the model, change a threshold, or
pick a value the user has not given. If a required argument is unknown, ask
the user (these are grill questions), or run `--demo` only when the user asks
for a demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

**Fire this skill when the user says any of:**
- "run MOFA on the pseudobulk RNA, ADT and Olink tables"
- "which views does each factor explain"
- "find a multi-omics factor that separates at-risk from healthy"
- "MOFA2 with 6 factors and scaled views"

**Do NOT fire when:**
- there is one view only (PCA or NMF does that);
- the integration is per cell (scVI, totalVI, WNN, Harmony);
- the user wants MEFISTO, smooth covariates, or the torch package mofaflex;
- the user wants factors regressed on group plus age (a model with covariates
  is a separate step on `factors.tsv`).

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--view name=path` | TSV, wide, repeated | first column sample ID; one numeric column per feature; empty or `NA` is missing | `examples/demo_rna.tsv` |
| `--metadata` | TSV | column `sample`, the `--group-col`; other columns are tested too | `examples/demo_meta.tsv` |

Feature names must be unique across views (prefix them, `rna_CD4`). Views are
joined on sample IDs: `skillkit.check_compatible` stops the run when no sample
is in every view, and `result.json` reports `n_shared_ids` (adopted 2026-10-28). IDs are
read as text, so `01` stays `01`. Gaussian views should be normalised and on a
log-like scale; the script does not transform them.

## Parameters (all required, no defaults)

Ask these through the analysis contract: `${CLAUDE_PLUGIN_ROOT}/skills/lucid-coder/references/contracts.md`.

Each row is a grill question. Ask the user; never fill one in.

| Argument | Grill question | Demo-only value |
|---|---|---|
| `--context` | Species, tissue or cell type, condition and dataset? It opens the report and the chat summary. | `synthetic demo: 60 samples, healthy vs disease, no real tissue` |
| `--view` | Which views, each as `name=path`? (repeat the flag) | `rna`, `adt`, `tf` from `examples/` |
| `--metadata` | Which TSV holds one row per sample? | `examples/demo_meta.tsv` |
| `--group-col` | Which metadata column holds the groups to colour and test? | `group` |
| `--likelihoods` | Per view, `gaussian` (continuous), `bernoulli` (0/1) or `poisson` (counts)? As `name=lik,...` | all `gaussian` |
| `--factors` | How many factors to train (>= 2)? | `6` |
| `--scale-views` | Scale each view to unit variance (`yes`), so large views do not dominate, or keep raw scales (`no`)? | `yes` |
| `--convergence` | ELBO tolerance: `fast`, `medium` or `slow`? | `fast` |
| `--min-var-explained` | Which R² (%) in a view counts as a factor explaining it? | `5` |
| `--n-top-features` | How many top weights per view go in each weight chart? | `10` |
| `--seed` | Which seed fixes the initialisation and the box jitter? | `0` |
| `--reference` | Which genome build and annotation release name the gene-level features (`GRCh38, GENCODE 44`)? Pass `none` when no view is gene-level. | `synthetic demo, no genome build` |
| `--outdir` | Which new directory gets the outputs? It must not exist. | given by the user |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them;
beside `--demo` only `--outdir` may be given.

`--reference` is required because the usual MOFA run has a gene-level RNA
view, and a gene name means nothing without its annotation release.

```bash
python scripts/mofa.py \
  --context "human PBMC pseudobulk, RNA + ADT + Olink, healthy vs SLE, cohort X" \
  --view rna=results/pb/rna_scaled.tsv --view adt=results/pb/adt.tsv \
  --view olink=results/olink/npx.tsv \
  --metadata results/meta/samples.tsv --group-col status \
  --likelihoods rna=gaussian,adt=gaussian,olink=gaussian \
  --factors 6 --scale-views yes --convergence slow \
  --min-var-explained 2 --n-top-features 10 --seed 1 \
  --reference 'GRCh38, GENCODE 44' \
  --outdir results/2026-10-27_mofa

# Demo (only when the user asks for one)
python scripts/mofa.py --demo --outdir demo_out
```

## Workflow (fixed order)

1. **Validate** arguments: ranges, `--view name=path` pairs, one likelihood per
   view.
2. **Read** each view: first column is the sample ID, every other cell numeric
   or empty (missing). Check values against the likelihood: bernoulli only 0/1,
   poisson only non-negative integers. `skillkit.check_unique` refuses a
   repeated feature name or sample ID within a view (read from the raw header,
   before pandas renames a repeat); a feature name in two views is refused too.
3. **Read** the metadata and keep the modelled samples, the sorted union over
   views. A sample missing from one view is left missing there; MOFA handles
   that. `check_unique` refuses a repeated metadata sample; a modelled sample
   absent from the metadata or without a group is refused. Metadata samples in
   no view are left out and counted per group (`dropped_per_group`).
   MOFA fits no design, so there is no design-rank check, and it joins no
   resource (gene sets, network), so there is no `n_matched` count.
4. **Train** with mofapy2's `entry_point`: views, samples and features sorted
   first, because the seeded initialisation maps onto input order;
   `set_data_options(scale_views=)`, `set_model_options(factors=)`,
   `set_train_options(convergence_mode=, seed=, dropR2=None, gpu_mode=False)`.
   BLAS threads are pinned to 1. `save()` writes `data/model.hdf5`, with
   factors sorted by variance explained.
5. **Extract** factor scores, weights and R² (percent, per factor and view)
   from the saved model.
6. **Flag** factors: a factor *explains* a view when its R² there is
   `>= --min-var-explained`; it is *active* when it explains at least one view.
   No factor is dropped.
7. **Test** every factor against every metadata column except `sample`:
   Spearman for a numeric column, Kruskal-Wallis otherwise; BH across all
   factor-covariate tests.
8. **Draw** the R² heatmap and the Factor1 vs Factor2 scatter, then a weight
   chart and a group box for each active factor.

## Output structure

- `variance_explained.tsv`: factor, view, r2 (%), explains, active; sorted by
  factor number, then view.
- `factors.tsv`: sample, Factor1..K; sorted by sample.
- `weights.tsv`: view, feature, factor, weight (unscaled); sorted by view,
  feature, factor number.
- `factor_vs_metadata.tsv`: factor, covariate, n, test (`spearman` or
  `kruskal-wallis`), statistic (rho or H), pvalue, padj (BH over all rows with
  a p-value); sorted by factor number, then covariate.

Weight and box figures exist for active factors only; the demo draws three of
each.

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # context, n_samples, views, n_factors, n_active, active_factors, r2_total, device, threads, dropped_per_group, n_shared_ids, params, seed, reference, chat_summary_lines, preferred_artifacts
├── tables/
│   ├── variance_explained.tsv
│   ├── factors.tsv
│   ├── weights.tsv
│   └── factor_vs_metadata.tsv
├── data/
│   └── model.hdf5             # mofapy2 model; load_model() in MOFA2 reads it
├── figures/
│   ├── variance_explained.pdf
│   ├── variance_explained.png
│   ├── variance_explained.csv
│   ├── factor1_vs_factor2.pdf
│   ├── factor1_vs_factor2.png
│   ├── factor1_vs_factor2.csv
│   ├── weights_Factor1.pdf
│   ├── weights_Factor1.png
│   ├── weights_Factor1.csv
│   ├── weights_Factor2.pdf
│   ├── weights_Factor2.png
│   ├── weights_Factor2.csv
│   ├── weights_Factor3.pdf
│   ├── weights_Factor3.png
│   ├── weights_Factor3.csv
│   ├── box_Factor1.pdf
│   ├── box_Factor1.png
│   ├── box_Factor1.csv
│   ├── box_Factor2.pdf
│   ├── box_Factor2.png
│   ├── box_Factor2.csv
│   ├── box_Factor3.pdf
│   ├── box_Factor3.png
│   └── box_Factor3.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python and package versions, device, thread pins
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

`data/model.hdf5` is not checksummed: HDF5 bytes are not stable.

## Gotchas

- **The model will want to drop factors below `--min-var-explained`** (MOFA2's
  `drop_factor_threshold`). Do not. Training keeps all `--factors`; the table
  flags them, and a dropped factor cannot be checked later.
- **The model will want to call Factor1 "the disease factor"** because it is
  first. Do not. Factors are ordered by variance explained, not by relevance;
  read `factor_vs_metadata.tsv`, and remember the sign of a factor is arbitrary
  (demo: disease scores lower on Factor1).
- **The model will want to change `--factors` or the seed** until a factor
  separates the groups. Do not. Both are fixed before the run; report what the
  fit gives.
- **The model will want to pass counts as `gaussian`.** Do not. Normalise and
  log first, or give `poisson`; ask the user which.
- **The model will want to read a weight as a marker** of the group. Do not. A
  weight ties a feature to a factor; only a factor's test ties it to a group.
- **The model will want to add age into the factor test** as in the corpus glm.
  Do not here. This skill tests one covariate at a time; a model with
  covariates is a separate step on `factors.tsv`.
- **The model will want to name a factor as a cell type or a mechanism**
  ("the T-cell activation factor"). Do not. A factor is a shared axis of
  variation; its weights and tests suggest a reading, and that reading is a
  hypothesis to state as one. A factor tied to the group is "associated with"
  it, never "driven by" it (house-style `claims.strength`).
- **The model will want to read a large R² as a significant factor.** Do not.
  Variance explained is not a test; only `factor_vs_metadata.tsv` padj is.
  A factor with high R² can tie to no covariate.
- **The model will want to read a sample missing from a view as zero** or as
  low. Do not. It was not measured; MOFA leaves it missing.
- **The model will want to guess `--reference`** from the gene names. Do not.
  Ask; pass `none` only when no view is gene-level.
- **The model will want to fill in a missing argument** ("MOFA's default is 15
  factors"). Do not. Ask, or run `--demo` only when the user asked for a demo.

## Determinism

mofapy2 seeds numpy with `--seed` when it initialises the model, and the script
sorts views, samples and features first, so row and column order of the input
do not change the fit (a test shuffles one view). BLAS threads are pinned to 1
(`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`) and
`gpu_mode=False`; both are recorded in `environment.txt` and `result.json`.
Reruns here were bit-identical, so the rerun test compares checksums. Another
BLAS build may move R² in the second decimal; the golden R² test allows ±0.5
percentage points for that.

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
