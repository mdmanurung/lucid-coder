---
name: longitudinal
description: >-
  Longitudinal mixed models per feature on repeated sample-level values (CLR
  cell-type frequencies, pseudobulk vst or log CPM, signature scores, Olink
  NPX): one lmerTest `lmer(value ~ time + covariates + (1|subject))` per
  feature, the named coefficient tested with Satterthwaite df, BH or Bonferroni
  across features, singular and unconverged fits flagged per feature, subject
  trajectory, forest and p-value figures, a report and a reproducibility bundle.
  Fire on "mixed model per gene over visits", "lmer per cell type with a random
  subject intercept", "does the CLR frequency change over time", "time effect
  across visits per feature", "longitudinal Olink with lme4". Do NOT fire for a
  two-visit paired test without covariates (group-tests), raw counts that need
  pseudobulk DE (pseudobulk-de, dreamlet), compositional models (crumblr,
  scCODA), random slopes or nested random effects, GLMMs, or survival models.
license: MIT
metadata:
  version: "0.1.0"
  domain: statistics
  tags:
    - lme4
    - lmerTest
    - mixed model
    - longitudinal
    - repeated measures
    - satterthwaite
  inputs:
    - name: context
      type: string
      description: species, tissue or cell type, condition and dataset
      required: true
    - name: input
      type: file
      format:
        - tsv
      description: long table, one value per feature, subject and time
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: parameters, fit counts and captioned figures
    - name: result
      type: file
      format:
        - json
      description: counts, significant, singular and unconverged features, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - numpy
      - pandas
      - scipy
      - matplotlib
      - plotnine-extra
    r:
      - lme4
      - lmerTest
      - reformulas
      - readr
      - dplyr
      - purrr
      - forcats
  demo_data:
    - path: examples/demo.tsv
      description: synthetic, 24 subjects x 3 visits, 30 features, 6 planted time effects, 1 singular fit
  endpoints:
    cli: python scripts/longitudinal.py --input {input} --formula ... --outdir {outdir}
  trigger_keywords:
    - mixed model per feature
    - lmer per gene
    - random intercept per subject
    - longitudinal visits
    - time effect over visits
---

# Longitudinal mixed models

You are the dispatcher for **longitudinal**. The scripts do the analysis; you
collect the arguments, run them, and explain the outputs.

## Agent boundary

Run `scripts/longitudinal.py`; it calls `scripts/lmer.R` through Rscript. Never
re-implement the model, change a threshold, change the formula or the reference
time, or pick a value the user has not given. If a required argument is unknown,
ask the user (these are grill questions), or run `--demo` only when the user asks
for a demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
Always say how many fits were singular or unconverged. If the script prints
`ERROR: ...` and exits 2, tell the user the message and the argument it names.
Do not work around it.

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--input` | TSV, long | columns `feature`, `subject`, `time`, `value` (numeric) | `examples/demo.tsv` |
| covariates | TSV columns | every other name in `--formula`; text columns become factors with levels in sorted order, numbers stay numeric | `age`, `sex` |

One row per feature, subject and time; a repeat is an error. A row with a
missing value in a column the formula uses is dropped and counted. Labels are
read as text, so `01` stays `01`. Values are already on the scale to model
(CLR, vst, log CPM, NPX); the script does not transform them.

Required checks (docs/pipeline-skills.md, adopted 2026-10-27):
- **Unique ids**: a repeated (feature, subject, time) row stops the run
  (`duplicate_rows`); in a long table that is the uniqueness check.
- **Design rank**: the formula's fixed terms (time and covariates, one row per
  subject and time) go through `skillkit.check_design_rank`; a covariate aliased
  with time, such as a batch that equals the visit, stops the run before R.
- **Donor versus sample**: the subject is the unit and `(1|subject)` is required,
  so repeated visits are never treated as independent.
- **No `--reference`**: features are generic (frequencies, scores, NPX), so no
  genome build is recorded. Gene-level input belongs in `pseudobulk-de`.
- **Drops per group**: rows dropped for a missing value are counted per time in
  `result.json` (`n_rows_dropped_by_time`).

## Parameters (all required, no defaults)

Each row is a grill question. Ask the user; never fill one in.

| Argument | Grill question | Demo-only value |
|---|---|---|
| `--context` | Species, tissue or cell type, condition and dataset? It opens the report and the chat summary. | `synthetic demo: 24 subjects, 3 visits, 30 features, no real tissue` |
| `--input` | Which TSV holds the values? | `examples/demo.tsv` |
| `--formula` | Which fixed covariates beside time? The form is `value ~ time + covariates + (1\|subject)`; `a:b` and `a*b` interactions allowed, no other random term. | `value ~ time + age + (1\|subject)` |
| `--time-type` | Is time a set of visits (`factor`, one coefficient per visit against the reference) or a number with a linear trend (`numeric`, one slope)? | `factor` |
| `--reference-time` | Which visit is the baseline? For `numeric`, time is shifted so this value is 0; the slope does not change. | `1` |
| `--value-unit` | What unit or scale is `value` on (e.g. `CLR`, `NPX, log2`)? It labels the y axis and the forest x axis. | `a.u.` |
| `--time-unit` | What unit is `time` in (e.g. `visit`, `week`, `day`)? It labels the x axis. | `visit` |
| `--term` | Which coefficient is tested? R names it `time<level>` for a visit (`time3` = visit 3 minus the reference) or `time` for a slope. | `time3` |
| `--adjust` | Correct across features with `bh` (FDR) or `bonferroni` (FWER)? | `bh` |
| `--padj-max` | Which adjusted-p cutoff marks a feature significant? | `0.05` |
| `--n-plot` | How many top features get a trajectory figure? | `4` |
| `--seed` | Which seed is recorded? lmer draws no random numbers, so it is set for provenance. | `0` |
| `--outdir` | Which new directory gets the outputs? It must not exist. | given by the user |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them;
beside `--demo` only `--outdir` may be given.

```bash
python scripts/longitudinal.py \
  --context "human PBMC, CD4 T cell frequencies, pre vs post vaccination, cohort X" \
  --input results/2026-10-27_freqs/clr_long.tsv \
  --formula "value ~ time + age + sex + (1|subject)" \
  --time-type factor --reference-time pre --term timepost \
  --value-unit CLR --time-unit visit \
  --adjust bh --padj-max 0.05 --n-plot 6 --seed 0 \
  --outdir results/2026-10-27_clr-lmer

# Demo (only when the user asks for one)
python scripts/longitudinal.py --demo --outdir demo_out
```

## Workflow (fixed order)

1. **Validate** arguments: value ranges; `--formula` has `value` on the left,
   `time` and exactly one `(1|subject)` on the right, and only column names,
   `:` or `*` otherwise. A numeric `--reference-time` for `--time-type numeric`.
2. **Read** the TSV: keep `feature`, `subject`, `time`, `value` and the formula's
   columns; drop rows with a missing value and count them per time; refuse a
   repeated (feature, subject, time), an absent factor reference time and a
   rank-deficient fixed design.
3. **Fit** each feature in `lmer.R`: `lmerTest::lmer(formula, data)` with REML
   (the lme4 default). Factor time has the reference as its first level; numeric
   time is `time - reference`. `--term` must be a fixed-effect column of the
   model matrix, or the run stops and lists the valid names. Each fit records
   estimate, SE, Satterthwaite df, t and p for `--term`, `isSingular()`, whether
   the optimizer converged, and every lme4 warning or message. A fit that errors
   keeps its row with NA statistics and the error text.
4. **Correct** across the features with a p-value (`--adjust`); NA stays NA.
   Singular and unconverged fits stay in the family.
5. **Call** significance: `padj < --padj-max`. Sort by padj, then feature, NA last.
6. **Draw** a trajectory figure for each of the `--n-plot` features with the
   lowest padj, a forest plot of the top 20, and the p-value histogram.

## Output structure

`model_results.tsv` has one row per feature: feature, estimate, se, df, t,
pvalue, padj, significant, singular, converged, n_obs, n_subjects, fit_message.
Figures are named by padj rank, `top01` to `top<n-plot>`; the demo draws four.

```
outdir/
├── report.md                  # params, fit counts, figures with padj captions
├── result.json                # context, term, n_features, n_fitted, n_significant, n_singular, n_not_converged, n_rows_dropped, n_rows_dropped_by_time, significant/singular/not_converged/failed_features, params, seed, chat_summary_lines, preferred_artifacts
├── tables/
│   └── model_results.tsv
├── figures/
│   ├── top01.pdf
│   ├── top01.png
│   ├── top01.csv
│   ├── top02.pdf
│   ├── top02.png
│   ├── top02.csv
│   ├── top03.pdf
│   ├── top03.png
│   ├── top03.csv
│   ├── top04.pdf
│   ├── top04.png
│   ├── top04.csv
│   ├── forest.pdf
│   ├── forest.png
│   ├── forest.csv
│   ├── pvalue_hist.pdf
│   ├── pvalue_hist.png
│   └── pvalue_hist.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python, R, lme4 and lmerTest versions
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

Figures are matplotlib in the vendored `theme.py`, saved by `save_panel` as PDF,
PNG and source CSV, with no title and no significance stars:
- `topNN` (60 x 60 mm): one grey line and points per subject over time, with the
  raw mean per time and its t-based 95% CI in red. Factor time sits on visit
  positions in level order. A line never bridges a missed visit: each subject is
  drawn on the full visit grid (every factor level, or every time seen for
  numeric time) with NaN at a gap, so the line breaks there. The caption counts
  subjects with a gap and subjects dropped for having no values. `topNN.csv`
  holds that drawn grid (subject, x, value); an empty value is a break.
- `forest` (89 mm wide): estimate and 95% CI (`estimate ± t(0.975, df) · se`)
  for the 20 features with the lowest padj, lowest at the top, dashed line at 0.
  Red rows pass `--padj-max`; hollow points are singular or unconverged fits.
- `pvalue_hist` (89 x 70 mm): raw p-values in 20 bins.

## Gotchas

- **The model will want to drop singular or unconverged features** to "clean"
  the table. Do not. They stay in the table, the BH family and the report
  counts; tell the user how many there are and let them decide.
- **The model will want to fit `lme4::lmer` or call `anova()`** for p-values.
  Do not. lme4 gives no df or p; the script uses lmerTest's Satterthwaite t-test
  on the one named `--term`.
- **The model will want to put subject in as a fixed effect** or run an
  unpaired test across visits. Do not. Repeated visits of one subject are not
  independent; `(1|subject)` is required.
- **The model will want to change the reference visit or `--term`** after seeing
  results, or try `numeric` when `factor` finds little. Do not. Both are the
  user's answers, fixed before the run.
- **The model will want to add random slopes or `(1|batch)`.** Do not. This
  skill fits a random subject intercept only; a richer model is a new question
  for the user.
- **The model will want to "fix" an unconverged fit** with another optimizer,
  `allFit()` or rescaled covariates inside this run. Do not. Report the flagged
  features; an optimizer check is a follow-up the user asks for.
- **The model will want to trust a random intercept over a handful of
  subjects.** Below about 6 subjects the subject variance is biased low and
  singular fits are common (glmm-practical-guide); say so when `n_subjects` is
  small.
- **The model will want to add significance stars** to a trajectory figure.
  Do not. padj and its family are in each caption.
- **The model will want to join a subject's points across a missed visit** to
  make the figure look complete. Do not. A bridged gap invents a trajectory;
  the line breaks and the caption counts the gaps.
- **The model will want to say time or treatment "drives" or "causes" a change.**
  Do not. The `--term` estimate is change over visits in this cohort, associated
  with time; without a control arm it is not a treatment effect (house-style
  `claims.strength`).
- **The model will want to read padj as the chance the feature changed.** Do not.
  A p-value is not a probability of the hypothesis; padj bounds an error rate.
- **The model will want to read the red mean line as the model fit.** It is the
  raw mean per time; the tested estimate adjusts for the covariates.

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
