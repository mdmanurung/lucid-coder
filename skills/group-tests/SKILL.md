---
name: group-tests
description: >-
  Group-comparison and association tests on sample-level values (cell-type
  frequencies or CLR values, signature scores, Olink NPX): one Wilcoxon
  rank-sum, paired Wilcoxon signed-rank, Spearman or Pearson test per feature,
  BH or Bonferroni across features, an effect size for every test, a table,
  box or scatter figures for the top features, a report and a reproducibility
  bundle. Fire on "Wilcoxon per cell type", "compare frequencies between
  groups", "paired test pre vs post per donor", "Olink NPX healthy vs disease",
  "correlate scores with age", "Spearman per feature with FDR". Do NOT fire for
  gene-level differential expression from counts (pseudobulk DE), cell-level
  tests (marker genes, rank_genes_groups), compositional models (scCODA,
  crumblr, propeller), regression with covariates, or combining p-values across
  cohorts.
license: MIT
metadata:
  version: "0.1.0"
  domain: statistics
  tags:
    - wilcoxon
    - mann-whitney
    - signed-rank
    - spearman
    - pearson
    - multiple testing
  inputs:
    - name: context
      type: string
      description: species, tissue or cell type, condition and dataset
      required: true
    - name: input
      type: file
      format:
        - tsv
      description: long table, one value per feature and subject
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
      - numpy
      - pandas
      - scipy
      - matplotlib
      - plotnine-extra
  demo_data:
    - path: examples/demo.tsv
      description: synthetic, 20 donors x pre/post, 30 features with planted effects
  endpoints:
    cli: python scripts/group_tests.py --input {input} --test ... --outdir {outdir}
  trigger_keywords:
    - wilcoxon per feature
    - paired wilcoxon
    - mann-whitney per cell type
    - frequency comparison
    - olink npx comparison
    - spearman correlation per feature
---

# Group tests

You are the dispatcher for **group-tests**. The script does the analysis; you
collect its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/group_tests.py`. Never re-implement the tests, change a threshold,
or pick a value the user has not given. If a required argument is unknown, ask
the user (these are grill questions), or run `--demo` only when the user asks
for a demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--input` | TSV, long | columns `feature`, `subject`, `value` (numeric) | `examples/demo.tsv` |
| group column | TSV column | `--group-col`, holding both `--levels` | `visit` with `pre`, `post` |
| pair column | TSV column | `--pair-col`, one sample per pair and level | `donor` |
| covariate | TSV column | `--covariate`, numeric | `age` |

One row per feature and subject. A subject is one sample: in paired data each
donor has two subjects (`D01_pre`, `D01_post`) that share a pair id. Labels are
read as text, so `01` stays `01`.

Required checks (docs/pipeline-skills.md, adopted 2026-10-27):
- **Unique ids**: a repeated (feature, subject) or (feature, group, pair) row stops
  the run (`duplicate_rows`); in a long table that is the uniqueness check.
- **Donor versus sample**: `paired-wilcoxon` pairs samples by `--pair-col`. An
  unpaired test assumes one subject per donor; ask before running it on repeats.
- **No design**: the tests fit no design matrix, so there is no rank check.
- **No `--reference`**: features are generic (frequencies, scores, NPX), not genes,
  so no genome build is recorded. Gene-level input belongs in `pseudobulk-de`.
- **Drops per group**: rows dropped for a missing value are counted per level in
  `result.json` (`n_missing_dropped_by_group`; `all` for correlation tests).

## Parameters (all required, no defaults)

Ask these through the analysis contract: `${CLAUDE_PLUGIN_ROOT}/skills/lucid-coder/references/contracts.md`.

Each row is a grill question. Ask the user; never fill one in.

| Argument | Grill question | Demo-only value |
|---|---|---|
| `--context` | Species, tissue or cell type, condition and dataset? It opens the report and the chat summary. | `synthetic demo: 20 donors, pre vs post, 30 features, no real tissue` |
| `--input` | Which TSV holds the values? | `examples/demo.tsv` |
| `--test` | Two groups (`wilcoxon`), paired samples (`paired-wilcoxon`), or a trend with a continuous covariate (`spearman` rank, `pearson` linear)? | `paired-wilcoxon` |
| `--group-col` | Which column holds the groups? (group tests only) | `visit` |
| `--levels` | Which two levels, reference first (`a,b`)? | `pre,post` |
| `--pair-col` | Which column pairs the samples, e.g. donor? (paired only) | `donor` |
| `--covariate` | Which numeric column do we correlate with? (correlation only) | absent |
| `--value-unit` | What unit or scale is `value` on (e.g. `CLR`, `NPX, log2`, `% of CD4 T cells`)? It labels the y axis. | `a.u.` |
| `--adjust` | Correct across features with `bh` (FDR) or `bonferroni` (FWER)? | `bh` |
| `--alpha` | Which adjusted-p cutoff marks a feature significant? | `0.05` |
| `--n-plot` | How many top features get a figure? | `4` |
| `--seed` | Which seed fixes the box-plot jitter? | `0` |
| `--outdir` | Which new directory gets the outputs? It must not exist. | given by the user |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them;
beside `--demo` only `--outdir` may be given.

## Workflow (fixed order)

1. **Validate** arguments: `--test` decides which of `--group-col`, `--levels`,
   `--pair-col`, `--covariate` it needs; an extra one is an error.
2. **Read** the TSV and keep the rows the test uses: both `--levels` of
   `--group-col`, or a numeric `--covariate`. Drop rows with a missing value and
   count them per level. Refuse a repeated (feature, subject) or (feature, group, pair).
3. **Test** each feature:
   - `wilcoxon`: two-sided Mann-Whitney U of b against a
     (`scipy.stats.mannwhitneyu`, `method="auto"`); effect = Hodges-Lehmann
     shift, the median of all b - a differences.
   - `paired-wilcoxon`: complete pairs only; signed-rank test of b - a
     (`scipy.stats.wilcoxon`, `zero_method="wilcox"`, `method="auto"`); effect =
     median of the Walsh averages of b - a.
   - `spearman` / `pearson`: `scipy.stats.spearmanr` / `pearsonr` of value on
     the covariate; effect = rho / r.
   A feature the test cannot run on (an empty group, no non-zero paired
   difference, fewer than 3 points or a constant column) gets an NA p-value.
4. **Correct** across the features with a p-value (`--adjust`); NA stays NA.
5. **Call** significance: `padj < --alpha`. Sort by padj, then feature, NA last.
6. **Draw** the `--n-plot` features with the lowest padj.

Levels are `a,b`: a is the reference, so every shift is b minus a.

## CLI reference

```bash
python scripts/group_tests.py \
  --context "human PBMC, CD4 T cell frequencies, pre vs post vaccination, cohort X" \
  --input results/2026-10-27_freqs/clr_long.tsv \
  --test paired-wilcoxon --group-col visit --levels pre,post --pair-col donor \
  --value-unit CLR --adjust bh --alpha 0.05 --n-plot 6 --seed 0 \
  --outdir results/2026-10-27_freq-tests

python scripts/group_tests.py --context "human serum, Olink NPX vs age, cohort X" --input scores.tsv --test spearman --covariate age --value-unit "score, a.u." \
  --adjust bh --alpha 0.05 --n-plot 6 --seed 0 --outdir results/score-age

# Demo (only when the user asks for one)
python scripts/group_tests.py --demo --outdir demo_out
```

## Output structure

`tests.tsv` has one row per feature, sorted by padj, then feature, NA last.
Its columns depend on `--test`:
- `wilcoxon`: feature, n_a, n_b, median_a, median_b, hodges_lehmann,
  statistic (U of b), pvalue, padj, significant;
- `paired-wilcoxon`: feature, n_pairs, median_a, median_b, hodges_lehmann,
  statistic (scipy's two-sided signed-rank statistic), pvalue, padj, significant;
- `spearman` / `pearson`: feature, n, rho / r, pvalue, padj, significant.

Here `a` and `b` are the level names (`median_pre`, `median_post`).
Figures are named by padj rank, `top01` to `top<n-plot>`; the demo draws four.

```
outdir/
├── report.md                  # params, summary, figures with padj captions
├── result.json                # context, test, n_features, n_tested, n_significant, n_rows_dropped, n_missing_dropped_by_group, significant_features, params, seed, chat_summary_lines, preferred_artifacts
├── tables/
│   └── tests.tsv
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
│   └── top04.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python and package versions
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

## Gotchas

- **The model will want to subtract baseline** and test the change scores between
  groups. Do not, in a parallel-group design: model the follow-up value with
  baseline as a covariate (Harrell, R Workflow "analysis": "it is not
  appropriate to compute change from baseline").
- **The model will want to use `wilcoxon` on paired samples** because the
  paired run "loses" incomplete donors. Do not. Paired samples are not
  independent; ask whether the design is paired and pass `--pair-col`.
- **The model will want to add significance stars** (`*`, `**`, `ns`) to a
  figure. Do not. Stars hide the correction family and are often drawn on raw p
  in the corpus; padj with its family is in the caption.
- **The model will want to relax `--alpha` or switch `--adjust`** after seeing
  few hits. Do not. Both are the user's answers, fixed before the run; report
  the count as it is.
- **The model will want to combine p-values across cohorts** (Fisher, `sumlog`)
  on adjusted p. Do not. That is a separate meta-analysis question; run one
  cohort per call and ask the user.
- **The model will want to read a correlation as an effect of the covariate.**
  Do not. rho and r are associations; a group difference can also come from a
  confounder, as F17 shows in the demo.
- **The model will want to say a group "drives" or "causes" a feature change.**
  Do not. A significant test is an association with the group label; say
  "associated with" (house-style `claims.strength`). Groups differ in more than
  the label.
- **The model will want to read an NA p-value as "no difference".** Do not. NA
  means the feature was not tested (empty group, no non-zero paired difference,
  under 3 points); say "not tested".
- **The model will want to read padj as the chance the feature changed.** Do not.
  A p-value is not a probability of the hypothesis; padj bounds an error rate.
- **The model will want to fill in a missing argument** with a "standard" value
  (0.05, BH, the first two levels). Do not. Ask, or run `--demo` only when the
  user asked for a demo.

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
