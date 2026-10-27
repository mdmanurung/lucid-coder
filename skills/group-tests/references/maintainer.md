# group-tests: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Trigger

**Fire this skill when the user says any of:**
- "test each cell-type frequency between healthy and disease"
- "paired Wilcoxon of the scores pre vs post, by donor"
- "which Olink proteins differ between groups, with FDR"
- "Spearman of each signature score against age"

**Do NOT fire when:**
- the values are gene counts that need pseudobulk DE (use `pseudobulk-de`);
- the test is between cells (markers, `rank_genes_groups`);
- the user wants a compositional model (scCODA, crumblr) or a regression with
  covariates;
- the user wants to combine p-values across cohorts.

## Why this exists

- **Without it**: each analysis writes its own loop. Paired data get unpaired
  tests, incomplete pairs slip in, the correction family changes, and stars on
  raw p reach the figure.
- **With it**: one tested script. Every choice is a named argument in
  `commands.sh`, every test has an effect size, and padj is in the table and the
  caption.

## Demo

`examples/make_demo.py` (seed 0) writes 20 donors (10 healthy, 10 disease), each
sampled `pre` and `post`, with 30 features. F01-F06 rise by 1 from pre to post,
F07-F12 are 2.5 higher in disease, F13-F18 rise 0.1 per year of age, and
F19-F30 are null. The `--demo` paired run tests 30 features and calls 7
significant: F01-F06 and F15 (an age feature, null for pre vs post), a
false discovery BH allows.

On the `pre` rows only (one value per donor; the tests build this file), BH at
0.05 gives:
- `wilcoxon` on `cohort` (`healthy,disease`): 7, F07-F12 plus F17. F17 is an age
  feature, and age differs by cohort by chance: a confounder, not an error.
- the same with `bonferroni`: 5, F07-F10 and F12.
- `spearman` and `pearson` on `age`: 6 each, exactly F13-F18.

## Fixed plot spec

Figures come from the vendored `scripts/plots.py`, saved by `save_panel` as
PDF, PNG and source CSV at 60 x 60 mm. No titles and **no significance stars**:
padj sits in `tests.tsv` and in each caption in `report.md`.

- Group tests: `plot_box`, box without outliers under one point per subject.
  Under 10 points in a group (`FEW_POINTS`), points and a median tick replace the
  box, and the caption says "Dot plot with medians".
  Unpaired points jitter by ±0.2 box widths (`--seed`); paired points sit on the
  box centre with a line per pair (complete pairs only). The reference level is
  blue (`#4477AA`), the other red (`#EE6677`).
- Correlation tests: `plot_scatter_fit`, one point per subject with a
  least-squares line and its 95% band, for both Spearman and Pearson.

The y axis reads `<feature> (<--value-unit>)`; tick labels carry `(n = ...)`. `scripts/plots.py`,
`theme.py`, `palettes.py` and `skillkit.py` are byte copies of `templates/`; a
test fails when they drift.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`.
Row IDs from `docs/clawbio-determinism/parts/*.md`.

| Row | Procedure | Source |
|---|---|---|
| B-an-16 | Wilcoxon on CLR frequencies per cell type, BH; median difference | `repos/IHA-Figure/Extended_Figure4/02_IMM_Age_Metric/01_Calculate_Scores.ipynb:143` |
| B-an-17 | Paired Wilcoxon by subject; raw p on plots | `repos/IHA-Figure/Extended_Figure9/01_FigA.ipynb:153` |
| B-an-19 | Fisher combination of BH p across cohorts (flagged smell; not done here) | `repos/IHA-Figure/Figure3/01_FigA.ipynb:121,134` |
| B-an-25 | Spearman of score against age per dataset and cell type | `repos/IHA-Figure/Figure2/07_FigI.ipynb:220` |
| C-an-20 | Paired / unpaired exact Wilcoxon, complete pairs only, FDR | `repos/multiple-myeloma/analyses-rna/04-pseudobulk/l3-pseudobulk/03-run-timepoint_freq_analysis.ipynb:58-87,110,125` |
| C-an-46 | Pearson and Spearman `cor.test` | `repos/IFN_manuscript/Figures/Supplementary/IFN_Sup5.ipynb:625-626,763-764` |
| C-an-47 | Wilcoxon of scores, paired pre/post | `repos/IFN_manuscript/Figures/Main/IFN_Fig4.ipynb:477,499,660,838` |
| C-an-48 | `p.adjust(method='fdr')`, q-value columns | `repos/multiple-myeloma/analyses-rna/04-pseudobulk/l3-pseudobulk/13-plots-deseq2-volcano_fixed.ipynb:79-80` |
| C-an-49 | Linear fit between two responses (here only as the drawn line) | `repos/IFN_manuscript/Figures/Main/IFN_Fig2.ipynb:1157-1160` |
| A-an-36 | `scipy.stats.mannwhitneyu(..., alternative="two-sided")` per gene on donor values | `repos_new/spatial_skin_atlas/figs/eFig8_tzone/e6f_e6g_ccl17boxplots.ipynb:318-335` |
| B-pl-12 | Box without outliers under jitter width 0.2 | `repos/IHA-Figure/Extended_Figure9/01_FigA.ipynb:150-153` |
| C-pl-05 | `geom_boxplot(outlier.shape=NA) + geom_jitter(position_jitter(0.2))` | `repos/ALTRA-manuscript/Figures/Figure1/Figure_1CD_S2DEF.ipynb:174-175` |
| C-pl-06 | Box with paired stats annotation (here padj in the caption) | `repos/IFN_manuscript/Figures/Main/IFN_Fig4.ipynb:660,838` |
| C-pl-09 | Scatter with `geom_smooth(method='lm')` | `repos/IFN_manuscript/Figures/Main/IFN_Fig5.ipynb:1011-1016` |

## ClawBio mechanisms used

IDs from `docs/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data | `examples/demo.tsv`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict and `skillkit`; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-11, F-13, F-14 | Fixed plot spec, style, Agg | vendored `plots.py`, `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` |
| F-23, F-26 | Thresholds documented, closed option sets | parameter table; `--test`, `--adjust` choices |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", "Trigger" |
| F-27 | One skill, one task | correlation and group tests only; no regression or meta-analysis |
| F-28 | Golden-value demo test | `tests/test_group_tests.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

From `skills/group-tests`:

```bash
UV_CACHE_DIR=/tmp/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with numpy,pandas,scipy,matplotlib,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

Regenerate the demo with `python examples/make_demo.py`; the golden values in
the tests change with it.

## Reference skills

Calls and gotchas agree with `S/glmm-practical-guide`,
`S/emmeans-statistical-testing` and `S/marginaleffects-inference` (`S` =
`~/.claude/skills/synced/f6b88190-5f7e-45de-8798-8d1f5b4ccea4_25160026-14a0-4c5c-af37-b6bb232e66d4`).
A model with covariates, random effects or contrasts belongs there or in
`longitudinal`, not here.

## Citations

- Hollander, Wolfe, Chicken (2013). Nonparametric Statistical Methods, 3rd ed.
  Wiley (Hodges-Lehmann estimators).
- Benjamini, Hochberg (1995). Controlling the false discovery rate. JRSS B
  57:289-300.
- Virtanen et al. (2020). SciPy 1.0. Nature Methods 17:261-272.
