# mofa: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`.
Row ID from `docs/clawbio-determinism/parts/C_disease.md:40`.

| Row | Procedure | Source |
|---|---|---|
| C-an-32 | `create_mofa` on a long table of four views (rna, adt, olink, tf) | `repos/ALTRA-manuscript/Analysis/TEA-Seq_analysis/RA_105_MOFA_analysis_cd4naive.ipynb:1504` |
| C-an-32 | `data_opts$scale_views <- TRUE` | `repos/ALTRA-manuscript/Analysis/TEA-Seq_analysis/RA_105_MOFA_analysis_cd4naive.ipynb:1568` |
| C-an-32 | `model_opts$num_factors <- 6`, default likelihoods (gaussian) | `repos/ALTRA-manuscript/Analysis/TEA-Seq_analysis/RA_105_MOFA_analysis_cd4naive.ipynb:1671` |
| C-an-32 | `run_mofa(..., use_basilisk = TRUE)`, default training options (no seed set) | `repos/ALTRA-manuscript/Analysis/TEA-Seq_analysis/RA_105_MOFA_analysis_cd4naive.ipynb:1814` |
| C-an-32 | `get_factors`, joined to sample metadata | `repos/ALTRA-manuscript/Figures/Figure5/Figure_5BCDE_S11BC.ipynb:377` |
| C-an-32 | `r2_per_factor` table and `plot_variance_explained(x = "view", y = "factor")` | `repos/ALTRA-manuscript/Figures/Figure5/Figure_5BCDE_S11BC.ipynb:491,542` |
| C-an-32 | `glm(value ~ status + age)` per factor, BH (here: one covariate per test) | `repos/ALTRA-manuscript/Figures/Figure5/Figure_5BCDE_S11BC.ipynb:643,646` |
| C-an-32 | `ggboxplot` of Factor1 by status with jitter | `repos/ALTRA-manuscript/Figures/Figure5/Figure_5BCDE_S11BC.ipynb:690` |
| C-an-32 | `plot_weights(nfeatures = 10, scale = TRUE)` per view | `repos/ALTRA-manuscript/Figures/Figure5/Figure_5BCDE_S11BC.ipynb:800-821` |

Checked against the installed reference skill `mofaflex`
(`~/.claude/skills/mofaflex/SKILL.md`), the torch successor of mofapy2. It
seeds training through `TrainingOptions(seed=)`; this skill keeps mofapy2 and
passes `--seed` to `set_train_options(seed=)` the same way.

## ClawBio mechanisms used

IDs from `docs/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data | `examples/demo_*.tsv`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict and `skillkit`; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-11, F-13, F-14 | Fixed plot spec, style, Agg | vendored `plots.py`, `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/`, with thread pins |
| F-23, F-26 | Thresholds documented, closed option sets | parameter table; `--scale-views`, `--convergence`, likelihood choices |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", "Trigger" |
| F-27 | One skill, one task | factor model and one-covariate tests; no covariate regression |
| F-28 | Golden-value demo test | `tests/test_mofa.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

From `skills/mofa`:

```bash
UV_CACHE_DIR=$TMPDIR/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with mofapy2,h5py,numpy,pandas,scipy,matplotlib,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

Tested with mofapy2 0.7.5. Regenerate the demo with
`python examples/make_demo.py`; the golden values in the tests change with it.

### Demo

`examples/make_demo.py` (seed 0) writes 60 samples (30 healthy, 30 disease)
and three views: `rna` (200 features), `adt` (40) and `tf` (80). Three factors
are planted, each loading on a random half of a view's features over N(0, 1)
noise: A is shared by all views and is 1.5 higher in disease, B is shared by
rna and adt, and C is specific to tf. `age` is null.

The `--demo` run trains 6 factors; 3 are active at 5%:
- Factor1 (A) explains adt, rna and tf (R² 27.1, 34.0, 33.5%);
- Factor2 (B) explains adt and rna (25.9, 13.9%);
- Factor3 (C) explains tf (25.3%).

Factor4-6 stay below 2% in every view and are flagged, not dropped. Factor1 is
the only factor tied to group (Kruskal-Wallis padj 0.0002); no factor is tied
to age.

## Why this exists

- **Without it**: each analysis hand-scales views, leaves the seed to the
  clock, drops factors silently by an R² cutoff, and reads Factor 1 as "the
  disease factor" from a plot.
- **With it**: one tested script. Every choice is a named argument in
  `commands.sh`, the seed and threads are recorded, low factors are flagged
  but kept, and the group test is in a table with padj.

## Fixed plot spec

Figures are saved by `save_panel` as PDF, PNG and source CSV, with no titles
and no significance stars.

- `variance_explained`: factor x view heatmap, 60 x 60 mm, drawn with
  matplotlib, not `plots.plot_heatmap`. R² is never negative, so it uses the
  sequential `grey_to_red` map from `plots.py`, from 0 to the largest R², in
  place of the house diverging scale. Rows stay in factor order and views in
  name order, not clustered; the caption says so.
- `factor1_vs_factor2`: `plots.plot_embedding` with points enlarged, ticks
  restored and free aspect, coloured by group, 80 x 60 mm; `axis_name="Factor"`
  gives axes "Factor 1" and "Factor 2". The caption names the factor count,
  view scaling and seed.
- `weights_<factor>`: Cleveland dot chart (plots.py has none) of the
  `--n-top-features` largest absolute weights per view, coloured by view, with
  a dashed line at 0. Weights are not rescaled to [-1, 1], unlike the corpus
  `plot_weights(scale = TRUE)`.
- `box_<factor>`: `plots.plot_box`, one jittered point (±0.2, `--seed`) per
  sample, y titled "<factor> score (a.u.)"; under 10 samples in a group it
  draws points and a median, no box. The caption carries the Kruskal-Wallis
  padj.

Group and view colours come from one fixed Tol list keyed by the sorted labels;
more than 7 groups or views is an error. `scripts/plots.py`, `theme.py`,
`palettes.py` and `skillkit.py` are byte copies of `templates/`; a test fails
when they drift.

## Citations

- Argelaguet et al. (2018). Multi-Omics Factor Analysis. Mol Syst Biol
  14:e8124.
- Argelaguet et al. (2020). MOFA+. Genome Biology 21:111.
- Benjamini, Hochberg (1995). Controlling the false discovery rate. JRSS B
  57:289-300.
