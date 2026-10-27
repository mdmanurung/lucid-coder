# longitudinal: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Demo

`examples/make_demo.py` (seed 0) writes 24 subjects measured at visits 1, 2 and
3 with `age` and `sex`, 30 features, a subject baseline with sd 1 and noise with
sd 0.5. F01-F06 rise 0.75 per visit; F07-F29 are null. F30 has no subject
baseline and its noise is centred within each subject, so its subject variance
fits at 0. S24 misses visit 2 in every feature, so its figure line breaks there.

The `--demo` run (`time3`, factor, reference 1) fits 30 features, calls exactly
F01-F06 significant (estimates 1.38-1.65, planted 1.5), flags F30 as the one
singular fit and has no unconverged or failed fit. With `--time-type numeric
--term time` the same six are significant, with slopes near 0.75.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`.
Row IDs from `docs/clawbio-determinism/parts/C_disease.md` (lines
24, 30, 73, 86); `CATALOGUE.md:62` groups C-pl-10 and C-pl-23 as line plots.

| Row | Procedure | Source |
|---|---|---|
| C-an-16 | `lmerTest::lmer(formula, data)` per gene on vst values, tidy per gene; q-values with `qvalue` (here BH or Bonferroni) | `repos/ALTRA-manuscript/Analysis/scRNA/ALTRA_scRNA_R_helper_functions.r:374-455` |
| C-an-22 | CLR frequencies, then `lmerTest::lmer(clr ~ days_to_conversion + age + bmi + (1\|subject))` per cell type | `repos/ALTRA-manuscript/Analysis/scRNA/ALTRA_scRNA_R_helper_functions.r:124-180`; `repos/ALTRA-manuscript/Analysis/Flow_Cytometry/RA_06-flow_statistical_modeling_cluster_based_Longitudinal_converters.Rmd:18,72,110` |
| C-pl-10 | Subject trajectories over visits, `geom_line(aes(group = subject))` + points | `repos/multiple-myeloma/manuscript-figures/figure-4/fig_4g.ipynb:80-102`; `repos/ALTRA-manuscript/Analysis/scRNA/ALTRA_scRNA_R_helper_functions.r:710-751` |
| C-pl-23 | CLR frequency over time, line plus smooth (here the mean with a 95% CI) | `repos/ALTRA-manuscript/Figures/Figure1/Figure_S1B.ipynb:210` |

The corpus fits did not record singular fits or convergence; this skill adds
both. `parts/C_disease.md` has no IFN longitudinal mixed-model row.

The lmer handling follows two installed reference skills (under
`~/.claude/skills/synced/.../<name>/SKILL.md`):
`glmm-practical-guide` ("Convergence Troubleshooting", "Singular Fit": a
singular fit can be fine and is checked with `isSingular()`; keep a random
intercept only; Satterthwaite df with lmerTest for LMM fixed effects; at least
about 6 levels for a random factor) and `emmeans-statistical-testing`
("Workflow 4": `lmer.df = "satterthwaite"`). A factor-time `--term` is the same
treatment contrast, visit minus reference, that emmeans would give with
Satterthwaite df when no covariate interacts with time.

## ClawBio mechanisms used

IDs from `docs/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data | `examples/demo.tsv`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict and `skillkit`; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-11, F-13, F-14 | Fixed plot spec, style, Agg | vendored `theme.py`; matplotlib forms in the house theme |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/`, with R package versions |
| F-23, F-26 | Thresholds documented, closed option sets | parameter table; `--time-type`, `--adjust` choices |
| F-25 | "The model will want to ... Do not" | "Gotchas" |
| F-27 | One skill, one task | random-intercept LMM per feature only |
| F-28 | Golden-value demo test | `tests/test_longitudinal.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

From `skills/longitudinal`:

```bash
RSCRIPT=/share/software/statistical/R/R-4.5.2/gcc-8.5.0/R-4.5.2/bin/Rscript \
R_LIBS_USER=/exports/para-lipg-hpc/mdmanurung/R/4.5 \
LD_LIBRARY_PATH=/home/mdmanurung/.local/lib/conda-gcc-runtime \
UV_CACHE_DIR=/tmp/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with numpy,pandas,scipy,matplotlib,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

Without Rscript the R-dependent tests skip. Regenerate the demo with
`python examples/make_demo.py`; the golden values in the tests change with it.
Format with `ruff format --config ../../templates/ruff.toml` and
`air format scripts/lmer.R`.

## Reference skills

Calls and gotchas agree with `S/glmm-practical-guide` (random intercepts, few
subjects, singular fits), `S/emmeans-statistical-testing` and
`S/marginaleffects-inference` (contrasts on a fitted model) (`S` =
`~/.claude/skills/synced/f6b88190-5f7e-45de-8798-8d1f5b4ccea4_25160026-14a0-4c5c-af37-b6bb232e66d4`).

## Citations

- Bates, Mächler, Bolker, Walker (2015). Fitting linear mixed-effects models
  using lme4. J Stat Softw 67(1).
- Kuznetsova, Brockhoff, Christensen (2017). lmerTest package: tests in linear
  mixed effects models. J Stat Softw 82(13).
- Benjamini, Hochberg (1995). Controlling the false discovery rate. JRSS B
  57:289-300.
