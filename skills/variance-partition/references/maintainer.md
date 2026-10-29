# variance-partition: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Fixed plot spec

Figures are matplotlib in the house theme, saved by `save_panel` as PDF, PNG
and source CSV. No titles; captions sit in `report.md`. Colours follow the
corpus (E-pl-01): `categorical_palette` for the terms in `sortCols` order (the corpus used Set1, which is not colour-blind safe), grey85
for Residuals.

- `variance_violin`: `plotVarPart` redrawn; one violin per term of percent
  variance explained over all genes, with a quartile box; y from 0 to 100.
- `percent_bars`: `plotPercentBars` redrawn; stacked horizontal bars for the top
  genes of each term, grouped by the term they top; a gene topping two terms is
  drawn once.
- `cancor_heatmap`: `plotCorrMatrix` redrawn with `plot_heatmap` (limit 1);
  canonical correlations lie in [0, 1], white to red. Rows and columns keep
  formula order (not clustered); the caption says so.

`scripts/skillkit.py`, `plots.py`, `theme.py` and `palettes.py` are byte copies
of `templates/`; a test fails when they drift.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`;
row IDs from `docs/clawbio-determinism/parts/*.md`.

| ID | What | Source |
|---|---|---|
| C-an-41 | covariate selection, then `fitExtractVarPartModel`, `plotVarPart` | `repos/single_cell_aging/04_selection_of_covariates/02_Selection_of_covariates_v2.ipynb:509,544` |
| E-an-06 | `fitVarPart(res.proc, form)`; form `~ (1|SubID) + (1|poolID) + (1|Sex) + scale(Age) + (1|Dx)` | `repos_new/dreamlet_analysis/PsychAD_r0/PsychAD_r0_analysis.Rmd:282-298` |
| E-an-16 | edgeR counts, `filterByExpr`, varpart per cell type | `repos_new/variancePartition_analysis/yazar_2022.Rmd:75-107` |
| E-an-07 | gene-level variance fractions as a table joined to gene features | `repos_new/dreamlet_analysis/PsychAD_r0/PsychAD_r0_analysis.Rmd:610-625` |
| D-an-29 | variance partition of program loadings (eta^2) | `I/script/FigureS1.R:130-275` |
| — | `canCorPairs(form, colData(pb))`, `plotCorrMatrix(C)` | `repos_new/dreamlet_analysis/Nathan_NatImm_2021/Nathan_NatImm_2021.Rmd:146-148` |
| — | `fitExtractVarPartModel(cobj, form, colData(pb))` | `repos_new/dreamlet_analysis/Nathan_NatImm_2021/Nathan_NatImm_2021.Rmd:495` |
| E-pl-01, C-pl-15 | `plotVarPart(sortCols(res.vp))`, Set1 + grey85 | `repos_new/dreamlet_analysis/PsychAD_r0/PsychAD_r0_analysis.Rmd:292-298` |
| E-pl-02 | `plotPercentBars(vp[i,])`, genes by `which.max` per column | `repos_new/variancePartition_analysis/yazar_2022.Rmd:136-140` |

API checked against the installed dreamlet skill
(`~/.claude/skills/dreamlet/references/core-workflow.md` "Variance
partitioning", `references/api-reference.md` "fitVarPart"): `fitVarPart` runs
`fitExtractVarPartModel` per cell type on `processAssays` output; this skill has
one matrix, so it calls `voomWithDreamWeights`, then `fitExtractVarPartModel`.

## ClawBio mechanisms used

IDs from `docs/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data | `examples/demo_*.tsv`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict, `skillkit.run`; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-11, F-13, F-14 | Fixed plot spec, thresholds drawn, style, Agg | vendored `plots.py`, `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json`, `set.seed` in R |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` via `skillkit` |
| F-23, F-26 | Thresholds documented, closed option sets | parameter table; formula grammar; no defaults |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", frontmatter |
| F-27 | One skill, one task | variance fractions only; DE hands off to pseudobulk-de |
| F-28 | Golden-value demo test | `tests/test_variance_partition.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

The demo data come from `examples/make_demo.py` (seed 0): 12 donors, `ctrl` and
`stim` each, 24 samples, 120 genes, negative binomial counts. `G000`-`G019` are
donor-driven, `G020`-`G039` condition-driven, `G040`-`G059` pure noise; the
rest carry a moderate donor effect. `rin` has no effect. The demo keeps all 24
samples and 120 genes. Donor has the largest median fraction (0.606), then
Residuals (0.200), rin, condition. The top 5 genes for donor, condition and
Residuals all come from their planted blocks. donor and rin correlate at 0.761,
under the demo `--max-cancor 0.9`, because rin varies per sample while each
donor has two samples.

The R engine needs (each R start loads variancePartition, about 2 minutes here):

```bash
export RSCRIPT=/share/software/statistical/R/R-4.5.2/gcc-8.5.0/R-4.5.2/bin/Rscript
export LD_LIBRARY_PATH=/home/mdmanurung/.local/lib/conda-gcc-runtime
export R_LIBS_USER=/exports/para-lipg-hpc/mdmanurung/R/4.5
```

From `skills/variance-partition`:

```bash
UV_CACHE_DIR=$TMPDIR/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with anndata,numpy,pandas,scipy,matplotlib,plotnine-extra,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

## Reference skills

Calls and gotchas agree with the `dreamlet` skill (`~/.claude/skills/dreamlet`:
`fitVarPart`, `plotVarPart`, `plotPercentBars`; dreamlet runs variancePartition
per cell type on voom weights, as this script does on one cell type).

## Biology-core checks (2026-10-28)

From `docs/pipeline-skills.md` "Biology-core checks". Applied:

- context: required `--context`, first line of the report and `chat_summary_lines`.
- status: `variance_fractions.tsv` keeps every input gene; genes under
  `--min-frac` get `not_detected` (zero everywhere) or `not_tested` with NA
  fractions. Medians, top genes and figures use tested genes only.
- gene ids: `--counts` must start with `gene_id`, `gene_symbol` (fixed ERROR);
  `demo_counts.tsv` was regenerated by `examples/make_demo.py` with the same
  counts. An all-empty `gene_id` column means no ids. `--adata` takes a required
  `--gene-id-col`. `check_gene_ids` replaces `check_unique` on gene names. The
  R side is unchanged: `run_r` still writes a symbol-keyed `counts.tsv`, and the
  `gene_id` column joins the fractions in Python.
- claim wording: one gotcha says a fraction is attribution, not a cause.

Not applied: compatible inputs (`--counts` and `--meta` are joined on sample
ids, already checked by `samples_mismatch`). The description keeps the user
trigger phrase "which covariates drive expression".

## Citations

- Hoffman, Schadt (2016). variancePartition. BMC Bioinformatics 17:483.
- Hoffman, Roussos (2021). dream. Bioinformatics 37:192-201.
- Law et al. (2014). voom. Genome Biology 15:R29.
