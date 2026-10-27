# mashr-sharing: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Fixed plot spec

Plain matplotlib in the vendored theme, saved by `save_panel` as PDF, PNG and
source CSV. No titles; captions sit in `report.md`.

- `sharing_heatmap`: K x K sharing, sequential white to red fixed at 0-1, values
  printed in each cell. Rows and columns keep `--de` order (not clustered); the
  caption says so.
- `significant_hist`: bars of genes per number of significant cell types, 1..K.
- `raw_vs_posterior`: the pair with the highest sharing (ties go to the first
  pair in name order). Left, raw effects; right, posterior means. Red genes pass
  `--lfsr-max` in either type; the dashed line is y = x.

`scripts/skillkit.py`, `plots.py`, `theme.py` and `palettes.py` are byte copies
of `templates/`; a test fails when they drift.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`;
row IDs from `mycelium-extra/plans/clawbio-determinism/parts/*.md`.

| ID | What | Source |
|---|---|---|
| C-an-15 | `run_mash(fit, coef='Age')`; `get_lfsr < 0.05` per gene and per cell type; `get_pm`, `get_psd`, `get_estimated_pi` | `repos/single_cell_aging/05_Differential_gene_expression_analysis/01_DGE_analysis_Age_psychDisease_mashR_v2.ipynb:3523-3543` |
| E-an-09 | `run_mash(res.dl, coef='TB_statusCASE')`; counts of genes with lfsr < 0.05 | `repos_new/dreamlet_analysis/Nathan_NatImm_2021/Nathan_NatImm_2021.Rmd:553-580` |
| C-an-40 | `cor.test(logFC_discovery, logFC_replication, method="spearman")` across cell types | `repos/single_cell_aging/05_Differential_gene_expression_analysis/06_Validation_Age_DE_across_cell_types_v2.ipynb:321-328` |

The installed dreamlet skill documents `run_mash`, the E-an-09 path:
`/home/mdmanurung/.claude/skills/dreamlet/SKILL.md` and
`/home/mdmanurung/.claude/skills/dreamlet/references/mashr.md`. This skill matches
its canonical covariances, `se = logFC / t` input and NA masking of unmeasured
cells. It differs in two ways. `V` stays at identity, where `run_mash` runs
`mash_estimate_corr_em`. `canonical+ed` adds ED covariances, which `run_mash`
does not fit.

## ClawBio mechanisms used

IDs from `mycelium-extra/plans/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data | `examples/demo_*.tsv`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict, `skillkit.run`; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-11, F-13, F-14 | Fixed plot spec, thresholds drawn, style, Agg | "Fixed plot spec", vendored `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` to `set.seed` and `mash`; `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` via `skillkit`; R versions appended |
| F-23, F-26 | Thresholds documented, closed option sets | parameter table; `--align`, `--covariances` choices |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", frontmatter |
| F-27 | One skill, one task | sharing only; DE comes from pseudobulk-de tables |
| F-28 | Golden-value demo test | `tests/test_mashr_sharing.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

The demo data come from `examples/make_demo.py` (seed 0): 4 cell types, 800
genes. `G0000`-`G0079` share one effect (±1.2) in all four; `G0080`-`G0159`
are specific, 20 per cell type in table order (±1.5); the rest are null. With
the demo values, 154 strong genes train ED. 153 genes pass lfsr < 0.05
somewhere: all 80 shared genes in all 4 types, and 73 specific genes in exactly
their own type. No null gene passes. Pairwise sharing is 0.67-0.70.

The script finds Rscript from `RSCRIPT`, else `Rscript` on `PATH`. On this HPC:

```bash
export RSCRIPT=/share/software/statistical/R/R-4.5.2/gcc-8.5.0/R-4.5.2/bin/Rscript
export LD_LIBRARY_PATH=/home/mdmanurung/.local/lib/conda-gcc-runtime
export R_LIBS_USER=/exports/para-lipg-hpc/mdmanurung/R/4.5   # mashr 0.2.79, ashr 2.2.63
```

From `skills/mashr-sharing` (the demo takes about a minute; R tests skip
without Rscript):

```bash
UV_CACHE_DIR=$TMPDIR/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with numpy,pandas,scipy,matplotlib,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

Regenerate the demo with `python examples/make_demo.py`; the golden values in
the tests change with it.

## Reference skills

Calls and gotchas agree with the `dreamlet` skill
(`~/.claude/skills/dreamlet/references/mashr.md`: `run_mash`, `get_pm`,
`get_lfsr`, `get_significant_results`).

## Citations

- Urbut, Wang, Carbonetto, Stephens (2019). Flexible statistical methods for
  estimating and testing effects in genomic studies with multiple conditions.
  Nature Genetics 51:187-195.
- Hoffman et al. (2023). dreamlet: efficient and scalable differential
  expression analysis of single-cell data. bioRxiv.
