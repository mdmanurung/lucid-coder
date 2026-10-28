# sc-normalize: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`.
Row IDs from `docs/clawbio-determinism/parts/*.md`
(`C_disease.md:12,42`, `E_methods.md:11,27,46`).

| Row | Step | Source |
|---|---|---|
| C-an-04 | SCTransform; Python side is an rpy2 wrapper of R sctransform | `repos/single_cell_aging/01_snRNA_seq_QC/sct_norm/sctransform.py:26-50`; `repos/ALTRA-manuscript/Analysis/TEA-Seq_analysis/RA_101_Initial_Object_Creation.ipynb:197-204` |
| C-an-34 | dsb on ADT via muon `pt.pp.dsb(..., isotype_controls=, random_state=)`, `dsb` layer into X | `repos/ALTRA-manuscript/Analysis/scRNA/Tfh17_validation_TEAseq/02_Teaseq_dsb_normalization.ipynb` (cells 81, 93, 99); `repos/ALTRA-manuscript/Figures/Figure4/Figure_S8CDEFG.ipynb:246-257` |
| E-an-03 | voom on pseudobulk (out of scope, pseudobulk-de) | `repos_new/dreamlet_analysis/PsychAD_r0/PsychAD_r0_analysis.Rmd:268-276` |
| E-an-19 | GTestimate (out of scope, R only) | `repos_new/GTestimate-Paper/notebooks/downstream_analysis_results.Rmd:45-71` |
| E-pl-04 | voom mean-variance trend plot | `repos_new/dreamlet_analysis/PsychAD_r0/PsychAD_r0_analysis.Rmd:279` |

Reference skills (docs/pipeline-skills.md table; `S` is
`~/.claude/skills/synced/f6b88190-5f7e-45de-8798-8d1f5b4ccea4_25160026-14a0-4c5c-af37-b6bb232e66d4`):
`S/single-cell-best-practices` (`SKILL.md:199-201`, shifted log via
`normalize_total` + `log1p`, and `sc.experimental.pp.normalize_pearson_residuals`,
the two RNA engines here); `single-cell-recipes` (`references/advanced_snippets.md`,
CITE GMM/DSB); `S/scvi-tools` for the later scVI engine, not offered yet.

## ClawBio mechanisms used

IDs from `docs/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" |
| F-04 | Bundled demo data | `examples/demo*.h5ad`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict and `skillkit`; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-13, F-14 | Fixed plot spec, style, Agg | vendored `plots.py`, `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` |
| F-23, F-26 | Parameters documented, closed option sets | parameter table; `--method` choices |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", "Agent boundary" |
| F-27 | One skill, one task; hand off through outputs | `data/normalized.h5ad` |
| F-28 | Golden-value demo test | `tests/test_sc_normalize.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

From `skills/sc-normalize`:

```bash
UV_CACHE_DIR=$TMPDIR/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with scanpy,muon,plotnine-extra,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

Demo golden values (`examples/make_demo.py`, seed 0):

- `--demo` (pearson): |Spearman rho| of PC1 with library size 0.2022, below
  the fixed 0.3; `--demo --method log1p` gives 0.9724, above 0.9.
- `--demo --method dsb`: CD3, CD19 and CD14 are above 3.5 in exactly their
  planted group (g1, g2, g3); CD16 and the isotypes stay below 3.5 everywhere.

Rebuild the demo with `python examples/make_demo.py` under the same
`uv run --with` line; the goldens in the test then need a re-check.

## Citations

- Lause, Berens, Kobak (2021). Analytic Pearson residuals for normalization of
  single-cell RNA-seq UMI data. Genome Biology 22:258.
- Choudhary, Satija (2022). Comparison and evaluation of statistical error models
  for scRNA-seq. Genome Biology 23:27.
- Hafemeister, Satija (2019). sctransform. Genome Biology 20:296.
- Mulè, Martins, Tsang (2022). Normalizing and denoising protein expression data
  from droplet-based single cell profiling. Nature Communications 13:2099.
- Bredikhin, Kats, Stegle (2022). MUON. Genome Biology 23:42.

## Biology-core checks (2026-10-28)

- `--context` is required and opens the report and `chat_summary_lines`.
- `gene_mean_variance.tsv` has a `status` column: an all-zero feature stays as
  `not_detected` with NA statistics; `result.json` has `status_counts`. Figures
  use tested features only.
- log1p/pearson take a required `--gene-id-col` (`none` writes NA `gene_id`);
  `check_gene_ids` runs before normalising. dsb does not take it: its features
  are antibodies, not genes. `examples/make_demo.py` now writes
  `var["gene_ids"]`; X, obs and var_names are unchanged.
- In an h5ad, `var_names` are unique, so the reachable ambiguous map is one id
  on two symbols (`duplicate gene ids`).
- The depth_pc1 caption says "tracks", not "drives"; a gotcha keeps |rho| at
  "associated with" (house-style `claims.strength`).
- `check_compatible` does not apply: the cells and empty droplets share
  features by a fixed-order check already.
