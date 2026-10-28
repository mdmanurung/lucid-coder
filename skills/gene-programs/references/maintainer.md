# gene-programs: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`;
row IDs from `docs/clawbio-determinism/parts/*.md`.

| ID | What | Source |
|---|---|---|
| C-an-29 | NMF gene programs, usage tables, starCAT and Spectra | `repos/ALTRA-manuscript/Analysis/scRNA/NMF_Analysis/RA_03_NMF_Analysis_Non_Naive_CD4_Tcells.ipynb:84,158-159`; `repos/multiple-myeloma/analyses-rna/03-subset-analyses/02-t_cells/03-pbmc-mm-starcat-tfh.ipynb:121,343` |
| B-an-28 | NMF projection of CD4 T programs, top-20 genes per factor | `repos/IHA-Figure/Figure5/01_NMF_CD4_T_cell_Projection/01_NMF_Score.ipynb:130,250` |
| D-an-20 | GP discovery (flashier EBMF), L cells × GP, F genes × GP | `repos/immgenT-GP-analysis/analysis/Methods_FlashierFit.Rmd:60-110` |
| D-an-22 | NMF/topic-model fit with a fixed seed, structure plots | `repos/immgenT-GP-analysis/analysis/old/t_cell_nmf_k20.Rmd:24,144-187` |
| D-an-23 | One-vs-rest AUC of GP loading per group | `repos/immgenT-GP-analysis/code/R/roc_auc.R:1-60`; `repos/immgenT-GP-analysis/code/pipeline/02_compute_auc.R:1-46` |
| D-an-24 | Active cell: loading / per-GP max > 0.1 | `repos/immgenT-GP-analysis/script/Figure2.R:61-62`; `repos/immgenT-GP-analysis/script/Figure1.R:115,125` |
| D-an-21, D-an-25 | Fixed-loading projection; activation shift (out of scope here) | `repos/immgenT-GP-analysis/code/pipeline/04_protein_projection.R:38-43`; `repos/immgenT-GP-analysis/script/Figure5.R:69-125` |
| D-pl-28 | GP × group mean loading heatmap | `repos/immgenT-GP-analysis/script/Figure1.R:296-310` |
| D-pl-29 | Structure plot of GP usage per cell by group | `repos/immgenT-GP-analysis/script/Figure3.R:132-170` |
| D-pl-30 | AUC swarm per GP by lineage | `repos/immgenT-GP-analysis/script/Figure3.R:82-130` |
| D-pl-34 | Row-centred mean-activity heatmap | `repos/immgenT-GP-analysis/code/R/centered_mean_heatmap.R:77-145` |

Checked against the installed reference skills `cnmf`
(`~/.claude/skills/cnmf/SKILL.md`) and `hdwgcna`
(`~/.claude/skills/hdwgcna/SKILL.md`). The engine follows cnmf's steps:
`n_iter` replicates (here `--n-replicates`; cnmf advises 100-200 for real
data), Frobenius loss (cnmf's default `beta_loss`), the local-density filter
(`density_threshold`) and the z-scored `gene_spectra_score`. hdWGCNA builds
co-expression modules on metacells, a different model; its kME and hub genes
are not computed here.

## ClawBio mechanisms used

IDs from `docs/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data | `examples/demo.h5ad`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract"; unique ids, `--reference` in `result.json` |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict, `skillkit.run`; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-11, F-13, F-14 | Fixed plot spec, thresholds drawn, style, Agg | vendored `plots.py`, `theme.py`; AUC 0.5 line |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` via `skillkit` |
| F-23, F-26 | Thresholds documented, no defaults | parameter table |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", frontmatter |
| F-27 | One skill, one task | programs at one k; no k sweep, no donor-level tests |
| F-28 | Golden-value demo test | `tests/test_gene_programs.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

The engine re-implements cNMF's consensus with sklearn instead of depending on
the `cnmf` package: cNMF writes a file tree per run and drives its own
multiprocessing, which fights the fixed output tree and single-thread pinning.
BLAS and OpenMP run on one thread (`threadpoolctl`), so reruns give identical
tables. Ceilings: the pairwise distance matrix is `(k·n-replicates)²` and NMF
cost grows with cells × HVGs × replicates; subsample cells upstream if a run is
too slow.

The demo data come from `examples/make_demo.py` (seed 0): groups `g1`-`g3`, 60
cells each, 300 genes. `G000`-`G009` are planted program P1 (mostly in `g1`),
`G010`-`G019` P2 (`g2`), `G020`-`G029` P3 (`g3`), `var["program_of"]`. With the
demo values all 60 replicate components pass the filter, 20 per program, and the
top 10 genes of `GP1`, `GP2`, `GP3` are exactly P2, P1, P3; each program's best
group has AUC 1.

From `skills/gene-programs`:

```bash
UV_CACHE_DIR=$TMPDIR/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with anndata,scanpy,scikit-learn,matplotlib,plotnine-extra,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

## Citations

- Kotliar et al. (2019). Identifying gene expression programs of cell-type
  identity and cellular activity with single-cell RNA-Seq. eLife 8:e43803.
- Wolf, Angerer, Theis (2018). SCANPY. Genome Biology 19:15.

## Biology-core checks (2026-10-28)

- `--context` is required and opens the report and `chat_summary_lines`.
- `gene_loadings.tsv` keeps every input gene: a gene with no counts gets one
  `not_detected` row with NA program and score; `status_counts` (per gene) in
  `result.json`. Figures and `top_genes.tsv` use scored genes only.
- Gene ids: required `--gene-id-col` (or `none`); `check_gene_ids` runs before
  any analysis and replaces the `var_names` uniqueness check, so a duplicated
  symbol with two ids stops as ambiguous. `gene_id` sits beside `gene` in
  `gene_loadings.tsv` and `top_genes.tsv`; the `gene` column keeps its name.
  The demo's `var["gene_ids"]` comes from `examples/make_demo.py`; X and the
  goldens are unchanged.
- A gotcha keeps claims at "associated with" (house-style `claims.strength`).
- `check_compatible` does not apply: one input.
