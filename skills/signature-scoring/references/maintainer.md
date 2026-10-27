# signature-scoring: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Corpus sources

Reference skills checked against (`docs/pipeline-skills.md`): `decoupler` and
`S/decoupler-enrichment` (`dc.mt.aucell`, its `score_aucell` obsm key).

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`;
row IDs from `mycelium-extra/plans/clawbio-determinism/parts/*.md`.

| ID | What | Source |
|---|---|---|
| A-an-14 | `sc.tl.score_genes(..., ctrl_size=1000, n_bins=25)`, then per-group means | `repos_new/skin_fibroblast_atlas/figures/fig4_by_scarring_risk.ipynb:904-905,976-979` (`parts/A_skin.md:22`) |
| A-an-15 | `score_genes` on OpenTargets GWAS sets, CONTROL_SIZE=1000, N_BINS=25 | `repos_new/spatial_skin_atlas/figs/revision/supp24_GWAS_all_celltyprs.ipynb:1236-1255` (`parts/A_skin.md:23`) |
| B-an-24 | composite score from DEGs, per-sample score table | `repos/IHA-Figure/Figure2/04_Composite_Score/02A_Aging_CompositeScore_SF4.ipynb:161,182` (`parts/B_pbmc.md:33`) |
| C-an-28 | `dc.mt.aucell(adata, net, tmin=3)`, scores from `obsm['score_aucell']` | `repos/multiple-myeloma/analyses-rna/03-subset-analyses/01-plasma_cells/09-py-literature_states.ipynb:542-547` (`parts/C_disease.md:36`) |
| C-an-30 | NNLS IFN scores per donor (out of scope here; per-donor summary kept) | `repos/IFN_manuscript/Scripts/IFN_Score_Generate.R:71-81` (`parts/C_disease.md:38`) |
| D-an-16 | UCell module scores, `maxRank=2000` (the rank cut-off `--n-up` mirrors) | `repos/reflex_manuscript_v2/Fig2g_h_VST3/06-R_VST3_Fig2h.ipynb:1532` (`parts/D_reflex.md:24`) |
| D-an-18 | curated gene list (pathway -> gene) averaged per group | `repos/reflex_manuscript_v2/Fig2b_c_d_VST2/05-R_Fig2c.ipynb:222-287` (`parts/D_reflex.md:26`) |
| A-pl-10 | score by group (violin there; box with sample points here, house-style `fig.wilke`) | `repos_new/spatial_skin_atlas/figs/revision/supp24_GWAS_all_celltyprs.ipynb:1324-1335` (`parts/A_skin.md:59`) |
| A-pl-11 | score dot-matrix by cell type | `repos_new/spatial_skin_atlas/figs/revision/supp24_GWAS_all_celltyprs.ipynb:1236-1242` (`parts/A_skin.md:60`) |
| B-pl-02 | UMAP feature plot, continuous colour | `repos/IHA-Figure/Extended_Figure9/02_TEAseq_CD8_T_cells_traj/01_Cleaning.ipynb:408` (`parts/B_pbmc.md:71`) |
| PIPELINES.md | "score vs random-gene-set control", the missing standard plot | `mycelium-extra/plans/clawbio-determinism/PIPELINES.md:119` |

## ClawBio mechanisms used

IDs from `mycelium-extra/plans/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data | `examples/demo.h5ad`, `examples/demo_gene_sets.tsv`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict, `skillkit.run`; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-13, F-14 | Fixed plot spec, style, Agg | vendored `plots.py`, `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` via `skillkit` |
| F-23, F-26 | Values documented, closed option sets | parameter table; `--method` choices; no defaults |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", frontmatter |
| F-27 | One skill, one task | scoring only; group tests hand off to group-tests |
| F-28 | Golden-value demo test | `tests/test_signature_scoring.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

The demo data come from `examples/make_demo.py` (seed 0): groups `healthy`
(H1-H4) and `disease` (D1-D4), 50 cells per sample, 300 genes. `sig_up`
(G000-G009) is 4-fold up in disease cells; `sig_null` (G010-G019) is flat and
lists `GMISSING`, absent from the data. The demo finds 10 of 10 and 10 of 11
genes, reports `GMISSING`, and gives `beats_random` sig_up true
(`p_empirical` 1/51) and sig_null false. Every disease sample's sig_up mean is
above every healthy one. AUCell with `--n-up 30` gives the same verdicts.

From `skills/signature-scoring`:

```bash
UV_CACHE_DIR=$TMPDIR/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with anndata,scanpy,decoupler,matplotlib,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

## Fixed plot spec

Figures are saved by `save_panel` as PDF, PNG and source CSV, 60 mm high. No
titles; captions sit in `report.md`.

- `score_by_group_<set>` (`plots.plot_box`): one point per sample over a box
  without outliers, jittered ±0.2 box widths (`--seed`). Under 10 samples in
  any group (`plots.FEW_POINTS`) the box gives way to a bar at the median, and
  the caption opens "Dot plot" instead of "Box plot". The y title names the set
  and the unit (`a.u.` for score_genes, `AUC` for AUCell). Groups in
  `--group-col` order, coloured Paul Tol bright in that order. Width 40 mm plus
  20 mm per group.
- `score_umap_<set>` (matplotlib in the house theme, since `plot_embedding`
  takes only a categorical palette): cells at `X_umap`, colour the per-cell
  score on `plots.grey_to_red`, highest drawn last, rasterised; axes "UMAP 1"
  and "UMAP 2". Only when `X_umap` exists. The caption says the UMAP is the
  input's, made upstream.
- `random_control`: per set, the random sets' spreads as grey jittered points
  and the signature's spread as a red diamond (`#B2182B`).

`scripts/skillkit.py`, `plots.py`, `theme.py` and `palettes.py` are byte copies
of `templates/`; a test fails when they drift.

## Citations

- Wolf, Angerer, Theis (2018). SCANPY. Genome Biology 19:15.
- Aibar et al. (2017). SCENIC (AUCell). Nature Methods 14:1083-1086.
- Badia-i-Mompel et al. (2022). decoupleR. Bioinformatics Advances 2:vbac016.
- Tirosh et al. (2016). Melanoma single-cell RNA-seq (the score_genes method). Science 352:189-196.
