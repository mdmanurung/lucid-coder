# ccc: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Corpus sources

Reference skills checked against (`docs/pipeline-skills.md`): `liana-ccc`
(`li.mt.rank_aggregate`, `magnitude_rank`, `specificity_rank`) and
`S/liana-cell-communication`.

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`;
row IDs from `mycelium-extra/plans/clawbio-determinism/parts/*.md`
(`CATALOGUE.md:47`).

| ID | What | Source |
|---|---|---|
| A-an-25 | LARIS spatial LR scoring, CellChatDB as a local CSV resource, `random_seed=27`, `expressed_pct=0.1` (`parts/A_skin.md:33`) | `repos_new/spatial_skin_atlas/figs/revision/e10f_LARIS_sebimmune.ipynb:478,710-725` |
| B-an-31 | CellPhoneDB statistical analysis per sample, database v5.0.0 downloaded once (`parts/B_pbmc.md:40`) | `repos/IHA-Figure/Extended_Figure8/05_CellPhoneDB/01_Calculate_Scores.ipynb:33-34`; `00_DownloadDB.ipynb:71` |
| C-an-31 | liana `cellphonedb(groupby=, use_raw=False, resource_name='consensus', expr_prop=0.1, n_perms=5000)`; `rank_aggregate` inspected (`parts/C_disease.md:39`) | `repos/ALTRA-manuscript/Analysis/scRNA/Tfh17_Bmem_interaction/02a_python_tfh17_Beff_liana_analysis.ipynb:822-828`; `03_python_tfh17_Beff_liana_results_viz.ipynb:71,133` |
| A-pl-22 | LR dotplot of the top 8-10 interactions (`parts/A_skin.md:71`) | `repos_new/spatial_skin_atlas/figs/revision/e10f_LARIS_sebimmune.ipynb:762-849` |
| B-pl-16 | circos / chord diagram, R `circlize` (`parts/B_pbmc.md:85`); not drawn here | `repos/IHA-Figure/Figure4/01_FigB.ipynb:111`; `repos/MOCHA_Manuscript/Fig5/Fig5_networkPlot.R:13-43` |

## ClawBio mechanisms used

IDs from `mycelium-extra/plans/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data | `examples/demo.h5ad`, `examples/demo_resource.tsv`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict, `skillkit.run`; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-11, F-13, F-14 | Fixed plot spec, style, Agg | vendored `plots.py`, `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` into liana `seed=`; `commands.sh`, `result.json`; rerun test compares bytes |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` via `skillkit` |
| F-23, F-26 | Thresholds documented, closed option sets | parameter table; `--methods` names; no defaults |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", frontmatter |
| F-27 | One skill, one task | steady-state CCC only; no spatial, factorisation or condition tests |
| F-28 | Golden-value demo test | `tests/test_ccc.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

The demo data come from `examples/make_demo.py` (seed 0): cell types `sender`,
`receiver`, `bystander`, 40 cells each, samples `s1`/`s2` alternating, 60 genes.
`LIG1` and `LIG2` are 8-fold up in `sender`, `REC1` and `REC2` 8-fold up in
`receiver`; the other 18 resource pairs use background genes. With the demo
values exactly 2 of 180 interactions pass: `sender → receiver LIG1-REC1` and
`LIG2-REC2`, the top two rows. `receiver → receiver LIG1-REC1` ranks third
(magnitude rank 0.0025) because `REC1` is high in `receiver`.

From `skills/ccc`:

```bash
UV_CACHE_DIR=$TMPDIR/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with "liana>=1.10",anndata,scanpy,matplotlib,threadpoolctl,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

## Fixed plot spec

Both figures are saved by `save_panel` as PDF, PNG and source CSV. No titles;
captions sit in `report.md`.

- `pair_counts_heatmap`: sources as rows, targets as columns, in `--groupby`
  order, each cell labelled with its count. Counts are unsigned, so it is drawn
  with matplotlib and the sequential `grey_to_red` map from `plots.py`, not the
  diverging `plot_heatmap`. With `--sample` the counts are summed over samples.
  Not clustered: the order is fixed by label, and the caption says so.
- `lr_dotplot` (`plot_dotplot`): ligand → receptor on y, source → target on x
  (prefixed by sample with `--sample`), in table order. Colour is
  `-log10(magnitude_rank)`; dot area is `1 - specificity_rank`. The diameter at 1
  is 4 mm (`DOT_MAX_SIZE`), as in liana's own dotplot (colour magnitude, size
  specificity). It stays a constant, not a grill question, because the panel
  widens 6 mm per pair: a 4 mm dot fills its column without touching the next
  one. It changes no result.
- No chord or circos plot: matplotlib has no chord form, and the corpus one
  (B-pl-16) is R `circlize`. The count heatmap carries the same source x target
  numbers.

`scripts/skillkit.py`, `plots.py`, `theme.py` and `palettes.py` are byte copies
of `templates/`; a test fails when they drift.

## Citations

- Dimitrov et al. (2022). Comparison of methods and resources for cell-cell
  communication inference from single-cell RNA-Seq data. Nature Communications
  13:3224.
- Dimitrov et al. (2024). LIANA+ provides an all-in-one framework for cell-cell
  communication inference. Nature Cell Biology 26:1613-1622.
