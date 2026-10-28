# sc-markers: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Corpus sources

Reference skills checked against (`docs/pipeline-skills.md`): `scanpy-plotting`
and `S/single-cell-best-practices` (`references/differential_analysis.md`).

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`;
row IDs from `docs/clawbio-determinism/parts/*.md`.

| ID | What | Source |
|---|---|---|
| A-an-12 | `rank_genes_groups(..., method='t-test')` per annotation | `repos_new/skin_fibroblast_atlas/figures/extended_fig5_scrna_disease.ipynb:758` |
| B-an-09 | `rank_genes_groups(adata_subset, res_name, method='wilcoxon')` | `repos/aifi-healthy-pbmc-reference/03-Subclustering/07-Python_B_cells.ipynb:1421` |
| C-an-19 | `rank_genes_groups` then `rank_genes_groups_dotplot` per Leiden resolution | `repos/ALTRA-manuscript/Figures/Figure3/Figure_3AB_S6BC.ipynb:696,739` |
| D-an-14 | `FindAllMarkers(..., group.by='CTaa')` (Seurat, Wilcoxon default) | `repos/reflex_manuscript_v2/Supp_Fig10_VST_5prime/02-R_Supp.ipynb:1916` |
| D-an-40 | Wilcoxon, `corr_method='benjamini-hochberg'`, `reference='rest'`, `use_raw=False` | `repos_new/cns-haematopoiesis/analysis/scRNA_seq/02_Integration/10_batch_correction_scvi_haem_best_integration.py:82` |
| A-pl-02 | dotplot, `dendrogram=False`, `categories_order`, colour "Mean expression (log1p norm)" | `repos_new/spatial_skin_atlas/figs/sFig1_to_11_dotplots_by_lineage/s1a_b.ipynb:430-441` |
| B-pl-03 | dotplot of markers by cluster | `repos/aifi-healthy-pbmc-reference/03-Subclustering/07-Python_B_cells.ipynb:261-265` |
| D-pl-19 | dotplot, `dendrogram=False`, `categories_order`, Reds | `repos_new/cns-haematopoiesis/figures/Chapter_3/Bone_Atlas_Figs.ipynb:880-915` |

## ClawBio mechanisms used

IDs from `docs/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data | `examples/demo.h5ad`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict, `skillkit.run`; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-11, F-13, F-14 | Fixed plot spec, thresholds drawn, style, Agg | vendored `plots.py`, `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` via `skillkit` |
| F-23, F-26 | Thresholds documented, closed option sets | parameter table; `--method` choices; no defaults |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", frontmatter |
| F-27 | One skill, one task | markers only; condition DE hands off to pseudobulk-de |
| F-28 | Golden-value demo test | `tests/test_sc_markers.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

The demo data come from `examples/make_demo.py` (seed 0): clusters `c1`-`c3`,
50 cells each, 200 genes. `G000`-`G004` are planted in `c1`, `G005`-`G009` in
`c2`, `G010`-`G014` in `c3` (`var["marker_of"]`, 8-fold up). `var["gene_ids"]`
holds Ensembl-like ids `ENSG90000000000.1` onward. With the demo
values both methods recover exactly the 15 planted markers, 5 per cluster, as
the only passing genes; `top_markers` per cluster is its planted set.

From `skills/sc-markers`:

```bash
UV_CACHE_DIR=$TMPDIR/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with anndata,scanpy,matplotlib,plotnine-extra,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

## Fixed plot spec

Both figures come from the vendored `scripts/plots.py` and are saved by
`save_panel` as PDF, PNG and source CSV. No titles; captions sit in `report.md`.
Genes run along x in top-marker order, grouped by group; groups run top-down in
`--groupby` order (`dendrogram=False`, `categories_order` in the corpus).

- `marker_dotplot` (`plot_dotplot`): colour is the group mean of the
  log-normalised values, unscaled (A-pl-02, B-pl-03; D-pl-19's `standard_scale`
  is not used); dot area is the fraction of the group's cells with a value > 0.
  The diameter at 100% is 4 mm (`DOT_MAX_SIZE`). It stays a constant, not a
  grill question, because the panel widens 4 mm per gene: a 4 mm dot at 100%
  fills its column without touching the next one. It changes no result.
- `marker_heatmap` (`plot_heatmap`): group means z-scored per gene across groups,
  diverging palette centred at 0, saturating at ±`--limit`. The order is fixed
  by label, not clustered: rows in `--groupby` order, columns in top-marker
  order. The caption says so.

Width grows with the gene count from 89 to 183 mm; height grows with the group
count from 85 mm. `scripts/skillkit.py`, `plots.py`, `theme.py` and
`palettes.py` are byte copies of `templates/`; a test fails when they drift.

## Biology-core checks (2026-10-28)

From `docs/pipeline-skills.md` "Biology-core checks". Applied:

- context: required `--context`, first line of the report and `chat_summary_lines`.
- status: `markers.tsv` has `status`; every gene is tested in every group (no
  filter), so every row is `tested` and `status_counts` is `{"tested": 600}`.
- gene ids: required `--gene-id-col` (demo `gene_ids`, added to `examples/make_demo.py`
  after the random draws, so the counts and goldens are unchanged). `check_gene_ids`
  replaces the `check_unique(var_names)` check, so a duplicated `var_name` with two
  ids is the ambiguous-map error. `gene` became `gene_symbol` in both tables.
- claim wording: a gotcha says a marker is associated with its group, not a
  driver or an identity.

Not applied: compatible inputs (one input file).

## Citations

- Wolf, Angerer, Theis (2018). SCANPY. Genome Biology 19:15.
