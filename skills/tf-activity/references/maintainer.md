# tf-activity: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Corpus sources

Reference skills checked against (`docs/pipeline-skills.md`): `decoupler`
(`dc.mt.ulm`, `dc.pp.get_obsm`, `dc.tl.rankby_group`) and `S/decoupler-enrichment`.

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`;
row IDs from `docs/clawbio-determinism/parts/*.md`. The corpus
uses decoupler 1.x names; the 2.x calls here are `dc.op.collectri`,
`dc.op.progeny`, `dc.mt.ulm`, `dc.mt.mlm`, `dc.pp.get_obsm` and
`dc.tl.rankby_group`.

| ID | What | Source |
|---|---|---|
| A-an-27 | `dc.get_collectri(organism='human', split_complexes=False)`, `dc.run_ulm(mat=adata, net, source, target, weight)`, `dc.get_acts` | `repos_new/skin_fibroblast_atlas/figures/extended_fig5_scrna_disease.ipynb:1275,1302,1350` |
| A-an-27 | `dc.rank_sources_groups(acts, groupby=..., reference='rest', method='t-test_overestim_var')` | `repos_new/skin_fibroblast_atlas/figures/extended_fig2_GSEA_TFs_healthy.ipynb:1027,1114,1191` (the row cites `supp_fig4_prurigo.ipynb`; the lines sit in this file) |
| A-an-28 | `dc.get_progeny(organism='human', top=500)`, `dc.run_mlm(...)` | `repos_new/skin_fibroblast_atlas/figures/fig3_part2_PROGENY.ipynb:218,265,384` |
| A-pl-25 | `sc.pl.violin(acts, keys=['KLF5'], groupby=...)`, commented out | `repos_new/skin_fibroblast_atlas/figures/extended_fig2_GSEA_TFs_healthy.ipynb:1162` |
| B-an-30 | ArchR `addMotifAnnotations(motifSet="JASPAR2020")` on ATAC (related, not done here) | `repos/IHA-Figure/Figure5/06_FigF.ipynb:2146` |
| B-an-51 | MOCHA motif set enrichment, `p.adjust(p_val, "fdr")` (related, not done here) | `repos/MOCHA_Manuscript/R/MotifSetEnrichmentAnalysis.R:209` |

## ClawBio mechanisms used

IDs from `docs/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data | `examples/demo.h5ad`, `demo_network.tsv`, `demo_de.tsv`, `make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict, `skillkit.run`; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree per mode |
| F-10, F-11, F-13, F-14 | Fixed plot spec, style, Agg | vendored `plots.py`, `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` via `skillkit` |
| F-23, F-26 | Thresholds documented, closed option sets | parameter table; `--method`, `--rank-method` choices; no defaults |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", frontmatter |
| F-27 | One skill, one task | activity scoring only; the network export and DE stay outside |
| F-28 | Golden-value demo test | `tests/test_tf_activity.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

The demo data come from `examples/make_demo.py` (seed 0): groups `g1`-`g3`, 60
cells each, 300 genes. TF1-TF6 have 15 targets each (12 up, 3 down); TF1 is
planted active in `g1` and TF2 in `g2`. TF7 (3 targets) and TF8 (0 in the data)
fall under `--tmin 5`. Golden values: in cells mode the top source of `g1` is
TF1 and of `g2` is TF2, and these are the only group-source pairs with
`padj < 0.05` and `stat > 0`; in contrast mode (`g1` against `g2`) the top-2
sources by absolute score are TF1 (positive) and TF2 (negative). Both methods
give the same top sources.

From `skills/tf-activity`:

```bash
UV_CACHE_DIR=$TMPDIR/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with "decoupler>=2",matplotlib,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

## Fixed plot spec

Figures are saved by `save_panel` as PDF, PNG and source CSV. No titles;
captions sit in `report.md`.

- `activity_heatmap` (`plots.plot_heatmap`): group-mean activity, sources as
  rows (top sources in group order), groups as columns, diverging palette
  centred at 0, saturating at ±`--limit`. Not clustered: rows in top-source
  order, columns in `--groupby` order; the caption says so. Colour, violin and
  bar axes carry the unit: decoupler's ULM and MLM scores are t-values.
- `activity_violin` (matplotlib, plots.py has no violin): per-cell activity of
  the same sources, one panel per source, groups on x with `(n = ...)`, median
  bar, grey fill, a grey line at 0. Panels keep their own y axis; the caption says so.
  Cells are many, so a violin replaces house-style's box with every point.
- `activity_bar` (matplotlib): signed score of the `--n-top` sources, bars from
  0, largest on top, one grey fill.

`scripts/skillkit.py`, `plots.py`, `theme.py` and `palettes.py` are byte copies
of `templates/`; a test fails when they drift.

## Citations

- Badia-i-Mompel et al. (2022). decoupleR: ensemble of computational methods to
  infer biological activities from omics data. Bioinformatics Advances 2:vbac016.
- Müller-Dott et al. (2023). Expanding the coverage of regulons from high-confidence
  prior knowledge for accurate estimation of transcription factor activities
  (CollecTRI). Nucleic Acids Research 51:10934.
- Schubert et al. (2018). Perturbation-response genes reveal signaling footprints
  in cancer gene expression (PROGENy). Nature Communications 9:20.
