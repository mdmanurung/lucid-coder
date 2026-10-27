# sc-annotate: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`
(row IDs from `docs/clawbio-determinism/CATALOGUE.md` and `parts/*.md`):

- B-an-03, C-an-07, D-an-01: CellTypist after `normalize_total(target_sum=1e4)`
  + `log1p`, `majority_voting=True` (AIFI) or off (reflex):
  `repos/aifi-healthy-pbmc-reference/02-Reference_Labeling/02-Python_Label_Predictions_Celltypist.ipynb:259-262`;
  `repos/multiple-myeloma/analyses-rna/01-preprocessing/02-doublets_and_label_transfer/02-predict_labels.py:42,83`;
  `repos/reflex_manuscript_v2/Fig1c_d_e_f_g_Benchmarking/01-Python_CellTypist_Scrublet.ipynb:563-580`.
- B-an-04, D-an-41: training a CellTypist model on a log1p-1e4 reference by a
  label column, then annotating with a majority-voting toggle:
  `repos_new/cns-haematopoiesis/analysis/scRNA_seq/03_annotation/cell_typist_model_training.py:36-47`;
  `repos_new/cns-haematopoiesis/analysis/scRNA_seq/03_annotation/celltypist_bone_atlas.py:13-33`;
  `repos/aifi-healthy-pbmc-reference/06-Modeling/30-Python_celltypist_L1_model.ipynb:358-364`.
- D-an-42: kNN label propagation with the winning-neighbour fraction as the
  score: `repos_new/cns-haematopoiesis/analysis/scRNA_seq/03_annotation/label_propagation_faiss.py:51-57`.
- B-an-08, D-an-43: marker dotplots to check labels:
  `repos/aifi-healthy-pbmc-reference/03-Subclustering/07-Python_B_cells.ipynb:70-76`;
  `repos_new/cns-haematopoiesis/analysis/scRNA_seq/03_annotation/annotation.py:45,126`.
- A-an-07, A-an-34, B-an-05, C-an-08 (scArches, scANVI/scPoli, Seurat
  `MapQuery`) and B-an-10 (L1/L2/L3 cleanup) are out of scope: see "Do NOT fire
  when".

Reference skills (docs/pipeline-skills.md table; `S` is
`~/.claude/skills/synced/f6b88190-5f7e-45de-8798-8d1f5b4ccea4_25160026-14a0-4c5c-af37-b6bb232e66d4`):
`S/celltype-annotation-workflow` (marker-based critique of draft labels: the
`--markers` dot plot is that check); `S/celltypist-immune-markers` (CellTypist
Immune_All_Low/High labels and their markers, for building `--markers` against
a CellTypist model); `S/azimuth-pbmc-markers` (Azimuth PBMC L1-L3 markers, for
`--markers` against a PBMC reference).

## ClawBio mechanisms used

IDs from `docs/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data, offline | `examples/demo_*`, `examples/make_demo.py`; the demo trains its own model |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-13, F-14 | Fixed plot spec, style, Agg | vendored `plots.py`, `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` via `skillkit.py` |
| F-23, F-26 | Settings documented, closed option sets | parameter table; `--engine`, `--majority-voting` choices |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", "Trigger" |
| F-27 | One skill, one task | clustering, UMAP and label cleanup stay upstream or downstream |
| F-28 | Golden-value demo test | `tests/test_sc_annotate.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

CellTypist 1.7.1 passes `multi_class` to scikit-learn's `LogisticRegression`,
which scikit-learn 1.8 removed; pin `scikit-learn<1.8`. From `skills/sc-annotate`:

```bash
UV_CACHE_DIR=$TMPDIR/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with celltypist,scanpy,"scikit-learn<1.8",pytest \
  python -m pytest -q -p no:cacheprovider tests
```

Threads are pinned to 1 (`threadpoolctl`) during the transfer, so reruns give
identical tables and figure PDFs.

## Trigger

**Fire this skill when the user says any of:**
- "annotate my cells with CellTypist using <model.pkl>"
- "train a CellTypist model on this reference and label the query"
- "transfer the labels from the reference with scanpy ingest"
- "which cell types are in my clusters, using the AIFI reference"

**Do NOT fire when:**
- there is no reference or model, only marker genes (manual annotation);
- the user wants scANVI, scArches or Seurat anchor mapping;
- the query still needs clustering or a UMAP (run that first);
- the task is merging or cleaning L1/L2/L3 label columns.

## Why this exists

- **Without it**: the agent writes a fresh annotation script each time. The
  normalisation, majority voting and the confidence definition drift, and a
  model may be downloaded mid-run.
- **With it**: one tested script with fixed steps, fixed file names and fixed
  figures. Every choice is a named argument in `commands.sh`, and nothing is
  downloaded.

## Demo

`examples/make_demo.py` (seed 0) simulates 4 types (B cell, Monocyte, NK cell,
T cell) with 2 markers each and 150 background genes. NK and T cells share most
of their programme. The reference has 75 cells per type. The query has 100 per
type at half the depth, `true_label`, 8 clusters (each type split in two, 5% of
cells moved to another type's cluster) and a 2D layout in `X_umap`.
`demo_markers.tsv` lists 8 markers; all three files total under 200 kB.

Golden values, as cells of 400 matching `true_label`:
- `--demo` (celltypist trained on the reference, voting on): 388; labels B 102,
  Monocyte 98, NK 102, T 98; 158 shared genes;
- `--demo --engine ingest` (10 PCs, k = 15): 379;
- the demo's saved model with `--majority-voting no`: 385.

## Fixed plot spec

Figures come from the vendored `scripts/plots.py` (`plot_embedding`,
`plot_dotplot`, `plot_box`, `plot_heatmap`), saved by `save_panel` as PDF, PNG
and source CSV. No titles; captions sit in `report.md`. Labels take `tab10`
colours in sorted order.

- `umap_label`: the query's `X_umap`, coloured by transferred label.
- `marker_dotplot`: marker genes in `--markers` order by label; colour is mean
  log-normalised expression scaled 0-1 per gene, size the fraction expressing
  (`DOT_MAX_SIZE` 4 mm at 100%).
- `confidence_box`: per-cell confidence by label, points jittered
  (`JITTER_WIDTH` 0.3, seeded by `--seed`); y title "CellTypist probability
  (0-1)" or "neighbour agreement (fraction)". A label with under 10 cells gets
  a median line instead of a box (`plots.FEW_POINTS`).
- `label_cluster_heatmap`: fraction of each cluster per label, columns sum to 1,
  scale 0-1 (`HEATMAP_LIMIT`). Order fixed by label (rows by sorted label,
  columns in cluster order), not clustered; the caption says so.
- The UMAP caption points to the upstream report for PCs, neighbours,
  `min_dist` and seed, because this skill only reuses `obsm["X_umap"]`.

`scripts/skillkit.py`, `plots.py`, `theme.py` and `palettes.py` are byte copies
of `templates/`; a test fails when they drift.

## Citations

- Domínguez Conde et al. (2022). CellTypist. Science 376:eabl5197.
- Wolf, Angerer, Theis (2018). SCANPY. Genome Biology 19:15.
