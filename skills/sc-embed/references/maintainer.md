# sc-embed: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`.
Row IDs from `mycelium-extra/plans/clawbio-determinism/parts/*.md`.

| Row | Step | Source |
|---|---|---|
| A-an-04 | `normalize_total(target_sum=1e4)`, `log1p`, HVG with `batch_key` | `repos_new/spatial_skin_atlas/processing_scripts/2_acr_map_for_milo.py:19-23,214-225` |
| A-an-05 | scVI integration (out of scope) | `repos_new/spatial_skin_atlas/processing_scripts/2_acr_map_for_milo.py:358-405` |
| A-an-10 | neighbours, `umap(min_dist=)`, Leiden resolution | `repos_new/spatial_skin_atlas/processing_scripts/3_subcluster_coarse.py:32-39` |
| B-an-06 | normalise, HVG, scale, arpack PCA, Harmony (`max_iter_harmony=20`), kNN, UMAP | `repos/aifi-healthy-pbmc-reference/03-Subclustering/06-Python_qc_and_clustering.ipynb:482-504,528,587-590,612-615` |
| B-an-07 | Leiden at a set resolution | `repos/aifi-healthy-pbmc-reference/03-Subclustering/07-Python_B_cells.ipynb:82-95` |
| C-an-05 | scale, PCA, neighbours, `umap(random_state=0)`, `leiden(flavor="igraph", n_iterations=2)` | `repos/multiple-myeloma/analyses-rna/00-utilities/functions/python/process_scrna_data.py:60-97` |
| C-an-06 | Harmony on PCs replaces `X_pca` | `repos/multiple-myeloma/analyses-rna/00-utilities/functions/python/process_scrna_data.py:19-20,83-85` |
| C-an-10 | Leiden resolution chosen per analysis, `n_iterations=2` | `repos/ALTRA-manuscript/Analysis/scRNA/NMF_Analysis/RA_03_NMF_Analysis_Non_Naive_CD4_Tcells.ipynb:216-217` |
| D-an-04 | CCA integration (Seurat, out of scope) | `repos/reflex_manuscript_v2/Fig1c_d_e_f_g_Benchmarking/05-R_Integration.ipynb:262-285` |
| D-an-05 | Harmony by donor after PCA | `repos/reflex_manuscript_v2/Fig2g_h_VST3/05-R_Peptide_Expansion.ipynb:442` |
| D-an-37 | HVG `n_top_genes` with `batch_key`; scVI grid (out of scope) | `repos_new/cns-haematopoiesis/analysis/scRNA_seq/02_Integration/05_batch_correction_scvi_job_array.py:23-66` |
| D-an-40 | neighbours on the integrated space, UMAP, Leiden resolution | `repos_new/cns-haematopoiesis/analysis/scRNA_seq/03_annotation/leiden_clustering_bone_atlas.py:11-13` |

Reference skills (docs/pipeline-skills.md table; `S` is
`~/.claude/skills/synced/f6b88190-5f7e-45de-8798-8d1f5b4ccea4_25160026-14a0-4c5c-af37-b6bb232e66d4`):
`S/single-cell-best-practices` (`SKILL.md:132-142`, `normalize_total`, `log1p`,
`neighbors`, `umap`; `SKILL.md:224`, batch-aware HVG). It integrates with
`sce.pp.harmony_integrate` (`SKILL.md:175`); this skill calls
`harmonypy.run_harmony` directly because the scanpy wrapper transposes
harmonypy 2.x output. `single-cell-recipes` (`references/advanced_snippets.md`,
Harmony variants); `S/scvi-tools` for the later `--integration scvi` engine.

## ClawBio mechanisms used

IDs from `mycelium-extra/plans/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data | `examples/demo.h5ad`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict and `skillkit`; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-13, F-14 | Fixed plot spec, style, Agg | vendored `plots.py`, `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` |
| F-23, F-26 | Parameters documented, closed option sets | parameter table; `--flavor`, `--integration` choices |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", "Trigger" |
| F-27 | One skill, one task; hand off through outputs | `data/embedded.h5ad` feeds annotation and markers |
| F-28 | Golden-value demo test | `tests/test_sc_embed.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

From `skills/sc-embed`:

```bash
UV_CACHE_DIR=$TMPDIR/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with scanpy,harmonypy,igraph,scikit-misc,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

Rebuild the demo with `python examples/make_demo.py` under the same
`uv run --with` line; the golden values in the test then need a re-check.

## Trigger

**Fire this skill when the user says any of:**
- "normalise, find HVGs, PCA, UMAP and Leiden on my QC'd h5ad"
- "integrate the batches with Harmony and cluster"
- "show the UMAP by batch and by cluster"
- "cluster the cells at resolution 0.8 with seed 0"

**Do NOT fire when:**
- the cells still need QC or doublet removal (use sc-qc first);
- the user wants scVI or scANVI integration (not built yet; see below);
- the question is cell-type labels, marker genes or DE between groups.

## Why this exists

- **Without it**: the agent writes a new scanpy block each time. HVG settings,
  the PC count, the seeds and the Leiden flavor drift, so clusters change
  between runs of the same data.
- **With it**: one tested script with a fixed step order and fixed seeds. Every
  choice is a named argument in `commands.sh`, and a rerun gives the same
  clusters.

## Demo

`examples/make_demo.py` (seed 0) simulates 2 batches x 3 groups x 100 cells
and 300 genes. Genes `G000`-`G059` mark the groups (`obs["true_group"]`).
Batch `b2` scales every gene by its own lognormal(0, 0.8) factor.

- `--demo` (Harmony): 3 clusters of 200 cells. Each cluster is one true group
  and half of each batch.
- `--demo --integration none`: 6 clusters of 94-106 cells, each over 93% one
  batch. The batch shift splits every group in two.

## Fixed plot spec

Figures come from the vendored `scripts/plots.py` and plain matplotlib in the
same theme, saved by `save_panel` as PDF, PNG and source CSV. No titles;
captions sit in `report.md`. Batch and cluster colours come from one fixed
10-colour list, assigned in sorted label order.

- `pca_variance`: variance explained (%) per PC, dashed line at `--n-pcs`.
- `umap_batch`, `umap_cluster`: `plots.plot_embedding(axis_name="UMAP")`,
  rasterised points. Each caption names the PCs and representation, neighbours,
  `min_dist`, seed and integration (house-style `fig.text`).
- `batch_composition`: `plots.plot_composition`, % of each cluster's cells per
  batch.

`scripts/skillkit.py`, `plots.py`, `theme.py` and `palettes.py` are byte copies
of `templates/`; a test fails when they drift.

## Citations

- Wolf, Angerer, Theis (2018). SCANPY. Genome Biology 19:15.
- Korsunsky et al. (2019). Harmony. Nature Methods 16:1289-1296.
- Traag, Waltman, van Eck (2019). From Louvain to Leiden. Scientific Reports 9:5233.
- McInnes, Healy, Melville (2018). UMAP. arXiv:1802.03426.
