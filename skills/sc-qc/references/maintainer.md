# sc-qc: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`.
Row IDs from `docs/clawbio-determinism/CATALOGUE.md`.

| Row | What it shows | Source |
|---|---|---|
| A-an-01 | `calculate_qc_metrics(qc_vars=["mt","ribo"])`, three threshold tiers | `repos_new/spatial_skin_atlas/processing_scripts/0_load_preprocess_scrna.py:504,546,604-645` |
| A-an-02 | CellBender before QC (out of scope here) | `repos_new/skin_fibroblast_atlas/integration/0_run_cellbender_gse.sh:1-20` |
| A-an-03 | Scrublet rate from the recovered-cell count | `repos_new/skin_fibroblast_atlas/integration/run_scrublet.ipynb:48-74,118-122` |
| B-an-01 | Scrublet per sample, `random_state = 3030` | `repos/aifi-healthy-pbmc-reference/02-Reference_Labeling/04-Python_Doublet_detection.ipynb:199-205` |
| B-an-02 | mito < 10, genes 200-5000, drop `predicted_doublet` | `repos/aifi-healthy-pbmc-reference/03-Subclustering/06-Python_qc_and_clustering.ipynb:209,272,330,380` |
| C-an-01 | `min_genes = 200`, `MT-` flag | `repos/ALTRA-manuscript/Analysis/scRNA/Preprocessing/ALTRA/04-Python_subset_qc_filtering.ipynb:60,146` |
| C-an-02 | Scrublet per 10x file, keep `== FALSE` | `repos/ALTRA-manuscript/Analysis/scRNA/Preprocessing/ALTRA/02-Python_doublet_detection.ipynb:220-226,401` |
| C-an-03 | DoubletDetection instead of Scrublet (not used) | `repos/single_cell_aging/01_snRNA_seq_QC/03_DoubletDetection/03_DoubletDetection_Downsampling_QC.ipynb:271-279` |
| D-an-02 | Scrublet per sample, drop called doublets | `repos/reflex_manuscript_v2/Fig1c_d_e_f_g_Benchmarking/01-Python_CellTypist_Scrublet.ipynb:584-595` |
| D-an-03 | Hard thresholds, drawn as lines on the QC plot | `repos/reflex_manuscript_v2/Fig2b_c_d_VST2/02-R_Clean_Object_Manual.ipynb:341,539-540` |
| D-an-34 | Per-sample MAD outliers (`upper_thresh`, 3 MADs) | `repos_new/cns-haematopoiesis/analysis/scRNA_seq/01_processing/05_basic_qc_filtering_add_metadata.py:51-52,90-117` |
| D-an-35 | Scrublet per sample in a helper | `repos_new/cns-haematopoiesis/analysis/scRNA_seq/01_processing/helper_functions.py:241-262` |

The MAD rule mirrors cns `is_outlier` / `is_upper_outlier`
(`repos_new/cns-haematopoiesis/analysis/scRNA_seq/01_processing/helper_functions.py:16-31,103-105`):
two-sided on `log1p` counts and genes (5 MADs), upper on mito (3 MADs). One
change: cns thresholds mito on `log1p_total_counts_mt`; this skill uses
`pct_counts_mito`, so the limit reads as a percentage. D-an-35's cluster-level
BH doublet test is not ported.

Reference skills (docs/pipeline-skills.md table):
`S/single-cell-rna-qc` (`scripts/qc_core.py:64` `detect_outliers_mad`, median ±
n MADs on counts, genes and % mito; it pools all cells and works on the raw
scale, where this skill runs per sample on `log1p` with upper-only mito) and
`S/single-cell-best-practices` (`references/preprocessing.md:149`,
`sc.pp.scrublet(adata, batch_key="sample")`, the per-sample call used here). `S` is
`~/.claude/skills/synced/f6b88190-5f7e-45de-8798-8d1f5b4ccea4_25160026-14a0-4c5c-af37-b6bb232e66d4`.

## ClawBio mechanisms used

IDs from `docs/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data | `examples/demo.h5ad`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-11, F-13, F-14 | Fixed plot spec, thresholds drawn, style, Agg | "Fixed plot spec"; vendored `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` |
| F-23, F-26 | Thresholds documented, closed option sets | parameter table; `--mode` choices; no defaults |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", "Trigger" |
| F-27 | One skill, one task | stops at the filtered h5ad; normalisation is the next skill |
| F-28 | Golden-value demo test | `tests/test_sc_qc.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

From `skills/sc-qc`:

```bash
UV_CACHE_DIR=$TMPDIR/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with scanpy==1.12.4,scikit-image,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

A demo run takes about a minute on this HPC, mostly numba compiling scanpy's
QC kernels. Scrublet output can shift across scanpy versions; the goldens are
for scanpy 1.12.4.

## Trigger

**Fire this skill when the user says any of:**
- "QC my scRNA-seq h5ad and remove doublets"
- "filter cells with more than 10% mito and fewer than 200 genes"
- "flag outlier cells at 3 MADs per sample"
- "run Scrublet per sample"

**Do NOT fire when:**
- the task is ambient RNA removal (CellBender, SoupX) or empty-droplet calling;
- the data are already filtered and the user wants normalisation, HVGs,
  embedding, clustering or annotation;
- the samples are bulk RNA-seq libraries.

## Why this exists

- **Without it**: the agent writes fresh QC code each time. Thresholds drift,
  Scrublet runs on the pooled object or without a seed, and one global cutoff
  hides a shallow sample.
- **With it**: one tested script with fixed steps, fixed file names and fixed
  figures. Every threshold is a named argument in `commands.sh`.

## Demo

`examples/make_demo.py` (seed 0) builds 3 samples (`S1`-`S3`, median 2500, 1500
and 4000 UMIs) of 345 cells and 533 genes (13 `MT-`, 20 `RPS`/`RPL`). Each sample
has 300 good cells from 5 populations, 20 low-quality cells (12% of the library,
40% mito) and 25 heterotypic doublets. `obs["planted"]` holds the truth.

With `--demo` (mad 5/5/3, rate 0.08, auto threshold) 877 of 1035 cells are kept
(S1 289, S2 293, S3 295). 86 cells fail QC and Scrublet calls 83 doublets. Every
planted low-quality cell and doublet is removed. 23 good cells also go: the
synthetic mito fraction is so tight that 3 MADs sits near 3.7%. With
`--mode fixed` (genes 200-2500, counts >= 500, mito <= 10%) all 900 good cells
are kept.

## Fixed plot spec

Plain matplotlib in the vendored `theme.py`, saved by `save_panel` as PDF, PNG
and source CSV; the mito colour scale is `plots.grey_to_red`. No titles; a
sample name sits in the corner of each facet; captions are in `report.md`.

- `qc_violins`: genes and counts (log y) and % mito, one violin per sample,
  before filtering; each sample's limits as dashed segments.
- `counts_vs_genes`: counts against genes, log axes, one facet per sample with
  shared axes, coloured by % mito, that sample's limits dashed.
- `doublet_scores`: 40 bins from 0 to the highest score, one facet per sample,
  that sample's Scrublet threshold dashed.
- `cells_per_sample`: stacked counts per sample by `fate`; bar height is the
  count before QC.

`scripts/skillkit.py`, `plots.py`, `theme.py` and `palettes.py` are byte copies
of `templates/`; a test fails when they drift.

## Citations

- Wolock, Lopez, Klein (2019). Scrublet. Cell Systems 8:281.
- Wolf, Angerer, Theis (2018). SCANPY. Genome Biology 19:15.
- Heumos et al. (2023). Best practices for single-cell analysis. Nat Rev Genet 24:550.
