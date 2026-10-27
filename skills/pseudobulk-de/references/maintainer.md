# pseudobulk-de: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`:

- Aggregation per sample and the 10% detection filter (`min_detection <- 0.1`):
  `repos/aifi-healthy-pbmc-reference/07-Pseudobulk/35-R_pseudobulk_seurat.ipynb:57`.
- Minimum cells and counts per pseudobulk sample (`min_cells=20, min_counts=1000`):
  `repos_new/spatial_skin_atlas/figs/3_Fig3/mintflow_tdc_degs.ipynb:516-523`.
- Additive DESeq2 design with covariates
  (`~ cohort.cohortGuid + subject.biologicalSex + subject.cmv`):
  `repos/aifi-healthy-pbmc-reference/08-Differential_expression/36-R_pseudobulk_deseq2.ipynb:200-213`.

## ClawBio mechanisms used

IDs from `mycelium-extra/plans/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data | `examples/demo.h5ad`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-11, F-13, F-14 | Fixed plot spec, thresholds drawn, style, Agg | vendored `plots.py`, `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` |
| F-23, F-26 | Thresholds documented, closed option sets | parameter table; `--engine` choices; no defaults |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", "Trigger" |
| F-28 | Golden-value demo test | `tests/test_pseudobulk_de.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

Tests, from `skills/pseudobulk-de` (the R tests skip without Rscript):

```bash
UV_CACHE_DIR=/tmp/uvcache PYTHONDONTWRITEBYTECODE=1 \
  RSCRIPT=/share/software/statistical/R/R-4.5.2/gcc-8.5.0/R-4.5.2/bin/Rscript \
  R_LIBS_USER=/exports/para-lipg-hpc/mdmanurung/R/4.5 \
  LD_LIBRARY_PATH=/home/mdmanurung/.local/lib/conda-gcc-runtime \
  uv run --python 3.12 --with anndata,numpy,pandas,scipy,matplotlib,pydeseq2,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

The script finds Rscript from `RSCRIPT`, else `Rscript` on `PATH`. On this HPC
the minimal working environment is:

```bash
export RSCRIPT=/share/software/statistical/R/R-4.5.2/gcc-8.5.0/R-4.5.2/bin/Rscript
export LD_LIBRARY_PATH=/home/mdmanurung/.local/lib/conda-gcc-runtime  # libstdc++ with GLIBCXX_3.4.29
export R_LIBS_USER=/exports/para-lipg-hpc/mdmanurung/R/4.5           # holds DESeq2 1.50.2, dreamlet 1.8.0
```

`R_HOME` is not needed. `module load statistical/R/4.5.2/gcc.8.5.0` sets the
first two; `~/.Renviron` sets `R_LIBS_USER` only when `HOME` is the user's.


## Trigger

**Fire this skill when the user says any of:**
- "pseudobulk DE for <cell type> between <conditions>"
- "run DESeq2 on pseudobulk counts per donor"
- "which genes change in monocytes after stimulation, paired by patient"
- "sum the counts per sample and test stim vs ctrl"

**Do NOT fire when:**
- the test is between cells (Wilcoxon, t-test, `rank_genes_groups`, markers);
- the data are bulk RNA-seq without cells (use a bulk DE skill);
- the question is about cell-type proportions or differential abundance.

## Why this exists

- **Without it**: the agent writes a fresh pseudobulk script each time. Filters,
  the design and the cutoffs drift from run to run, and cell-level tests creep in.
- **With it**: one tested script with fixed steps, fixed file names and fixed
  figures. Every analysis choice is a named argument in `commands.sh`.

## Demo

The demo data come from `examples/make_demo.py` (seed 0): 7 patients, `ctrl` and
`stim` samples for P1-P6, one `ctrl` sample for P7 with 5 T cells, 300 genes.
Genes `G000`-`G029` are truly DE in T cells (`var["true_de"]`, log2FC ±1.5); the
last 40 genes are near-silent. With `--engine python` the demo keeps 12 samples
from 6 donors (P7_ctrl dropped: `dropped_by_level` is
`{"stim": [], "ctrl": ["P7_ctrl"]}`), tests 250 genes and calls 28 significant, all truly DE.
`--engine r` gives the same 28. `--engine dreamlet` (design
`~ condition + (1|patient)`) keeps the same 12 samples, tests 250 genes and calls
29 significant, all 29 among the 30 truly DE genes. Its R step takes several
minutes here, most of it loading the packages.

## Fixed plot spec

Figures come from the vendored `scripts/plots.py` (`plot_volcano`,
`plot_heatmap`) and plain matplotlib in the same theme, saved by `save_panel` as
PDF, PNG and source CSV. No titles; captions sit in `report.md`.

- `sample_pca`: log2 CPM, top 500 variable genes, numerator red, denominator blue.
- `pvalue_hist`: raw p-values, 20 bins of width 0.05.
- `ma`: log2FC against mean normalised count (log x), dashed at ±`--lfc-min`.
  With `--engine dreamlet` the x axis is `AveExpr`, mean log2 CPM (linear x).
- `volcano`: genes with a `padj`, dashed at both cutoffs.
- `top_genes_heatmap`: 20 genes with the lowest `padj`, row z-scored log2 CPM,
  not clustered (rows by adjusted-p rank, columns by contrast label),
  colours saturate at ±2 (`HEATMAP_LIMIT`, a display constant, not a threshold).

`scripts/plots.py`, `theme.py`, `palettes.py` and `skillkit.py` are byte copies of
`templates/`; a test fails when they drift.

## Reference skills

`dreamlet` (`~/.claude/skills/dreamlet/SKILL.md`) documents the dreamlet calls
and the `~ condition + (1|donor_id)` form used by `--engine dreamlet`; this
skill's calls agree with it.

## Citations

- Love, Huber, Anders (2014). DESeq2. Genome Biology 15:550.
- Muzellec et al. (2023). PyDESeq2. Bioinformatics 39:btad547.
- Hoffman GE et al. (2023). Efficient differential expression analysis of
  large-scale single cell transcriptomics data using dreamlet. bioRxiv.
  doi:10.1101/2023.03.17.533005. Vignette:
  <https://diseaseneurogenomics.github.io/dreamlet/articles/dreamlet.html>.
