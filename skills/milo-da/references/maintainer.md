# milo-da: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`.
Row IDs from `mycelium-extra/plans/clawbio-determinism/parts/*.md`.

| Row | Step | Source |
|---|---|---|
| A-an-16 (`parts/A_skin.md:24`) | milopy `make_nhoods(prop=)`, `count_nhoods(sample_col=)`, `annotate_nhoods`, `DA_nhoods(design=)`, `build_nhood_graph` | `repos_new/spatial_skin_atlas/figs/1_Fig1/1d_1e_milo.ipynb:258-286,296-312` |
| B-an-27 (`parts/B_pbmc.md:36`) | `DA_nhoods(design=, model_contrasts='levelA-levelB')`, SpatialFDR 0.05 | `repos/IHA-Figure/Extended_Figure2/03_FigC.ipynb:123,178`; `repos/IHA-Figure/Figure4/04_FigE.ipynb:576` |
| A-pl-09 (`parts/A_skin.md:58`) | nhood graph on the UMAP, SpatialFDR 0.1 | `repos_new/spatial_skin_atlas/figs/1_Fig1/1d_1e_milo.ipynb:309-312` |
| B-pl-27 (`parts/B_pbmc.md:96`) | nhood graph, `RdBu_r`, `alpha=0.1` | `repos/IHA-Figure/helper_function/helper_function_IHA.py:36`; `repos/IHA-Figure/Figure4/04_FigE.ipynb:71` |

Reference skills (the `docs/pipeline-skills.md` table for milo-da: `milo-r`,
`pertpy-tools`, `S/pertpy-perturbation-analysis`); the calls and gotchas here
agree with them:

- `~/.claude/skills/milo-r/SKILL.md`: miloR method; `references/core-workflow.md:77,83`
  for `prop` 0.05-0.2 and nhood sizes peaking at 50-100 cells;
  `references/gastrulation-example.md:25-28` for the `Mixed` label convention.
- `~/.claude/skills/pertpy-tools/SKILL.md` and `references/milo.md`: pertpy
  `make_nhoods`, `count_nhoods`, `da_nhoods(solver=)`, `annotate_nhoods` and the
  `Mixed` relabel.
- `~/.claude/skills/synced/f6b88190-5f7e-45de-8798-8d1f5b4ccea4_25160026-14a0-4c5c-af37-b6bb232e66d4/pertpy-perturbation-analysis/SKILL.md`: pertpy
  install (`pertpy[de]` for pydeseq2) and when to choose Milo over scCODA.

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
| F-23, F-26 | Parameters documented, closed option sets | parameter table; solver fixed to pydeseq2 |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", "Trigger" |
| F-27 | One skill, one task; hand off through outputs | reads sc-embed `data/embedded.h5ad`; writes `tables/nhoods.tsv` |
| F-28 | Golden-value demo test | `tests/test_milo_da.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

From `skills/milo-da`:

```bash
UV_CACHE_DIR=$TMPDIR/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with pertpy,pydeseq2,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

Rebuild the demo with `python examples/make_demo.py` under the same
`uv run --with` line; the golden values in the test then need a re-check.

## Trigger

**Fire this skill when the user says any of:**
- "run Milo on my embedded h5ad, condition B versus A"
- "which neighbourhoods are differentially abundant between disease and control"
- "plot the Milo nhood graph and the DA beeswarm by cell type"
- "port this milopy notebook to a script"

**Do NOT fire when:**
- the cells have no embedding or kNN graph yet (use sc-embed first);
- the user wants proportions of labelled cell types tested (scCODA, crumblr,
  CLR plus a linear model);
- the question is genes (DE, markers) rather than cell abundance.

## Why this exists

- **Without it**: the agent writes a milopy or pertpy block each time. The
  index-cell sample, k, the solver and the SpatialFDR cutoff drift, and an
  unseeded `make_nhoods` gives different neighbourhoods on every rerun.
- **With it**: one tested script with a fixed step order. Every choice is a
  named argument in `commands.sh`, and a rerun gives the same table.

## Demo

`examples/make_demo.py` (seed 0) simulates 6 samples of 200 cells: A1-A3 in
condition A, B1-B3 in B. Four cell types (T, B, NK, Mono) are Gaussian blobs in
a 10-D `X_pca`; per-sample shares jitter by Dirichlet(200 p). Mono is planted at
3x its condition-A share in B. Each sample has its own donor (D1-D6). The
file stores a k = 30 graph and a UMAP, and no
genes.

`--demo` gives 75 neighbourhoods, 18 at SpatialFDR < 0.1: 17 up, all Mono, and
1 down in B (the compositional dip the Mono expansion causes). The dominant
annotation is Mono. `--demo --n-neighbors 20` rebuilds the graph.

## Fixed plot spec

Plain matplotlib in the vendored `theme.py` style, saved by `save_panel` as PDF,
PNG and source CSV; the diverging map comes from `plots.py`. No titles; captions
sit in `report.md`. Non-significant neighbourhoods are light grey; significant
ones are coloured by log2 fold change on a symmetric scale.

- `nhood_graph`: one dot per neighbourhood at its index cell's UMAP position,
  area by size, grey edges for shared cells (milopy `plot_nhood_graph`).
- `da_beeswarm`: log2 fold change per annotation, seeded vertical jitter,
  annotations ordered by median logFC, dashed line at 0; y axis "neighbourhood
  label", group sizes in the tick labels.
- `nhood_size`: histogram of cells per neighbourhood (numpy `auto` bins).

`scripts/skillkit.py`, `plots.py`, `theme.py` and `palettes.py` are byte copies
of `templates/`; a test fails when they drift.

## Citations

- Dann, Henderson, Teichmann, Morgan, Marioni (2022). Differential abundance
  testing on single-cell data using k-nearest neighbor graphs. Nature
  Biotechnology 40:245-253.
- Heumos et al. (2025). Pertpy: an end-to-end framework for perturbation
  analysis. Nature Methods.
- Muzellec et al. (2023). PyDESeq2. Bioinformatics 39:btad547.
