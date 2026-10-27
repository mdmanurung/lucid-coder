# trajectory: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Trigger

**Fire this skill when the user says any of:**
- "order these cells in pseudotime from the progenitors"
- "run PAGA on the clusters and show the graph"
- "diffusion pseudotime rooted in <group>"
- "how do <genes> change along the trajectory"

**Do NOT fire when:**
- the user wants RNA velocity or fate probabilities (scVelo, CellRank: a
  separate skill);
- the user names Slingshot or Monocle3 (not engines here yet; see "Out of scope");
- the object has no neighbour graph or UMAP yet (run sc-embed first);
- the question is about cell-type proportions or differential abundance.

## Why this exists

- **Without it**: the agent writes a fresh PAGA/DPT notebook each time. The
  root cell is picked by eye, the PAGA threshold drifts per plot, and the gene
  trends get a smoother whose span nobody recorded.
- **With it**: one tested script with a fixed root rule, fixed file names and
  fixed figures. Every analysis choice is a named argument in `commands.sh`.

## Demo

`examples/make_demo.py` (seed 0) simulates one linear trajectory: 600 cells with
a known time `obs["true_time"]` in [0, 1], 4 stages cut at its quartiles
(`progenitor`, `early`, `mid`, `late`), and 40 genes that rise, fall or peak
along it. It log-normalises, then runs PCA, neighbours and UMAP, so the file
looks like sc-embed output. Obs order is shuffled.

`--demo` roots at `cell0005`, keeps the 3 chain edges
(progenitor-early, early-mid, mid-late) at threshold 0.1, and finds 0
unreachable cells. The Spearman correlation of DPT with `true_time` is 0.98
(2 decimals); the test computes it from `tables/pseudotime.tsv`, because real
data have no true time.

## Fixed plot spec

Figures come from the vendored `scripts/plots.py` (`plot_box`) and plain
matplotlib in the same theme, saved by `save_panel` as PDF, PNG and source CSV.
No titles; captions sit in `report.md`.

- `paga_graph`: groups at their median UMAP position (no layout RNG); lines
  join pairs at or above `--paga-threshold`, width 2 x connectivity; node size
  by group size.
- `umap_pseudotime`: cells coloured by DPT (viridis); unreachable cells grey.
- `pseudotime_by_group`: `plot_box`, one jittered point per cell (jitter 0.3,
  `--seed`), group colours from tab20 in group order. Under 10 cells in any
  group it draws points and a median, no box; the caption says which.
- Captions name how the view was made: the stored UMAP's neighbour count,
  metric, PCs and seed (from `uns["neighbors"]["params"]`; `min_dist` is not
  stored), and the DPT diffusion components and seed. Axis titles say the scale:
  "diffusion pseudotime (0-1, unitless)".
- `gene_trends`: one panel per gene, up to 3 per row, bin mean as a line and
  the 95% CI of the bin mean as a band (step 8).

`scripts/skillkit.py`, `plots.py`, `theme.py` and `palettes.py` are byte copies
of `templates/`; a test fails when they drift.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`;
row IDs from `mycelium-extra/plans/clawbio-determinism/parts/*.md`.

- A-an-11, PAGA on annotation groups then PAGA-initialised UMAP:
  `repos_new/spatial_skin_atlas/processing_scripts/3_subcluster_coarse.py:44-51`.
- A-an-31 and A-pl-18, PAGA with a hard-coded `THRESHOLD=0.1`,
  `edge_width_scale=2, max_edge_width=2`:
  `repos_new/skin_fibroblast_atlas/figures/fig5_sfig7_wounds.ipynb:1347,1430,1456`.
- C-an-35 and C-pl-21, `sc.tl.paga(groups='leiden')`, `sc.pl.paga(plot=False)`:
  `repos/multiple-myeloma/analyses-rna/03-subset-analyses/01-plasma_cells/01-py-paga_analysis.ipynb:134-145`.
- B-an-26, seeded diffusion map (`set.seed(123)`, `destiny::DiffusionMap`):
  `repos/IHA-Figure/Extended_Figure4/02_IMM_Age_Metric/01_Calculate_Scores.ipynb:60,473`.
- B-pl-20, distribution of pseudotime per group (ridges there, a box plot here):
  `repos/IHA-Figure/Extended_Figure9/04_FigD.ipynb:503`.
- B-pl-23, trend along pseudotime (LOESS with `se = FALSE` there; binned means
  with a stated 95% CI here): `repos/IHA-Figure/Extended_Figure9/07_FigG.ipynb:120`.

Checked against the installed reference skill `cellrank`
(`~/.claude/skills/cellrank/references/kernel-pseudotime.md:49-83`): the root
goes in `adata.uns["iroot"]` before `sc.tl.dpt`, as here. cellrank picks the
root from a diffusion-component extremum; this skill uses the first
`--root-group` cell instead, because a DC sign is arbitrary (step 5). Use
cellrank's `PseudotimeKernel` on `dpt_pseudotime` for fate probabilities.

## ClawBio mechanisms used

IDs from `mycelium-extra/plans/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data | `examples/demo.h5ad`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-11, F-13, F-14 | Fixed plot spec, threshold drawn, style, Agg | vendored `plots.py`, `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` via `skillkit.py` |
| F-23, F-26 | Thresholds documented, fixed root rule | parameter table; step 5; no defaults |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", "Trigger" |
| F-27 | One task | PAGA + DPT only; velocity and other engines elsewhere |
| F-28 | Golden-value demo test | `tests/test_trajectory.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

From `skills/trajectory`:

```bash
UV_CACHE_DIR=$TMPDIR/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with scanpy,igraph,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

The PyPI name is `igraph`; `sc.tl.paga` imports it.

## Citations

- Wolf et al. (2019). PAGA. Genome Biology 20:59.
- Haghverdi et al. (2016). Diffusion pseudotime. Nature Methods 13:845-848.
- Wolf, Angerer, Theis (2018). Scanpy. Genome Biology 19:15.
