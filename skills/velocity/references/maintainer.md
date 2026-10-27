# velocity: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`
(`FB` = `repos_new/skin_fibroblast_atlas`); row IDs from
`mycelium-extra/plans/clawbio-determinism/parts/A_skin.md:37,38,68`.

- A-an-29, scVelo merge, moments, velocity, velocity pseudotime:
  `FB/integration/prep_velocity_a.ipynb:922,1342,1429,1439`;
  `FB/figures/fig5_sfig7_velocity.ipynb:175,242`.
- A-an-30, CellRank `VelocityKernel(adata)`, `compute_transition_matrix()`:
  `FB/figures/fig5_sfig7_velocity.ipynb:496-497,529`.
- A-pl-19, velocity stream and pseudotime scatter on the UMAP:
  `FB/figures/fig5_sfig7_velocity.ipynb:37,177-186,348,395-408`.

Checked against the installed reference skill `cellrank`
(`~/.claude/skills/cellrank/SKILL.md`, `references/kernel-rna-velocity.md`,
`references/estimator-initial-terminal.md`,
`references/estimator-fate-probabilities.md`): `VelocityKernel` then `GPCCA`,
`compute_macrostates`, `predict_terminal_states` or `set_terminal_states`, and
`compute_fate_probabilities`, as in steps 7-8.

## ClawBio mechanisms used

IDs from `mycelium-extra/plans/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data | `examples/demo.h5ad`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-13, F-14 | Fixed plot spec, style, Agg | vendored `plots.py`, `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` via `skillkit.py` |
| F-23 | Thresholds documented, no defaults | parameter table; step 7 names the fixed 0.96 rule |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", "Agent boundary" |
| F-27 | One task | velocity + fates only; DPT/PAGA in trajectory |
| F-28 | Golden-value demo test | `tests/test_velocity.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

From `skills/velocity`:

```bash
UV_CACHE_DIR=$TMPDIR/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with 'scvelo,cellrank,numpy<2.4,pandas<3,pytest' \
  python -m pytest -q -p no:cacheprovider tests
```

The pins are required. scvelo 0.3.4 stochastic velocity fails on numpy >= 2.4
(`leastsq_generalized` assigns a 1-element array to a scalar), and dynamical
mode fails on pandas 3 (`make_unique_list` passes a list to `pd.unique`).
Tested with scvelo 0.3.4, cellrank 2.3.3, scanpy 1.12.4, numpy 2.3, pandas 2.3.
The suite takes a few minutes: one demo run, one `auto` run, and a dynamical
smoke run on a 150-cell subset.

## Why this exists

- **Without it**: the agent writes a fresh scVelo notebook each time. The gene
  filter, HVG count and neighbour graph drift per run, terminal states get
  picked from a plot, and scvelo crashes on current numpy/pandas.
- **With it**: one tested script with pinned versions, a fixed step order and
  fixed figures. Every analysis choice is a named argument in `commands.sh`.

## Demo

`examples/make_demo.py` (seed 0) runs `scvelo.datasets.simulation` offline:
500 cells on one time axis `obs["true_t"]`, 50 genes, cut into 4 stages
(`s1`..`s4`) at quartiles of `true_t`. X and the layers hold the simulated
abundances; `X_umap` comes from a log-normalised copy. The file is 318 kB.

`--demo` keeps 30 genes, names `s4` as the terminal state, and gives every
stage a mean fate probability of 1.00 towards it (one lineage). The Spearman
correlation of velocity pseudotime with `true_t` is 0.996 (3 decimals); the test
computes it from `tables/velocity_pseudotime.tsv`, because real data have no
true time. Top velocity genes: G08, G20, G28, G44. With `--terminal-states auto
--n-states 3` the demo returns `s1_2` and `s4`: GPCCA also marks a stable
early macrostate, which is why `auto` results need checking against the
pseudotime.

## Fixed plot spec

Figures come from the vendored `scripts/plots.py` (`plot_composition`),
`scv.pl.velocity_embedding_stream`, and plain matplotlib in the same theme,
saved by `save_panel` as PDF, PNG and source CSV. No titles; captions sit in
`report.md`. The UMAP captions name the PCs, highly variable genes,
neighbours and seed of the velocity graph, and say the UMAP came from the input
(no `min_dist` here). Group colours come from tab20 in group order.

- `velocity_stream`: cells coloured by group, velocity streamlines on the UMAP;
  the CSV holds UMAP and projected velocity per cell.
- `umap_velocity_pseudotime`: cells coloured by velocity pseudotime (viridis).
- `fate_by_group`: stacked bar of mean fate probability per group, one segment
  per terminal state (Set2). Probabilities lie in [0, 1], so the diverging
  `plot_heatmap` scale does not fit.
- `phase_portraits`: Ms against Mu for the 4 top genes, coloured by group;
  the dashed line is the steady-state ratio (`velocity_gamma`, or
  `fit_scaling * fit_gamma / fit_beta` for dynamical).

`scripts/skillkit.py`, `plots.py`, `theme.py` and `palettes.py` are byte copies
of `templates/`; a test fails when they drift.

## Citations

- Bergen et al. (2020). scVelo. Nature Biotechnology 38:1408-1414.
- Weiler et al. (2024). CellRank 2. Nature Methods 21:1196-1205.
- La Manno et al. (2018). RNA velocity of single cells. Nature 560:494-498.
