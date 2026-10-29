---
name: velocity
description: >-
  RNA velocity and fate probabilities from a single-cell h5ad with spliced and
  unspliced layers: scVelo velocity (stochastic or dynamical), velocity
  pseudotime, and CellRank 2 VelocityKernel + GPCCA fate probabilities towards
  named or inferred terminal states, with tables, four fixed figures, a report
  and a reproducibility bundle. Fire on "RNA velocity", "scVelo", "CellRank",
  "velocity stream", "velocity pseudotime", "fate probabilities", "terminal
  states", "phase portraits", "spliced unspliced".
  Do NOT fire for diffusion pseudotime or PAGA without velocity (trajectory),
  data with no spliced/unspliced layers, clustering or embedding (sc-embed),
  or CellRank kernels other than velocity (pseudotime, CytoTRACE, real time).
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - single-cell
    - RNA velocity
    - scVelo
    - CellRank
    - fate probabilities
  inputs:
    - name: adata
      type: file
      format:
        - h5ad
      description: AnnData with spliced/unspliced layers, X_umap and a group column
      required: true
    - name: gene_id_col
      type: string
      description: var column of Ensembl gene ids (var_names are symbols), or none
      required: true
    - name: context
      type: string
      description: species, tissue or cell type, condition and dataset
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: parameters, summary and captioned figures
    - name: result
      type: file
      format:
        - json
      description: terminal states, top velocity genes, median pseudotime per group, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - scvelo
      - cellrank
      - anndata
      - scanpy
      - numpy<2.4
      - pandas<3
      - scipy
      - matplotlib
      - plotnine-extra
  demo_data:
    - path: examples/demo.h5ad
      description: scvelo.datasets.simulation, 500 cells, 50 genes, 4 stages, known time
  endpoints:
    cli: python scripts/velocity.py --adata {adata} ... --outdir {outdir}
  trigger_keywords:
    - rna velocity
    - scvelo
    - cellrank
    - velocity pseudotime
    - fate probabilities
    - terminal states
---

# Velocity

You are the dispatcher for **velocity**. The script does the analysis; you
collect its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/velocity.py`. Never re-implement the analysis, change a
threshold, or pick a value the user has not given. If a required argument is
unknown, ask the user (these are grill questions), or run `--demo` only when the
user asks for a demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

**Fire this skill when the user says any of:**
- "run scVelo on this object and show the velocity stream"
- "which fate do the <group> cells head towards"
- "CellRank fate probabilities towards <group> and <group>"
- "velocity pseudotime", "phase portraits of the top velocity genes"

**Do NOT fire when:**
- the object has no `spliced`/`unspliced` layers (run velocyto, STARsolo
  `--soloFeatures Velocyto` or alevin-fry first; that is not this skill);
- the user wants diffusion pseudotime or PAGA without velocity (trajectory);
- the user wants a pseudotime, CytoTRACE or real-time CellRank kernel;
- the object has no UMAP yet (run sc-embed first).

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--adata` | `.h5ad` | one object, cells x genes | `data/velocity_input.h5ad`, `examples/demo.h5ad` |
| `layers["spliced"]`, `layers["unspliced"]` | counts | from velocyto, STARsolo or alevin-fry, same genes as `X` | merged with `scv.utils.merge` |
| `obsm["X_umap"]` | 2-D UMAP | plotting only | sc-embed |
| `--groupby` | obs column | cluster or annotation labels, at least 2 | `leiden`, `stage` |
| `--gene-id-col` | var column | Ensembl gene ids, one per symbol; `var_names` are symbols | `gene_ids` |

`X` may hold raw or log-normalised values: scVelo normalises the layers itself
and skips `log1p` when `uns["log1p"]` is present.

## Parameters (all required, no defaults)

Ask these through the analysis contract: `${CLAUDE_PLUGIN_ROOT}/skills/lucid-coder/references/contracts.md`.

Each row is a grill question. Ask the user; never fill one in.

| Argument | Grill question | Demo-only value |
|---|---|---|
| `--adata` | Which h5ad holds the spliced and unspliced layers? | `examples/demo.h5ad` |
| `--gene-id-col` | Which var column holds the Ensembl gene ids (`none` if the source has none)? | `gene_ids` |
| `--mode` | Stochastic (steady-state, fast) or dynamical (full kinetics, slow) velocity? | `stochastic` |
| `--min-shared-counts` | How many counts must a gene have in spliced and unspliced together to stay? | `10` |
| `--n-top-genes` | How many highly variable genes go into velocity? | `30` |
| `--n-pcs` | How many PCs for the neighbour graph? | `10` |
| `--n-neighbors` | How many neighbours for the graph and the moments? | `30` |
| `--groupby` | Which obs column holds the groups? | `stage` |
| `--terminal-states` | Which groups are end points (comma list), or `auto` to let GPCCA infer them? | `s4` |
| `--n-states` | Only with `auto`: how many GPCCA macrostates? | not used |
| `--context` | Species, tissue or cell type, condition and dataset? | `synthetic demo: scVelo simulation, 500 cells, 4 stages, no real tissue` |
| `--reference` | Which genome build and annotation release named the genes (for example `GRCh38, GENCODE 44`)? | `synthetic demo, no genome build` |
| `--seed` | Which seed (PCA, neighbours)? | `0` |
| `--outdir` | Which new directory gets the outputs? It must not exist. | `velocity_demo` |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them.
With `--demo`, only `--outdir` may also be given; any other analysis argument is
an error.

## Workflow (fixed order)

1. **Validate** arguments and input: `--mode`, the integer ranges, the
   `--terminal-states`/`--n-states` pairing, the file, both layers,
   `obsm["X_umap"]`, `--groupby` in obs and every named terminal group in it.
   Gene names and cell ids must be unique (`skillkit.check_unique`).
   `skillkit.check_gene_ids` stops on a duplicated gene id; `--gene-id-col none`
   writes `gene_id` as NA, never a guessed id.
   Bad input gives `ERROR: <message>`, exit 2.
2. **Groups** become categorical; an existing category order is kept, else
   labels sort alphabetically. Unused levels are dropped; at least 2 groups.
3. **Genes**: `scv.pp.filter_genes(min_shared_counts=)`, then
   `scv.pp.normalize_per_cell`, `sc.pp.log1p` (skipped when `uns["log1p"]`
   exists, so sc-embed output is not logged twice), and
   `sc.pp.highly_variable_genes(n_top_genes=, subset=True)`. This is the
   CellRank tutorial order; scvelo 0.3.4's `filter_and_normalize` no longer
   takes `n_top_genes`. `--n-top-genes` above the genes left, or `--n-pcs` at or
   above min(cells, `--n-top-genes`), stops the run.
4. **Graph**: `sc.pp.pca(n_comps=--n-pcs, random_state=--seed)`,
   `sc.pp.neighbors(n_pcs=, n_neighbors=, random_state=--seed)` on the HVG
   subset, then `scv.pp.moments(n_pcs=None, n_neighbors=None)` on that graph.
   Any neighbour graph in the input is replaced; the input `X_umap` is used only
   for plotting.
5. **Velocity**: `--mode dynamical` first runs
   `scv.tl.recover_dynamics(n_jobs=1)` (slow: about 1 min per 500 cells x 30
   genes here); then `scv.tl.velocity(mode=)`.
6. **Velocity graph and pseudotime**: `scv.tl.velocity_graph(n_jobs=1)`,
   `scv.tl.velocity_pseudotime()` (scVelo picks root and end cells from the
   velocity graph), `scv.tl.velocity_embedding(basis="umap")`.
7. **Kernel and terminal states**: `cr.kernels.VelocityKernel` with
   `compute_transition_matrix(model="deterministic")`, then `GPCCA`.
   - A comma list: `set_terminal_states(states={group: all its cells})`, so
     every cell of each named group is absorbing.
   - `auto`: `compute_macrostates(n_states=--n-states, cluster_key=--groupby)`,
     then `predict_terminal_states(method="stability", stability_threshold=0.96)`,
     CellRank's documented rule. Macrostates are named after their majority
     group (`s1_2` when one group holds several).
8. **Fate probabilities**: `compute_fate_probabilities(solver="gmres",
   use_petsc=False)`, one thread, tolerance 1e-6. The mean per group gives
   `fate_by_group.tsv`.
9. **Top genes** for the phase portraits: the 4 highest `velocity_r2`
   (stochastic) or `fit_likelihood` (dynamical), ties by gene name.

Required checks that do not apply: no design is fitted (no rank check), no
resource is joined (no `n_matched`), and no cell is dropped. Step 3 filters genes
over all cells, not per group, so `result.json` reports `n_genes_input` beside
`n_genes`.

## Output structure

`velocity_pseudotime.tsv` and `fate_probabilities.tsv` are sorted by cell;
`velocity_genes.tsv` by `gene_symbol`; `fate_by_group.tsv` by group, then terminal state.

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # context, n_cells, n_genes, n_genes_input, status_counts, gene_id_col, mode, terminal_states, top_genes, median_velocity_pseudotime, params, seed, reference, chat_summary_lines, preferred_artifacts
├── data/
│   └── velocity.h5ad          # HVG subset with Ms/Mu/velocity layers, velocity_graph, velocity_pseudotime, obsm lineages_fwd (plain array; names in uns lineages_fwd_names)
├── tables/
│   ├── velocity_genes.tsv     # gene_id, gene_symbol, scVelo velocity_* (stochastic) or fit_* and velocity_genes (dynamical) columns, status
│   ├── velocity_pseudotime.tsv  # cell, <groupby>, velocity_pseudotime
│   ├── fate_probabilities.tsv   # cell, <groupby>, one column per terminal state
│   └── fate_by_group.tsv      # <groupby>, terminal_state, mean_fate_probability, n_cells
├── figures/
│   ├── velocity_stream.pdf
│   ├── velocity_stream.png
│   ├── velocity_stream.csv
│   ├── umap_velocity_pseudotime.pdf
│   ├── umap_velocity_pseudotime.png
│   ├── umap_velocity_pseudotime.csv
│   ├── fate_by_group.pdf
│   ├── fate_by_group.png
│   ├── fate_by_group.csv
│   ├── phase_portraits.pdf
│   ├── phase_portraits.png
│   └── phase_portraits.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python and package versions
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

`data/velocity.h5ad` is not checksummed: h5ad bytes are not stable.

`velocity_genes.tsv` has one row per input gene. `status` is `tested` for the
genes the model fit (the HVG subset), `not_detected` for genes with zero
spliced and unspliced counts, and `not_tested` for genes removed by
`--min-shared-counts` or the HVG cut; untested rows have NA statistics.
`result.json` `status_counts` counts each status.

## Gotchas

- **The model will want to pick terminal states from the stream plot.** Do
  not. They are the user's answer: named groups, or `auto` with the user's
  `--n-states`. Report which rule ran.
- **The model will want to trust `auto` terminal states as-is.** Do not. GPCCA
  can mark a stable early state as terminal (the demo's `s1_2`). Show the user
  the states next to the velocity pseudotime and ask before relying on them.
- **The model will want to lower `--min-shared-counts` or raise `--n-pcs`
  until the run passes.** Do not. Relay the `ERROR:` line and ask.
- **The model will want to switch to `--mode dynamical` "for accuracy"** on a
  big object. Do not without asking: it is many times slower. It is the user's
  choice.
- **The model will want to `pip install -U scvelo numpy pandas`** when scvelo
  crashes. Do not. Run with the pinned versions in the uv line below.
- **The model will want to read velocity pseudotime as real time** or fate
  probability as a lineage proof. Do not. Both come from a Markov chain built
  on the velocity graph; say so.
- **The model will want to describe the stream arrows as cells moving.** Do
  not. Arrows are a model projection of the velocity graph onto an input UMAP,
  not measured movement; the UMAP was not recomputed here.
- **The model will want to read a fate probability or `velocity_r2` as the
  chance a cell becomes that state.** Do not. Both are model scores, not
  calibrated probabilities of an outcome.
- **The model will want to call a top velocity gene a driver of the fate.** Do
  not. `velocity_r2` and `fit_likelihood` score how well the splicing model fits
  a gene; say "scored" or "predicted", never "drives" or "regulates"
  (house-style `claims.strength`).
- **The model will want to call a filtered gene "unchanged" or "not
  dynamic".** Do not. A gene removed by `--min-shared-counts` or outside the
  HVG subset was not tested.
- **The model will want to guess `--reference`** from the gene names. Do not.
  The genome build and annotation release are a grill question.
- **The model will want to re-run `scv.pl.*` in-session** to tweak a figure.
  Do not. The tables are the result; a second run has no provenance.

## CLI reference

```bash
python scripts/velocity.py \
  --adata data/velocity_input.h5ad --mode stochastic \
  --min-shared-counts 20 --n-top-genes 2000 --n-pcs 30 --n-neighbors 30 \
  --groupby leiden --terminal-states Alpha,Beta --gene-id-col gene_ids \
  --context "mouse pancreas, E15.5 endocrinogenesis, Bastidas-Ponce 2019" \
  --reference "GRCm39, GENCODE M33" \
  --seed 0 --outdir results/2026-10-27_velocity

# Inferred terminal states
python scripts/velocity.py ... --terminal-states auto --n-states 6 ...

# Demo (only when the user asks for one)
python scripts/velocity.py --demo --outdir demo_out
```

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
