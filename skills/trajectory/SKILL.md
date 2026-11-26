---
name: trajectory
description: >-
  Trajectory and pseudotime from an embedded single-cell h5ad: PAGA group
  connectivity and diffusion pseudotime (scanpy DPT) from a fixed root cell,
  with tables, four fixed figures, a report and a reproducibility bundle. Fire
  on "pseudotime", "trajectory", "PAGA", "diffusion pseudotime", "DPT",
  "order cells along differentiation", "gene trends along pseudotime".
  Do NOT fire for RNA velocity (scVelo, CellRank), Slingshot or Monocle3 runs,
  clustering or embedding (sc-embed), or differential abundance.
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - single-cell
    - trajectory
    - pseudotime
    - PAGA
    - diffusion pseudotime
  inputs:
    - name: adata
      type: file
      format:
        - h5ad
      description: embedded AnnData with neighbours, X_umap and a group column
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
      description: root cell, PAGA edges, median pseudotime per group, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - anndata
      - scanpy
      - igraph
      - numpy
      - pandas
      - scipy
      - matplotlib
      - plotnine-extra
  demo_data:
    - path: examples/demo.h5ad
      description: synthetic linear trajectory, 600 cells, 4 stages, 40 genes, known time
  endpoints:
    cli: python scripts/trajectory.py --adata {adata} ... --outdir {outdir}
  trigger_keywords:
    - pseudotime
    - trajectory
    - paga
    - diffusion pseudotime
    - dpt
    - gene trends along pseudotime
---

# Trajectory

You are the dispatcher for **trajectory**. The script does the analysis; you
collect its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/trajectory.py`. Never re-implement the analysis, change a
threshold, or pick a value the user has not given. If a required argument is
unknown, ask the user (these are grill questions), or run `--demo` only when the
user asks for a demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--adata` | `.h5ad` | log-normalised values in `X` (sc-embed output, not raw counts) | `data/embedded.h5ad`, `examples/demo.h5ad` |
| neighbour graph | `uns["neighbors"]`, `obsp["connectivities"]` | from `sc.pp.neighbors` | sc-embed |
| embedding | `obsm["X_umap"]` | 2-D UMAP | sc-embed |
| `--groups` | obs column | cluster or annotation labels, at least 2 | `leiden`, `stage` |

The root cell does not depend on the cell order in `obs`: it is a named cell
(`--root-cell`) or comes from a rule on the diffusion map (`--root-rule`, step 5).

## Parameters (all required, no defaults)

Ask these through the analysis contract: `${CLAUDE_PLUGIN_ROOT}/skills/lucid-coder/references/contracts.md`.

Each row is a grill question. Ask the user; never fill one in.

| Argument | Grill question | Demo-only value |
|---|---|---|
| `--adata` | Which embedded h5ad holds the cells? | `examples/demo.h5ad` |
| `--groups` | Which obs column holds the groups for PAGA? | `stage` |
| `--root-cell` (or `--root-rule`) | Which cell (an `obs_names` entry) is the root? Give exactly one of `--root-cell` and `--root-rule`; never both, never neither. | `cell0005` |
| `--root-rule` (or `--root-cell`) | Pick the root by rule instead: `diffusion-extreme` takes the cell of `--root-group` furthest out on diffusion component 1 (step 5). Needs `--root-group`. | not used |
| `--root-group` | With `--root-rule` only: which group holds the start of the trajectory? Refused with `--root-cell`. | not used |
| `--n-dcs` | How many diffusion components (diffmap `n_comps` and DPT `n_dcs`)? | `10` |
| `--paga-threshold` | Which PAGA connectivity must a group pair reach to count as an edge? | `0.1` |
| `--genes` | Which genes get a trend plot? A comma list, or a file with one gene per line. | `G00,G01,G02` |
| `--reference` | Which genome build and annotation release name the genes (for example `GRCh38, GENCODE 44`)? Recorded in `result.json`, never guessed. | `synthetic demo, no genome build` |
| `--context` | Species, tissue or cell type, condition and dataset? Printed first in the report and the chat summary. | `synthetic demo: 600 cells, 4 stages along one simulated path, no real tissue` |
| `--seed` | Which seed (diffusion map `random_state`, box-plot jitter)? | `0` |
| `--outdir` | Which new directory gets the outputs? It must not exist. | `trajectory_demo` |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them.
With `--demo`, only `--outdir` may also be given; any other analysis argument is
an error.

## CLI reference

```bash
python scripts/trajectory.py \
  --adata data/embedded.h5ad --groups leiden \
  --root-rule diffusion-extreme --root-group 3 \
  --n-dcs 10 --paga-threshold 0.1 --genes genes.txt \
  --reference "GRCh38, GENCODE 44" \
  --context "human bone marrow CD34+, healthy donors, Setty 2019" --seed 0 --outdir results/2026-10-27_trajectory

# Demo (only when the user asks for one)
python scripts/trajectory.py --demo --outdir demo_out
```

## Layer contract

- Input: `X` holds log-normalised values (sc-embed's `embedded.h5ad`), used only
  for the gene trends; the neighbour graph and `X_umap` come from sc-embed.
  `layers["counts"]` is not read.
- Output `data/trajectory.h5ad`: the input with `X` and layers unchanged, plus
  `obs["dpt_pseudotime"]`, `uns["paga"]`, `uns["iroot"]` and `obsm["X_diffmap"]`.

## Fixed method settings

Not asked, recorded here. `sc.tl.paga` (scanpy 1.12.4): `model="v1.2"` default.
`sc.tl.diffmap`: package defaults apart from `n_comps` and `random_state`.
`sc.tl.dpt(n_branchings=0)`. Gene trends: 20 equal-width pseudotime bins with a
t-based 95% CI (plot only).

## Workflow (fixed order)

1. **Validate** arguments and input: the file, `--groups` in obs, a neighbour
   graph (`uns["neighbors"]`, `obsp["connectivities"]`), `obsm["X_umap"]`,
   every `--genes` entry in `var_names`, and `--n-dcs` below the cell count.
   Gene names and cell ids must be unique (`skillkit.check_unique`).
   Bad input gives `ERROR: <message>`, exit 2.
2. **Groups** become categorical; an existing category order is kept, else
   labels sort alphabetically. Unused levels are dropped. At least 2 groups, and
   `--root-group` (with `--root-rule`) must be one of them, and `--root-cell`
   must be an `obs_names` entry.
3. **PAGA**: `sc.tl.paga(adata, groups=--groups)` on the stored neighbour graph.
4. **Diffusion map**: `sc.tl.diffmap(n_comps=--n-dcs, random_state=--seed)`.
5. **Root cell**: `--root-cell` sets it directly. `--root-rule
   diffusion-extreme` uses DC1 = `obsm["X_diffmap"][:, 1]` (column 0 is the
   trivial constant component). DC1's sign is arbitrary, so the side comes from
   the data: if the median DC1 of `--root-group` is at or above the median DC1 of
   all cells, the side is `+1`, else `-1`. The root is the `--root-group` cell
   with the largest `side * DC1`; an exact tie goes to the smallest `obs_names`
   entry. Flipping the eigenvector sign picks the same cell, and reordering
   `obs` does not change it. The root is stored in `uns["iroot"]`.
6. **DPT**: `sc.tl.dpt(n_dcs=--n-dcs)` from that root, no branching. Cells not
   reachable from the root get `inf`: they stay `inf` in the table, are counted
   in the summary and per group (`dropped_per_group`), and are left out of the
   box plot and gene trends.
7. **Edges**: a group pair is kept when its PAGA connectivity is
   `>= --paga-threshold`.
8. **Gene trends**: for each `--genes` entry, the mean of `X` in 20 equal-width
   bins from 0 to the largest finite pseudotime, with a t-based 95% confidence
   interval of each bin mean (`n - 1` df; a bin with one cell has no interval).

## Output structure

`pseudotime.tsv` is sorted by cell; `paga_connectivities.tsv` holds every group
pair once (upper triangle), in group order.

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # n_cells, n_groups, root_cell, root_group, root_rule, n_unreachable, dropped_per_group, n_edges_above_threshold, edges_above_threshold, median_pseudotime, genes, params, seed, context, reference, chat_summary_lines, preferred_artifacts
├── data/
│   └── trajectory.h5ad        # input plus obs dpt_pseudotime, uns paga and iroot, obsm X_diffmap
├── tables/
│   ├── pseudotime.tsv         # cell, <groups>, dpt_pseudotime, is_root
│   └── paga_connectivities.tsv  # group_a, group_b, connectivity, above_threshold
├── figures/
│   ├── paga_graph.pdf
│   ├── paga_graph.png
│   ├── paga_graph.csv
│   ├── umap_pseudotime.pdf
│   ├── umap_pseudotime.png
│   ├── umap_pseudotime.csv
│   ├── pseudotime_by_group.pdf
│   ├── pseudotime_by_group.png
│   ├── pseudotime_by_group.csv
│   ├── gene_trends.pdf
│   ├── gene_trends.png
│   └── gene_trends.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python and package versions
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

`data/trajectory.h5ad` is not checksummed: h5ad bytes are not stable.

## Out of scope

- **Slingshot and Monocle3** are engines to add later, as R scripts behind an
  `--engine` flag. Until then, say this skill runs scanpy PAGA and DPT only.
- **RNA velocity** (scVelo, CellRank) is a separate skill.
- **Branch detection** (`n_branchings > 0`) is off; DPT returns one ordering.

## Gotchas

- **The model will want to pick the root cell by eye** or by the lowest DC1.
  Do not. Ask the user for `--root-cell`, or for `--root-rule` with a root
  group; DC1's sign is arbitrary, which is why the rule sets the side from the
  group's median. Report `root_cell` and `root_rule`.
- **The model will want to lower `--paga-threshold`** until the graph looks
  connected. Do not. The threshold is the user's answer; report the edges found.
- **The model will want to re-run `sc.tl.dpt` or `sc.pl.paga` in-session** to
  tweak a figure. Do not. The tables are the result; a second run has no
  provenance.
- **The model will want to read pseudotime as real time** or as a lineage
  claim. Do not. DPT orders cells by diffusion distance from one root; say so,
  and name the root.
- **The model will want to say a gene whose trend rises along pseudotime
  "drives" the trajectory** or that the root group "gives rise to" the late
  groups. Do not. Pseudotime and PAGA edges are model predictions from one
  snapshot; a trend is "associated with" the predicted ordering
  (house-style `claims.strength`).
- **The model will want to drop `inf` pseudotime cells silently** or fill them.
  Do not. They are counted in `n_unreachable`; tell the user.
- **The model will want to read PAGA connectivity as a transition
  probability** or a lineage claim. Do not. It is a score of how much two
  groups' neighbour graphs connect, not a probability or a direction.
- **The model will want to call a gene "off" where its trend is near 0.** Do
  not. A low bin mean is low measured expression; non-detection is not absence.
- **The model will want to fill in `--reference`** from the gene names. Do not.
  The build and release are the user's answer.
- **The model will want to feed raw counts** in `X`. Do not. Trends are means of
  `X`, which must hold log-normalised values (sc-embed output).

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
