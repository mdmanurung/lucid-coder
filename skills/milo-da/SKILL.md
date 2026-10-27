---
name: milo-da
description: >-
  Milo differential abundance on a kNN graph with pertpy (pt.tl.Milo, pydeseq2
  solver): seeded index-cell sampling, counts per sample per neighbourhood, NB
  GLM per neighbourhood, SpatialFDR, annotation with a Mixed label. Writes
  data/milo.h5mu, a neighbourhood table, significant counts per annotation,
  three fixed figures, a report and a reproducibility bundle. Fire on "Milo",
  "differential abundance on the kNN graph", "neighbourhood DA",
  "which cell states expand in disease", "milopy", "nhood graph", "DA
  beeswarm". Do NOT fire for cell-type proportion tests on labels (scCODA,
  crumblr, CLR tests), differential expression, or when cells still need
  embedding (sc-embed first).
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - single-cell
    - differential abundance
    - Milo
    - kNN graph
    - SpatialFDR
  inputs:
    - name: adata
      type: file
      format:
        - h5ad
      description: embedded cells with obsm embedding, X_umap, kNN graph and sample, design and annotation columns
      required: true
  outputs:
    - name: milo
      type: file
      format:
        - h5mu
      description: MuData with "rna" (cells, nhoods in obsm) and "milo" (samples x nhoods counts, DA results in var)
    - name: report
      type: file
      format:
        - md
      description: params, summary and captioned figures
    - name: result
      type: file
      format:
        - json
      description: nhood and significant counts, dominant annotation, params, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - pertpy
      - pydeseq2
      - scanpy
      - anndata
      - mudata
      - numpy
      - pandas
      - scipy
      - matplotlib
      - threadpoolctl
  demo_data:
    - path: examples/demo.h5ad
      description: synthetic, 6 samples x 200 cells in 2 conditions, 4 cell types, Mono 3x in B
  endpoints:
    cli: python scripts/milo_da.py --adata {adata} ... --outdir {outdir}
  trigger_keywords:
    - Milo
    - differential abundance
    - neighbourhood DA
    - SpatialFDR
    - nhood graph
---

# milo-da

You are the dispatcher for **milo-da**. The script does the analysis; you
collect its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/milo_da.py`. Never re-implement the analysis, change a threshold,
or pick a value the user has not given. If a required argument is unknown, ask
the user (these are grill questions), or run `--demo` only when the user asks
for a demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--adata` | `.h5ad` | unique cell ids; `obsm[--use-rep]`, `obsm["X_umap"]`; a kNN graph in `obsp` (rebuilt if absent or different) | sc-embed `data/embedded.h5ad`, `examples/demo.h5ad` |
| sample column | obs column | `--sample-col`, one level per biological replicate, no missing values | `sample` with `A1`-`B3` |
| donor column | obs column | `--donor-col`, the donor of each sample, no missing values | `donor` with `D1`-`D6` |
| design columns | obs columns | every column in `--design`, constant within a sample, no missing values | `condition` with `A`, `B` |
| annotation column | obs column | `--annotation-col`, cell labels, no missing values | `cell_type` |

## Parameters (all required, no defaults)

Each row is a grill question. Ask the user; never fill one in.

| Argument | Grill question | Demo-only value |
|---|---|---|
| `--adata` | Which embedded h5ad holds the cells and graph? | `examples/demo.h5ad` |
| `--sample-col` | Which obs column is the replicate (sample, donor)? | `sample` |
| `--donor-col` | Which obs column is the donor? A donor with several samples needs it in `--design`. | `donor` |
| `--design` | Which design formula, with the test variable last (`~batch + condition`)? | `~condition` |
| `--model-contrasts` | Which contrast, naming both levels (`conditionB-conditionA`)? | `conditionB-conditionA` |
| `--prop` | What fraction of cells to sample as index cells (0.05-0.2 is usual)? | `0.1` |
| `--use-rep` | Which obsm embedding defines neighbourhoods (`X_pca`, `X_pca_harmony`, `X_scVI`)? | `X_pca` |
| `--n-neighbors` | Which k for the graph? A different k from the input graph rebuilds it. | `30` |
| `--spatial-fdr-max` | Which SpatialFDR cutoff calls a neighbourhood significant? | `0.1` |
| `--annotation-col` | Which obs column labels the neighbourhoods? | `cell_type` |
| `--min-annotation-frac` | Under what top-label fraction is a neighbourhood `Mixed`? | `0.6` |
| `--seed` | Which seed for index-cell sampling, a rebuilt graph and the plot jitter? | `0` |
| `--outdir` | Which new directory gets the outputs? It must not exist. | `milo_da_demo` |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them.
With `--demo`, only `--n-neighbors` and `--outdir` may also be given; any other
analysis argument is an error.

## Workflow (fixed order)

1. **Validate** arguments and input; bad input gives `ERROR: <message>`, exit 2.
   Cell ids must be unique. Every design column must take one value per
   `--sample-col` level. On the per-sample rows, the design terms must not be
   aliased (`check_design_rank`), and a donor with several samples needs the
   donor column as a design term (`check_replication`).
2. **kNN graph**: keep the input graph when `uns["neighbors"]["params"]` has
   the same `use_rep` and `n_neighbors`; otherwise rebuild it with
   `sc.pp.neighbors(n_neighbors=--n-neighbors, use_rep=--use-rep)` on all
   columns of the representation, seeded. The summary says which.
3. **Neighbourhoods**: `make_nhoods(prop=--prop, seed=--seed)` samples
   `round(prop * n_cells)` index cells, then refines each to the cell nearest
   its neighbours' median in `--use-rep`.
4. **Count** cells per sample per neighbourhood: `count_nhoods(sample_col=)`.
5. **Test**: `da_nhoods(design=--design, model_contrasts=--model-contrasts,
   solver="pydeseq2")`. With a contrast the design drops its intercept, so the
   contrast names both levels (`conditionB-conditionA`). SpatialFDR is the
   Benjamini-Hochberg adjustment weighted by 1 / distance to the k-th neighbour.
6. **Annotate**: `annotate_nhoods(anno_col=--annotation-col)` gives each
   neighbourhood its most frequent label; a neighbourhood whose top label is
   under `--min-annotation-frac` of its cells becomes `Mixed`. Then
   `build_nhood_graph(basis="X_umap")`.

BLAS runs on one thread (`threadpoolctl`). pydeseq2 fans its dispersion fits out
over all CPUs; each neighbourhood is fitted alone, so results do not change.

## Output structure

`nhoods.tsv` is sorted by SpatialFDR, then PValue, then nhood;
`da_by_annotation.tsv` by annotation. The data file is **h5mu** (MuData).

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # n_cells, n_samples, n_donors, n_nhoods, n_significant, n_up, n_down, dominant_annotation, significant_by_annotation, graph_rebuilt, solver, seed, params, chat_summary_lines, preferred_artifacts
├── data/
│   └── milo.h5mu              # "rna": cells, obsm nhoods, obs nhood_ixs_*; "milo": samples x nhoods counts, var logFC, PValue, SpatialFDR, nhood_annotation(_frac), Nhood_size
├── tables/
│   ├── nhoods.tsv             # nhood, index_cell, size, logFC, PValue, SpatialFDR, annotation, annotation_frac
│   └── da_by_annotation.tsv   # annotation, n_nhoods, n_up, n_down (SpatialFDR < cutoff)
├── figures/
│   ├── nhood_graph.pdf
│   ├── nhood_graph.png
│   ├── nhood_graph.csv
│   ├── da_beeswarm.pdf
│   ├── da_beeswarm.png
│   ├── da_beeswarm.csv
│   ├── nhood_size.pdf
│   ├── nhood_size.png
│   └── nhood_size.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python and package versions
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

## Gotchas

- **The model will want to drop `--seed`** or call `make_nhoods` again in-session
  to "try another sample". Do not. The index cells are a seeded random draw;
  rerun the script into a new outdir so each draw has provenance.
- **The model will want to treat each cell or each condition as a replicate.**
  Do not. `--sample-col` must be the biological replicate; with one sample per
  condition the dispersion is not estimable and pertpy refuses the model.
- **The model will want to treat two samples of one donor as independent.** Do
  not. A donor with several samples stops the run until `--design` has the
  donor column as a fixed term (the pydeseq2 solver has no `(1|donor)`). With
  one sample per donor, leave the donor out: `~donor + condition` is then
  rank-deficient and stops the run.
- **The model will want to read a significant Mono neighbourhood as "Mono
  expand"** and a single down neighbourhood elsewhere as a real loss. Do not
  over-read either. A proportion shift is relative (compositional): counts are
  normalised per sample, so one population's expansion pushes others down. The
  report names the contrast direction; an absolute change needs cell counts
  per volume.
- **The model will want to loosen `--spatial-fdr-max`** or switch to raw
  PValue when nothing is significant. Do not. The cutoff is the user's call;
  report zero hits as zero.
- **The model will want to lower `--prop` or `--n-neighbors`** to speed things
  up. Do not change them silently: they set the neighbourhood size
  (`nhood_size.png`), which the miloR workflow says should peak around 50-100
  cells. Ask the user.
- **The model will want to call `Mixed` neighbourhoods a new cell type.** Do
  not. `Mixed` means no label reaches `--min-annotation-frac`.

## Engine scope

pertpy `pt.tl.Milo` with `solver="pydeseq2"` (pertpy 1.4.0, pydeseq2 0.5.4):
pure Python, no R. It normalises with TMM like miloR, fits dispersions with
pydeseq2 and tests the contrast with a likelihood ratio test. This comes close
to miloR's edgeR quasi-likelihood F-test but does not reproduce it exactly
(pertpy docstring). The edgeR solver needs rpy2 and R edgeR, and R miloR does
not build on this host (beachmat needs a newer C++ compiler), so both are out
of scope. Random-effect designs (`(1|donor)`) and `de_nhoods` (miloDE) are
also out of scope.

## CLI reference

```bash
python scripts/milo_da.py \
  --adata results/2026-10-27_embed/data/embedded.h5ad \
  --sample-col sample --donor-col donor --design "~donor + condition" \
  --model-contrasts conditiondisease-conditionhealthy \
  --prop 0.1 --use-rep X_pca_harmony --n-neighbors 30 \
  --spatial-fdr-max 0.1 --annotation-col cell_type --min-annotation-frac 0.6 \
  --seed 0 --outdir results/2026-10-27_milo

# Demo (only when the user asks for one); add --n-neighbors 20 to see a rebuild
python scripts/milo_da.py --demo --outdir demo_out
```

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
