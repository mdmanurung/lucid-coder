---
name: ccc
description: >-
  Ligand-receptor cell-cell communication between cell types from a
  log-normalised single-cell h5ad with LIANA+ rank_aggregate (CellPhoneDB,
  Connectome, log2FC, NATMI, SingleCellSignalR, or a subset), a bundled or
  exported ligand-receptor resource TSV, a fixed pass rule on magnitude and
  specificity ranks, a source x target count heatmap, an LR dot plot, a report
  and a reproducibility bundle. Optional per-sample runs. Fire on "cell-cell
  communication", "ligand-receptor analysis", "which cell types talk to which",
  "LIANA", "CellPhoneDB", "rank_aggregate", "CCC between clusters".
  Do NOT fire for spatial CCC (bivariate, MISTy), for factorising CCC across many
  samples (Tensor-Cell2cell, MOFAtalk), for DE-driven CCC between conditions, for
  non-human data without a translated resource, or for raw counts.
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - single-cell
    - cell-cell communication
    - ligand-receptor
    - liana
    - cellphonedb
  inputs:
    - name: context
      type: string
      description: species, tissue or cell type, condition and dataset
      required: true
    - name: gene_id_col
      type: string
      description: var column of Ensembl gene ids, or none when the source has none
      required: true
    - name: adata
      type: file
      format:
        - h5ad
      description: cell-level AnnData, log-normalised values in X, a cell-type obs column
      required: true
    - name: resource
      type: file
      format:
        - tsv
      description: ligand-receptor pairs with columns ligand, receptor
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
      description: counts, top interactions, methods, params, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - "liana>=1.10"  # liana 1.8.0 crashes on anndata 0.13.4
      - anndata
      - numpy
      - pandas
      - scipy
      - matplotlib
      - plotnine-extra
      - threadpoolctl
  demo_data:
    - path: examples/demo.h5ad
      description: synthetic, 3 cell types x 40 cells, 2 samples, 60 genes
    - path: examples/demo_resource.tsv
      description: 20 ligand-receptor pairs, 2 planted sender -> receiver
  endpoints:
    cli: python scripts/ccc.py --adata {adata} --resource {resource} ... --outdir {outdir}
  trigger_keywords:
    - cell-cell communication
    - ligand-receptor
    - liana
    - cellphonedb
    - rank_aggregate
---

# Cell-cell communication

You are the dispatcher for **ccc**. The script does the analysis; you collect
its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/ccc.py`. Never re-implement the analysis, change a threshold, or
pick a value the user has not given. If a required argument is unknown, ask the
user (these are grill questions), or run `--demo` only when the user asks for a
demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--adata` | `.h5ad` | log-normalised values (`normalize_total` + `log1p`) in `X`; not raw counts, not scaled; `var_names` are gene symbols, `--gene-id-col` names the Ensembl id column | sc-annotate's output h5ad; `examples/demo.h5ad` |
| `--groupby` | obs column | at least 2 cell-type labels | `cell_type`; demo `cell_type` |
| `--sample` (optional) | obs column | sample or donor ID; one LIANA run per sample | demo `sample` (not used by `--demo`) |
| `--resource` | TSV | columns `ligand`, `receptor`; complex subunits joined by `_`; same identifiers as `var_names` | `consensus_resource.tsv`; `examples/demo_resource.tsv` |

The script refuses a matrix with negative values (scaled) or only integers (raw
counts), and a resource with no pair fully present in `var_names` (wrong species
or identifiers). `.raw` is never read (`use_raw=False`). A categorical
`--groupby` keeps its category order; any other column is cast to strings and
sorted. Gene names and cell ids must be unique (`skillkit.check_unique`).

The resource is a file so the run never downloads and never depends on the
liana version's bundled copy. liana ships its consensus resource offline
(`liana/resource/omni_resource.csv`); export it once and pass the file:

```python
import liana as li

li.rs.select_resource("consensus").to_csv(
    "consensus_resource.tsv", sep="\t", index=False
)  # human symbols, 4620 pairs in liana 1.10.0; "mouseconsensus" for mouse
```

## Parameters (all required, no defaults)

Each row is a grill question. Ask the user; never fill one in.

| Argument | Meaning | Grill question | Demo-only value |
|---|---|---|---|
| `--adata` | input h5ad | Which h5ad holds the log-normalised cells? | `examples/demo.h5ad` |
| `--context` | species, tissue or cell type, condition, dataset | Species, tissue or cell type, condition and dataset? | `synthetic demo: 3 cell types x 40 cells, 2 samples, planted LIG1-REC1 and LIG2-REC2, no real tissue` |
| `--gene-id-col` | `var` column of Ensembl gene ids, or `none` | Which `var` column holds the Ensembl gene ids? (`none` if the source has symbols only; `gene_id` is then NA, never guessed.) | `gene_ids` |
| `--resource` | ligand-receptor TSV | Which resource file (consensus export, or a curated one)? | `examples/demo_resource.tsv` |
| `--groupby` | obs column of cell types | Which obs column holds the cell-type labels? | `cell_type` |
| `--sample` (optional) | obs column of samples | Run once per sample? If so, which column? Absent means one pooled run. | absent |
| `--methods` | `rank_aggregate` (all 5) or a comma list of `cellphonedb`, `connectome`, `logfc`, `natmi`, `singlecellsignalr` | Which methods go into the rank aggregate? | `rank_aggregate` |
| `--expr-prop` | min fraction of a group expressing each subunit | In what fraction of a cell type must the ligand and receptor be detected? | `0.1` |
| `--min-cells` | min cells per group | Below how many cells is a cell type dropped? | `5` |
| `--n-perms` | label permutations | How many permutations for the specificity p-values? | `100` |
| `--magnitude-max` | magnitude rank cutoff | Which magnitude rank (0-1, lower is stronger) calls an interaction? | `0.002` |
| `--specificity-max` | specificity rank cutoff | Which specificity rank (0-1, lower is more specific) calls an interaction? | `0.01` |
| `--sources` | sender labels, comma list or `all` | Which cell types send? | `all` |
| `--targets` | receiver labels, comma list or `all` | Which cell types receive? | `all` |
| `--n-top` | interactions in the dot plot | How many top passing interactions go in the dot plot? | `10` |
| `--reference` | genome build and annotation release of the gene names | Which genome build and annotation were the genes counted against (e.g. `GRCh38, GENCODE 44`)? | `synthetic demo genes, no genome build` |
| `--seed` | permutation seed | Which seed? | `0` |
| `--outdir` | new output directory | Which new directory gets the outputs? It must not exist. | `ccc_demo` |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them.
With `--demo`, only `--outdir` may also be given; any other argument is an error.

```bash
python scripts/ccc.py \
  --adata results/2026-10-27_annotate/data/annotated.h5ad \
  --resource resources/consensus_resource.tsv --groupby cell_type \
  --methods rank_aggregate --expr-prop 0.1 --min-cells 5 --n-perms 1000 \
  --magnitude-max 0.05 --specificity-max 0.05 --sources all --targets all \
  --n-top 20 --reference "GRCh38, GENCODE 44" --seed 0 \
  --outdir results/2026-10-27_ccc

# Demo (only when the user asks for one)
python scripts/ccc.py --demo --outdir demo_out
```

## Workflow (fixed order)

1. **Validate** arguments, `--methods` and input; bad input gives
   `ERROR: <message>`, exit 2. The method set needs at least one magnitude score
   (all but `logfc`) and one specificity score (all but `singlecellsignalr`).
2. **Check** gene names and cell ids are unique, `X` is log-normalised, the
   resource overlaps `var_names`, and every `--sources`/`--targets` label exists.
   `result.json` reports `n_matched` distinct ligand and receptor genes
   (complex subunits split on `_`) in `var_names` out of `n_in_resource`
   (demo: 40 of 40).
3. **Infer** with liana `AggregateClass(aggregate_meta, methods=...)`, the object
   behind `li.mt.rank_aggregate`: `resource=` the TSV, `expr_prop`, `min_cells`,
   `n_perms`, `seed`, `n_jobs=1`, `use_raw=False`, under one BLAS thread.
   Sources or targets other than `all` go in as `groupby_pairs`, which shrinks
   the set the ranks are computed over. `--sample` runs `.by_sample`.
4. **Call interactions**: `passes` when `magnitude_rank <= --magnitude-max` and
   `specificity_rank <= --specificity-max`. Zero passing stops the run (`no_pass`).
5. **Count** passing interactions per source x target (per sample with
   `--sample`); every pair of labels gets a row, zeros kept.
6. **Draw** the count heatmap and the dot plot of the first `--n-top` passing rows.
7. **Write** report, result and reproducibility bundle.

## Output structure

`interactions.tsv` holds every tested interaction, sorted by sample (if any),
then `magnitude_rank`, `specificity_rank`, source, target, ligand, receptor.
Then one row per resource pair with no tested interaction, with NA source,
target and ranks: `status` is `not_measured` when a subunit is missing from
`var_names`, `not_tested` when every subunit is measured but no source-target
row passed `--expr-prop` / `--min-cells`. Tested rows have `status` `tested`;
`status_counts` in `result.json` counts rows. `ligand_gene_id` and
`receptor_gene_id` map each subunit through `--gene-id-col` (subunits joined by
`_`; NA when a subunit has no id). Figures use tested rows only.
`ligand` and `receptor` are liana's `ligand_complex` and `receptor_complex`.
`pair_counts.tsv` is sorted by sample (if any), source, target in `--groupby`
order.

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # n_cells, n_groups, n_resource_pairs, n_matched, n_in_resource, n_interactions, n_pass, status_counts, gene_id_col, top_interactions, methods, seed, params, context, reference, chat_summary_lines (first: "Context: ..."), preferred_artifacts
├── tables/
│   ├── interactions.tsv       # [sample], source, target, ligand, receptor, ligand_gene_id, receptor_gene_id, magnitude_rank, specificity_rank, passes, status
│   └── pair_counts.tsv        # [sample], source, target, n_pass
├── figures/
│   ├── pair_counts_heatmap.pdf
│   ├── pair_counts_heatmap.png
│   ├── pair_counts_heatmap.csv
│   ├── lr_dotplot.pdf
│   ├── lr_dotplot.png
│   └── lr_dotplot.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python and package versions
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

## Gotchas

- **The model will want to call `li.mt.rank_aggregate(resource_name="consensus")`**
  and let liana pick the resource. Do not. Pass the exported TSV, so the resource
  is a recorded input and not a property of the installed liana version.
- **The model will want to relax `--magnitude-max` or `--specificity-max`**
  after a run with few or no passing interactions. Do not. Report what passed;
  the cutoffs are the user's answers, fixed before the run.
- **The model will want to read the ranks as p-values or probabilities** or say
  an interaction "is significant". Do not. A magnitude or specificity rank is a
  score aggregated over methods that orders interactions within this run; say
  "passes the rank cutoffs".
- **The model will want to call a pair with a gene outside `var_names` "not
  communicating".** Do not. Non-detection is not absence: the pair was not
  tested; its `status` in `interactions.tsv` is `not_measured`, and
  `n_matched` says how many resource genes were.
- **The model will want to say a sender drives or regulates the receiver** via a
  passing pair. Do not. A pair is predicted from co-expression against a prior
  resource; say "predicted interaction". Signalling needs protein or
  perturbation data (house-style `claims.strength`).
- **The model will want to run the human resource on mouse genes**, or to
  uppercase `var_names` to get past `no resource pair has its ligand and
  receptor in var_names`. Do not. Ask for the right species' resource.
- **The model will want to pass `--sources`/`--targets` to "just filter" the
  plot.** They change which interactions are ranked. Filter `interactions.tsv`
  afterwards when the user only wants a view.
- **The model will want to compare conditions by subtracting two runs' counts.**
  Do not. Per-sample runs give one table per sample; condition testing is a
  separate task.
- **The model will want to fill in a missing argument** with liana's defaults
  (`expr_prop=0.1`, `n_perms=1000`, `min_cells=5`). Do not. Ask, or run
  `--demo` only when the user asked for one.

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
