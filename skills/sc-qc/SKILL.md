---
name: sc-qc
description: >-
  scRNA-seq cell QC and doublet removal from a raw-count h5ad with a sample
  column: computes per-cell QC metrics (genes, counts, % mito, % ribo), scores
  doublets with Scrublet per sample, applies fixed or per-sample MAD thresholds,
  drops failed cells and predicted doublets, and writes a filtered h5ad, tables,
  four fixed figures, a report and a reproducibility bundle. Fire on "QC my
  single-cell data", "filter low-quality cells", "remove doublets", "run
  Scrublet", "mito / gene count thresholds", "MAD-based QC per sample".
  Do NOT fire for ambient RNA removal (CellBender, SoupX), normalisation,
  clustering or annotation, gene filtering, or QC of bulk RNA-seq samples.
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - single-cell
    - quality control
    - doublets
    - Scrublet
    - scanpy
  inputs:
    - name: context
      type: string
      description: species, tissue or cell type, condition and dataset
      required: true
    - name: adata
      type: file
      format:
        - h5ad
      description: cell-level AnnData with raw integer counts in X or a layer
      required: true
  outputs:
    - name: filtered
      type: file
      format:
        - h5ad
      description: kept cells, raw counts unchanged, QC columns in obs
    - name: report
      type: file
      format:
        - md
      description: parameters, summary and captioned figures
    - name: result
      type: file
      format:
        - json
      description: counts before and after, thresholds, params, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - anndata
      - scanpy==1.12.4
      - scikit-image
      - numpy
      - pandas
      - scipy
      - matplotlib
      - plotnine-extra
  demo_data:
    - path: examples/demo.h5ad
      description: synthetic, 3 samples x 345 cells, planted low-quality cells and doublets
  endpoints:
    cli: python scripts/sc_qc.py --adata {adata} ... --outdir {outdir}
  trigger_keywords:
    - single-cell QC
    - cell filtering
    - doublet removal
    - Scrublet
    - mitochondrial threshold
    - MAD outliers
---

# scRNA QC and doublet removal

You are the dispatcher for **sc-qc**. The script does the analysis; you collect
its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/sc_qc.py`. Never re-implement the analysis, change a threshold, or
pick a value the user has not given. If a required argument is unknown, ask the
user (these are grill questions), or run `--demo` only when the user asks for a
demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--adata` | `.h5ad` | raw counts (non-negative integers) in `X`, or in `--layer` | `examples/demo.h5ad` |
| cell-level `obs` | column | `--sample-col`, one value per library (Scrublet runs per value) | `sample` |
| `var_names` | gene symbols | mitochondrial genes start with `--mito-prefix` | `MT-CO1` |

Each sample needs more than 30 cells (one per Scrublet principal component).

## Parameters (all required, no defaults)

Each row is a grill question. Ask the user; never fill one in.

| Argument | Grill question | Demo-only value |
|---|---|---|
| `--adata` | Which h5ad holds the unfiltered cells? | `examples/demo.h5ad` |
| `--layer` (optional) | Which layer holds the raw counts? Absent means `X`. | absent (`X`) |
| `--sample-col` | Which obs column names the sample (one 10x library)? | `sample` |
| `--mito-prefix` | Which prefix marks mitochondrial genes (`MT-` human, `mt-` mouse)? | `MT-` |
| `--mode` | Fixed thresholds for every sample (`fixed`), or outliers per sample by MAD (`mad`)? | `mad` |
| `--min-genes` (fixed) | Fewest genes detected for a kept cell? | — |
| `--max-genes` (fixed) | Most genes detected for a kept cell? | — |
| `--min-counts` (fixed) | Fewest UMI counts for a kept cell? | — |
| `--max-pct-mito` (fixed) | Highest % mitochondrial counts for a kept cell? | — |
| `--mads-counts` (mad) | How many MADs from the median of log1p counts, each side? | `5` |
| `--mads-genes` (mad) | How many MADs from the median of log1p genes, each side? | `5` |
| `--mads-mito` (mad) | How many MADs above the median % mito? | `3` |
| `--expected-doublet-rate` | Which doublet rate does the 10x loading table give for these libraries' recovered cells? | `0.08` |
| `--doublet-threshold` | Let Scrublet pick a cutoff per sample (`auto`), or one fixed score in (0, 1)? | `auto` |
| `--seed` | Which seed for Scrublet? | `0` |
| `--reference` | Which genome build and annotation release named the genes (e.g. `GRCh38, GENCODE 44`)? Recorded in `result.json`. | `synthetic demo genes, no genome build` |
| `--context` | Species, tissue or cell type, condition and dataset? | `synthetic demo: 3 samples of 5 populations, planted low-quality cells and doublets, no real tissue` |
| `--outdir` | Which new directory gets the outputs? It must not exist. | `sc_qc_demo` |

The four fixed-mode flags are required with `--mode fixed` and refused with
`--mode mad`; the three MAD flags the other way round. The last column is
**demo-only, not defaults**. `--demo` alone pins all of them; with `--demo`, only
`--outdir` may also be given.

## Workflow (fixed order)

1. **Validate** arguments and input; bad input gives `ERROR: <message>`, exit 2.
   Gene names and cell ids must be unique (`skillkit.check_unique`).
2. **QC metrics**: `sc.pp.calculate_qc_metrics` on the raw counts, with
   `mito` = genes starting with `--mito-prefix` and `ribo` = genes starting
   with `RPS` or `RPL` (any case). Ribo is reported, never thresholded.
3. **Scrublet** per sample: `sc.pp.scrublet(batch_key=<sample>,
   expected_doublet_rate=..., random_state=--seed, n_prin_comps=30,
   use_approx_neighbors=False)` on all cells. A cell is a predicted doublet when
   its score exceeds its sample's threshold (`--doublet-threshold auto` lets
   Scrublet pick one per sample).
4. **Thresholds**, per sample. A cell passes when it sits inside every limit.
   - `--mode fixed`: `--min-genes <= genes <= --max-genes`,
     `counts >= --min-counts`, `% mito <= --max-pct-mito`; the same in every sample.
   - `--mode mad`: per sample, median ± n × MAD on `log1p` counts
     (`--mads-counts`) and `log1p` genes (`--mads-genes`), both sides; median +
     n × MAD on % mito (`--mads-mito`), upper side only. MAD is unscaled
     (`scipy.stats.median_abs_deviation`, scale 1).
   - `qc_status` names the first failing rule in the order `low_counts`,
     `high_counts`, `low_genes`, `high_genes`, `high_mito`; else `pass`.
5. **Drop** cells that fail QC or are predicted doublets; write the filtered h5ad.

## Output structure

`qc_per_cell.tsv` is sorted by sample, then cell; `cells_per_sample.tsv` by sample.

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # n_cells_before, n_cells_after, n_failed_qc, n_doublets, cells_after_per_sample, failed_qc_per_sample, doublets_per_sample, doublet_thresholds, mode, seed, params, reference, context, chat_summary_lines (first: "Context: ..."), preferred_artifacts
├── data/
│   └── qc_filtered.h5ad       # kept cells; X and layers as input; obs gains the QC columns
├── tables/
│   ├── qc_per_cell.tsv        # cell, sample, total_counts, n_genes_by_counts, pct_counts_mito, pct_counts_ribo, doublet_score, predicted_doublet, qc_status, fate, kept
│   └── cells_per_sample.tsv   # sample, n_cells_before, n_failed_qc, n_doublets, n_cells_after, min/max counts, min/max genes, max_pct_mito, doublet_threshold
├── figures/
│   ├── qc_violins.pdf
│   ├── qc_violins.png
│   ├── qc_violins.csv
│   ├── counts_vs_genes.pdf
│   ├── counts_vs_genes.png
│   ├── counts_vs_genes.csv
│   ├── doublet_scores.pdf
│   ├── doublet_scores.png
│   ├── doublet_scores.csv
│   ├── cells_per_sample.pdf
│   ├── cells_per_sample.png
│   └── cells_per_sample.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python and package versions
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

`n_failed_qc` and `n_doublets` overlap: a doublet that also fails QC counts in
both. `fate` is exclusive: `kept`, `failed QC`, or `doublet` (passed QC only).
Limits in `cells_per_sample.tsv` are on the raw scale; `max_counts` is `inf` in
fixed mode.

Drops are reported per sample, the group QC runs in. sc-qc takes no condition
column on purpose: QC often runs before the sample sheet is joined, and samples
nest in conditions. The skill that aggregates samples next (pseudobulk-de,
composition) sums `failed_qc_per_sample` and `doublets_per_sample` per condition.

## Gotchas

- **The model will want to run Scrublet on the pooled object** to save time.
  Do not. Doublets form within one library; pooled scores mix samples with
  different depths and doublet rates. The script runs it per sample.
- **The model will want one global MAD threshold** across samples. Do not. A
  shallow sample then loses most of its cells; `--mode mad` is per sample by design.
- **The model will want to tune thresholds after seeing the violins** or the
  cell counts. Do not. Thresholds are the user's answers, fixed before the run;
  report the counts as they are and ask before a rerun.
- **The model will want to guess `--mito-prefix`** or fall back to `MT-` for
  mouse. Do not. Ask; the `no gene name starts with` error exists because a wrong
  prefix silently gives 0% mito.
- **The model will want to drop a whole sample** that loses many cells, or to
  set `--doublet-threshold` by eye. Do not. Report it and ask the user.
- **The model will want to read the doublet score as a probability** ("82%
  likely a doublet"). Do not. A Scrublet score is a kNN fraction of simulated
  doublets; only the per-sample threshold turns it into a call.
- **The model will want to call a `predicted_doublet` cell a doublet** ("83
  doublets in S1"). Do not. It is a Scrublet call on a score, so write
  "predicted doublets" (house-style `claims.strength`).
- **The model will want to guess `--reference`** from the gene names. Do not.
  Ask for the build and annotation release; it is recorded, never inferred.
- **The model will want to feed normalised or log values** to clear a
  `counts ... must be non-negative integers` error. Do not. Point `--layer` at the
  raw counts or ask where they are.

## CLI reference

```bash
python scripts/sc_qc.py \
  --adata data/raw_cells.h5ad --layer counts --sample-col sample \
  --mito-prefix MT- --mode mad --mads-counts 5 --mads-genes 5 --mads-mito 3 \
  --expected-doublet-rate 0.06 --doublet-threshold auto \
  --seed 0 --reference 'GRCh38, GENCODE 44' \
  --context 'human PBMC, IFN-beta vs control, Kang 2018' --outdir results/2026-10-27_qc

python scripts/sc_qc.py \
  --adata data/raw_cells.h5ad --sample-col sample --mito-prefix MT- \
  --mode fixed --min-genes 200 --max-genes 5000 --min-counts 500 --max-pct-mito 10 \
  --expected-doublet-rate 0.06 --doublet-threshold auto \
  --seed 0 --reference 'GRCh38, GENCODE 44' \
  --context 'human PBMC, IFN-beta vs control, Kang 2018' --outdir results/2026-10-27_qc-fixed

# Demo (only when the user asks for one)
python scripts/sc_qc.py --demo --outdir demo_out
```

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
