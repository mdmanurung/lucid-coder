---
name: tf-activity
description: >-
  TF or pathway activity from a prior-knowledge network (CollecTRI, PROGENy or
  any source-target-weight TSV) with decoupler 2.x ULM or MLM. Two inputs: a
  log-normalised single-cell h5ad (score each cell, then each group against the
  rest) or a DE table (score one contrast's gene statistic). Writes an activity
  table, the sources dropped under tmin, a heatmap and violin plot (cells) or a
  bar plot (contrast), a report and a reproducibility bundle. Fire on "TF
  activity", "transcription factor activity", "CollecTRI", "PROGENy", "pathway
  activity", "decoupler ulm", "decoupler mlm", "which TFs drive this contrast".
  Do NOT fire for gene-set over-representation or GSEA (fgsea, ORA), for ATAC
  motif deviation or motif enrichment (chromVAR, ArchR, MOCHA), for building or
  downloading the network itself, or for raw counts that are not normalised.
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - transcription factor
    - pathway activity
    - decoupler
    - CollecTRI
    - PROGENy
  inputs:
    - name: context
      type: string
      description: species, tissue or cell type, condition and dataset
      required: true
    - name: network
      type: file
      format:
        - tsv
      description: source, target, weight edges, exported once with decoupler
      required: true
    - name: adata
      type: file
      format:
        - h5ad
      description: cell-level AnnData, log-normalised X (--input cells)
      required: false
    - name: de_table
      type: file
      format:
        - tsv
      description: one contrast, gene plus a statistic column (--input contrast)
      required: false
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
      description: counts, top sources, params, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - decoupler>=2.0
      - anndata
      - numpy
      - pandas
      - scipy
      - matplotlib
      - plotnine-extra
  demo_data:
    - path: examples/demo.h5ad
      description: synthetic, 3 groups x 60 cells, 300 genes; TF1 active in g1, TF2 in g2
    - path: examples/demo_network.tsv
      description: 8 sources; TF7 and TF8 fall under tmin 5
    - path: examples/demo_de.tsv
      description: Welch t per gene, g1 against g2
  endpoints:
    cli: python scripts/tf_activity.py --input {cells|contrast} --network {network} ... --outdir {outdir}
  trigger_keywords:
    - TF activity
    - CollecTRI
    - PROGENy
    - pathway activity
    - decoupler
---

# TF and pathway activity

You are the dispatcher for **tf-activity**. The script does the analysis; you
collect its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/tf_activity.py`. Never re-implement the analysis, change a
threshold, or pick a value the user has not given. If a required argument is
unknown, ask the user (these are grill questions), or run `--demo` only when the
user asks for a demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--network` | TSV | columns `source`, `target`, numeric `weight`; no duplicate source-target pair; other columns ignored | `collectri_human.tsv`; `examples/demo_network.tsv` |
| `--adata` (cells) | `.h5ad` | log-normalised `X` (`normalize_total` + `log1p`); not raw counts, not scaled | sc-embed's `data/embedded.h5ad`; `examples/demo.h5ad` |
| `--groupby` (cells) | obs column | at least 2 groups, each with at least 2 cells | `leiden`, `cell_type`; demo `group` |
| `--de-table` (contrast) | TSV | column `gene`, unique, and the `--stat-col` column; one contrast per file | pseudobulk-de's `tables/de_results.tsv`; `examples/demo_de.tsv` |

Gene names in the data must match the network's `target` names (symbols for
CollecTRI and PROGENy). Rows of the DE table with a missing or non-finite stat
are dropped and counted in the summary. Gene names, cell ids and DE-table genes
must be unique (`skillkit.check_unique`).

The skill does not download. Export each network once with decoupler and pass
the file (needs network access; CollecTRI's licence is academic by default):

```python
import decoupler as dc

dc.op.collectri(organism="human", remove_complexes=False, license="academic").to_csv(
    "collectri_human.tsv", sep="\t", index=False
)
dc.op.progeny(organism="human", top=500, thr_padj=0.05, license="academic").to_csv(
    "progeny_human_top500.tsv", sep="\t", index=False
)
```

## Parameters (all required, no defaults)

Ask these through the analysis contract: `${CLAUDE_PLUGIN_ROOT}/skills/lucid-coder/references/contracts.md`.

Each row is a grill question. Ask the user; never fill one in.

| Argument | Meaning | Grill question | Demo-only value |
|---|---|---|---|
| `--input` | `cells` or `contrast` | Score each cell of an h5ad, or one DE contrast? | `cells` |
| `--network` | network TSV | Which exported network file (CollecTRI for TFs, PROGENy for pathways)? | `examples/demo_network.tsv` |
| `--method` | `ulm` or `mlm` | Univariate (ULM, one source at a time) or multivariate (MLM, sources fitted jointly)? | `ulm` |
| `--tmin` | minimum targets per source | How many of a source's targets must be in the data for it to be scored? | `5` |
| `--n-top` | top sources in the figures | How many top sources per group (cells) or in total (contrast) go in the figures? | `2` |
| `--reference` | genome build and annotation release of the gene names | Which genome build and annotation were the genes counted against (e.g. `GRCh38, GENCODE 44`)? | `synthetic demo genes, no genome build` |
| `--context` | biological context of the data | Species, tissue or cell type, condition and dataset? | `synthetic demo: 3 groups x 60 cells, no real tissue` |
| `--seed` | seed | Which seed? (Recorded; no step draws random numbers.) | `0` |
| `--adata` (cells) | input h5ad | Which h5ad holds the log-normalised cells? | `examples/demo.h5ad` |
| `--groupby` (cells) | obs column of groups | Which obs column holds the clusters or labels? | `group` |
| `--rank-method` (cells) | `t-test_overestim_var`, `t-test` or `wilcoxon` | Which test compares each group's cell scores with the rest? | `t-test_overestim_var` |
| `--limit` (cells) | heatmap colour saturation | At which absolute mean activity should the heatmap colours saturate? | `3` |
| `--de-table` (contrast) | DE table | Which DE table holds the contrast? | `examples/demo_de.tsv` |
| `--stat-col` (contrast) | statistic column | Which column to score: the Wald or t `stat`, or `log2_fc`? | `stat` |
| `--contrast` (contrast) | contrast label | What label names this contrast in the tables? | `g1_vs_g2` |
| `--outdir` | new output directory | Which new directory gets the outputs? It must not exist. | `tf_activity_demo` |

The last column is **demo-only, not defaults**. `--demo` alone pins the cells
values; `--demo --input contrast` pins the contrast values. With `--demo`, only
`--outdir` and `--input` may also be given. Arguments marked (cells) or
(contrast) are required in that mode and refused in the other.

```bash
python scripts/tf_activity.py --input cells \
  --adata results/2026-10-27_embed/data/embedded.h5ad --groupby leiden \
  --network ref/collectri_human.tsv --method ulm --tmin 5 \
  --rank-method t-test_overestim_var --n-top 5 --limit 3 \
  --reference "GRCh38, GENCODE 44" --seed 0 \
  --context "human PBMC, IFN-beta vs control, Kang 2018" \
  --outdir results/2026-10-27_tf_cells

python scripts/tf_activity.py --input contrast \
  --de-table results/2026-10-27_pbde/tables/de_results.tsv --stat-col stat \
  --contrast disease_vs_healthy --network ref/progeny_human_top500.tsv \
  --method mlm --tmin 5 --n-top 10 --reference "GRCh38, GENCODE 44" --seed 0 \
  --context "human PBMC, IFN-beta vs control, Kang 2018" \
  --outdir results/2026-10-27_tf_de

# Demo (only when the user asks for one)
python scripts/tf_activity.py --demo --outdir demo_out
python scripts/tf_activity.py --demo --input contrast --outdir demo_contrast
```

## Workflow (fixed order)

1. **Validate** arguments, the mode's arguments and the network; bad input gives
   `ERROR: <message>`, exit 2.
2. **Load** the cells (log-normalised check, at least 2 groups of 2 cells) or the
   DE table (unique genes, finite stat).
3. **Prune** the network: count each source's targets present in the data
   (`n_targets`; in cells mode a gene counts only when some cell has a nonzero
   value, because decoupler drops all-zero genes); sources under `--tmin` go to
   `sources_dropped.tsv`. None left stops the run (`no_sources`). `result.json`
   reports `n_matched` distinct network targets in the data out of
   `n_in_resource` distinct targets (demo: 93 of 99).
4. **Score** with `dc.mt.ulm` or `dc.mt.mlm` (`tmin=--tmin`).
   - cells: one score per cell and source, then `dc.tl.rankby_group(...,
     reference="rest", method=--rank-method)` per group. `score` is the group
     mean; `stat`, `pvalue`, `padj` are the group-vs-rest test, BH within group.
   - contrast: one score per source for the `--stat-col` vector
     (`empty=False`). BH across sources. decoupler returns BH p for ULM and raw p
     for MLM; the script rebuilds the other one (ULM raw p from the t value with
     genes - 2 df), and a test checks `BH(pvalue) == padj` for both.
5. **Pick top sources**: cells, the first `--n-top` per group by `stat`
   (most active against the rest); contrast, the `--n-top` by absolute score.
6. **Draw** the heatmap and violin plot (cells) or the bar plot (contrast).
7. **Write** report, result and reproducibility bundle.

## Output structure

`activity.tsv` columns: cells `source, group, score, stat, pvalue, padj,
status, n_targets`, sorted by group (category order), then `stat` descending,
then source; contrast `source, contrast, score, pvalue, padj, status,
n_targets`, sorted by score descending, then source. Every network source has a
row: scored sources have `status` `tested`; sources under `--tmin` follow as one
`not_tested` row each, sorted by source, with NA group, statistics and p-values
(never 0). Figures use only `tested` rows. `result.json` `status_counts` counts
sources per status. `sources_dropped.tsv` is `source, n_targets`, sorted by
source.

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # input, n_cells or n_genes, n_sources_scored, n_sources_dropped, status_counts, n_matched, n_in_resource, top_sources, method, seed, params, context, reference, chat_summary_lines, preferred_artifacts
├── tables/
│   ├── activity.tsv
│   └── sources_dropped.tsv
├── figures/
│   ├── activity_heatmap.pdf   # (cells)
│   ├── activity_heatmap.png   # (cells)
│   ├── activity_heatmap.csv   # (cells)
│   ├── activity_violin.pdf    # (cells)
│   ├── activity_violin.png    # (cells)
│   ├── activity_violin.csv    # (cells)
│   ├── activity_bar.pdf       # (contrast)
│   ├── activity_bar.png       # (contrast)
│   └── activity_bar.csv       # (contrast)
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python and package versions
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

## Gotchas

- **The model will want to download CollecTRI or PROGENy inside the run.** Do
  not. The run is offline; export the network once with the calls above and
  pass the file, so the run records which network it used.
- **The model will want to lower `--tmin`** when a favourite TF lands in
  `sources_dropped.tsv`. Do not. A source with few targets in the data gives an
  unstable score; report it as dropped, and ask the user before any rerun.
- **The model will want to read per-cell p-values as group significance.** Do
  not. Cells are the unit only when the condition was assigned to each cell (a
  per-cell perturbation); otherwise the donor is the unit, and the cells-mode
  p-values are descriptive, ranking sources within one dataset. Condition
  effects go through pseudobulk-de, then `--input contrast`.
- **The model will want to read a ULM or MLM score as a probability** that the
  TF is active. Do not. It is a t-value of the target regression, signed by the
  network weights; an active source is not a change in every target.
- **The model will want to call a target outside the data "not regulated".** Do
  not. Non-detection is not absence: the gene was filtered or never measured, so
  it is "not tested"; `n_matched` says how many targets were.
- **The model will want to score a DE table's `padj` or `-log10(p)` column.**
  Do not. Score a signed statistic (`stat` or `log2_fc`), as the user chose;
  an unsigned column cannot say activated or repressed.
- **The model will want to say a top TF "drives", "causes" or "regulates" the
  contrast or cluster.** Do not. TF activity is inferred from its targets'
  expression and a prior network, not measured; say the TF's activity is
  "predicted higher" or "scored higher" and "associated with" the group
  (house-style `claims.strength`). A causal claim needs a perturbation.
- **The model will want to call a positive score "TF expression is up".** Do
  not. The score is from the targets' values and the network's signed weights;
  the TF's own transcript is not used.
- **The model will want to feed scaled or raw-count data** past the
  `does not look log-normalised` error. Do not. Ask the user where the
  log-normalised matrix is.
- **The model will want to fill in a missing argument** with a "standard" value
  (ulm, tmin 5, the corpus `t-test_overestim_var`). Do not. Ask, or run `--demo`
  only when the user asked for one.

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
