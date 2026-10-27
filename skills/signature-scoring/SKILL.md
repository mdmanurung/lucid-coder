---
name: signature-scoring
description: >-
  Gene-signature scores per cell from a log-normalised single-cell h5ad (scanpy
  score_genes or decoupler AUCell), the mean per sample and group in a long
  table ready for group-tests, a size-matched random-gene-set control per
  signature, box, UMAP and control figures, a report and a reproducibility
  bundle. Fire on "score this gene set", "signature score per cell",
  "score_genes", "module score", "AddModuleScore in scanpy", "AUCell score",
  "IFN signature per sample", "does my signature beat random genes".
  Do NOT fire for marker discovery (use sc-markers), for DE between conditions
  (use pseudobulk-de), for TF or pathway activity from a prior-knowledge network
  (ULM/MLM, PROGENy, CollecTRI), for UCell or NNLS reference scoring, or for raw
  counts that are not yet normalised.
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - single-cell
    - gene signature
    - score_genes
    - AUCell
    - decoupler
  inputs:
    - name: adata
      type: file
      format:
        - h5ad
      description: cell-level AnnData, log-normalised values in X, sample and group in obs
      required: true
    - name: gene_sets
      type: file
      format:
        - tsv
      description: long table with columns set, gene
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
      description: counts, genes found, random-control verdicts, params, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - anndata
      - scanpy
      - decoupler>=2
      - numpy
      - pandas
      - scipy
      - matplotlib
  demo_data:
    - path: examples/demo.h5ad
      description: synthetic, 2 groups x 4 samples x 50 cells, 300 genes, X_umap
    - path: examples/demo_gene_sets.tsv
      description: sig_up (planted, disease-up), sig_null (flat, one missing gene)
  endpoints:
    cli: python scripts/signature_scoring.py --adata {adata} --gene-sets {gene_sets} ... --outdir {outdir}
  trigger_keywords:
    - gene signature score
    - score_genes
    - module score
    - AUCell
    - signature per sample
---

# Gene-signature scoring

You are the dispatcher for **signature-scoring**. The script does the analysis;
you collect its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/signature_scoring.py`. Never re-implement the scoring, change a
value, or pick a value the user has not given. If a required argument is
unknown, ask the user (these are grill questions), or run `--demo` only when the
user asks for a demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it. To test whether a score differs
between groups, hand `tables/sample_scores.tsv` to the group-tests skill.

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--adata` | `.h5ad` | log-normalised values (`normalize_total` + `log1p`) in `X`; not raw counts, not scaled | sc-embed's `data/embedded.h5ad`; `examples/demo.h5ad` |
| `--sample-col` | obs column | one sample (donor, library) per cell; each sample in one group | `sample` |
| `--group-col` | obs column | 2 to 6 groups | `condition`; demo `group` |
| `--gene-sets` | TSV, long, header | columns `set` and `gene`, one row per gene; set names use only letters, digits, `_`, `.`, `-` | `examples/demo_gene_sets.tsv` |
| `obsm["X_umap"]` | optional | 2-D embedding; without it the UMAP figures are skipped | demo has one |

Gene symbols must match `var_names` exactly. Duplicate rows are dropped. The
script refuses a matrix with negative values (scaled) or only integers (raw
counts). A categorical `--group-col` keeps its category order; any other column
is cast to strings and sorted. An h5ad `.raw` is ignored (`use_raw=False`).
Gene names and cell ids must be unique (`skillkit.check_unique`).

## Parameters (all required, no defaults)

Each row is a grill question. Ask the user; never fill one in.

| Argument | Meaning | Grill question | Demo-only value |
|---|---|---|---|
| `--adata` | input h5ad | Which h5ad holds the log-normalised cells? | `examples/demo.h5ad` |
| `--gene-sets` | set, gene TSV | Which file lists the signature genes? | `examples/demo_gene_sets.tsv` |
| `--sample-col` | obs column of samples | Which obs column names the sample (the replicate)? | `sample` |
| `--group-col` | obs column of groups | Which obs column holds the groups to compare? | `group` |
| `--method` | `score_genes` or `aucell` | Mean-minus-control score (scanpy `score_genes`) or rank-based AUCell (decoupler)? | `score_genes` |
| `--ctrl-size` | score_genes only: control genes drawn per bin | How many control genes per expression bin? (corpus: 1000) | `50` |
| `--n-bins` | score_genes only: expression bins | How many expression bins for control genes? (corpus: 25) | `25` |
| `--n-up` | aucell only: top-ranked genes per cell in the AUC | How many top-ranked genes per cell count for the AUC? (decoupler's hidden default is 5% of genes) | not used |
| `--min-genes` | genes a set must have in the data | How many of a set's genes must be in the data to score it? | `5` |
| `--n-random` | random sets per signature | How many size-matched random gene sets per signature? | `50` |
| `--reference` | genome build and annotation release of the gene names | Which genome build and annotation were the genes counted against (e.g. `GRCh38, GENCODE 44`)? | `synthetic demo genes, no genome build` |
| `--seed` | seed | Which seed? Fixes score_genes control genes, random sets and jitter. | `0` |
| `--outdir` | new output directory | Which new directory gets the outputs? It must not exist. | given by the user |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them;
beside `--demo` only `--outdir` may be given. `--ctrl-size` and `--n-bins` are
required with `score_genes` and refused with `aucell`; `--n-up` the reverse.

```bash
python scripts/signature_scoring.py \
  --adata results/2026-10-27_embed/data/embedded.h5ad \
  --gene-sets data/hallmark_ifn.tsv --sample-col donor --group-col condition \
  --method score_genes --ctrl-size 1000 --n-bins 25 --min-genes 5 \
  --n-random 100 --reference "GRCh38, GENCODE 44" --seed 0 \
  --outdir results/2026-10-27_ifn-score

python scripts/signature_scoring.py --adata embedded.h5ad --gene-sets sets.tsv \
  --sample-col donor --group-col condition --method aucell --n-up 500 \
  --min-genes 5 --n-random 100 --reference "GRCh38, GENCODE 44" --seed 0 \
  --outdir results/ifn-aucell

# Demo (only when the user asks for one)
python scripts/signature_scoring.py --demo --outdir demo_out
```

## Workflow (fixed order)

1. **Validate** arguments; bad input gives `ERROR: <message>`, exit 2. Gene
   names and cell ids must be unique.
2. **Check** the matrix is log-normalised, each sample sits in one group, and
   there are 2 to 6 groups. Match each set's genes to `var_names`; record
   missing genes per set; stop if any set has fewer than `--min-genes` found.
   `result.json` reports `n_matched` distinct set genes found out of
   `n_in_resource` distinct genes listed (demo: 20 of 21).
3. **Score** each set per cell on found genes only.
   `score_genes`: `sc.tl.score_genes(ctrl_size=, n_bins=, random_state=seed,
   use_raw=False)`, one call per set. `aucell`: one `dc.mt.aucell(net,
   tmin=--min-genes, n_up=, empty=False)` call for all sets.
4. **Summarise** per sample: the mean score over the sample's cells, with its
   group and cell count.
5. **Control**: for each set draw `--n-random` random gene sets of the same
   size from the genes not in the set (`numpy.random.default_rng(seed)`), score
   them with the same engine and arguments, and summarise them the same way.
   The statistic is the spread of group means of sample means (max minus min
   group mean). A set **beats random** when its spread exceeds every random
   set's; `p_empirical = (1 + n random >= observed) / (n_random + 1)`.
6. **Write** `data/scored.h5ad` (obs `<method>_<set>`), the tables, figures,
   report, result and reproducibility bundle.

## Output structure

`cell_scores.tsv` is sorted by cell; `sample_scores.tsv` by feature, then
subject; `genes_found.tsv` and `random_control.tsv` by set. `sample_scores.tsv`
uses group-tests' input names: `feature` is the set, `subject` the sample, and
the group column keeps its obs name. Figures are named per set; the tree lists
the demo's sets `sig_null` and `sig_up`.

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # n_cells, n_samples, n_sets, n_matched, n_in_resource, n_found, missing_genes, beats_random, p_empirical, method, seed, params, reference, chat_summary_lines, preferred_artifacts
├── data/
│   └── scored.h5ad            # input plus obs <method>_<set>; not checksummed
├── tables/
│   ├── cell_scores.tsv        # cell, sample col, group col, one column per set
│   ├── sample_scores.tsv      # feature, subject, group col, value, n_cells
│   ├── genes_found.tsv        # set, n_listed, n_found, missing
│   └── random_control.tsv     # set, n_genes, n_random, observed_spread, random_mean, random_max, p_empirical, beats_random
├── figures/
│   ├── score_by_group_sig_null.pdf
│   ├── score_by_group_sig_null.png
│   ├── score_by_group_sig_null.csv
│   ├── score_by_group_sig_up.pdf
│   ├── score_by_group_sig_up.png
│   ├── score_by_group_sig_up.csv
│   ├── score_umap_sig_null.pdf
│   ├── score_umap_sig_null.png
│   ├── score_umap_sig_null.csv
│   ├── score_umap_sig_up.pdf
│   ├── score_umap_sig_up.png
│   ├── score_umap_sig_up.csv
│   ├── random_control.pdf
│   ├── random_control.png
│   └── random_control.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python and package versions
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

## Gotchas

- **The model will want to test groups on the per-cell scores.** Do not. Cells
  are not replicates; the test unit is the sample. Hand `sample_scores.tsv` to
  group-tests.
- **The model will want to call a score a probability, or say every gene in the
  set moved.** Do not. A score_genes value is a log-expression difference and
  an AUCell value is a rank-based AUC; neither is the chance a cell "has" the
  programme. A set that separates groups is not a change in every member gene.
- **The model will want to call a gene in `missing` "not expressed".** Do not.
  Non-detection is not absence: a gene not in `var_names` was filtered or named
  differently, so it is "not tested".
- **The model will want to treat each sample as an independent unit** when one
  donor gave several samples. Do not. The sample means and `p_empirical` assume
  independent samples; the donor term belongs in group-tests' design.
- **The model will want to read `beats_random` as a group test.** Do not. It
  says the signature separates groups more than random genes of the same size;
  direction and significance between groups come from group-tests.
- **The model will want to fill in `--ctrl-size 1000 --n-bins 25`** because the
  corpus uses them, or let AUCell pick its 5% `n_up`. Do not. Ask; the corpus
  values go in the grill question, not the command.
- **The model will want to lower `--min-genes`** or drop the set when a set
  fails it. Do not. Report the error and the missing genes; the user fixes the
  symbols (case, aliases, Ensembl IDs) or chooses the threshold.
- **The model will want to compare score_genes and AUCell values** across runs
  or sets as if on one scale. Do not. score_genes is a log-expression
  difference, AUCell an AUC in [0, 1]; compare within one method and run.
- **The model will want to score on raw counts or scaled data** to get past the
  `does not look log-normalised` error. Do not. Ask where the log-normalised
  values are.

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
