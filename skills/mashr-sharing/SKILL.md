---
name: mashr-sharing
description: >-
  Cross-cell-type effect sharing with mashr: takes one DE table per cell type
  (gene_id, gene_symbol, logFC, SE), checks the ids are compatible, aligns genes,
  fits mashr with canonical or canonical plus data-driven (ED) covariances, and
  writes posterior means, posterior SDs and lfsr per gene and cell type, counts
  of genes significant in 1..K cell types, mashr pairwise sharing, raw-logFC
  Spearman per pair, three figures, a report and a reproducibility bundle.
  Fire on "mashr across cell types", "run_mash", "which DE genes are shared
  between cell types", "cell-type-specific vs shared effects", "effect sharing
  heatmap", "replicate logFC across cell types", "borrow strength across cell
  types". Do NOT fire for running the per-cell-type DE itself (pseudobulk-de),
  meta-analysis of the same cell type across cohorts, gene-set enrichment on
  the posterior (zenith/enrichment skill), or cell-type composition questions.
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - mashr
    - multivariate adaptive shrinkage
    - effect sharing
    - differential expression
    - cell types
  inputs:
    - name: context
      type: string
      description: species, tissue or cell type, condition and dataset
      required: true
    - name: de
      type: file
      format:
        - tsv
      description: one DE table per cell type, given as --de name=path (repeatable)
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
      description: shared and specific counts, params, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - numpy
      - pandas
      - scipy
      - matplotlib
      - plotnine-extra
    r:
      - mashr
      - ashr
  demo_data:
    - path: examples/demo_B.tsv
      description: synthetic, 4 cell types x 800 genes, 80 shared, 20 specific per type
  endpoints:
    cli: python scripts/mashr_sharing.py --de {name}={path} ... --outdir {outdir}
  trigger_keywords:
    - mashr
    - run_mash
    - effect sharing
    - shared DE genes
    - cell-type-specific DE
    - lfsr
---

# mashr sharing

You are the dispatcher for **mashr-sharing**. The script does the analysis; you
collect its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/mashr_sharing.py`. Never re-implement the analysis, change a
threshold, or pick a value the user has not given. If a required argument is
unknown, ask the user (these are grill questions), or run `--demo` only when the
user asks for a demo.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--de name=path` | TSV, one per cell type, at least 2 | one row per gene; `name` labels the cell type | `CD4T=out_cd4/tables/de_results.tsv` |
| gene columns | `gene_id`, `gene_symbol` (fixed names) | `gene_id` unique within a table (Ensembl, as the source gives it), or NA in every row when the source has no ids | `ENSG00000111640.15`, `GAPDH` |
| effect column | `--effect-col` | numeric, no NA (logFC, same contrast in every table) | `log2_fc` |
| SE column | `--se-col` | numeric, no NA, > 0 | `lfc_se` |

A pseudobulk-de `tables/de_results.tsv` carries `gene_id`, `gene_symbol`,
`log2_fc`, `lfc_se`. Drop its NA rows first (genes the engine could not fit). A limma or
dreamlet table has no SE column: add `se = logFC / t`, as dreamlet's `run_mash`
does.

Required checks (docs/pipeline-skills.md, adopted 2026-10-27):
- **Unique ids**: a repeated gene in a table stops the run (`duplicate_gene`), as
  does a repeated cell-type name (`duplicate_name`).
- **Gene ids** (adopted 2026-10-28): genes are joined on `gene_id`. A table with
  `gene_id` NA in some rows stops (`partial_ids`); NA in every row (a source with
  symbols only) joins on `gene_symbol`, writes `gene_id` as NA and records
  `joined_on: gene_symbol`. A symbol with several ids stops the run
  (`skillkit.check_gene_ids`).
- **Compatible inputs** (adopted 2026-10-28): `skillkit.check_compatible` stops
  when the tables use different id types (versioned vs unversioned Ensembl, ids vs
  symbols) or share no gene; `result.json` reports `n_shared_ids`.
- **Reference record**: `--reference` (genome build and annotation release) is
  required and written to `result.json`; every table must use the same one.
- **Drops per group**: genes each cell type loses to `--align intersect` are
  reported per cell type (`genes_dropped_by_align`; all 0 for `union`).
- **No design**: mashr fits no design matrix here, so there is no rank or
  donor-versus-sample check; replication belongs to the upstream DE tables.

## Parameters (all required, no defaults)

Ask these through the analysis contract: `${CLAUDE_PLUGIN_ROOT}/skills/lucid-coder/references/contracts.md`.

Each row is a grill question. Ask the user; never fill one in.

| Argument | Grill question | Demo-only value |
|---|---|---|
| `--context` | Species, tissue or cell type, condition and dataset? It opens the report and the chat summary. | `synthetic demo: 4 cell types x 800 genes, no real tissue` |
| `--de` (repeat) | Which DE table belongs to which cell type? Do all tables test the same contrast? | `B=`, `CD4T=`, `CD8T=`, `Mono=` `examples/demo_*.tsv` |
| `--effect-col` | Which column holds the effect (logFC)? | `log2_fc` |
| `--se-col` | Which column holds its standard error? | `lfc_se` |
| `--align` | Keep only genes tested in every cell type (`intersect`), or every gene with the missing cells treated as unmeasured (`union`)? | `intersect` |
| `--covariances` | Canonical covariances only (`canonical`, as dreamlet's `run_mash`), or canonical plus data-driven ED (`canonical+ed`)? | `canonical+ed` |
| `--null-corr` | Do these DE tables come from the same donors? Then `em` (null correlation `V` by EM, as dreamlet's `run_mash`); `identity` only for independent cells and samples. | `identity` |
| `--strong-lfsr` | Which lfsr in a one-by-one ash fit picks the strong signals that train the ED covariances? Used only with `canonical+ed`. | `0.05` |
| `--lfsr-max` | Which lfsr marks a gene significant in a cell type (also the sharing threshold)? | `0.05` |
| `--sharing-factor` | Within what factor must two posterior means agree to count as shared (mashr's `factor`)? | `0.5` |
| `--reference` | Which genome build and annotation release are the genes from (e.g. `GRCh38, GENCODE 44`)? | `synthetic demo genes, no genome build` |
| `--seed` | Which seed? (Passed to `set.seed` and `mash(seed=)`.) | `0` |
| `--outdir` | Which new directory gets the outputs? It must not exist. | `mashr_sharing_demo` |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them;
only `--outdir` may be given beside it.

```bash
python scripts/mashr_sharing.py \
  --context "human PBMC, IFN-beta vs control, Kang 2018" \
  --de B=b/tables/de_results.tsv --de CD4T=cd4/tables/de_results.tsv \
  --de Mono=mono/tables/de_results.tsv \
  --effect-col log2_fc --se-col lfc_se --align intersect \
  --covariances canonical+ed --strong-lfsr 0.05 --lfsr-max 0.05 \
  --sharing-factor 0.5 --reference "GRCh38, GENCODE 44" --seed 0 \
  --outdir results/2026-10-27_mash

# Demo (only when the user asks for one)
python scripts/mashr_sharing.py --demo --outdir demo_out
```

## Workflow (fixed order)

1. **Validate** arguments and every table; bad input gives `ERROR: <message>`, exit 2.
2. **Align** genes across tables by `--align`: `intersect` keeps genes in every
   table; `union` keeps all, with NA where a table lacks the gene. Genes and cell
   types are sorted by name.
3. **`mash_set_data(Bhat, Shat)`**. It sets a missing cell to `Bhat = 0`,
   `Shat = 1e6`, so it carries no information. With `--null-corr em`, the null
   correlation `V` is estimated by `mash_estimate_corr_em` over the canonical
   covariances and set with `mash_update_data`, before any ED step; with
   `identity` it stays at identity (see "Null correlation V" below).
4. **Covariances:** `cov_canonical`. With `canonical+ed`, also: strong genes are
   those with lfsr < `--strong-lfsr` in any cell type in `mash_1by1`;
   `cov_pca(npc = K)` on them seeds `cov_ed`. K is the number of cell types, so
   no PC is dropped and npc is no choice. Fewer strong genes than K is an error.
5. **`mash(data, U, seed = --seed)`**, then `get_pm`, `get_psd`, `get_lfsr`. As in
   dreamlet's `run_mash`, a cell the gene was not measured in is set to NA.
6. **Pairwise sharing:** for each pair, mashr's `get_pairwise_sharing(factor =
   --sharing-factor, lfsr_thresh = --lfsr-max)` on the genes measured in both.
7. **Raw Spearman** of the input effects per pair, over genes measured in both
   (C-an-40).
8. **Count** genes by the number of cell types with lfsr < `--lfsr-max`.
9. **Write** tables, figures, report and result. Every refusal comes before the
   outdir is made, so a refused run leaves no outdir.

## Output structure

`posterior.tsv` is long, sorted by gene symbol then cell type. Every gene of
every input table has a row: `status` is `tested`, or `not_tested` with NA
statistics when `--align intersect` dropped it. A gene a table never held has
no row for that cell type. Counts, figures and sharing use `tested` rows only;
`result.json` counts rows per status (`status_counts`). The pair tables are sorted by `cell_type_a`, `cell_type_b`.

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # context, n_genes, cell_types, genes_dropped_by_align, n_shared_ids, joined_on, status_counts, n_strong, null_corr, n_significant_any, n_shared_all, n_specific, significant_per_cell_type, scatter_pair, params, seed, reference, chat_summary_lines, preferred_artifacts
├── tables/
│   ├── posterior.tsv          # gene_id, gene_symbol, cell_type, raw_effect, raw_se, pm, psd, lfsr, significant, status
│   ├── significant_counts.tsv # n_cell_types (1..K), n_genes
│   ├── pairwise_sharing.tsv   # cell_type_a, cell_type_b, sharing
│   ├── raw_spearman.tsv       # cell_type_a, cell_type_b, n_genes, rho, pvalue
│   └── null_corr.tsv          # cell_type_a, cell_type_b, v (0 everywhere with identity)
├── figures/
│   ├── sharing_heatmap.pdf
│   ├── sharing_heatmap.png
│   ├── sharing_heatmap.csv
│   ├── significant_hist.pdf
│   ├── significant_hist.png
│   ├── significant_hist.csv
│   ├── raw_vs_posterior.pdf
│   ├── raw_vs_posterior.png
│   └── raw_vs_posterior.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every --de
    ├── environment.txt        # Python packages, Rscript path, R, mashr, ashr versions
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

## Gotchas

- **The model will want to call `get_pairwise_sharing(m)` with its defaults.**
  Do not. Its `lfsr_thresh = 0.05` and `factor = 0.5` are hidden thresholds; the
  script passes `--lfsr-max` and `--sharing-factor`.
- **The model will want to fill a missing gene with logFC 0** to get a full
  matrix. Do not. A zero with a real SE is evidence of no effect; `--align
  union` marks it unmeasured (SE 1e6) and reports NA there.
- **The model will want to read lfsr as an FDR** or adjust it again. Do not.
  lfsr is a per-gene, per-cell-type local false sign rate, not an FDR over a
  list and not the probability that the effect is zero; report the counts as
  they are.
- **The model will want to read a sharing fraction as P(shared).** Do not. It is
  the fraction of genes significant in either type whose posterior means agree
  in sign and within `--sharing-factor`; it is not a probability that any gene
  is shared.
- **The model will want to call a gene "specific" because it passes in one
  cell type.** Do not say more than the table shows. One pass and three misses
  can still be a shared effect that is underpowered elsewhere; point to the
  sharing heatmap and the posterior means.
- **The model will want to read shared effects as a shared mechanism** or say
  one cell type "drives" another. Do not. Sharing is agreement in sign and size
  of posterior effects; say "shared effect" or "associated with" (house-style
  `claims.strength`).
- **The model will want to read `not_tested` as no effect.** Do not. The gene
  was missing from another table; it was not tested in this run.
- **The model will want to switch to `canonical`** when `canonical+ed` stops on
  too few strong signals, or to relax `--strong-lfsr`. Do not. Ask the user;
  both change the prior.
- **The model will want to mix contrasts or effect scales** across tables (log2
  in one, natural log in another, or flipped numerator). Do not. mashr assumes
  the same contrast on the same scale; ask the user.
- **The model will want to add a residual-correlation step** because dreamlet's
  `run_mash` estimates `V` with `mash_estimate_corr_em`. Do not patch it in.
  This skill keeps `V` at identity; if the cell types share donors, say so in
  the report and raise it with the user (see "Null correlation V").
- **The model will want to guess `--reference`** from gene symbols. Do not. Ask
  the user for the genome build and annotation release.

## Null correlation V

`mash.R` keeps the null correlation `V` at identity (`mash_set_data(Bhat, Shat)`
with no `V`). Identity assumes the per-cell-type DE tables come from independent
cells and samples, so their null errors do not correlate. dreamlet's `run_mash`
instead estimates `V` with `mash_estimate_corr_em(data, U.c)` (checked in
dreamlet 1.8.0 source); cell types measured in the
same donors usually have correlated errors, and identity then overstates
sharing evidence. `--null-corr` is required and has no default: `em` when the
tables share donors (dreamlet's choice), `identity` only when they come from
independent cells and samples. Its grill question: do these DE tables come from
the same donors?

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
