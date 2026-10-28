---
name: enrichment
description: >-
  Gene-set enrichment on one or more DE tables (contrasts or cell types):
  preranked fgsea in R, or over-representation (one-sided Fisher /
  hypergeometric) of the significant genes against a stated universe, BH within
  each contrast, with an NES dot plot, an NES heatmap, running enrichment-score
  curves, a report and a reproducibility bundle. Fire on "GSEA on the DE
  results", "fgsea Hallmark per cell type", "pathway enrichment of the DE
  genes", "over-representation analysis", "ORA with a background", "which
  pathways are up in stim", "NES dotplot".
  Do NOT fire for scoring a gene set per cell or sample (signature scoring),
  for TF or pathway activity inference (decoupler ULM/MLM, PROGENy), for
  enrichment that needs the dreamlet fit object (zenith), or when there is no
  DE table yet (run pseudobulk-de first).
license: MIT
metadata:
  version: "0.1.0"
  domain: transcriptomics
  tags:
    - enrichment
    - GSEA
    - fgsea
    - ORA
    - MSigDB
  inputs:
    - name: context
      type: string
      description: species, tissue or cell type, condition and dataset
      required: true
    - name: de
      type: file
      format:
        - tsv
      description: one DE table per contrast, a gene_symbol column plus statistics
      required: true
    - name: gene_sets
      type: file
      format:
        - gmt
        - tsv
      description: gene sets, GMT or long TSV (set, gene)
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
      description: per-contrast counts and significant sets, params, chat summary
  dependencies:
    python: ">=3.11"
    packages:
      - numpy
      - pandas
      - scipy
      - matplotlib
      - plotnine-extra
    r:
      - fgsea
      - BiocParallel
  demo_data:
    - path: examples/demo_T.tsv
      description: synthetic DE table, 1000 genes, one planted up and one down set
    - path: examples/demo_Mono.tsv
      description: synthetic null DE table, 1000 genes
    - path: examples/demo.gmt
      description: 10 gene sets
  endpoints:
    cli: python scripts/enrichment.py --de {name}={de} --gene-sets {gene_sets} ... --outdir {outdir}
  trigger_keywords:
    - GSEA
    - fgsea
    - gene set enrichment
    - pathway enrichment
    - over-representation
    - Hallmark
---

# Gene-set enrichment

You are the dispatcher for **enrichment**. The script does the analysis; you
collect its arguments, run it, and explain its outputs.

## Agent boundary

Run `scripts/enrichment.py`. Never re-implement the analysis, change a
threshold, or pick a value the user has not given. If a required argument is
unknown, ask the user (these are grill questions), or run `--demo` only when the
user asks for a demo. Never download gene sets: the user supplies the file.

Then read `result.json` and relay `chat_summary_lines` and `preferred_artifacts`.
If the script prints `ERROR: ...` and exits 2, tell the user the message and the
argument it names. Do not work around it.

**Fire when the user says:** "run GSEA on the pseudobulk DE per cell type",
"fgsea with Hallmark on the stat column", "which GO terms are over-represented
in the up genes", "NES dotplot across contrasts".

**Do NOT fire when:** the user wants a per-cell or per-sample set score
(signature scoring), TF/pathway activity (decoupler), zenith on a dreamlet fit,
or there is no DE table yet.

## Input contract

| Input | Format | Required content | Example |
|---|---|---|---|
| `--de name=path` (repeat per contrast) | TSV, one row per gene | `gene_symbol` (joined to the gene sets); fgsea: `--rank-col` (numeric); ora: `padj`, `log2_fc` (numeric) | pseudobulk-de `tables/de_results.tsv` |
| `--gene-sets` | `.gmt` (name, description, genes) or `.tsv` | GMT: unique set names; TSV: columns `set`, `gene` (long) | `examples/demo.gmt` |

Gene IDs must be the same type in both files (symbols, or Ensembl). A contrast
name holds only letters, digits, `.`, `_` and `-`. Rows with an NA in a used
column leave that contrast (DESeq2 independent filtering leaves NA `padj`).

Gene sets are never downloaded at run time. Export MSigDB Hallmark once with
msigdbr (data ship with the package) and pass the file, recording the version:

```r
h <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
sets <- unique(data.frame(set = h$gs_name, gene = h$gene_symbol))
write.table(sets, "hallmark_msigdb_2026.1.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

## Parameters (all required, no defaults)

Each row is a grill question. Ask the user; never fill one in.

| Argument | Grill question | Demo-only value |
|---|---|---|
| `--de` | Which DE tables, and what is each contrast called (`name=path`, repeat)? | `T=examples/demo_T.tsv`, `Mono=examples/demo_Mono.tsv` |
| `--gene-sets` | Which gene-set file, from which source and version? | `examples/demo.gmt` |
| `--method` | Preranked GSEA (`fgsea`) or over-representation (`ora`)? | `fgsea` |
| `--min-size` | Smallest set size (genes in the set and in the contrast) to test? | `10` |
| `--max-size` | Largest set size to test? | `500` |
| `--padj-max-sets` | Which BH-adjusted p marks a set significant? | `0.05` |
| `--limit` | At which absolute NES (or log2 odds ratio) do the plot colours saturate? | `2.5` |
| `--seed` | Which seed? fgsea's permutations use it; ORA draws none. | `0` |
| `--context` | Species, tissue or cell type, condition and dataset? | `synthetic demo: 1000 genes, planted T contrast and null Mono contrast, no real tissue` |
| `--reference` | Which genome build and annotation release gave the gene IDs (e.g. `GRCh38, GENCODE 44`)? Recorded in `result.json`. | `synthetic demo, no genome build` |
| `--rank-col` (fgsea) | Which DE column ranks the genes (Wald `stat`, signed -log10 p, logFC)? | `stat` |
| `--n-perm-simple` (fgsea) | How many permutations for fgsea's simple step? | `10000` |
| `--universe` (ora) | Background: every tested gene (`tested`) or tested genes in any set (`tested-in-sets`)? | `tested` (with `--demo --method ora`) |
| `--padj-max` (ora) | Which gene padj marks a DE gene? | `0.05` (ora demo) |
| `--lfc-min` (ora) | Which absolute log2 fold change must a DE gene exceed? | `0.5` (ora demo) |
| `--direction` (ora) | Foreground of `up`, `down` or `both` DE genes? | `both` (ora demo) |

The last column is **demo-only, not defaults**. `--demo` alone pins all of them.
With `--demo`, only `--method` and `--outdir` may also be given. A method
argument given to the other method is an error.

```bash
python scripts/enrichment.py \
  --de T=results/de_T/tables/de_results.tsv --de Mono=results/de_Mono/tables/de_results.tsv \
  --gene-sets refs/hallmark_msigdb_2026.1.tsv --method fgsea --rank-col stat \
  --min-size 10 --max-size 500 --n-perm-simple 10000 --padj-max-sets 0.05 \
  --limit 2.5 --seed 42 --reference "GRCh38, GENCODE 44" \
  --context "human PBMC, IFN-beta vs control, Kang 2018" --outdir results/2026-10-27_enrichment

# Demo (only when the user asks for one)
python scripts/enrichment.py --demo --outdir demo_out
```

## Workflow (fixed order)

1. **Validate** arguments: `--method` decides which method arguments it needs;
   `--de` items are unique `name=path` pairs.
2. **Read** each DE table: `gene` plus the used columns, numeric and finite;
   refuse a repeated gene (`skillkit.check_unique`); drop rows with an NA in a
   used column and count the drops per contrast.
3. **Read** the gene sets (GMT or long TSV); refuse a repeated GMT set name
   (`skillkit.check_unique`). Count the distinct set genes found in each DE
   table (`n_matched` of `n_in_resource`).
4. **Size-filter** per contrast: a set's size is its genes found in the
   contrast (fgsea: ranked genes; ora: the universe). Keep
   `--min-size <= size <= --max-size`. No overlap at all, or no set left, is an
   error.
5. **Test** per contrast:
   - `fgsea`: `scripts/fgsea.R` runs `fgsea::fgsea(minSize, maxSize,
     nPermSimple, BPPARAM = SerialParam())` with `set.seed(--seed)` before each
     contrast; effect = NES; leading edge from fgsea.
   - `ora`: foreground = genes with `padj < --padj-max` and `log2_fc` beyond
     `--lfc-min` in `--direction`, inside the universe. One-sided Fisher test
     (`scipy.stats.fisher_exact(alternative="greater")`, the hypergeometric
     upper tail, as R `phyper(k - 1, ...)`); effect = sample odds ratio;
     leading edge = the overlap genes. A contrast with no foreground gets p = 1
     for every set; all contrasts empty is an error.
6. **Adjust**: BH over the sets tested within each contrast (the family is one
   contrast). `significant` = `padj < --padj-max-sets`.
7. **Sort** by contrast (in `--de` order), `padj`, set; NA last. Every set in
   the file gets a row per contrast (`skillkit.with_status`): a set outside the
   size range is `not_tested`, with NA statistics.
8. **Draw** the dot plot and heatmap (both methods) and the running-ES curves
   (fgsea only).

There is no design matrix here, so no design-rank check runs; the DE skill
that made each table owns that check.

## Output structure

Tree for `--method fgsea` (the demo). `--method ora` writes no `running_es.*`,
and its table has `overlap` and `odds_ratio` in place of `nes`.

```
outdir/
├── report.md                  # params, summary, figures with captions
├── result.json                # method, n_contrasts, n_sets_in_file, contrasts (per contrast: n_genes, n_sets_tested, n_matched, n_in_resource, n_significant, significant_sets; ora: n_universe, n_foreground), dropped_per_group (NA rows per contrast), status_counts, params, seed, reference, context, chat_summary_lines (first: "Context: ..."), preferred_artifacts
├── tables/
│   └── enrichment.tsv         # contrast, set, set_size (genes in the file), size (found in the contrast), nes | overlap, odds_ratio, pvalue, padj, significant, leading_edge, status
├── figures/
│   ├── dotplot.pdf
│   ├── dotplot.png
│   ├── dotplot.csv
│   ├── heatmap.pdf
│   ├── heatmap.png
│   ├── heatmap.csv
│   ├── running_es.pdf
│   ├── running_es.png
│   └── running_es.csv
└── reproducibility/
    ├── commands.sh            # argv as run, and the resolved command with every argument
    ├── environment.txt        # Python packages; Rscript, R and fgsea versions, SerialParam
    └── checksums.sha256       # tables and figure PDFs; sha256sum -c from outdir
```

`leading_edge` is a comma-separated gene list. `status` is `tested`, or
`not_tested` for a set outside `--min-size`/`--max-size` in that contrast (too
few of its genes matched); those rows keep NA statistics, and `size` says how
many genes matched. Figures draw tested rows only. Figures come from the vendored
`plots.py` (`plot_heatmap`) and matplotlib in the same theme, saved by
`save_panel`, with no titles.

- `dotplot`: the 20 sets (display constant) with the lowest padj in any
  contrast × contrasts; area = -log10 padj, colour = NES (ora: log2 odds
  ratio) clipped at ±`--limit`, black outline = significant.
- `heatmap`: same sets and colours, `*` = significant, blank = not tested.
  Both captions state the fixed order: sets by lowest padj, contrasts in
  `--de` order, never clustered.
- `running_es`: one panel per contrast for its top table row; GSEA running sum
  with weight 1 (as `fgsea::plotEnrichment`), ticks at the set's genes.

## Gotchas

- **The model will want to download gene sets** (msigdbr, Enrichr, a URL)
  inside the run. Do not. The user names a file and its version; the export
  snippet above runs once, outside this skill.
- **The model will want to rank by `log2_fc`** because it "looks like GSEA".
  Do not choose. The rank column is a grill question; the Wald `stat` and
  signed -log10 p rank differently from a raw fold change.
- **The model will want to use the whole genome as the ORA universe.** Do not.
  The universe is the genes that were tested (`--universe`); genes never
  measured cannot be in the foreground and inflate every odds ratio.
- **The model will want to relax `--padj-max-sets`** or switch to raw p after
  seeing no significant set. Do not. Report the counts as they are.
- **The model will want to drop the null contrast** from the plot or table.
  Do not. A contrast with no significant set is a result.
- **The model will want to pool contrasts into one BH family** or merge
  several DE tables into one ranking. Do not. Each contrast is its own ranking
  and its own family.
- **The model will want to say "the genes in the set change".** Do not. Say
  "the set is enriched among the changed genes" (ora) or "the set sits at the
  top or bottom of the ranking" (fgsea). The leading edge is a subset, and
  many members do not move.
- **The model will want to call a `not_tested` set "not enriched".** Do not.
  Too few of its genes were in the DE table to test it.
- **The model will want to say the pathway "drives" the response or is
  "activated".** Do not. Enrichment is an association between a set and the
  ranking; say "enriched among" or "associated with" (house-style
  `claims.strength`). Activity or mechanism needs a perturbation.
- **The model will want to read an NES or a padj as a probability** ("95%
  likely active"). Do not. NES is a normalised enrichment score; padj is a
  BH-adjusted p for a set-level test.
- **The model will want to call a set gene that is absent from the DE table
  "unchanged".** Do not. It was not tested: filtered, not detected, or a
  different gene ID. `n_matched` of `n_in_resource` says how many were found.
- **The model will want to guess `--reference`** from the gene IDs. Do not.
  Ask for the genome build and annotation release.
- **The model will want to run zenith here.** zenith needs the dreamlet fit,
  not a DE table; it is a later engine.

## Maintainer notes

Corpus sources, ClawBio mechanisms, tests, citations and the plot spec are in `references/maintainer.md`. The agent running the script does not need them.
