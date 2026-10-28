# enrichment: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`;
rows from `docs/clawbio-determinism/parts/`.

| Row | What it gave | Source |
|---|---|---|
| B-an-22 (`B_pbmc.md:31`) | fgsea on Hallmark, seeded, `padj < 0.05` | `repos/IHA-Figure/Figure4/08_DEGs_and_EnrichmentAnalysis/02_GSEA.ipynb:131-139` |
| C-an-25 (`C_disease.md:33`) | `minSize=10, maxSize=500, nPermSimple=10000, BPPARAM` | `repos/ALTRA-manuscript/Analysis/scRNA/ALTRA_scRNA_R_helper_functions.r:492-570`; `repos/multiple-myeloma/analyses-rna/04-pseudobulk/l3-pseudobulk/15-all-fgsea-no-filter.R:109-125` |
| C-an-27 (`C_disease.md:35`) | ORA with BH, MSigDB GMTs | `repos/IFN_manuscript/Figures/Supplementary/IFN_Sup3.ipynb:264,380-383,516` |
| C-an-38 (`C_disease.md:46`) | overlap test: odds ratio + p (GeneOverlap) | `repos/single_cell_aging/09_Overlap_Age_DE_psychopathology_DE/01_Overlap_Age_DE_&_psychiatric_Disease_DE_v2.ipynb:61,331-332` |
| D-an-15 (`D_reflex.md:23`) | fgsea Hallmark + GO per treatment | `repos/reflex_manuscript_v2/Fig2g_h_VST3/06-R_VST3_Fig2h.ipynb:260-300,374-384` |
| E-an-08 (`E_methods.md:16`) | zenith (later engine) | `repos_new/dreamlet_analysis/PsychAD_r0/PsychAD_r0_analysis.Rmd:743-753` |
| C-pl-12 (`C_disease.md:75`) | NES dotplot | `repos/multiple-myeloma/analyses-rna/04-pseudobulk/l3-pseudobulk/16-plot-fgsea-no-filter.ipynb:209,266` |
| D-pl-37 (`D_reflex.md:93`) | GSEA dot plot, size and colour channels | `repos/immgenT-GP-analysis/script/FigureS3.R:107-125` |
| D-pl-10 (`D_reflex.md:66`) | NES heatmap with significance marks | `repos/reflex_manuscript_v2/Fig2g_h_VST3/06-R_VST3_Fig2h.ipynb:1312-1460` |
| B-pl-36 (`B_pbmc.md:105`) | top-pathway bar (not drawn; the dot plot carries it) | `repos/MOCHA_Manuscript/Fig4/Figure4.R:74-78` |

The running-ES curve is the PIPELINES.md "std" plot (no corpus row).

Reference skills: `decoupler` and `S/decoupler-enrichment` (see
`docs/pipeline-skills.md`, "Reference skills to check against").
Checked against the installed skills `decoupler`
(`~/.claude/skills/decoupler/SKILL.md`) and `decoupler-enrichment`
(`decoupler-enrichment/SKILL.md:69-96,160-165`). Both agree on three points.
A DEG list goes to ORA and a ranked statistic to GSEA. The Wald `stat` from DE
is fed in directly, without normalisation. Sets below a minimum size are
dropped (decoupler `tmin`, default 5). This skill makes that minimum the
required `--min-size`, with no default. decoupler's ORA background is the
genes in the input, which matches `--universe tested`. Use decoupler for
ULM/MLM activity; this skill does not run it.

## ClawBio mechanisms used

IDs from `docs/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" |
| F-04 | Bundled demo data | `examples/demo_*.tsv`, `examples/demo.gmt`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-11, F-13, F-14 | Fixed plot spec, thresholds drawn, style, Agg | vendored `plots.py`, `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json`; `set.seed` per contrast in `fgsea.R` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` |
| F-23, F-26 | Thresholds documented, closed option sets | parameter table; `--method`, `--universe`, `--direction` choices |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", "Agent boundary" |
| F-27 | One skill, one task; hand-off by table | reads pseudobulk-de's `de_results.tsv` |
| F-28 | Golden-value demo test | `tests/test_enrichment.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

## Running and testing

The demo (`examples/make_demo.py`, seed 0) has 1000 genes and two contrasts. In
`T`, `SET_UP` (40 genes) has its Wald stat shifted +3 and `SET_DOWN` -3; `Mono`
is null. `SET_TINY` (8 genes) falls below `--min-size`. With `--demo`, fgsea
tests 9 sets per contrast and calls exactly `SET_UP` (NES +3.37) and `SET_DOWN`
(NES -3.03) in `T`, none in `Mono`. `--demo --method ora` gives a foreground of
31 genes in `T` (0 in `Mono`) and the same two significant sets.

fgsea needs Rscript from `RSCRIPT` or `PATH`. On this HPC:

```bash
export RSCRIPT=/share/software/statistical/R/R-4.5.2/gcc-8.5.0/R-4.5.2/bin/Rscript
export LD_LIBRARY_PATH=/home/mdmanurung/.local/lib/conda-gcc-runtime
export R_LIBS_USER=/exports/para-lipg-hpc/mdmanurung/R/4.5   # fgsea 1.36.2
```

From `skills/enrichment` (fgsea tests skip without Rscript):

```bash
UV_CACHE_DIR=$TMPDIR/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with numpy,pandas,scipy,matplotlib,plotnine-extra,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

## Citations

- Korotkevich et al. (2021). Fast gene set enrichment analysis. bioRxiv
  060012 (fgsea).
- Subramanian et al. (2005). Gene set enrichment analysis. PNAS 102:15545.
- Liberzon et al. (2015). The MSigDB hallmark gene set collection. Cell
  Systems 1:417.
- Benjamini, Hochberg (1995). Controlling the false discovery rate. JRSS B
  57:289-300.

## Biology-core checks (2026-10-28)

- `--context` is required and opens the report and `chat_summary_lines`.
- `enrichment.tsv` keeps every set in the file per contrast with a `status`
  column: a set outside `--min-size`/`--max-size` is `not_tested` with NA
  statistics; `result.json` has `status_counts`. Figures use tested rows only.
- Gotchas keep claims at "enriched among" or "associated with" (house-style
  `claims.strength`) and say a `not_tested` set is not "not enriched".
- Gene ids do not apply: the output is keyed by gene set, and sets join the DE
  table on the `gene` column as given. If pseudobulk-de renames `gene` to
  `gene_symbol`, `read_de` and the input contract must follow.
- `check_compatible` does not apply: each `--de` table is tested on its own.
