# sc-de: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`.
Row IDs from `docs/clawbio-determinism/parts/{B_pbmc,C_disease}.md`.

| Row | Procedure | Source |
|---|---|---|
| C-an-17 | MAST per donor with a CDR latent variable, seeded downsampling to 1000 cells per group | `repos/IFN_manuscript/Scripts/IFN_N1_Downsample_MAST.ipynb:147-165,263,299` |
| C-an-18 | MAST on an ISG-restricted gene set, `min.pct = 0.01` | `repos/IFN_manuscript/Scripts/IFN_Scores_Run_MAST_NMF.ipynb:188-192,568-569` |
| B-an-20 | Olink `lm(NPX ~ cohort + CMV + sex)` per protein | `repos/IHA-Figure/Figure1/01_Olink_Analysis/01_BRI_Olink_Y1D0.ipynb:129` |
| C-an-21 | Olink `lm(NPX ~ status + age + sex + BMI)` per assay, then q-values | `repos/ALTRA-manuscript/Analysis/OLINK_analysis/02-R_aim1_ARI_vs_CON_modeling.ipynb:295` |

Reference skill: the installed `muscat` skill (`~/.claude/skills/muscat/SKILL.md`, and
`references/ds-analysis.md` on `mmDS`) covers the donor-aware alternatives:
`pbDS` pseudobulk and `mmDS` cell-level mixed models with a sample random
effect. It has no MAST path; its mixed-model track agrees with this skill's
`(1|donor)` option and its pseudobulk track with `pseudobulk-de`.

## ClawBio mechanisms used

IDs from `docs/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" |
| F-04 | Bundled demo data | `examples/demo.h5ad`, `examples/demo_npx.tsv`, `make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-11, F-13, F-14 | Fixed plot spec, thresholds drawn, style, Agg | vendored `plots.py`, `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` |
| F-23, F-26 | Thresholds documented, closed option sets | parameter table; `--mode`, `--cdr`, `--log-base` choices |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", "Agent boundary" |
| F-27 | One skill, one task | per-feature DE only; hand-off through `de_results.tsv` |
| F-28 | Golden-value demo test | `tests/test_sc_de.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

## Running and testing

```bash
python scripts/sc_de.py --mode mast --adata data/cells.h5ad \
  --cell-type-col cell_type --cell-type "CD14 Mono" --donor-col donor \
  --formula "~ condition" --contrast-term condition,stim,ctrl \
  --n-cells-per-donor 1000 --min-frac 0.01 --cdr yes --log-base e \
  --padj-max 0.01 --lfc-min 0.2 --seed 0 --reference "GRCh38, GENCODE 44" \
  --gene-id-col none --context "human PBMC, IFN-beta vs control, Kang 2018" \
  --outdir results/2026-10-27_mono-mast

python scripts/sc_de.py --mode lm --input olink_long.tsv \
  --formula "~ status + age + sex + bmi" --term status,ARI,CON \
  --padj-max 0.05 --lfc-min 0 --seed 0 --reference "Olink Explore 3072" \
  --context "human plasma, ARI vs control, Olink cohort" \
  --outdir results/2026-10-27_olink-lm

python scripts/sc_de.py --demo --outdir demo_mast           # only on request
python scripts/sc_de.py --demo --mode lm --outdir demo_lm
```

The demo data come from `examples/make_demo.py` (seed 0). Mast demo: 600 T cells
from 6 donors (50 per donor and condition), 200 genes tested, 16 significant:
the 15 planted genes G000-G014 and G141, a cell-level false positive. Lm demo:
40 samples, 50 proteins, 8 significant, exactly the planted P01-P08; with
`--term age --lfc-min 0` exactly P09-P12. The mast demo takes about a minute.

R comes from `$RSCRIPT`, else `Rscript` on PATH. On this HPC, from
`skills/sc-de`:

```bash
export RSCRIPT=/share/software/statistical/R/R-4.5.2/gcc-8.5.0/R-4.5.2/bin/Rscript
export R_LIBS_USER=/exports/para-lipg-hpc/mdmanurung/R/4.5
export LD_LIBRARY_PATH=/home/mdmanurung/.local/lib/conda-gcc-runtime
UV_CACHE_DIR=$TMPDIR/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with anndata,numpy,pandas,scipy,matplotlib,plotnine-extra,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

Tests that need R skip when no Rscript is found.

## Biology-core checks (2026-10-28)

From `docs/pipeline-skills.md` "Biology-core checks". Applied:

- context: required `--context` in both modes, first line of the report and
  `chat_summary_lines` (the cell-unit warning moved to line 2).
- status: `de_results.tsv` keeps every input feature with `status`; mast marks
  genes never detected in the kept cells `not_detected` and genes under
  `--min-frac` `not_tested`; lm marks features with no estimable term or no rows
  left `not_tested`. `status_counts` in `result.json`; the volcano uses tested rows.
- gene ids, mast only: required `--gene-id-col` (demo `gene_ids`, added to
  `examples/make_demo.py` after the random draws; `demo_npx.tsv` is unchanged).
  `check_gene_ids` replaces `check_unique(var_names)`. Mast writes `gene_symbol`
  and `gene_id` (the names mashr-sharing reads); lm keeps `feature`.
- claim wording: one gotcha says a DE call is an association, not a mechanism.

Not applied: gene ids in lm mode (its features are proteins or assays named by
the input table, not genes), compatible inputs (one input file).

## Citations

- Finak et al. (2015). MAST. Genome Biology 16:278.
- Benjamini, Hochberg (1995). Controlling the false discovery rate. JRSS B 57:289-300.
- Crowell et al. (2020). muscat. Nature Communications 11:6077.
