# composition: maintainer notes

Moved from SKILL.md so the agent reads only what it needs to run the script.

## Corpus sources

Paths relative to
`/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`.
Row IDs from `mycelium-extra/plans/clawbio-determinism/parts/*.md`.

| Row | Procedure | Source |
|---|---|---|
| E-an-10 (`E_methods.md:18`) | `crumblr(cellCounts(pb))`, `dream(cobj, form, colData(pb))`, `eBayes`, `topTable(coef=)` with `(1\|donor)` | `repos_new/dreamlet_analysis/Nathan_NatImm_2021/Nathan_NatImm_2021.Rmd:487-520`; `repos_new/dreamlet_analysis/archived/Mathys_Nature_2019/Mathys_Nature_2019_crumblr.Rmd:78-110` |
| B-an-15 (`B_pbmc.md:24`) | Cell counts per sample x label, CLR with a pseudocount | `IHA/Extended_Figure9/05_FigE.ipynb:77` |
| B-an-16 (`B_pbmc.md:25`) | Wilcoxon on CLR + BH (the test this skill replaces with a model) | `IHA/Extended_Figure4/02_IMM_Age_Metric/01_Calculate_Scores.ipynb:143` |
| C-an-20 (`C_disease.md:28`) | Paired Wilcoxon on frequencies, FDR | `repos/multiple-myeloma/analyses-rna/04-pseudobulk/l3-pseudobulk/03-run-timepoint_freq_analysis.ipynb:58-87,110,125` |
| C-an-22 (`C_disease.md:30`) | CLR + `lmer(clr ~ time + ... + (1\|subject))` | `repos/ALTRA-manuscript/Analysis/scRNA/ALTRA_scRNA_R_helper_functions.r:124-180` |
| C-an-23 (`C_disease.md:31`) | scCODA (later engine) | `repos/ALTRA-manuscript/Analysis/scRNA/Differential_abundance/02a_scCODA_reanalysis_DA_AIM1.ipynb:963-966,999` |
| C-an-24 (`C_disease.md:32`) | Proportions from `xtabs(~ sample + celltype)`, `lm` per cell type | `repos/single_cell_aging/10_Cell_type_composition_analysis/01_Cell_type_composition_analysis.ipynb:273,378` |
| D-an-09 (`D_reflex.md:17`) | Count and proportion per group and cluster | `R2/Fig2b_c_d_VST2/04-R_Fig2b.ipynb:597-650` |
| A-pl-05 (`A_skin.md:54`) | Stacked proportion bar from a row-normalised crosstab | `SA/figs/3_Fig3/3a_bar_plots.ipynb:217-252` |
| C-pl-08 (`C_disease.md:71`) | Stacked composition bar | `repos/ALTRA-manuscript/Figures/Figure1/Figure_1CD_S2DEF.ipynb:1948,2019` |
| D-pl-08 (`D_reflex.md:64`) | Stacked bar of category proportion per group | `R2/Fig1b_Jurkat_Raji/03-R_Fig1b.ipynb:984-1030` |
| E-pl-15 (`E_methods.md:57`) | Per-sample fraction box by group | `repos_new/dreamlet_analysis/Kfoury_CancerCell_2021/Kfoury_CancerCell_2021.Rmd:447-460` |
| C-pl-19 (`C_disease.md:82`) | Effect plot of DA estimates (here a CI forest) | `repos/ALTRA-manuscript/Figures/Figure1/Figure_S3E.ipynb:443,971,1156-1157` |

## ClawBio mechanisms used

IDs from `mycelium-extra/plans/clawbio-determinism/parts/F_clawbio_patterns.md`.

| ID | Mechanism | Here |
|---|---|---|
| F-01, F-02 | Agent boundary, mandatory run path | "Agent boundary" above |
| F-04 | Bundled demo data | `examples/demo.h5ad`, `examples/make_demo.py` |
| F-05 | Input contract table | "Input contract" |
| F-06, F-07 | Fixed fail strings, `ERROR:` + exit 2, no traceback | `ERRORS` dict and `skillkit`; one test per string |
| F-08, F-09 | Fixed output tree, output-contract test | "Output structure"; the test parses this tree |
| F-10, F-11, F-13, F-14 | Fixed plot spec, style, Agg | vendored `plots.py`, `theme.py` |
| F-16, F-17 | Seed as a flag, echoed into provenance | `--seed` in `commands.sh`, `result.json` |
| F-18, F-19, F-20 | `commands.sh`, environment, checksums | `reproducibility/` |
| F-23, F-26 | Thresholds documented, no defaults | parameter table |
| F-25, F-30 | "The model will want to ... Do not", "Do NOT fire when" | "Gotchas", "Agent boundary" |
| F-27 | One skill, one task | composition only; no tree test, no scCODA |
| F-28 | Golden-value demo test | `tests/test_composition.py` |
| F-33 | Structured result handoff | `chat_summary_lines`, `preferred_artifacts` |

F-03 (demo when no input) is replaced on purpose: the demo runs only when the
user asks for it.

## Running and testing

The R step needs Rscript from `RSCRIPT` or `PATH`. On this HPC:

```bash
export RSCRIPT=/share/software/statistical/R/R-4.5.2/gcc-8.5.0/R-4.5.2/bin/Rscript
export LD_LIBRARY_PATH=/home/mdmanurung/.local/lib/conda-gcc-runtime  # libstdc++ with GLIBCXX_3.4.29
export R_LIBS_USER=/exports/para-lipg-hpc/mdmanurung/R/4.5           # crumblr 1.2.0, variancePartition 1.40.2
```

```bash
python scripts/composition.py \
  --adata data/cells.h5ad --sample-col sample --label-col cell_type --donor-col donor \
  --formula "~ condition + (1|donor)" --coef conditionstim --group-col condition \
  --tau 1 --min-count 50 --padj-max 0.05 --seed 0 --outdir results/2026-10-27_da

python scripts/composition.py --counts counts.tsv --meta meta.tsv ...  # counts mode
python scripts/composition.py --demo --outdir demo_out                 # only on request
```

**Demo.** `examples/make_demo.py` (seed 0) writes 8 donors x ctrl/stim, 300-500
cells per sample, 7 labels; Mono is planted 3-fold up and B 3-fold down in stim,
and pDC has about 2 cells per sample. With the demo values: 16 samples from 8 donors, pDC
dropped, 6 labels tested, significant **Mono up** (logFC 1.08) and **B down**
(logFC -0.89); CD4 T, CD8 T, DC and NK are not (padj >= 0.16). Each R call takes
about 3 minutes here, most of it loading packages.

From `skills/composition`:

```bash
UV_CACHE_DIR=/tmp/uvcache PYTHONDONTWRITEBYTECODE=1 \
  uv run --python 3.12 --with anndata,matplotlib,scipy,pytest \
  python -m pytest -q -p no:cacheprovider tests
```

Regenerate the demo with `python examples/make_demo.py`; the golden values in
the tests change with it.

## Citations

- Hoffman GE, Roussos P (2021). dream: powerful differential expression analysis
  for repeated measures designs. Bioinformatics 37:192-201.
- crumblr: <https://gabrielhoffman.github.io/crumblr/> (Hoffman GE et al.,
  count ratio uncertainty modeling for cell-type composition).
