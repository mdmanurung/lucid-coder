# Sources for what was built

Written 2026-10-27. This file names the repo and script behind each pattern in what was built.

## Citation conventions

- **Corpus paths** are relative to `$B=/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings`.
- **ClawBio paths** are relative to `C=~/.claude/plugins/cache/clawbio/clawbio/0.3.0`.
- **Notebook line numbers** count raw `.ipynb` JSON lines.
- **The "Row" column** points to the full entry in `CATALOGUE.md` and `parts/`. Those entries list up to 4 sources each.
- **Lines marked ✓** were opened and checked by the parent session on 2026-10-27.
- **Unmarked lines** were checked by the catalogue agents.
- **The "new" label** marks a choice made for this module that no corpus repo makes.

## `lucid-coder/templates/plots.R` and `plots.py`

### Shared mechanics

| Pattern | Taken from | Source |
|---|---|---|
| One tested script holds the plot choices. The agent calls it and does not redraw. | ClawBio | `C/templates/SKILL-TEMPLATE.md:243-247` (F-01) |
| `check_cols()` stops with one fixed message | ClawBio `rnaseq-de` input validation | `C/skills/rnaseq-de/rnaseq_de.py:53-138` (F-06) |
| An unknown palette label stops the script (`palette_for`, already in `palettes.R`) | immgenT-GP-analysis | `repos/immgenT-GP-analysis/code/R/level2_group_palette.R:40,46` ✓ |
| No hidden thresholds; every cut-off is a required argument | ClawBio states the rule (`C/skills/sample-qc-triage/SKILL.md:125-133,168`, F-23). The user decided on 2026-10-27 that cut-offs carry no default and become grill questions. | user decision; the ClawBio rule |
| Cut-off lines drawn on the plot | ClawBio `rnaseq-de` (F-11), which bakes in 0.05 and 1.0. Here the values are arguments instead. | `C/skills/rnaseq-de/rnaseq_de.py:352-375` |
| Headless backend (`Agg`), and `pdf(NULL)` in the R self-check | ClawBio | `C/skills/diff-visualizer/diff_visualizer.py:23` (F-14) |
| Self-check on toy data in the `__main__` / `sys.nframe()` block | existing lucid-coder `palettes.R`/`palettes.py`. ClawBio does the same with demo tests (F-28). | `lucid-coder/templates/palettes.R:30-33` |

### Per form

| Function | Pattern | Repo and script | Row |
|---|---|---|---|
| `plot_volcano` | dashed cut-off lines at ±lfc and -log10(padj) | multiple-myeloma `analyses-rna/04-pseudobulk/l3-pseudobulk/13-plots-deseq2-volcano_fixed.ipynb:94-95` ✓ | C-pl-04 |
| | `grey` / `darkred` for not-significant / significant | dreamlet_analysis `PsychAD_r0/PsychAD_r0_compare_methods.Rmd:366` ✓; variancePartition_analysis `plotVolcano.R:39` ✓ | E-pl-06 |
| | significant genes labelled | IHA-Figure `Figure4/05_FigF.ipynb:103-116`; MOCHA `Fig5/Fig5_ZeroInflated_promoters.R:252` | B-pl-11 |
| `plot_box` | `geom_boxplot(outlier.shape = NA)` + jitter | ALTRA-manuscript `Figures/Figure1/Figure_1CD_S2DEF.ipynb:174,2534` ✓ | C-pl-05 |
| | jitter with a fixed seed | rule from lucid-coder `rules/r.md` (`seed =` on `position_jitter()`). Most corpus jitter is unseeded. | `parts/D_reflex.md` Recurring defaults |
| | paired lines across groups (`pair`) | multiple-myeloma `manuscript-figures/figure-4/fig_4g.ipynb:80-102`; ALTRA `Analysis/scRNA/ALTRA_scRNA_R_helper_functions.r:710-751` | C-pl-10 |
| | `(n = k)` in tick labels | house-style row `fig.labels` (spatial_skin_atlas, variancePartition) | `skills/lucid-coder/references/house-style.md:125` |
| `plot_embedding` | `scattermore::geom_scattermore()` for dense clouds | immgenT-GP-analysis `script/Figure1.R:34,81` ✓ | D-pl-27 |
| | rasterised points in Python | sound-life `02-reference_labeling/07a-R_filter_review_plots.ipynb:79-164` | B-pl-10 |
| | UMAP coloured by label, no ticks | skin_fibroblast_atlas `figures/extended_fig1.ipynb:109,158` (`frameon=False`) | A-pl-01 |
| `plot_heatmap` | diverging scale centred at 0, symmetric `limits = c(-cap, cap)`, `oob = scales::squish` | immgenT-GP-analysis `code/R/cross_gp_helpers.R:16,22,27` ✓ | house-style `fig.heatmap` |
| | blue-white-red diverging fill | dreamlet_analysis `COVID_combat/COVID_combat_major.Rmd:246-278` (`scale_fill_gradient2(low="blue", mid="white", high="red")`). The hex values `#2166AC`/`#B2182B` are RdBu ends, as in ClawBio `RdBu_r` (F-13). **New:** the exact hex. | E-pl-17 |
| | no clustering; the caller's factor order | **new**, for determinism. Corpus heatmaps cluster with ComplexHeatmap or pheatmap defaults (C-pl-03). | — |
| `plot_dotplot` | dot size = fraction expressing, colour = mean | spatial_skin_atlas `figs/sFig1_to_11_dotplots_by_lineage/s1a_b.ipynb:430-441`; aifi `visualizations/00-Python_b_cell_markers.ipynb:327-331`; reflex v2 Seurat `DotPlot` `Fig1c_d_e_f_g_Benchmarking/02-R_Build_REFLEX_Object.ipynb:3110` | A-pl-02, B-pl-03, D-pl-07 |
| | `grey90` → red colour ramp | cns-haematopoiesis `color_map="Reds"` (`figures/Chapter_3/Bone_Atlas_Figs.ipynb:880-915`). **New:** the exact ramp. | D-pl-19 |
| `plot_composition` | stacked bar of proportions per sample | spatial_skin_atlas `figs/3_Fig3/3a_bar_plots.ipynb:217-252`; reflex v2 `Fig1c_d_e_f_g_Benchmarking/06-R_Fig1c.ipynb:332-371`; ALTRA `Figures/Figure1/Figure_1CD_S2DEF.ipynb:1948,2019` | A-pl-05, D-pl-08, C-pl-08 |
| `plot_scatter_fit` | `geom_smooth(method = "lm", se = TRUE)` | immgenT-GP-analysis `script/Figure2.R:140` ✓; GTestimate `notebooks/cta_seq_results.Rmd:405-440`; IFN `Figures/Main/IFN_Fig5.ipynb:1011-1016` | D-pl-26, E-pl-32, C-pl-09 |
| | Python fit colour `#3366FF` | ggplot2's `geom_smooth` default, chosen to match the R twin. **New** in Python. | — |

### Shared theme and export

Neither file below was built this session.

- **`theme.R` / `theme.py`:** `theme_bw(base_size = 8)`, mm sizes, `cairo_pdf`, `pdf.fonttype = 42`. The house-style rows `fig.theme`, `fig.size` and `fig.export` cite their sources (`house-style.md:118-120`).
- **`palettes.R` / `palettes.py`:** keyed palettes. The house-style row `fig.palette` cites its source (`house-style.md:117`).

## `mycelium-extra/skills/grill/references/analysis-decisions.md`, row "Plot parameters"

| Pattern | Source |
|---|---|
| Plot cut-offs become one user-owned grill question instead of defaults | user instruction, 2026-10-27 |
| The parameters listed are the required arguments of `plots.R`/`plots.py` | `lucid-coder/templates/plots.R:7-13` (header) |
| Why not copy defaults across projects: the corpus disagrees | `CATALOGUE.md`, section "Where repos disagree" |
| Corpus examples of the disagreement | IHA `padj < 0.05`/`abs(log2FC) > 0.1` (B-an-13); multiple-myeloma `Q < 0.1`/`abs(log2FC) > 0.5` (C-pl-04, `13-plots-deseq2-volcano_fixed.ipynb:79-80` ✓); dreamlet FDR 0.05 (E-an-04) |

## `CATALOGUE.md` and `PIPELINES.md`

Every row ID resolves to a `path:line` in `parts/A_skin.md` to `parts/E_methods.md`. ClawBio mechanisms F-01 to F-36 resolve in `parts/F_clawbio_patterns.md`.

Some rows are only partly verified. `CATALOGUE.md` lists them under "Unverified rows".
