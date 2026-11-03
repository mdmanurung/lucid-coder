# kai findings for five pipeline skills

Issue #20, part of #19. Every count here is **tut**: seen in the kai tutorial
corpus (Zenodo 17660667, version `251121`, 2,845 notebooks, code cells only,
full-line comments stripped). Tutorials show what people demo, not published
practice. No value here is a default. Each proposal below is **PROPOSED** and
lands in a `SKILL.md` only after the user approves it.

Counts are notebooks of 2,845 unless a smaller denominator is named. Where one
repository supplies most hits, the repo count follows. Counts come from regexes
over cell text, so they are heuristics: they miss calls built from variables
and can catch look-alike strings.

## How to rerun

The corpus and extracts stay outside the repo. The scripts sit in
`/tmp/claude-149488/-exports-para-lipg-hpc-mdmanurung-lucid-coder/80951150-83a9-4eed-9340-2d41f52de923/scratchpad/kai/work20/`:

- `dumpc.py` writes `cells.jsonl`, one code cell per line (run with `python3 -I`).
- `qc.py cells.jsonl REGEX` counts notebooks, and notebooks per captured value.
- `co.py` counts notebooks with pattern A, then those also holding B.
- `alias.py` finds `v.layers["counts"] = v.X` followed by in-place `normalize_total(v)` or `log1p(v)`.
- `view.py` finds a bare subset `v = w[...]` followed by a write into `v`.
- `scaled_rgg.py` finds `sc.pp.scale(v)` followed by `rank_genes_groups(v)` without `use_raw` or `layer`.
- `repos.py` gives notebook and repo counts for one regex.
- `alias_test.py` (run with `pixi run python`) shows the counts-layer alias and sparse PCA on the repo's scanpy 1.12.4, anndata 0.13.4 and scikit-learn 1.7.2.

## Carried-over question: should sc-normalize keep log1p output sparse?

**Yes for `log1p`; the question does not arise for `pearson` or `dsb`.** scanpy
already returns a sparse matrix. `normalize_total` keeps CSR and only casts an
integer matrix to float32. `log1p_sparse` rewrites `.data` in place
(`scanpy/preprocessing/_simple.py`, scanpy 1.12.4). The one dense step is in
`skills/sc-normalize/scripts/sc_normalize.py:252`, `adata.X = dense(adata.X)`,
which also fills `layers[--layer-out]` with a dense copy. Pearson residuals are
dense by construction, because a zero count gives a nonzero residual. dsb output
is dense too.

The dense `X` exists for `tables()`. Two calls there need a dense array:
`adata.X.var(axis=0, ddof=1)`, which a scipy sparse matrix lacks, and
`adata.X.sum(axis=1)`, which returns a matrix rather than a vector.
`PCA(svd_solver="arpack").fit_transform` takes sparse input in scikit-learn ≥ 1.4.
`alias_test.py` confirms this on the installed 1.7.2.

Corpus evidence: 36 of 530 `log1p` notebooks assign a dense `X`. 66
notebooks densify `X` for any reason. In the 30 contexts read, the dense array
feeds a tool that needs one: torch, cupy, TFvelo, model training. Several
convert straight back. `X` is assigned a fresh `csr_matrix(...)` in 22
notebooks, for example `adata.X = adata.X.toarray()` around
`calculate_qc_metrics`, then `adata.X = csr_matrix(adata.X)`
(`BayraktarLab/cell2location/bayraktarlab_cell2location_cell2location_short_demo`).
So kai densifies at the point of use and keeps the stored matrix sparse.

**PROPOSED change (needs approval, not made):** for `--method log1p`, keep `X`
and `layers[--layer-out]` sparse in `normalized.h5ad`. Densify only inside
`tables()`, or compute the per-feature variance from sparse sums. Leave
`pearson` and `dsb` dense, and say so in the output tree. The workflow line
"Copy the dense normalised X" then names the method. Cost: a test that pins the
storage type per method.

## sc-normalize

### Findings

- **target_sum.** `normalize_total` is called in 376 notebooks. 128 leave `target_sum` unset, so scanpy scales to the median library size. Of those that set it: `1e4` in 195 (plus `10000` in 3), `1e6` in 19, `10e4` in 7, other values in 1-2 each. Legacy `normalize_per_cell(counts_per_cell_after=)` adds `1e4` ×47, `1e6` ×43 and `10000` ×12 in 102 notebooks.
- **Pearson residuals.** 2 notebooks in 2 repos call `normalize_pearson_residuals`. Neither sets `theta` or `clip`, so both take scanpy's defaults: `normalize_pearson_residuals(adata)` (`scverse/scanpy-tutorials/scverse_scanpy_tutorials_tutorial_pearson_residuals`).
- **dsb.** 2 notebooks, both muon tutorials. One sets `empty_counts_range=(1.5, 2.8)` with no isotypes (`scverse/muon-tutorials/scverse_muon_tutorials_1_cite_seq_pbmc_5k`). The other passes `isotype_controls=isotypes, random_state=1` (`..._1_tea_seq_pbmc`).
- **Scaling after normalising.** `sc.pp.scale` appears in 87 of 530 `log1p` notebooks. `max_value=10` appears in 69 notebooks across 31 repos.
- **`.raw` and `.X` shuffles.** `.raw =` in 256 notebooks (53 repos; PanGIAtlas 47). `X = <obj>.raw.X` restores in 17, `raw.to_adata()` in 29.
- **Counts rebuilt from log values.** 2 notebooks write inverted log values where counts belong. `sc_adata.X = (np.around(np.exp(sc_adata.X.toarray())) + 1).astype(int)` adds 1 where it should subtract 1, and it inverts normalised values, not counts (`scverse/spatialdata-notebooks/scverse_spatialdata_notebooks_03_annotate_visium_cell2location`). `nichecompass-reproducibility/..._data_preparation` writes `np.round(np.exp(X) - 1)` into `layers["counts"]` (3 cells). The one correct rebuild multiplies `expm1(X)` by the stored per-cell UMI totals over 10,000 (`theislab/Covid_meta_analysis/theislab_covid_meta_analysis_figure1cd`), so it needs the size factors.
- **Aliased counts layer.** `v.layers["counts"] = v.X` with no `.copy()` appears in 78 notebooks. 53 then run in-place `normalize_total(v)` or `log1p(v)` on the same object before reassigning it. That is 9 repos, 39 notebooks from `teichlab/PanGIAtlas`. With a float `X`, the layer is normalised and logged too. With an integer `X`, `normalize_total` casts to a new float array, so the layer survives. `alias_test.py` shows both outcomes. `layers["counts"] = X.copy()` appears in 171 notebooks.

### Proposed gotchas (PROPOSED)

- **The model will want to stash counts with `layers["counts"] = adata.X`** before normalising an input by hand. Do not. Without `.copy()`, the layer and `X` are one object. `normalize_total` and `log1p` then rewrite both when `X` is float, which is what `read_10x_h5` gives. 53 kai notebooks do this. The script keeps its own copy, so pass the raw file as it is. (tut, `alias.py`, `alias_test.py`)
- **The model will want to rebuild counts with `np.expm1` or `np.exp(X) - 1`** when only a log-normalised `X` is at hand. Do not. That gives size-normalised values, not counts; one kai tutorial even adds 1 instead of subtracting it. Only stored library sizes make the rebuild exact. Ask for the raw counts. This extends the existing `.raw`/`.X` gotcha with the corpus case.
- **The model will want to add `sc.pp.scale(max_value=10)`** as part of "normalising". Do not. Scaling is a separate, dense step for PCA; 87 kai log1p notebooks scale, but this skill's output is normalised, not scaled, and DE skills refuse scaled data.

### Proposed grill-question citations (PROPOSED)

- `--target-sum`: "To what total is each cell scaled first? (kai tutorials: `1e4` in 195 of 376 `normalize_total` calls, `1e6` in 19, unset (median) in 128.)"
- `--theta` and `--clip`: add "kai's 2 Pearson tutorials leave both at scanpy's defaults; that choice is still yours." It supports the existing gotcha against picking them silently.
- `--empty-droplets`: the existing gotcha's `empty_counts_range=(1.2, 2.8)` is from the obs corpus. kai's muon tutorial uses `(1.5, 2.8)`. Two sources with two ranges is a reason to ask.

## sc-qc

### Findings

- **Mito prefix.** `startswith('MT-')` in 80 notebooks, `'mt-'` in 40, both together in 7.
- **Fixed thresholds.** `filter_cells(min_genes=)` is set in 106 notebooks: `200` ×39, `100` ×11, `400` ×9, `500` ×6, then `700`, `600`, `350`, `250` and others. `filter_genes(min_cells=)` is set in 182: `3` ×60, `10` ×48, `20` ×29, `5` ×15, `1` ×14. `max_genes` is rare (18). Mito cutoffs appear in 59 notebooks on two scales. Fractions: `< 0.2` ×23, `0.15` ×6, `0.1` ×6. Percents: `15` ×5, `5` ×5, `10` ×4, `20` ×2, `50` ×3.
- **MAD outliers.** 15 notebooks, 12 of them in `lueckenlab/masterpraktikum_fibrosis_atlas` (4 repos). The pattern is the sc-best-practices one: `is_outlier(adata, "log1p_total_counts", 5)` and `is_outlier(adata, "pct_counts_mt", 3)` (`lueckenlab/.../masterpraktikum_fibrosis_atlas_1_preprocessing`). One notebook uses 7 MADs for counts.
- **Scrublet.** 28 notebooks, 11 repos. 5 pass `batch_key` and 6 call Scrublet inside a loop or function, probably per sample. 17 call it once on the whole object; the code alone cannot say whether that object holds one library. `expected_doublet_rate` is set in 9: `0.1` ×4, `0.02` ×4, `0.06` ×1. No notebook sets `random_state` or a fixed `threshold`. Both scanpy's and Scrublet's `random_state` default to 0, so these runs are seeded by default, not unseeded.
- **Subsets without `.copy()`.** Writes into a bare subset are covered under `py.views` in `house-style.md`. A broader regex (`view.py`) counts 332 notebooks across 90 repos that write `.obs`, `.var`, `.X`, `.layers`, `.uns` or call an in-place `sc.pp`/`sc.tl` on such a subset. anndata warns (`ImplicitModificationWarning`) and turns the view into a copy, so the write lands in the subset and the parent is untouched. About 70 of the 332 also mute all warnings, and `ImplicitModificationWarning` is named in 0 notebooks. A QC example: `gran_adata = gran_adata[gran_adata.obs['percent_mito'] < ...]`, then `sc.pp.normalize_per_cell(gran_adata, ...)` (`ShalekLab/mTB_CD8_Depletion/shaleklab_mtb_cd8_depletion_0_preprocessing`).
- **Muted warnings.** `filterwarnings("ignore")` with no category in 263 notebooks (64 repos; 264 in `house-style.md`, counted before comment stripping), plus `simplefilter("ignore")` in 54. `house-style.md` `code.warnings` already cites the first count.

### Proposed gotchas (PROPOSED)

- **The model will want to copy a tutorial threshold** ("min_genes 200, mito 20%, as in the scanpy tutorial"). Do not. kai tutorials use 200 in 39 of 106 notebooks, then 100 to 700. The mito cutoff ranges from 5% to 50%, and fractions and percents are mixed (`mt_frac < 0.2` and `pct_counts_mt < 20` are the same cut). Ask, and state the unit `--max-pct-mito` takes (percent).
- **The model will want to run Scrublet once on a merged object** because tutorials do. Do not; 17 of 28 kai notebooks call it once with no `batch_key`. This adds the corpus count to the existing pooled-object gotcha.

### Proposed grill-question citations (PROPOSED)

- `--mito-prefix`: "(kai: `MT-` 80 notebooks, `mt-` 40, both 7.)"
- `--min-genes`: "(kai: 200 in 39 of 106 notebooks; 100-700 in the rest.)"
- `--max-pct-mito`: "In percent. (kai cutoffs run 5-50%; 23 notebooks write it as the fraction 0.2.)"
- `--mads-counts`, `--mads-genes`, `--mads-mito`: "(kai: 5, 5 and 3 in the sc-best-practices pattern, 15 notebooks, mostly one repo.)"
- `--expected-doublet-rate`: "(kai: 0.1, 0.02 or 0.06, never tied to the 10x loading table.)" It supports asking for the table value.

## sc-de

### Findings

- **Cell-level tests.** `rank_genes_groups` appears in 217 notebooks (65 repos). The method is set in 151: `wilcoxon` 118, `t-test` 33, `t-test_overestim_var` 3, `logreg` 2. `corr_method` is set once.
- **Condition groups at cell level.** 26 notebooks group by a condition-like column (disease, condition, treatment, perturbation, label). One example is `rank_genes_groups(fibroblast, groupby="disease", method="wilcoxon", layer="normalized_counts")` (`YosefLab/scvi-hub-reproducibility/yoseflab_scvi_hub_reproducibility_fig3_hlca_emphysema`). Per-cell perturbations, such as Norman's `perturbation_name`, are the legitimate case.
- **The silent `.raw` switch.** `rank_genes_groups` tests `.raw` when it exists and `use_raw` is unset. Only 17 notebooks set `use_raw` (15 `False`, 2 `True`), and 6 set `layer`. 13 notebooks scale an object and then call `rank_genes_groups` on it without `use_raw` or `layer`; 11 of them get log values only because `.raw` was set earlier (`scaled_rgg.py`).
- **Counts in a log test.** 2 notebooks pass `layer="counts"` to `rank_genes_groups`, whose fold change assumes log1p values.
- **MAST.** `zlm(` appears in 3 notebooks, all in `theislab/organoid_dexamethasone_stimulation`: `zlm(formula = ~treatment + ngeneson, sca=sca_sub)` with `doLRT='treatmentVeh'`. The CDR covariate is in all 3, and none has a random effect. Their gene filter is `rowSums(assay(sca_sub)) != 0`, which drops only all-zero genes.
- **Downsampling seeds.** `sc.pp.subsample` appears in 72 notebooks; 57 leave `random_state` at its default 0, and 15 set it (42 ×6, 0 ×6). A data frame or `.obs` is sampled with `.sample()` and no `random_state` in 15 notebooks, against 2 with it. `np.random.choice` appears in 88 notebooks, 58 with no seed anywhere in the notebook. Those two are truly unseeded.

### Proposed gotchas (PROPOSED)

- **The model will want to point `--layer` at whatever `X` holds** because scanpy DE "just works" on it. Do not. scanpy silently tests `.raw` when it exists. In 11 kai notebooks, a scaled `X` gave sensible results only because `.raw` held log values. This script reads exactly `--layer` (or `X`) and never `.raw`, so ask which matrix holds the log-normalised values.
- **The model will want to run mast mode by condition because tutorials call `rank_genes_groups(groupby="disease")`.** Do not; 26 kai notebooks test a condition at cell level. This adds the corpus count to the existing donor-level gotcha.

### Proposed grill-question citations (PROPOSED)

- `--cdr`: "(kai's only MAST tutorials, 3 notebooks in one repo, include it: `~treatment + ngeneson`.)"
- `--min-frac`: "(kai's MAST tutorials drop only all-zero genes; the fraction is yours.)"
- `--n-cells-per-donor` and `--seed`: "(kai: 15 notebooks sample a data frame with no seed, and 58 call `np.random.choice` with no seed in the notebook; this script always seeds.)" It supports keeping `--seed` required.

## pseudobulk-de

### Findings

- **Aggregation.** `get_pseudobulk` appears in 30 notebooks (9 repos; `theislab/spatial_pancreas` 9). Of those setting each argument: `min_cells=10` ×14 of 16, `min_counts=1000` ×13 of 16 (`10000` ×2), and `mode="sum"` in 19 of 20. `layer="counts"` is set in 20 of 27, `raw_counts` ×3, `soupX_counts` ×3.
- **Normalised values summed.** One notebook sums a normalised layer into a pseudobulk model, then filters on `min_counts=1000` of those sums: `layer="normalised_counts"` (`scverse/multi-condition-comparisions/scverse_multi_condition_comparisions_breast_cancer`).
- **Filtering before DESeq2.** All 16 `DeseqDataSet` notebooks (7 repos) set `refit_cooks=True`. 14 of them run decoupler's `filter_by_expr` first, typically `min_count=10, min_total_count=15, large_n=10, min_prop=0.7`. None sets `alpha` or `independent_filter`, so PyDESeq2's independent filtering runs on top. `decoupler_rna_psbk` adds `filter_by_prop(min_prop=0.1, min_smpls=2)` as a third filter. The skill's own order is the same two-stage one: `--min-frac`, then independent filtering.
- **Paired design without the patient term.** Two liana notebooks run PyDESeq2 on Kang 2018 with `design_factors=condition_key` alone, though their markdown says the data are "8 pooled patient lupus samples, each before and after IFN-beta stimulation". The DESeq2 cell also starts with `%%capture` and `quiet = True`, so fit warnings are hidden (`saezlab/liana-py/saezlab_liana_py_targeted`, and its copy in `saezlab/lianaplus_manuscript`). Designs elsewhere: `["disease", "sex"]` ×2, `['disease', 'sex', 'assay']` ×1; the rest are single factors.
- **R engines.** `DESeqDataSetFromMatrix` 2, `DGEList` 6 (mostly differential abundance), `voom` 4, `dreamlet` 0.

### Proposed gotchas (PROPOSED)

- **The model will want to pre-filter the pseudobulk with `filter_by_expr`** (or `filter_genes`) before the script, as decoupler's tutorials do. Do not. The script filters by `--min-frac`, and DESeq2 filters again independently. A third filter changes which genes are tested and goes unrecorded. 14 of 16 kai DESeq2 notebooks stack two filters, and one stacks three.
- **The model will want to sum a normalised layer** because it is called `normalised_counts`. Do not; DESeq2 and dreamlet need raw counts. This adds the corpus case to the existing `counts ... must be non-negative integers` gotcha.
- **The model will want to drop the patient term on Kang-style data because the tutorial did.** Do not. The liana tutorial fits `condition` alone on 8 patients measured before and after stimulation. This adds a source to the existing gotcha.
- **The model will want to wrap a run in `%%capture` or `quiet=True`** to tidy a notebook. Do not; that hides convergence and Cook's warnings. Read `result.json` and the log instead.

### Proposed grill-question citations (PROPOSED)

- `--min-cells`: "(kai: 10 in 14 of 16 `get_pseudobulk` calls that set it.)"
- A per-sample count floor such as `min_counts=1000` (13 of 16 kai calls) has no flag in this skill. **PROPOSED for discussion only:** ask whether the user wants one, or keep it out and say why in `maintainer.md`. Do not add it without approval.
- `--design`: "Does a subject give more than one sample? (kai: the liana Kang tutorial omits the patient term for paired data.)"

## tf-activity

### Findings

- **Method.** ULM in 41 notebooks (11 repos) and MLM in 12 (7 repos). VIPER is absent; AUCell and wsum appear once each.
- **Networks.** CollecTRI in 34 notebooks and PROGENy in 32; DoRothEA appears in only 3. PROGENy `top=` is `5000` ×8, `500` ×6, `2500`, `300` and `100` ×1 each. `organism` is `human` ×32, `mouse` ×3.
- **`tmin`.** Set in 2 notebooks, both `tmin=3` on a marker set (`scverse/decoupler/scverse_decoupler_rna_sc`). Everywhere else it is left at decoupler's default.
- **Cells mode.** `rankby_group` or `rank_sources_groups` sets `method="t-test_overestim_var"` in 11 of 11 notebooks. 19 notebooks run ULM or MLM on an AnnData, and 0 of those notebooks scale that object. decoupler calls pass `use_raw=False` in 12 notebooks. Old decoupler defaulted to `use_raw=True` and read `.raw` silently.
- **Contrast mode.** The 11 notebooks that name the DE column score `stat`, for example `DeseqStats(...)` then `results_df[["stat"]].T` (`saezlab/corneto/saezlab_corneto_advanced_carnival`). Two more score log fold changes (`run_mlm(mat=logFCs, ...)`, `saezlab/marco_kidney_singlecell/..._ccrcc_ass1_final`). None scores `padj` or `-log10(p)`.

### Proposed gotchas (PROPOSED)

- **The model will want to pick PROGENy `top=500` or `5000` from a tutorial** when exporting the network. Do not. kai uses both (500 ×6, 5000 ×8), and the choice changes every score. Ask, and record the value with the network file.
- **The model will want to feed an h5ad whose `.raw` holds the log values and whose `X` is scaled.** Do not. Old decoupler read `.raw` silently, which is why 12 kai notebooks pass `use_raw=False` to a decoupler call. This script reads `X` only and refuses scaled data, so ask where the log-normalised matrix is.

### Proposed grill-question citations (PROPOSED)

- `--method`: "(kai: ULM 41 notebooks, MLM 12.)"
- `--tmin`: "(kai leaves decoupler's default in all but 2 notebooks.)"
- `--rank-method`: "(kai: `t-test_overestim_var` in 11 of 11.)" The existing gotcha already names it as a tutorial value, not a default.
- `--stat-col`: "(kai: `stat` in all 11 that name a column; log fold change in 2.)"

## Searched, not found

Zero hits are evidence too:

- Pearson `theta=` or `clip=`;
- Scrublet `random_state=` or `threshold=`, and `call_doublets(threshold=)`;
- `DeseqStats(alpha=)` and `independent_filter`;
- `sc.pp.log1p(base=)`;
- `ImplicitModificationWarning`;
- `recipe_pearson_residuals`;
- `dreamlet`.
