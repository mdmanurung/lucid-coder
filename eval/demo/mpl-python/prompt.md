Write a Python script `src/mpl_trend_test.py` for the analysis below. Do not run anything; the data are not here. Use pandas, numpy, statsmodels and scipy.

Data:
- `data/MPL_data.parquet`: one row per cell from MPL-stimulated samples. Columns: `patient_id`, `batch`, `group3`, `subclusterPG_merged` (subcluster name, formatted `<lineage>_<annotation>`, e.g. `CD4T_Th1`), and marker intensity columns from `CD45` through `IL10` (contiguous).
- `data/MPL_metadata.parquet`: one row per patient. Columns: `patient_id`, `group2` (categorical, three ordered levels, rural to urban), `age`, `sex`.
- `src/weighted_contrasts.py` defines `poly_weighted(levels, counts)`, which returns polynomial contrasts weighted by group size as a DataFrame indexed by level with columns `L`, `Q`.

Analysis:
1. Median intensity of every marker per subcluster. Save to `data/MPL_cluster_median_markers.parquet`.
2. Per patient: count cells per subcluster and per lineage. Frequency of a subcluster = its cells / cells in its lineage x 100. Transform with arcsin(sqrt(proportion)). Attach `group2`, `age`, `sex` from the metadata. Drop rows with missing values and exclude patient `LD413`.
3. Subcluster x patient table of the transformed frequencies, 0 where a patient has no cells in a subcluster.
4. For each subcluster, fit OLS of the transformed frequency on the weighted linear and quadratic `group2` contrasts (from `poly_weighted`, with counts = patients per level after filtering) plus `age`, `sex` and `batch`. Use robust (HC3) standard errors.
5. One table with lineage, cluster, the linear-term estimate and p-value (`pval_L`) and the quadratic-term p-value (`pval_Q`), with Benjamini-Hochberg FDR across subclusters for each term (`fdr_L`, `fdr_Q`). Save to `outputs/MPL_trend_test.parquet`.
