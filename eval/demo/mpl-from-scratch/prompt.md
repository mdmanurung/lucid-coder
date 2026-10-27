Write an R script `src/MPL_trend_test.R` for the analysis below. Do not run anything; the data are not here.

Data:
- `data/MPL_data.qs` (qs format): one row per cell from MPL-stimulated samples. Columns: `patient_id`, `batch`, `group3`, `subclusterPG_merged` (subcluster name, formatted `<lineage>_<annotation>`, e.g. `CD4T_Th1`), and marker intensity columns from `CD45` through `IL10` (contiguous).
- `data/MPL_metadata.rds`: one row per patient. Columns: `patient_id`, `group2` (factor, three levels ordered rural to urban), `age`, `sex`.
- `src/weightedcontrasts.R` defines `contr.poly.weighted(f)`, which returns polynomial contrasts weighted by group size for a factor `f` (columns named `.L`, `.Q`, ...).

Analysis:
1. Median intensity of every marker per subcluster. Save to `data/MPL_cluster_median_markers.qs`.
2. Per patient: count cells per subcluster and per lineage. Frequency of a subcluster = its cells / cells in its lineage x 100. Transform with asin(sqrt(proportion)). Attach `group2`, `age`, `sex` from the metadata. Drop rows with missing values and exclude patient `LD413`.
3. Subcluster x patient matrix of the transformed frequencies, 0 where a patient has no cells in a subcluster, including when the patient has no cells in that whole lineage.
4. limma: `lmFit(method = "robust", maxit = 100)` on `~ group2 + age + sex + batch`, with `contrasts(group2) <- contr.poly.weighted(group2)` so the group terms are the linear and quadratic trends; then `eBayes(robust = TRUE)`.
5. `topTable` for the linear and the quadratic term (all subclusters). Rename `P.Value`/`adj.P.Val` to `pval_L`/`fdr_L` and `pval_Q`/`fdr_Q`. One table with `lineage` and `cluster`, where `cluster` is the full subcluster name (e.g. `CD4T_Th1`, as downstream figure scripts join on it), plus the linear-term statistics and `pval_Q`, `fdr_Q`. Save to `outputs/MPL_limma_trend_test.rds`; create `outputs/` if it does not exist.
