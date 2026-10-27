"""Rural-to-urban trend test of MPL subcluster frequencies.

Purpose: per subcluster, OLS of arcsin-sqrt frequency (within lineage) on weighted
    linear and quadratic group2 contrasts, adjusted for age, sex and batch (HC3 SEs).
Inputs: data/MPL_data.parquet, data/MPL_metadata.parquet [gap: producers]
Outputs: data/MPL_cluster_median_markers.parquet, outputs/MPL_trend_test.parquet
Run: python src/mpl_trend_test.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

from weighted_contrasts import poly_weighted

ROOT = Path(__file__).resolve().parents[1]
out_dir = ROOT / "outputs"
EXCLUDED_PATIENTS = ["LD413"]
FORMULA = "asin_freq ~ L + Q + age + C(sex) + C(batch)"

cells = pd.read_parquet(ROOT / "data/MPL_data.parquet")
metadata = pd.read_parquet(ROOT / "data/MPL_metadata.parquet")

cluster_medians = (
    cells.loc[:, "CD45":"IL10"]
    .groupby(cells["subclusterPG_merged"], observed=True)
    .median()
    .reset_index()
)
cluster_medians.to_parquet(ROOT / "data/MPL_cluster_median_markers.parquet")

# validate="1:1" fails if a patient spans more than one batch
patients = (
    cells[["patient_id", "batch"]]
    .drop_duplicates()
    .merge(
        metadata[["patient_id", "group2", "age", "sex"]],
        on="patient_id",
        validate="1:1",
    )
)
n_start = len(patients)
patients = patients.dropna()
print(f"Dropped {n_start - len(patients)} patients with missing batch/group2/age/sex")
n_start = len(patients)
patients = patients[~patients["patient_id"].isin(EXCLUDED_PATIENTS)]
print(f"Dropped {n_start - len(patients)} excluded patients: {EXCLUDED_PATIENTS}")

# Subcluster x patient counts; reindexing adds 0 for patients without cells
kept_cells = cells[cells["patient_id"].isin(patients["patient_id"])]
cluster_counts = pd.crosstab(
    kept_cells["subclusterPG_merged"], kept_cells["patient_id"]
).reindex(columns=patients["patient_id"], fill_value=0)
lineage = cluster_counts.index.str.split("_", n=1).str[0]
lineage_counts = cluster_counts.groupby(lineage).transform("sum")
# A patient with no cells in a lineage gets frequency 0 (0/0 -> NaN -> 0), by design
freq_pct = (cluster_counts / lineage_counts * 100).fillna(0)
asin_freq = np.arcsin(np.sqrt(freq_pct / 100))

levels = metadata["group2"].cat.categories
level_counts = patients["group2"].value_counts().reindex(levels, fill_value=0)
contrasts = poly_weighted(levels, level_counts.to_numpy())

model_data = (
    asin_freq.rename_axis(index="cluster")
    .stack()
    .rename("asin_freq")
    .reset_index()
    .merge(patients, on="patient_id", validate="m:1")
    .merge(contrasts, left_on="group2", right_index=True, validate="m:1")
)

fits = {
    cluster: smf.ols(FORMULA, data=cluster_data).fit(cov_type="HC3")
    for cluster, cluster_data in model_data.groupby("cluster")
}
trend_test = pd.DataFrame(
    {
        "lineage": [cluster.split("_", 1)[0] for cluster in fits],
        "cluster": list(fits),
        "estimate_L": [fit.params["L"] for fit in fits.values()],
        "pval_L": [fit.pvalues["L"] for fit in fits.values()],
        "pval_Q": [fit.pvalues["Q"] for fit in fits.values()],
    }
)
# BH family: all subclusters, separately for each contrast term
trend_test["fdr_L"] = multipletests(trend_test["pval_L"], method="fdr_bh")[1]
trend_test["fdr_Q"] = multipletests(trend_test["pval_Q"], method="fdr_bh")[1]

out_dir.mkdir(exist_ok=True)
trend_test.to_parquet(out_dir / "MPL_trend_test.parquet", index=False)
