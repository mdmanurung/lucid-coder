"""Trend test of MPL subcluster frequencies across group2 (rural to urban).

Purpose: per-subcluster OLS of arcsine-sqrt lineage frequency on weighted linear and
    quadratic group2 contrasts, adjusted for age, sex and batch (HC3 SEs).
Inputs: data/MPL_data.parquet [gap], data/MPL_metadata.parquet [gap],
    src/weighted_contrasts.py
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
data_dir = ROOT / "data"
out_dir = ROOT / "outputs"
excluded_patients = ["LD413"]
formula = "asin_freq ~ L + Q + age + C(sex) + C(batch)"

cells = pd.read_parquet(data_dir / "MPL_data.parquet")
metadata = pd.read_parquet(data_dir / "MPL_metadata.parquet")

marker_cols = cells.loc[:, "CD45":"IL10"].columns
cluster_medians = cells.groupby("subclusterPG_merged", observed=True)[
    marker_cols
].median()
cluster_medians.to_parquet(data_dir / "MPL_cluster_median_markers.parquet")

# Batch is recorded per cell; validate="1:1" fails if a patient spans two batches.
patient_batch = cells[["patient_id", "batch"]].drop_duplicates()
patients = metadata[["patient_id", "group2", "age", "sex"]].merge(
    patient_batch, on="patient_id", validate="1:1"
)
n_before = len(patients)
patients = patients.dropna()
print(f"Dropped {n_before - len(patients)} patients with missing group2/age/sex/batch")
n_before = len(patients)
patients = patients[~patients["patient_id"].isin(excluded_patients)]
print(f"Dropped {n_before - len(patients)} excluded patients: {excluded_patients}")

group2_levels = metadata["group2"].cat.categories
patients_per_level = patients["group2"].value_counts().reindex(group2_levels)
contrasts = poly_weighted(group2_levels, patients_per_level.to_numpy())
patients = patients.merge(
    contrasts, left_on="group2", right_index=True, validate="m:1"
).set_index("patient_id")

# Full patient x subcluster grid so absent subclusters count as 0 cells; a patient
# with no cells in the whole lineage gets NaN (0/0) and drops out of that fit.
cluster_counts = (
    cells[cells["patient_id"].isin(patients.index)]
    .groupby(["patient_id", "subclusterPG_merged"], observed=True)
    .size()
    .unstack(fill_value=0)
    .stack()
    .rename("n_cells")
    .reset_index()
    .rename(columns={"subclusterPG_merged": "cluster"})
)
cluster_counts["lineage"] = cluster_counts["cluster"].str.split("_", n=1).str[0]
lineage_cells = cluster_counts.groupby(["patient_id", "lineage"])["n_cells"].transform(
    "sum"
)
cluster_counts["freq"] = cluster_counts["n_cells"] / lineage_cells * 100
cluster_counts["asin_freq"] = np.arcsin(np.sqrt(cluster_counts["freq"] / 100))

freq_table = cluster_counts.pivot(
    index="cluster", columns="patient_id", values="asin_freq"
)

trend_rows = []
for cluster, asin_freq in freq_table.iterrows():
    fit = smf.ols(formula, data=patients.assign(asin_freq=asin_freq)).fit(
        cov_type="HC3"
    )
    trend_rows.append(
        {
            "lineage": cluster.split("_", 1)[0],
            "cluster": cluster,
            "estimate_L": fit.params["L"],
            "pval_L": fit.pvalues["L"],
            "pval_Q": fit.pvalues["Q"],
        }
    )
trend = pd.DataFrame(trend_rows)

# BH family: all subclusters, corrected separately for the L and Q terms.
trend["fdr_L"] = multipletests(trend["pval_L"], method="fdr_bh")[1]
trend["fdr_Q"] = multipletests(trend["pval_Q"], method="fdr_bh")[1]

out_dir.mkdir(exist_ok=True)
trend.to_parquet(out_dir / "MPL_trend_test.parquet", index=False)
