"""Trend test of MPL subcluster frequencies across the rural-to-urban gradient.

Purpose: per-subcluster OLS of arcsin-sqrt frequency on weighted polynomial group2
    contrasts, adjusted for age, sex and batch (HC3 SEs, BH FDR per term).
Inputs: data/MPL_data.parquet, data/MPL_metadata.parquet [gap: producer],
    src/weighted_contrasts.py
Outputs: data/MPL_cluster_median_markers.parquet, outputs/MPL_trend_test.parquet
Run: python src/mpl_trend_test.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import false_discovery_control

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from weighted_contrasts import poly_weighted  # noqa: E402

EXCLUDED_PATIENTS = ["LD413"]
FORMULA = "asin_freq ~ L + Q + age + C(sex) + C(batch)"

cells = pd.read_parquet(ROOT / "data/MPL_data.parquet")
metadata = pd.read_parquet(ROOT / "data/MPL_metadata.parquet")

# 1. Median marker intensity per subcluster
cluster_medians = (
    cells.loc[:, "CD45":"IL10"].groupby(cells["subclusterPG_merged"]).median()
)
cluster_medians.to_parquet(ROOT / "data/MPL_cluster_median_markers.parquet")

# 2. Subcluster frequency within its lineage, per patient
cells["lineage"] = cells["subclusterPG_merged"].str.split("_", n=1).str[0]
freq = (
    cells.groupby(["patient_id", "lineage", "subclusterPG_merged"], observed=True)
    .size()
    .rename("n_cells")
    .reset_index()
)
freq["freq"] = (
    freq["n_cells"]
    / freq.groupby(["patient_id", "lineage"])["n_cells"].transform("sum")
    * 100
)
freq["asin_freq"] = np.arcsin(np.sqrt(freq["freq"] / 100))

# Batch is recorded per cell; the model needs it per patient.
patient_batch = cells[["patient_id", "batch"]].drop_duplicates()
if patient_batch["patient_id"].duplicated().any():
    raise ValueError("some patients span more than one batch")
covariates = metadata[["patient_id", "group2", "age", "sex"]].merge(
    patient_batch, on="patient_id", how="left", validate="1:1"
)
freq = freq.merge(covariates, on="patient_id", how="left", validate="m:1")

n_before = freq["patient_id"].nunique()
freq = freq.dropna()
print(f"missing values: {n_before - freq['patient_id'].nunique()} patients dropped")
n_before = freq["patient_id"].nunique()
freq = freq[~freq["patient_id"].isin(EXCLUDED_PATIENTS)]
print(f"exclusion list: {n_before - freq['patient_id'].nunique()} patients dropped")

# 3. Subcluster x patient table; absent subcluster = frequency 0
freq_wide = freq.pivot_table(
    index="subclusterPG_merged",
    columns="patient_id",
    values="asin_freq",
    fill_value=0,
    observed=True,
)

# 4. Weighted polynomial contrasts, weights = patients per group2 level after filtering
patients = covariates[covariates["patient_id"].isin(freq_wide.columns)]
group2_levels = list(metadata["group2"].cat.categories)
group2_counts = patients["group2"].value_counts().reindex(group2_levels)
contrasts = poly_weighted(group2_levels, group2_counts.to_numpy())
patients = patients.merge(
    contrasts[["L", "Q"]], left_on="group2", right_index=True, validate="m:1"
)

model_data = (
    freq_wide.rename_axis(columns="patient_id")
    .stack()
    .rename("asin_freq")
    .reset_index()
    .merge(patients, on="patient_id", validate="m:1")
)


def fit_trend(cluster_data):
    fit = smf.ols(FORMULA, data=cluster_data).fit(cov_type="HC3")
    return pd.Series(
        {
            "estimate_L": fit.params["L"],
            "pval_L": fit.pvalues["L"],
            "pval_Q": fit.pvalues["Q"],
        }
    )


# 5. One row per subcluster; BH family = all subclusters, separately per term
trend = (
    model_data.groupby("subclusterPG_merged")
    .apply(fit_trend, include_groups=False)
    .rename_axis("cluster")
    .reset_index()
)
trend.insert(0, "lineage", trend["cluster"].str.split("_", n=1).str[0])
trend["fdr_L"] = false_discovery_control(trend["pval_L"], method="bh")
trend["fdr_Q"] = false_discovery_control(trend["pval_Q"], method="bh")

(ROOT / "outputs").mkdir(exist_ok=True)
trend.to_parquet(ROOT / "outputs/MPL_trend_test.parquet", index=False)
