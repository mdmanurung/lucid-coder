"""
Purpose: Test MPL-stimulated subcluster frequencies for a rural-to-urban trend.
Inputs: data/MPL_data.parquet [gap], data/MPL_metadata.parquet [gap]
Outputs: data/MPL_cluster_median_markers.parquet, outputs/MPL_trend_test.parquet
Run: python src/mpl_trend_test.py
Note: a patient with no cells in a whole lineage has an undefined frequency (0/0);
it is NA and left out of that subcluster's fit, not set to 0.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

from weighted_contrasts import poly_weighted

ROOT = Path(__file__).resolve().parents[1]
EXCLUDED_PATIENTS = ["LD413"]
FORMULA = "asin_freq ~ L + Q + age + C(sex) + C(batch)"

cells = pd.read_parquet(ROOT / "data/MPL_data.parquet")
metadata = pd.read_parquet(ROOT / "data/MPL_metadata.parquet")

# Medians over all cells, before any patient filtering
markers = cells.loc[:, "CD45":"IL10"].columns
cells.groupby("subclusterPG_merged", observed=True)[markers].median().to_parquet(
    ROOT / "data/MPL_cluster_median_markers.parquet"
)

# validate="1:1" also fails if a patient spans more than one batch
patients = (
    cells[["patient_id", "batch"]]
    .drop_duplicates()
    .merge(
        metadata[["patient_id", "group2", "age", "sex"]],
        on="patient_id",
        how="left",
        validate="1:1",
    )
)
n_patients = len(patients)
patients = patients.dropna()
print(f"Removed {n_patients - len(patients)} patients with missing covariates")
n_patients = len(patients)
patients = patients[~patients["patient_id"].isin(EXCLUDED_PATIENTS)]
print(f"Removed {n_patients - len(patients)} excluded patients: {EXCLUDED_PATIENTS}")

# Category order of group2 is the trend direction (rural to urban)
level_counts = patients["group2"].value_counts(sort=False)
if (level_counts == 0).any():
    raise ValueError(f"group2 level without patients: {level_counts.to_dict()}")
contrasts = poly_weighted(level_counts.index, level_counts.to_numpy())
design = patients.merge(
    contrasts, left_on="group2", right_index=True, how="left", validate="m:1"
).set_index("patient_id")

# Subcluster x patient cell counts; 0 where a patient has no cells in a subcluster
counts = pd.crosstab(cells["subclusterPG_merged"], cells["patient_id"])
counts = counts.loc[:, design.index]
counts = counts[counts.sum(axis=1) > 0]
lineage = counts.index.astype(str).str.split("_", n=1).str[0]
asin_freq = np.arcsin(np.sqrt(counts / counts.groupby(lineage).transform("sum")))

fits = {
    subcluster: smf.ols(
        FORMULA, data=design.assign(asin_freq=asin_freq.loc[subcluster])
    ).fit(cov_type="HC3")
    for subcluster in asin_freq.index
}
trend = pd.DataFrame(
    {
        "lineage": lineage,
        "cluster": asin_freq.index.astype(str),
        "estimate_L": [fit.params["L"] for fit in fits.values()],
        "pval_L": [fit.pvalues["L"] for fit in fits.values()],
        "pval_Q": [fit.pvalues["Q"] for fit in fits.values()],
    }
)

# BH family: all subclusters, per term. HC3 gives NaN p-values when a point has
# leverage 1; multipletests mis-sorts NaN, so those stay NaN and are left out.
for term in ["L", "Q"]:
    tested = trend[f"pval_{term}"].notna()
    trend[f"fdr_{term}"] = np.nan
    trend.loc[tested, f"fdr_{term}"] = multipletests(
        trend.loc[tested, f"pval_{term}"], method="fdr_bh"
    )[1]

(ROOT / "outputs").mkdir(exist_ok=True)
trend.to_parquet(ROOT / "outputs/MPL_trend_test.parquet", index=False)
