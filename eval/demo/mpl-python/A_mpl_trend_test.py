"""Trend test of MPL subcluster frequencies across the rural-urban gradient.

For each subcluster, the arcsin-sqrt transformed frequency (within lineage) is
regressed on weighted linear and quadratic group2 contrasts, adjusting for age,
sex and batch, with HC3 robust standard errors.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

from weighted_contrasts import poly_weighted

ROOT = Path(__file__).resolve().parent.parent
DATA_PATH = ROOT / "data" / "MPL_data.parquet"
META_PATH = ROOT / "data" / "MPL_metadata.parquet"
MEDIAN_OUT = ROOT / "data" / "MPL_cluster_median_markers.parquet"
TREND_OUT = ROOT / "outputs" / "MPL_trend_test.parquet"

CLUSTER_COL = "subclusterPG_merged"
EXCLUDED_PATIENTS = ["LD413"]
FORMULA = "transformed ~ L + Q + age + C(sex) + C(batch)"


def lineage_of(cluster: pd.Series) -> pd.Series:
    """Lineage is the part of `<lineage>_<annotation>` before the first underscore."""
    return cluster.astype(str).str.split("_", n=1).str[0]


def cluster_medians(cells: pd.DataFrame) -> pd.DataFrame:
    marker_cols = cells.loc[:, "CD45":"IL10"].columns
    medians = cells.groupby(CLUSTER_COL, observed=True)[marker_cols].median()
    return medians.rename_axis("cluster").reset_index()


def patient_covariates(cells: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """One row per analysable patient: group2, age, sex, batch."""
    n_batches = cells.groupby("patient_id", observed=True)["batch"].nunique()
    mixed = n_batches[n_batches > 1]
    if not mixed.empty:
        raise ValueError(f"patients in multiple batches: {list(mixed.index)}")
    batch = (
        cells.groupby("patient_id", observed=True)["batch"].first().reset_index()
    )
    batch["patient_id"] = batch["patient_id"].astype(str)

    meta = meta[["patient_id", "group2", "age", "sex"]].copy()
    meta["patient_id"] = meta["patient_id"].astype(str)
    if not (isinstance(meta["group2"].dtype, pd.CategoricalDtype)
            and meta["group2"].cat.ordered):
        raise TypeError("group2 must be an ordered categorical (rural to urban)")

    covars = batch.merge(meta, on="patient_id", how="inner").dropna()
    covars = covars[~covars["patient_id"].isin(EXCLUDED_PATIENTS)]
    covars["group2"] = covars["group2"].cat.remove_unused_categories()
    return covars.reset_index(drop=True)


def transformed_frequencies(cells: pd.DataFrame, patients: pd.Index) -> pd.DataFrame:
    """Long table of arcsin-sqrt within-lineage subcluster frequencies."""
    counts = (
        cells.assign(patient_id=cells["patient_id"].astype(str))
        .groupby(["patient_id", CLUSTER_COL], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .rename(columns={CLUSTER_COL: "cluster"})
    )
    counts["cluster"] = counts["cluster"].astype(str)
    counts["lineage"] = lineage_of(counts["cluster"])
    counts["n_lineage"] = counts.groupby(["patient_id", "lineage"])[
        "n_cells"
    ].transform("sum")
    prop = counts["n_cells"] / counts["n_lineage"]
    counts["freq"] = prop * 100
    counts["transformed"] = np.arcsin(np.sqrt(prop))
    return counts[counts["patient_id"].isin(patients)]


def frequency_table(freqs: pd.DataFrame, patients: pd.Index) -> pd.DataFrame:
    """Subcluster x patient table; 0 where a patient has no cells in a subcluster."""
    table = freqs.pivot(index="cluster", columns="patient_id", values="transformed")
    return table.reindex(columns=patients).fillna(0)


def fit_trends(table: pd.DataFrame, covars: pd.DataFrame) -> pd.DataFrame:
    levels = covars["group2"].cat.categories
    n_per_level = covars["group2"].value_counts(sort=False).reindex(levels)
    contr = poly_weighted(list(levels), n_per_level.to_numpy())

    design = covars.copy()
    group = design["group2"].astype(object)
    design["L"] = group.map(contr["L"]).astype(float)
    design["Q"] = group.map(contr["Q"]).astype(float)

    long = (
        table.rename_axis(index="cluster", columns=None)
        .reset_index()
        .melt(id_vars="cluster", var_name="patient_id", value_name="transformed")
        .merge(design, on="patient_id", how="left")
    )

    rows = []
    for cluster, sub in long.groupby("cluster", sort=True):
        fit = smf.ols(FORMULA, data=sub).fit(cov_type="HC3")
        rows.append(
            {
                "cluster": cluster,
                "estimate_L": fit.params["L"],
                "pval_L": fit.pvalues["L"],
                "pval_Q": fit.pvalues["Q"],
            }
        )
    res = pd.DataFrame(rows)
    res.insert(0, "lineage", lineage_of(res["cluster"]))
    return res


def bh_fdr(pvals: pd.Series) -> pd.Series:
    """Benjamini-Hochberg adjustment, leaving NaN p-values as NaN."""
    out = pd.Series(np.nan, index=pvals.index)
    ok = pvals.notna()
    if ok.any():
        out[ok] = multipletests(pvals[ok], method="fdr_bh")[1]
    return out


def main() -> None:
    cells = pd.read_parquet(DATA_PATH)
    meta = pd.read_parquet(META_PATH)

    # 1. median marker intensity per subcluster
    MEDIAN_OUT.parent.mkdir(parents=True, exist_ok=True)
    cluster_medians(cells).to_parquet(MEDIAN_OUT, index=False)

    # 2. per-patient within-lineage frequencies, arcsin-sqrt transformed
    covars = patient_covariates(cells, meta)
    patients = pd.Index(covars["patient_id"], name="patient_id")
    freqs = transformed_frequencies(cells, patients)

    # 3. subcluster x patient table, zero-filled
    table = frequency_table(freqs, patients)

    # 4. per-subcluster OLS with weighted polynomial group2 contrasts, HC3 SEs
    res = fit_trends(table, covars)

    # 5. BH FDR across subclusters for each term
    res["fdr_L"] = bh_fdr(res["pval_L"])
    res["fdr_Q"] = bh_fdr(res["pval_Q"])
    res = res[
        ["lineage", "cluster", "estimate_L", "pval_L", "fdr_L", "pval_Q", "fdr_Q"]
    ]

    TREND_OUT.parent.mkdir(parents=True, exist_ok=True)
    res.to_parquet(TREND_OUT, index=False)


if __name__ == "__main__":
    main()
