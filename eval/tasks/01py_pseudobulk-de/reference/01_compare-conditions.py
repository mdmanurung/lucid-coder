"""Purpose: test each gene for a disease-vs-healthy difference; the donor is the unit.
Inputs:  data/counts.tsv   [gap: provided with the task]
         data/cells.tsv    [gap: provided with the task]
Outputs: <run_dir>/de.tsv, <run_dir>/session_info.txt
Run:     uv run --python 3.12 --with pandas,scipy,matplotlib python \
         scripts/01_compare-conditions.py results/run1
Notes:   [Stochastic] none. Inference unit: donor. Cells from one donor are not
         independent, so counts are summed per donor (pseudobulk) before testing
         3 vs 3 donors.
"""

# %% Setup
import argparse
import sys
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, ttest_ind

# Reference = healthy, so log2_fc reads as disease minus healthy.
CONDITION_LEVELS = ["healthy", "disease"]
# Pseudocount 1 on CPM keeps log2 finite for a gene with zero counts in a donor.
CPM_PSEUDOCOUNT = 1
REQUIRED_CELL_COLS = ["cell_id", "donor_id", "condition"]

parser = argparse.ArgumentParser()
parser.add_argument("run_dir", type=Path, help="new directory under results/")
run_dir = parser.parse_args().run_dir
if run_dir.exists():
    raise FileExistsError(
        f"run_dir {run_dir} already exists; pick a new name so nothing is overwritten."
    )
run_dir.mkdir(parents=True)

# %% Load inputs
project_root = Path(__file__).resolve().parents[1]
counts = pd.read_csv(project_root / "data" / "counts.tsv", sep="\t", index_col="gene")
cells = pd.read_csv(project_root / "data" / "cells.tsv", sep="\t")
missing_cols = set(REQUIRED_CELL_COLS) - set(cells.columns)
if missing_cols:
    raise KeyError(f"cells.tsv lacks columns {sorted(missing_cols)}.")
# Cell order in cells.tsv is not guaranteed to match the count columns.
cells = cells.set_index("cell_id").reindex(counts.columns)
if cells["donor_id"].isna().any():
    raise ValueError("some count columns have no row in cells.tsv; fix the cell table.")

# %% Sum counts per donor
donors = cells.drop_duplicates().reset_index(drop=True)
if donors["donor_id"].duplicated().any():
    raise ValueError("a donor has cells in more than one condition; fix cells.tsv.")
donors["condition"] = pd.Categorical(donors["condition"], categories=CONDITION_LEVELS)
if donors["condition"].isna().any():
    raise ValueError(
        "cells.tsv has conditions outside CONDITION_LEVELS; recode them at load."
    )
donor_counts = counts.T.groupby(cells["donor_id"].to_numpy()).sum().T
donor_log_cpm = np.log2(donor_counts / donor_counts.sum() * 1e6 + CPM_PSEUDOCOUNT)
condition_of = donors.set_index("donor_id")["condition"]
healthy_donors = condition_of.index[condition_of == "healthy"]
disease_donors = condition_of.index[condition_of == "disease"]


# %% Test each gene
# Pooled-variance t-test on donor log-CPM: one value per donor, so df = donors - 2.
def fit_gene(log_cpm):
    healthy = log_cpm[healthy_donors]
    disease = log_cpm[disease_donors]
    test = ttest_ind(disease, healthy, equal_var=True)
    return {"log2_fc": disease.mean() - healthy.mean(), "p_value": test.pvalue}


de_results = pd.DataFrame(
    [{"gene": gene, **fit_gene(log_cpm)} for gene, log_cpm in donor_log_cpm.iterrows()]
)
# Multiple-testing family: every gene in counts.tsv, adjusted together.
de_results["q_value"] = false_discovery_control(de_results["p_value"], method="bh")
de_results = de_results.sort_values("p_value")

# %% Write outputs
de_results.to_csv(run_dir / "de.tsv", sep="\t", index=False)

# %% Session info
session_info = [f"python {sys.version}"] + [
    f"{pkg} {version(pkg)}" for pkg in ["numpy", "pandas", "scipy"]
]
(run_dir / "session_info.txt").write_text("\n".join(session_info) + "\n")
