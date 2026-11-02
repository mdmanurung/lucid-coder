"""Purpose: score TF activity per cell from a signed regulatory network and test
         which TFs differ between groups A and B.
Inputs:  data/expression.tsv  [gap: log-normalised expression, provided with the task]
         data/cells.tsv       [gap: cell groups, provided with the task]
         data/network.tsv     [gap: signed TF-target network, provided with the task]
Outputs: <run_dir>/tf_scores.tsv        ulm score per cell and TF, long format
         <run_dir>/tf_group_test.tsv    per TF: difference (B - A) and paired t p-value
         <run_dir>/session_info.txt
Run:     uv run --python 3.12 --with decoupler==2.3.1,pandas,matplotlib python \
         scripts/15_score-regulators.py results/run1
Notes:   [Stochastic] none. decoupler reads edge signs only from a `weight` column
         and scores every edge +1 without one, so `mor` is renamed to `weight`. Every
         sample holds cells of both groups, so samples are the replicates: the group
         test is a paired t-test on the per-sample group means (n = 3), not over cells.
"""

# %% Setup
import argparse
import sys
from importlib.metadata import version
from pathlib import Path

import decoupler as dc
import pandas as pd
from scipy import stats

# Each TF has 3-4 targets, below decoupler's default tmin of 5.
TMIN = 3

parser = argparse.ArgumentParser()
parser.add_argument("run_dir", type=Path, help="new directory under results/")
run_dir = parser.parse_args().run_dir
if run_dir.exists():
    raise FileExistsError(f"run_dir {run_dir} already exists; pick a new name.")
run_dir.mkdir(parents=True)

# %% Load inputs
project_root = Path(__file__).resolve().parents[1]
expression = pd.read_csv(
    project_root / "data" / "expression.tsv", sep="\t", index_col="cell_id"
)
cells = pd.read_csv(project_root / "data" / "cells.tsv", sep="\t", index_col="cell_id")
network = pd.read_csv(project_root / "data" / "network.tsv", sep="\t").rename(
    columns={"mor": "weight"}
)

# %% Score TFs per cell
scores, _ = dc.mt.ulm(data=expression, net=network, tmin=TMIN)
tf_scores = (
    scores.rename_axis("cell_id")
    .reset_index()
    .melt(id_vars="cell_id", var_name="source", value_name="score")
)

# %% Test groups A vs B with samples as replicates
sample_means = (
    tf_scores.merge(cells, left_on="cell_id", right_index=True, validate="m:1")
    .groupby(["source", "sample", "group"])["score"]
    .mean()
    .unstack("group")
)
tf_group_test = (
    sample_means.groupby(level="source")
    .apply(
        lambda m: pd.Series(
            {
                "difference": (m["B"] - m["A"]).mean(),
                "p_value": stats.ttest_rel(m["B"], m["A"]).pvalue,
            }
        )
    )
    .reset_index()
)

# %% Write outputs
tf_scores.to_csv(run_dir / "tf_scores.tsv", sep="\t", index=False)
tf_group_test.to_csv(run_dir / "tf_group_test.tsv", sep="\t", index=False)

# %% Session info
session_info = [f"python {sys.version}"] + [
    f"{pkg} {version(pkg)}" for pkg in ["decoupler", "pandas", "scipy"]
]
(run_dir / "session_info.txt").write_text("\n".join(session_info) + "\n")
