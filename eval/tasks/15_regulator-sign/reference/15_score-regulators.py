"""Purpose: score TF activity per cell from a signed regulatory network.
Inputs:  data/expression.tsv  [gap: log-normalised expression, provided with the task]
         data/cells.tsv       [gap: cell groups, provided with the task]
         data/network.tsv     [gap: signed TF-target network, provided with the task]
Outputs: <run_dir>/tf_scores.tsv, <run_dir>/session_info.txt
Run:     uv run --python 3.12 --with decoupler==2.3.1,pandas,matplotlib python \
         scripts/15_score-regulators.py results/run1
Notes:   [Stochastic] none. decoupler reads edge signs only from a `weight` column and
         scores every edge +1 without one, so `mor` is renamed to `weight`.
"""

# %% Setup
import argparse
import sys
from importlib.metadata import version
from pathlib import Path

import decoupler as dc
import pandas as pd

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

# %% Write outputs
tf_scores.to_csv(run_dir / "tf_scores.tsv", sep="\t", index=False)

# %% Session info
session_info = [f"python {sys.version}"] + [
    f"{pkg} {version(pkg)}" for pkg in ["decoupler", "pandas"]
]
(run_dir / "session_info.txt").write_text("\n".join(session_info) + "\n")
