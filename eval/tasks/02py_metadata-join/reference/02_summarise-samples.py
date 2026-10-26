"""Purpose: attach condition and batch to each sample and report its library size.
Inputs:  data/counts.tsv     [gap: provided with the task]
         data/metadata.tsv   [gap: sample sheet from the sequencing core]
Outputs: <run_dir>/sample_summary.tsv, <run_dir>/session_info.txt
Run:     uv run --python 3.12 --with pandas,scipy,matplotlib python \
         scripts/02_summarise-samples.py results/run1
Notes:   [Stochastic] none. Count columns and sheet rows are in different orders,
         so samples are joined by sample_id, never by position.
"""

# %% Setup
import argparse
import sys
from importlib.metadata import version
from pathlib import Path

import pandas as pd

REQUIRED_METADATA_COLS = ["sample_id", "condition", "batch"]

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
metadata = pd.read_csv(project_root / "data" / "metadata.tsv", sep="\t")
missing_cols = set(REQUIRED_METADATA_COLS) - set(metadata.columns)
if missing_cols:
    raise KeyError(f"metadata.tsv lacks columns {sorted(missing_cols)}.")

# %% Join metadata by sample_id
library_sizes = (
    counts.sum().rename("library_size").rename_axis("sample_id").reset_index()
)
# validate: a duplicated sample_id on either side stops the join.
sample_summary = library_sizes.merge(
    metadata, on="sample_id", how="left", validate="one_to_one"
)[["sample_id", "condition", "batch", "library_size"]]
unmatched = sample_summary.loc[sample_summary["condition"].isna(), "sample_id"]
if not unmatched.empty:
    raise ValueError(
        f"no metadata row for samples {unmatched.tolist()}; fix the sample sheet."
    )

# %% Write outputs
sample_summary.to_csv(run_dir / "sample_summary.tsv", sep="\t", index=False)

# %% Session info
session_info = [f"python {sys.version}"] + [
    f"{pkg} {version(pkg)}" for pkg in ["pandas"]
]
(run_dir / "session_info.txt").write_text("\n".join(session_info) + "\n")
