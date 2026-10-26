"""Purpose: plot cell-type proportion by condition, one point per donor, lab palette.
Inputs:  data/proportions.tsv   [gap: provided with the task]
         data/palette.tsv       [gap: lab palette, keyed by label]
Outputs: <run_dir>/proportions_by-condition.{pdf,png,csv}, <run_dir>/session_info.txt
Run:     uv run --python 3.12 --with pandas,scipy,matplotlib python \
         scripts/03_plot-proportions.py results/run1
Notes:   [Stochastic] none. Unit: donor. Denominator: all of a donor's cells.
"""

# %% Setup
import argparse
import sys
from importlib.metadata import version
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

# Healthy first, so it sits left of disease in every cell-type group.
CONDITION_LEVELS = ["healthy", "disease"]
PANEL = "proportions_by-condition"
# Single-column journal width.
PANEL_WIDTH_MM = 89
PANEL_HEIGHT_MM = 60
# Dodge keeps the two conditions apart within a cell-type group.
DODGE = 0.15
REQUIRED_COLS = ["donor_id", "condition", "cell_type", "proportion"]

# Type 42 embeds editable text; 8 pt matches the R theme.
matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 8})

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
proportions = pd.read_csv(project_root / "data" / "proportions.tsv", sep="\t")
missing_cols = set(REQUIRED_COLS) - set(proportions.columns)
if missing_cols:
    raise KeyError(f"proportions.tsv lacks columns {sorted(missing_cols)}.")
palette_table = pd.read_csv(project_root / "data" / "palette.tsv", sep="\t")
palette = dict(zip(palette_table["label"], palette_table["hex"], strict=True))
# Colours by label, never by position: an unknown condition must stop the script.
unknown_labels = set(proportions["condition"]) - palette.keys()
if unknown_labels:
    raise KeyError(f"palette.tsv has no colour for {sorted(unknown_labels)}.")
proportions["condition"] = pd.Categorical(
    proportions["condition"], categories=CONDITION_LEVELS
)
if proportions["condition"].isna().any():
    raise ValueError(
        "proportions.tsv has conditions outside CONDITION_LEVELS; recode them at load."
    )

# %% Plot
cell_types = sorted(proportions["cell_type"].unique())
fig, ax = plt.subplots(
    figsize=(PANEL_WIDTH_MM / 25.4, PANEL_HEIGHT_MM / 25.4), layout="constrained"
)
for offset, condition in zip([-DODGE, DODGE], CONDITION_LEVELS, strict=True):
    rows = proportions[proportions["condition"] == condition]
    ax.scatter(
        rows["cell_type"].map(cell_types.index) + offset,
        rows["proportion"],
        s=8,
        color=palette[condition],
        label=condition,
    )
ax.set_xticks(range(len(cell_types)), cell_types, rotation=30, ha="right")
ax.set_ylabel("Proportion of the donor's cells")
ax.legend(frameon=False)
fig.supxlabel("One point per donor; denominator: all cells of that donor.", fontsize=7)

# %% Write outputs
fig.savefig(run_dir / f"{PANEL}.pdf", bbox_inches="tight")
fig.savefig(run_dir / f"{PANEL}.png", dpi=300, bbox_inches="tight")
plt.close(fig)
# Source data: the same frame the panel was drawn from.
proportions.to_csv(run_dir / f"{PANEL}.csv", index=False)

# %% Session info
session_info = [f"python {sys.version}"] + [
    f"{pkg} {version(pkg)}" for pkg in ["matplotlib", "pandas"]
]
(run_dir / "session_info.txt").write_text("\n".join(session_info) + "\n")
