"""Purpose: plot cell-type proportion by condition, one point per donor, lab palette.
Inputs:  data/proportions.tsv   [gap: provided with the task]
         data/palette.tsv       [gap: lab palette, keyed by label]
Outputs: <run_dir>/proportions_by-condition.{pdf,png,csv}, <run_dir>/session_info.txt
Run:     uv run --python 3.12 --with pandas,scipy,matplotlib,plotnine python \
         scripts/03_plot-proportions.py results/run1
Notes:   [Stochastic] jitter, seeded. Unit: donor. Denominator: all of a donor's cells.
"""

import sys
from importlib.metadata import version
from pathlib import Path

import pandas as pd
from plotnine import (
    aes,
    geom_point,
    ggplot,
    labs,
    position_jitterdodge,
    scale_colour_manual,
)

from theme import save_panel

ROOT = Path(__file__).resolve().parents[1]
# Healthy first, so it sits left of disease in every cell-type group.
CONDITION_LEVELS = ["healthy", "disease"]
PANEL = "proportions_by-condition"
# Single-column journal width.
PANEL_WIDTH_MM = 89
PANEL_HEIGHT_MM = 60
SEED = 1


def main():
    # Setup -------------------------------------------------------------------
    run_dir = Path(sys.argv[1])
    if run_dir.exists():
        raise FileExistsError(f"{run_dir} exists; pick a new run directory.")
    run_dir.mkdir(parents=True)

    # Load inputs -------------------------------------------------------------
    proportions = pd.read_csv(ROOT / "data" / "proportions.tsv", sep="\t")
    palette_table = pd.read_csv(ROOT / "data" / "palette.tsv", sep="\t")
    palette = dict(zip(palette_table["label"], palette_table["hex"], strict=True))

    # Build the plot frame ----------------------------------------------------
    # Colours by label, never by position: an unknown condition stops the script.
    unknown = set(proportions["condition"]) - palette.keys()
    if unknown:
        raise KeyError(f"palette.tsv has no colour for {sorted(unknown)}")
    frame = proportions.assign(
        condition=pd.Categorical(proportions["condition"], CONDITION_LEVELS)
    )

    # Plot --------------------------------------------------------------------
    plot = (
        ggplot(frame, aes("cell_type", "proportion", colour="condition"))
        + geom_point(
            size=1,
            position=position_jitterdodge(
                jitter_width=0.1, dodge_width=0.5, random_state=SEED
            ),
        )
        + scale_colour_manual(values={c: palette[c] for c in CONDITION_LEVELS})
        + labs(x="cell type", y="proportion of the donor's cells")
    )

    # Write outputs -----------------------------------------------------------
    save_panel(plot, frame, run_dir, PANEL, PANEL_WIDTH_MM, PANEL_HEIGHT_MM)
    versions = [f"python {sys.version}"] + [
        f"{pkg} {version(pkg)}" for pkg in ["pandas", "plotnine"]
    ]
    (run_dir / "session_info.txt").write_text("\n".join(versions) + "\n")


if __name__ == "__main__":
    main()
