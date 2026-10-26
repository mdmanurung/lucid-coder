# Naive answer: makes the figure but forgets the plotted-data CSV.
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)

df = pd.read_csv("data/proportions.tsv", sep="\t")
pal = pd.read_csv("data/palette.tsv", sep="\t")
cols = dict(zip(pal["label"], pal["hex"]))

fig, ax = plt.subplots(figsize=(6, 4))
for cond, rows in df.groupby("condition"):
    ax.scatter(rows["cell_type"], rows["proportion"], color=cols[cond], label=cond)
ax.legend()
fig.savefig(out / "proportions_by-condition.pdf")
fig.savefig(out / "proportions_by-condition.png", dpi=150)
