"""Purpose: build the frozen inputs and truth of task 15 once; the task never fetches.
Inputs:  decoupler 2.3.1 dc.ds.toy(seed=42): synthetic, generated offline (BSD-3-Clause)
Outputs: inputs/data/expression.tsv  30 cells x 20 genes, cell_id then one column per gene
         inputs/data/cells.tsv       cell_id, group (A/B), sample
         inputs/data/network.tsv     source, target, mor (the toy net's weight, renamed)
         truth/group_means.tsv       mean ULM score per group and TF, scored with the signs
Run:     uv run --python 3.12 --with-requirements eval/envs/python-tasks.txt \
         python -I eval/tasks/15_regulator-sign/make-inputs.py
Notes:   The planted trap is the column name: decoupler v2 reads only `weight` and scores
         every edge +1 when it is absent, so repressor T2 turns positive in group B.
"""

from pathlib import Path

import decoupler as dc

TASK = Path(__file__).resolve().parent
DATA = TASK / "inputs" / "data"

adata, net = dc.ds.toy(seed=42)
expression = adata.to_df().round(4).rename_axis("cell_id")
cells = adata.obs[["group", "sample"]].rename_axis("cell_id")

scores, _ = dc.mt.ulm(data=expression, net=net, tmin=0)
truth = (
    scores.groupby(cells["group"].to_numpy())
    .mean()
    .rename_axis("group")
    .reset_index()
    .melt(id_vars="group", var_name="source", value_name="mean_score")
)

expression.to_csv(DATA / "expression.tsv", sep="\t")
cells.to_csv(DATA / "cells.tsv", sep="\t")
net.rename(columns={"weight": "mor"}).to_csv(
    DATA / "network.tsv", sep="\t", index=False
)
truth.to_csv(TASK / "truth" / "group_means.tsv", sep="\t", index=False)
print(truth.pivot(index="group", columns="source", values="mean_score").round(2))
