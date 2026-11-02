"""Purpose: build the frozen inputs and truth of task 15b once; the task never fetches.
Inputs:  decoupler 2.3.1 dc.ds.toy(seed=42): synthetic, generated offline (BSD-3-Clause)
Outputs: inputs/data/expression.tsv  30 cells x 20 genes, cell_id then one column per gene
         inputs/data/cells.tsv       cell_id, group (A/B), sample
         inputs/data/network.tsv     source, target, mor (the toy net's weight, renamed)
         truth/group_means.tsv       mean ULM score per group and TF, scored with the signs
         truth/tf_group_test.tsv     source, difference (B - A), p_value: paired t on the
                                     per-sample group means, n = 3 samples
Run:     uv run --python 3.12 --with-requirements eval/envs/python-tasks.txt \
         python -I eval/tasks/15b_regulator-sign-groups/make-inputs.py
Notes:   The planted trap is the column name: decoupler v2 reads only `weight` and scores
         every edge +1 when it is absent, so repressor T2 turns positive in group B.
         Every sample holds cells of both groups, so the group test pairs the 3 sample
         means; cells are not replicates.
"""

from pathlib import Path

import decoupler as dc
import pandas as pd
from scipy import stats

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

# Paired t-test per TF on the per-sample group means.
sample_means = (
    scores.join(cells)
    .melt(id_vars=["group", "sample"], var_name="source", value_name="score")
    .groupby(["source", "sample", "group"], observed=True)["score"]
    .mean()
    .unstack("group")
)
group_test = sample_means.groupby(level="source").apply(
    lambda m: pd.Series(
        {
            "difference": (m["B"] - m["A"]).mean(),
            "p_value": stats.ttest_rel(m["B"], m["A"]).pvalue,
        }
    )
)

expression.to_csv(DATA / "expression.tsv", sep="\t")
cells.to_csv(DATA / "cells.tsv", sep="\t")
net.rename(columns={"weight": "mor"}).to_csv(
    DATA / "network.tsv", sep="\t", index=False
)
truth.to_csv(TASK / "truth" / "group_means.tsv", sep="\t", index=False)
group_test.to_csv(TASK / "truth" / "tf_group_test.tsv", sep="\t")
print(truth.pivot(index="group", columns="source", values="mean_score").round(2))
print(group_test)
