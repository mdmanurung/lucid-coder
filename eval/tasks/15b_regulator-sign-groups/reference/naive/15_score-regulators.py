# Naive answer: signs reach the scores (`mor` renamed to `weight`), but the A vs B test
# treats the 30 cells as replicates with a Welch t-test; the 3 samples are the replicates.
import sys
from pathlib import Path

import decoupler as dc
import pandas as pd
from scipy import stats

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)

expression = pd.read_csv("data/expression.tsv", sep="\t", index_col="cell_id")
cells = pd.read_csv("data/cells.tsv", sep="\t", index_col="cell_id")
network = pd.read_csv("data/network.tsv", sep="\t").rename(columns={"mor": "weight"})

scores, _ = dc.mt.ulm(data=expression, net=network, tmin=0)
res = (
    scores.rename_axis("cell_id")
    .reset_index()
    .melt(id_vars="cell_id", var_name="source", value_name="score")
)
res.to_csv(out / "tf_scores.tsv", sep="\t", index=False)

is_b = (cells.loc[scores.index, "group"] == "B").to_numpy()
test = pd.DataFrame(
    {
        "source": scores.columns,
        "difference": scores[is_b].mean() - scores[~is_b].mean(),
        "p_value": stats.ttest_ind(scores[is_b], scores[~is_b], equal_var=False).pvalue,
    }
)
test.to_csv(out / "tf_group_test.tsv", sep="\t", index=False)
