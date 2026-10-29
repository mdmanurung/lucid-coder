# Naive answer: passes the network as given, so decoupler ignores `mor` and scores every edge +1.
import sys
from pathlib import Path

import decoupler as dc
import pandas as pd

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)

expression = pd.read_csv("data/expression.tsv", sep="\t", index_col="cell_id")
network = pd.read_csv("data/network.tsv", sep="\t")

scores, _ = dc.mt.ulm(data=expression, net=network, tmin=0)
res = (
    scores.rename_axis("cell_id")
    .reset_index()
    .melt(id_vars="cell_id", var_name="source", value_name="score")
)
res.to_csv(out / "tf_scores.tsv", sep="\t", index=False)
