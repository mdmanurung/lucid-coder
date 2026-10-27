"""Purpose: write the enrichment demo, two synthetic DE tables and a small GMT.
Inputs:  none
Outputs: examples/demo_T.tsv, examples/demo_Mono.tsv (pseudobulk-de columns),
         examples/demo.gmt (10 gene sets)
Run:     python examples/make_demo.py
Notes:   1000 genes G0000-G0999. In T, SET_UP (G0000-G0039) has its Wald stat
         shifted +3 and SET_DOWN (G0040-G0079) -3; in Mono every gene is null.
         Seven null sets of 20-60 genes, plus SET_TINY (8 genes) that the demo
         --min-size 10 removes. 20 genes per table have NA padj, as DESeq2's
         independent filtering leaves them.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SEED = 0
N_GENES, N_PLANTED, SHIFT = 1000, 40, 3.0
OUT = Path(__file__).resolve().parent

rng = np.random.default_rng(SEED)
genes = np.array([f"G{i:04d}" for i in range(N_GENES)])
for contrast, shift in [("T", SHIFT), ("Mono", 0.0)]:
    stat = rng.normal(0, 1, N_GENES)
    stat[:N_PLANTED] += shift
    stat[N_PLANTED : 2 * N_PLANTED] -= shift
    pvalue = 2 * stats.norm.sf(np.abs(stat))
    padj = stats.false_discovery_control(pvalue, method="bh")
    padj[rng.choice(np.arange(2 * N_PLANTED, N_GENES), 20, replace=False)] = np.nan
    lfc_se = rng.uniform(0.3, 0.6, N_GENES)
    de = pd.DataFrame(
        {
            "gene": genes,
            "base_mean": rng.lognormal(5, 1, N_GENES),
            "log2_fc": stat * lfc_se,
            "lfc_se": lfc_se,
            "stat": stat,
            "pvalue": pvalue,
            "padj": padj,
        }
    )
    de.sort_values(["padj", "gene"], na_position="last").to_csv(
        OUT / f"demo_{contrast}.tsv", sep="\t", index=False, float_format="%.6g"
    )

null_genes = genes[2 * N_PLANTED :]
sets = {
    "SET_UP": genes[:N_PLANTED],
    "SET_DOWN": genes[N_PLANTED : 2 * N_PLANTED],
    "SET_TINY": rng.choice(null_genes, 8, replace=False),
}
for i, size in enumerate([20, 25, 30, 35, 40, 50, 60], start=1):
    sets[f"SET_NULL{i}"] = rng.choice(null_genes, size, replace=False)
(OUT / "demo.gmt").write_text(
    "".join(f"{name}\tdemo\t" + "\t".join(sorted(g)) + "\n" for name, g in sets.items())
)
