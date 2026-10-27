"""Purpose: fit the disease-vs-healthy contrast per gene on donor-level expression.
Inputs:  data/donor_expression.tsv   [gap: toy data built inline until it exists]
Outputs: <run_dir>/contrasts.tsv, <run_dir>/top-gene_by-condition.{pdf,png,csv},
         <run_dir>/session_info.txt
Run:     python scripts/01_fit-contrasts.py results/2026-10-26_contrasts
Notes:   [Stochastic] none. Inference unit: donor. One row per donor x gene.
"""

# %% Setup
import argparse
import sys
from importlib.metadata import version
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import false_discovery_control, ttest_ind

# align copies palettes.py and theme.py to the project root, one level above scripts/.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from palettes import pal_condition, palette_for
from theme import save_panel

# Reference = healthy, so log2_fc reads as disease minus healthy.
CONDITION_LEVELS = ["healthy", "disease"]

parser = argparse.ArgumentParser()
parser.add_argument("run_dir", type=Path, help="new directory under results/")
run_dir = parser.parse_args().run_dir
if run_dir.exists():
    raise FileExistsError(f"run_dir {run_dir} already exists; pick a new name")
run_dir.mkdir(parents=True)

# %% Load inputs
# Toy data stands in for data/donor_expression.tsv; replace with the read call.
expression = pd.DataFrame(
    {
        "donor_id": [f"D{i:02d}" for i in range(1, 9)] * 3,
        "condition": (["healthy"] * 4 + ["disease"] * 4) * 3,
        "gene": ["IFI27"] * 8 + ["ACTB"] * 8 + ["CD3E"] * 8,
        "log2_expr": [
            *[5.1, 4.8, 5.4, 5.0, 8.2, 7.9, 8.5, 8.1],
            *[10.1, 9.9, 10.2, 10.0, 10.0, 10.3, 9.8, 10.1],
            *[6.2, 6.0, 6.4, 6.1, 5.3, 5.6, 5.1, 5.4],
        ],
    }
)
expression["condition"] = pd.Categorical(
    expression["condition"], categories=CONDITION_LEVELS
)
# pd.Categorical turns an unknown label into NaN silently, so check it.
if expression["condition"].isna().any():
    raise ValueError("condition has labels outside CONDITION_LEVELS")


# %% Fit contrasts
# Pooled-variance t-test: the same model as lm(log2_expr ~ condition) in the R example.
fits = []
for gene, gene_expression in expression.groupby("gene"):
    healthy = gene_expression.loc[
        gene_expression["condition"] == "healthy", "log2_expr"
    ]
    disease = gene_expression.loc[
        gene_expression["condition"] == "disease", "log2_expr"
    ]
    test = ttest_ind(disease, healthy, equal_var=True)
    fits.append(
        {
            "gene": gene,
            "log2_fc": disease.mean() - healthy.mean(),
            "p_value": test.pvalue,
        }
    )
de_results = pd.DataFrame(fits)
# Multiple-testing family: every gene tested in this script, adjusted together.
de_results["q_value"] = false_discovery_control(de_results["p_value"], method="bh")
de_results = de_results.sort_values("q_value")

# %% Plot the top gene
# Donor points, because the inference unit is the donor.
top_gene = de_results["gene"].iloc[0]
top_gene_expression = expression[expression["gene"] == top_gene]
fig, ax = plt.subplots()
ax.scatter(
    top_gene_expression["condition"].cat.codes,
    top_gene_expression["log2_expr"],
    c=palette_for(top_gene_expression["condition"], pal_condition),
)
ax.set_xticks(range(len(CONDITION_LEVELS)), CONDITION_LEVELS)
ax.set_title(top_gene)
ax.set_ylabel("log2 expression (one point per donor)")

# %% Write outputs
de_results.to_csv(run_dir / "contrasts.tsv", sep="\t", index=False)
save_panel(fig, top_gene_expression, run_dir, "top-gene_by-condition", 45, 45)
plt.close(fig)

# %% Session info
session_info = [f"python {sys.version}"] + [
    f"{pkg} {version(pkg)}" for pkg in ["matplotlib", "pandas", "scipy"]
]
(run_dir / "session_info.txt").write_text("\n".join(session_info) + "\n")
