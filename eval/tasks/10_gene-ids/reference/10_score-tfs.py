"""Purpose: score TF activity for TGF-beta vs control in hepatic stellate cells.
Inputs:  data/de.tsv.gz         [gap: DESeq2 table, treated vs control]
         data/genes.tsv.gz      [gap: gene annotation]
         data/collectri.tsv.gz  [gap: CollecTRI network]
Outputs: <run_dir>/tf_activity.tsv, <run_dir>/session_info.txt
Run:     uv run --python 3.12 --with decoupler==2.3.1,pandas,matplotlib python \
         scripts/10_score-tfs.py results/run1
Notes:   [Stochastic] none. de.tsv mixes Ensembl IDs and symbols; CollecTRI uses
         symbols, so every ID is mapped to its gene_name before scoring. Genes that
         map to one symbol keep the row with the largest |stat|.
"""

# %% Setup
import argparse
import sys
from importlib.metadata import version
from pathlib import Path

import decoupler as dc
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("run_dir", type=Path, help="new directory under results/")
run_dir = parser.parse_args().run_dir
if run_dir.exists():
    raise FileExistsError(f"run_dir {run_dir} already exists; pick a new name.")
run_dir.mkdir(parents=True)

# %% Load inputs
project_root = Path(__file__).resolve().parents[1]
de = pd.read_csv(project_root / "data" / "de.tsv.gz", sep="\t")
genes = pd.read_csv(project_root / "data" / "genes.tsv.gz", sep="\t")
collectri = pd.read_csv(project_root / "data" / "collectri.tsv.gz", sep="\t")

# %% Map gene IDs to symbols
symbol_of = genes.dropna(subset="gene_name").set_index("gene_id")["gene_name"]
de["symbol"] = de["gene"].map(symbol_of).fillna(de["gene"])
in_net = de["symbol"].isin(collectri["target"])
print(
    f"network targets matched: {de['gene'].isin(collectri['target']).sum()} as given,"
    f" {in_net.sum()} after mapping"
)
de = de.loc[de["stat"].abs().groupby(de["symbol"]).idxmax()]

# %% Score TFs on the Wald stat
score, padj = dc.mt.ulm(data=de.set_index("symbol")[["stat"]].T, net=collectri)
tf_activity = pd.DataFrame(
    {"source": score.columns, "score": score.iloc[0], "padj": padj.iloc[0]}
).sort_values("score", ascending=False)

# %% Write outputs
tf_activity.to_csv(run_dir / "tf_activity.tsv", sep="\t", index=False)

# %% Session info
session_info = [f"python {sys.version}"] + [
    f"{pkg} {version(pkg)}" for pkg in ["decoupler", "pandas"]
]
(run_dir / "session_info.txt").write_text("\n".join(session_info) + "\n")
