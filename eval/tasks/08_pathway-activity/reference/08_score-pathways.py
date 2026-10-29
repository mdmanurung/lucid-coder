"""Purpose: score PROGENy pathway activity for IFN-beta stimulation in CD14+ monocytes.
Inputs:  data/de.tsv        [gap: DESeq2 table, stim vs ctrl, provided with the task]
         data/progeny.tsv   [gap: PROGENy top 100 footprints, provided with the task]
Outputs: <run_dir>/pathways.tsv, <run_dir>/session_info.txt
Run:     uv run --python 3.12 --with decoupler==2.3.1,pandas,matplotlib python \
         scripts/08_score-pathways.py results/run1
Notes:   [Stochastic] none. The input vector is the signed Wald `stat`: it carries the
         direction and the precision of each change. log2FoldChange ignores precision;
         p-values lose the sign.
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
de = pd.read_csv(project_root / "data" / "de.tsv", sep="\t", index_col="gene")
progeny = pd.read_csv(project_root / "data" / "progeny.tsv", sep="\t")

# %% Score pathways on the Wald stat
wald = de[["stat"]].T.rename(index={"stat": "stim_vs_ctrl"})
score, padj = dc.mt.ulm(data=wald, net=progeny)
pathways = pd.DataFrame(
    {"pathway": score.columns, "score": score.iloc[0], "padj": padj.iloc[0]}
).sort_values("score", ascending=False)

# %% Write outputs
pathways.to_csv(run_dir / "pathways.tsv", sep="\t", index=False)

# %% Session info
session_info = [f"python {sys.version}"] + [
    f"{pkg} {version(pkg)}" for pkg in ["decoupler", "pandas"]
]
(run_dir / "session_info.txt").write_text("\n".join(session_info) + "\n")
