"""Purpose: freeze the task 13 inputs and truth from two nf-core scdownstream fixtures.
Inputs:  <download_dir>/counts_in_X.h5ad      [nf-core/test-datasets, branch scdownstream,
         <download_dir>/counts_in_layer.h5ad    anndata-variations/; MIT]
Outputs: eval/tasks/13_counts-layer/inputs/data/lab_a.h5ad   (counts_in_X.h5ad, verbatim)
         eval/tasks/13_counts-layer/inputs/data/lab_b.h5ad   (counts_in_layer.h5ad, verbatim)
         eval/tasks/13_counts-layer/truth/cells.tsv   (lab, cell_id, library_size, n_genes)
Run:     uv run --python 3.12 --with anndata python -I \
         eval/tasks/13_counts-layer/make-inputs.py <download_dir>
Notes:   [Stochastic] none. Not wired into eval/fixtures/make-fixtures.R; rerun by hand.
         lab_a holds integer counts in X. lab_b holds counts in layers["counts"] and a
         transformed float (0.25 * count + 0.125, never zero) in X. Truth comes from the
         counts slot of each file.
"""

import shutil
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

download_dir = Path(sys.argv[1])
task_dir = Path(__file__).resolve().parent
SOURCES = {"lab_a": "counts_in_X.h5ad", "lab_b": "counts_in_layer.h5ad"}

(task_dir / "inputs" / "data").mkdir(parents=True, exist_ok=True)
(task_dir / "truth").mkdir(exist_ok=True)
truth = []
for lab, source in SOURCES.items():
    target = task_dir / "inputs" / "data" / f"{lab}.h5ad"
    shutil.copyfile(download_dir / source, target)
    adata = ad.read_h5ad(target)
    counts = np.asarray(adata.layers.get("counts", adata.X).toarray())
    if not np.array_equal(counts, np.round(counts)):
        raise ValueError(f"{source}: counts slot is not integer")
    truth.append(
        pd.DataFrame(
            {
                "lab": lab,
                "cell_id": adata.obs_names,
                "library_size": counts.sum(axis=1).astype(int),
                "n_genes": (counts > 0).sum(axis=1),
            }
        )
    )
pd.concat(truth).to_csv(task_dir / "truth" / "cells.tsv", sep="\t", index=False)
