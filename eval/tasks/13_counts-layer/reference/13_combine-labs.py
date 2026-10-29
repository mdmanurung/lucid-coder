"""Purpose: combine two labs' h5ad files on raw counts and log-normalise them.
Inputs:  data/lab_a.h5ad   [gap: lab A, provided with the task]
         data/lab_b.h5ad   [gap: lab B, provided with the task]
Outputs: <run_dir>/combined.h5ad   (X log1p CP10k, layers["counts"] raw)
         <run_dir>/cells.tsv       (lab, cell_id, library_size, n_genes)
         <run_dir>/session_info.txt
Run:     uv run --python 3.12 --with scanpy,anndata python \
         scripts/13_combine-labs.py results/run1
Notes:   [Stochastic] none. The labs store counts in different slots: lab_a in X,
         lab_b in layers["counts"] with transformed floats in X. Each file's counts
         are taken from whichever slot holds integers and must hold them.
"""

# %% Setup
import io
import sys
from contextlib import redirect_stdout
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
LABS = ["lab_a", "lab_b"]

run_dir = Path(sys.argv[1])
run_dir.mkdir(parents=True)  # fails if it exists: a run dir is never overwritten

# %% Load counts per lab
# layers["counts"] wins when present; X is used only when there is no counts layer.
adatas = {}
for lab in LABS:
    adata = ad.read_h5ad(ROOT / "data" / f"{lab}.h5ad")
    counts = adata.layers.get("counts", adata.X)
    values = counts.data if hasattr(counts, "data") else np.asarray(counts)
    if (values < 0).any() or not np.array_equal(values, np.round(values)):
        raise ValueError(f"{lab}.h5ad: no integer counts in layers['counts'] or X")
    adatas[lab] = ad.AnnData(
        X=counts.copy(), obs=adata.obs.copy(), var=adata.var.copy()
    )
    adatas[lab].obs["cell_id"] = adata.obs_names

# %% Combine and normalise
combined = ad.concat(adatas, label="lab", index_unique="-")
combined.layers["counts"] = combined.X.copy()
sc.pp.calculate_qc_metrics(combined, percent_top=None, log1p=False, inplace=True)
sc.pp.normalize_total(combined, target_sum=1e4)
sc.pp.log1p(combined)

# %% Write outputs
combined.write_h5ad(run_dir / "combined.h5ad")
cells = combined.obs[["lab", "cell_id", "total_counts", "n_genes_by_counts"]].rename(
    columns={"total_counts": "library_size", "n_genes_by_counts": "n_genes"}
)
cells.to_csv(run_dir / "cells.tsv", sep="\t", index=False)

# %% Session info
with redirect_stdout(io.StringIO()) as versions:
    sc.logging.print_versions()
(run_dir / "session_info.txt").write_text(versions.getvalue())
