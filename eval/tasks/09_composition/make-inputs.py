# Purpose: freeze the task 09 inputs from Haber et al. 2017 (pertpy haber_2017_regions), once.
# Run:     uv run --python 3.12 --with anndata python -I \
#              eval/tasks/09_composition/make-inputs.py <haber_2017_regions.h5ad>
# Inputs:  haber_2017_regions.h5ad  [https://exampledata.scverse.org/pertpy/haber_2017_regions.h5ad,
#                                    35,784,747 B; Haber et al. Nature 2017, GEO GSE92332;
#                                    no licence stated by pertpy]
# Outputs: eval/tasks/09_composition/inputs/data/cells.tsv  one row per cell, Control and
#                                                          Salmonella samples only
# Notes:   Only .obs is read. The Hpoly samples are dropped so the comparison is 4 Control vs
#          2 Salmonella samples, the scCODA tutorial contrast (only Enterocyte credibly up).
from pathlib import Path
import sys

import anndata as ad

CONDITIONS = ["Control", "Salmonella"]

out_path = Path(__file__).resolve().parent / "inputs" / "data" / "cells.tsv"
obs = ad.read_h5ad(sys.argv[1], backed="r").obs
cells = obs.loc[
    obs["condition"].isin(CONDITIONS), ["barcode", "batch", "condition", "cell_label"]
]
cells.columns = ["cell_id", "sample_id", "condition", "cell_type"]
out_path.parent.mkdir(parents=True, exist_ok=True)
cells.to_csv(out_path, sep="\t", index=False)
