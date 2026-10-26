I have single-cell RNA-seq data from six donors: three healthy and three with the disease, about 200 cells per donor, 50 genes.

- `data/counts.tsv`: raw UMI counts, genes × cells. The first column is `gene`; every other column is one cell, headed by its cell ID.
- `data/cells.tsv`: one row per cell, with `cell_id` (matches the count column names), `donor_id` and `condition` (`healthy` or `disease`).

Can you test which genes differ between disease and healthy? Write `de.tsv` (tab-separated) into the run directory with one row per gene and the columns `gene`, `log2_fc` (disease vs healthy), `p_value` and `q_value` (Benjamini-Hochberg over all genes).

Write one Python script at `scripts/01_compare-conditions.py` that takes the run directory as its only argument and writes its outputs there. Run it once with `uv run --python 3.12 --with pandas,scipy,matplotlib python scripts/01_compare-conditions.py results/run1`. Use only the standard library and: pandas, numpy, scipy, matplotlib.
