Two labs sent me single-cell RNA-seq of the same 128 genes, 128 cells each, as AnnData files.

- `data/lab_a.h5ad`: from lab A.
- `data/lab_b.h5ad`: from lab B.

Cell IDs restart in each file, so a cell is identified by its lab and its `cell_id` (the obs name in its file).

Can you combine the two into one object, normalised as log1p of counts per 10,000 (scanpy `normalize_total(target_sum=1e4)` then `log1p`), with the raw counts kept in `layers["counts"]`? Save it as `combined.h5ad` in the run directory. Also write `cells.tsv` (tab-separated) there with one row per cell and the columns `lab` (`lab_a` or `lab_b`), `cell_id`, `library_size` (total raw UMI count) and `n_genes` (genes with at least one count).

Write one Python script at `scripts/13_combine-labs.py` that takes the run directory as its only argument and writes its outputs there. Run it once with `uv run --python 3.12 --with-requirements requirements.txt python scripts/13_combine-labs.py results/run1`. Use only the standard library and: pandas, numpy, scipy, matplotlib, anndata, scanpy.
