I have expression for 30 cells in two groups (A and B) and a small signed regulatory network: 5 transcription factors, 3 or 4 targets each.

- `data/expression.tsv`: log-normalised expression, one row per cell; the first column is `cell_id`, every other column is one of 20 genes.
- `data/cells.tsv`: one row per cell, with `cell_id`, `group` (`A` or `B`) and `sample`.
- `data/network.tsv`: one row per TF-target edge, with `source` (the TF), `target` (the gene) and `mor`, the mode of regulation: its sign says whether the TF activates (+) or represses (-) the target, its size how strongly.

Which TFs are active in each group? Score every cell for TF activity with decoupler's univariate linear model (`ulm`) and write `tf_scores.tsv` (tab-separated) into the run directory in long format, one row per cell and TF, with the columns `cell_id`, `source` and `score`.

Write one Python script at `scripts/15_score-regulators.py` that takes the run directory as its only argument and writes its outputs there. Run it once with `uv run --python 3.12 --with-requirements requirements.txt python scripts/15_score-regulators.py results/run1`. Use only the standard library and: pandas, numpy, scipy, matplotlib, decoupler.
