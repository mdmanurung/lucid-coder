I have pseudobulk counts for CD14+ monocytes from eight donors, unstimulated and stimulated with IFN-beta (Kang et al. 2018).

- `data/counts.tsv`: raw counts summed per sample, one row per gene; the first column is `gene`, every other column is one of 16 samples.
- `data/samples.tsv`: one row per sample, with `sample`, `replicate` (the donor) and `label` (`ctrl` or `stim`).
- `data/progeny.tsv`: the PROGENy pathway footprints (top 100 genes per pathway), with `source` (pathway), `target` (gene), `weight` and `padj`.

Which signalling pathways does the stimulation activate or repress? Score each PROGENy pathway with decoupler's univariate linear model (`ulm`) and write `pathways.tsv` (tab-separated) into the run directory, one row per pathway, with the columns `pathway`, `score` (positive means more active after stimulation) and `padj`.

Write one Python script at `scripts/08b_score-pathways.py` that takes the run directory as its only argument and writes its outputs there. Run it once with `uv run --python 3.12 --with-requirements requirements.txt python scripts/08b_score-pathways.py results/run1`. Use only the standard library and: pandas, numpy, scipy, matplotlib, decoupler, pydeseq2.
