# Notice: task 07 data

The files under `inputs/data/` and `truth/` are a subset of the MOCHA manuscript repository:
https://github.com/aifimmunology/MOCHA_Manuscript, commit 48c75cd78075bc59a022d68f4090b0e831cd6a11,
licensed GPL-3.0.

What was taken, from `SuppFig17_ZeroInflatedGLMMs/`:

- `metadata.csv`: all 69 samples from 18 donors, reduced to sample, donor, time, time squared, age
  and sex (`inputs/data/samples.tsv`).
- `tf_matrix.csv`: the chromVAR z-scores of 30 of the 870 TFs, values copied as text
  (`inputs/data/tf_scores.tsv`).
- `tf_usage_model_results.csv`: the repo's time slope and p-value for the same 30 TFs
  (`truth/tf_time.tsv`).

`make-inputs.R` in this folder does the subsetting and records how the 30 TFs were picked.
These files are redistributed under GPL-3.0, as the source requires. The full text is in
`LICENSE-GPL-3.0.txt` next to this file. Only this folder (`eval/tasks/07_repeated-visits/`) is
GPL-3.0; the licence of the rest of the repository is undecided.
