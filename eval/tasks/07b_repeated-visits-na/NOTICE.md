# Notice: task 07b data

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

Task 07 (`eval/tasks/07_repeated-visits/`) made this subset; `make-inputs.R` there records how the
30 TFs were picked. Modified for task 07b: `make-inputs.R` in this folder sets 20 TF scores to
missing (HOXB2 in 17 samples, RUNX1 in 2, TET1 in 1) and refits the repo's model on each TF's
available samples for `truth/tf_time.tsv`.

These files are redistributed under GPL-3.0, as the source requires. The full text is in
`LICENSE-GPL-3.0.txt` next to this file. Only this folder (`eval/tasks/07b_repeated-visits-na/`)
and task 07's folder are GPL-3.0; the licence of the rest of the repository is undecided.
