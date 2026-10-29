I have chromVAR transcription-factor scores for 69 blood samples from 18 COVID-19 patients, each sampled at 3 to 5 visits after symptom onset.

- `data/samples.tsv`: one row per sample, with `sample_id`, `donor_id`, `time` (days since symptom onset, already centred), `time_sqrd` (`time` squared), `age` and `sex`.
- `data/tf_scores.tsv`: one row per sample, with `sample_id` and one column per TF (30 TFs), the chromVAR deviation z-score.

Which TFs change over the course of infection? For each TF, model the score on age, sex, `time` and `time_sqrd`, and test the `time` slope. Write `tf_time.tsv` (tab-separated) into the run directory with one row per TF and the columns `tf`, `estimate` (the `time` slope), `p_value` and `q_value` (Benjamini-Hochberg over the 30 TFs).

Write one R script at `scripts/07_model-tf-time.R` that takes the run directory as its only argument and writes its outputs there. Run it once with `Rscript scripts/07_model-tf-time.R results/run1`. Use only base R and: dplyr, tidyr, readr, purrr, stringr, forcats, magrittr, ggplot2, here, limma, edgeR, lme4, lmerTest.
