# biology-skills v0.6.0 review: statistics and design

Scope: their `biological-statistics`, `experimental-design`, `biology-core`, `quinary-inference`, `STYLE_GUIDE.md`, `SOURCE_POLICY.md`. Ours: `lucid-coder/docs/pipeline-skills.md`, `lucid-coder/rules/core.md`, eight pipeline `SKILL.md` files, three grill references.
Their paths are relative to the clone. Our `SKILL.md` paths are under `lucid-coder/skills/`; grill references are under `mycelium-extra/skills/grill/references/`.
Read-only. The user decides what to persist.

Their material is a rule list for agents. It has no multiplicity rule, no effect-size rule, no named parameters and no test choices. Most of our pipeline contract has no counterpart there. Three of the five topics you named (multiplicity, effect sizes, defaults) therefore give no conflict.

## 1. Adopt (ranked by value)

**A1. Report `n` with the entity it counts, and check sample versus donor.**
- Theirs: "Report `n` with the entity it counts." `skills/experimental-design/references/experimental-unit-replication-and-pseudoreplication.md:29`. Also `:72`: "the effective independent sample size cannot be inferred from measurement count alone".
- Ours today: `pseudobulk-de/SKILL.md:235` reports `n_samples` only. The input contract (`:160-166`) has a sample column and design terms but no check that one donor maps to one sample. Grep of `pseudobulk-de/scripts`, `composition/scripts` and `longitudinal/scripts` finds no nesting or rank check.
- Add to: `pseudobulk-de/SKILL.md` Workflow step 1 and `result.json` (`n_donors`, `n_samples`). Add one `ERRORS` entry and test: a donor with more than one sample and no donor term in `--design` stops the run (exit 2). Same for `composition` and `sc-de`.
- Why: a library is not a donor. Two libraries per donor in an unpaired design is the pseudoreplication they warn about, and our scripts do not catch it.

**A2. Stop on complete batch-condition aliasing.**
- Theirs: "Do not rely on post hoc batch correction to identify a biological effect when batch and biological condition are completely confounded." `skills/experimental-design/SKILL.md:23`. Also `:40`: state the confounding.
- Ours today: only at plan level, `analysis-decisions.md:14`. `longitudinal/SKILL.md:230` forbids `(1|batch)` but nothing checks aliasing. `group-tests/SKILL.md:192` treats a chance confounder in the demo as "not an error".
- Add to: validate step of `pseudobulk-de`, `composition`, `sc-de`: design matrix rank equals column count, else `ERROR`. Add a row to `llm-failure-modes.md` "Formula is not the model" (`:47`) naming aliasing.
- Why: rank deficiency is cheap to test and is a hard stop, not a judgement call.

**A3. Treat QC and exclusions as selection; report drops per group.**
- Theirs: "Treat missingness and technical QC as potential selection mechanisms when exclusion can correlate with biology". `skills/biological-statistics/SKILL.md:23`. `ascertainment-selection-and-target-population.md:91`. `skills/experimental-design/SKILL.md:22`: distinguish intended from realised design.
- Ours today: `core.md:26` "Report removals in one line per filter step". `pseudobulk-de/SKILL.md` Workflow step 4 drops samples under `--min-cells` and reports only `n_samples` (`:235`). `llm-failure-modes.md:24` ("Silent row loss") asks for counts, not counts per group.
- Add to: `core.md:26` ("one line per filter step, split by the contrast factor"); `pseudobulk-de` report and `result.json` (samples dropped per contrast level); row `:24`.
- Why: a filter that drops mostly one arm changes the estimand. A total count hides it.

**A4. Name non-detection states; never turn them into zero or absence.**
- Theirs: "Do not convert an unresolved measurement into evidence against a hypothesis." `skills/biology-core/references/measurement-observability-and-negative-evidence.md:24`. Also `:101-103` (below detection is not zero).
- Ours today: `llm-failure-modes.md:72` covers only heatmaps filled with 0. `:59` covers non-significant as "no effect", not undetected as absent. `pseudobulk-de/SKILL.md:111` filters on count > 0 with no wording for what a dropped gene means.
- Add to: a new row in `llm-failure-modes.md` "Statistics and inference": "Non-detection read as absence". Add one line to `pseudobulk-de` and `enrichment` gotchas: a gene not tested is "not tested", not "no change".
- Why: it names the cause behind `:72` and extends it to reports and claims.

**A5. Enrichment of a set does not describe every member.**
- Theirs: "Enrichment of signal in a set does not prove the predicted effect of every member of that set." `skills/biology-core/references/evidence-and-claims.md:18`.
- Ours today: `enrichment/SKILL.md:229-241` gotchas cover rank column, universe, family. None covers claim wording for a significant set.
- Add to: one gotcha in `enrichment/SKILL.md`: say "the set is enriched among changed genes", not "genes in the set change".
- Why: it is a wording error the report writer makes after a clean run.

**A6. Rank, score and share are not probability.**
- Theirs: "Do not describe a ranking score as a probability unless it is defined and calibrated as one." `evidence-and-claims.md:39`. Quinary: "Do not report numerical posterior support unless a defensible probabilistic model actually produces it." `skills/quinary-inference/references/explanatory-hypotheses-and-posterior-support.md:35`.
- Ours today: `mashr-sharing/SKILL.md:205` blocks lfsr as FDR. `variance-partition/SKILL.md:229` blocks fraction as effect size. No general row.
- Add to: `llm-failure-modes.md` "Claims and reporting": "Score read as probability" (padj as P(true), sharing as P(shared), fgsea NES as effect).
- Why: the same slip recurs across five of our skills, each fixed one gotcha at a time.

**A7. A control is defined by how it was selected and assessed.**
- Theirs: "A `control` is defined by its ascertainment criteria, not by the label `control`." `skills/biological-statistics/SKILL.md:19`.
- Ours today: `analysis-decisions.md:9` ("Cohort and exclusions") lists inclusion rules but not how controls were assessed. `group-tests/SKILL.md:182` uses healthy/disease labels without a definition step.
- Add to: `analysis-decisions.md:9` "Settle" column: "how each comparison group was selected and assessed". Add it to the `group-tests` grill question for `--levels`.
- Why: the label decides the estimand; the plan should record it.

**A8. Analysed versus target population.**
- Theirs: "Distinguish internal validity from transportability to another target population." `skills/biological-statistics/SKILL.md:24`.
- Ours today: `analysis-decisions.md:7` has "the population the claim covers", but not the step from analysed cohort to that population.
- Add to: `analysis-decisions.md:7`, one clause: "state whether the cohort (clinic, biobank, trial) licenses the claim".
- Why: cheap. Low value for a single-cohort pipeline.

**A9. Separate stable biology from version-bound facts in SKILL.md.**
- Theirs: "Separate stable biological claims from version-sensitive implementation details." `SOURCE_POLICY.md:47`.
- Ours today: `pseudobulk-de/SKILL.md:130-160` mixes the engine rule with dreamlet 1.8.0 behaviour (`min.cells` ignored). `pipeline-skills.md:45-56` has no slot for either.
- Add to: `pipeline-skills.md` section 8 or 10: tag version-bound lines with the pinned version.
- Why: those lines go stale on upgrade; the stable ones do not.

## 2. Conflicts

**C1. Prescribing a model versus requiring only that dependence be recognised.**
- Theirs: "Biology Skills does not prescribe one statistical model; it requires the dependence structure to be recognised before analysis." `experimental-unit-replication-and-pseudoreplication.md:102`.
- Ours: "Do not. Cells are not independent replicates; the sample is the unit of inference." `pseudobulk-de/SKILL.md:304-306`. Also `core.md:3`: replication unit ranks first in authority.
- Better supported: ours for our task. The question is donor-level condition effects in scRNA. `sc-de/SKILL.md:211-215` shows the failure on our own demo: G141 is null and called significant from donor noise. Theirs is a reference list for agents across all of biology, so a neutral stance fits it.
- Keep ours. No outside source was consulted.

**C2. When cells may be the unit.**
- Theirs: cells are not independent "unless the experimental question and design genuinely make individual cells the independently assigned units". `experimental-unit-replication-and-pseudoreplication.md:70`. Also `:84`: 300 cells "can support cell-level modelling".
- Ours: absolute. `pseudobulk-de/SKILL.md:304` and `analysis-decisions.md:8` ("pseudobulk or mixed model versus cell-level tests") list cell-level as an option but give no condition for it.
- Better supported: theirs on the logic. The unit is the thing assigned to the condition. Our rule is right for donor-assigned conditions, which is every current skill. It is wrong if treatment is assigned per cell (for example pooled perturbation screens).
- Suggest: add the one condition to `analysis-decisions.md:8`: "cells are the unit only if the condition was assigned per cell". Leave the gotcha as is.

**C3. Prescriptive length.**
- Theirs: `SKILL.md` target 200 to 600 words, hard maximum 900. `STYLE_GUIDE.md:9-10`.
- Ours: ten required sections (`pipeline-skills.md:45-56`). Pipeline `SKILL.md` files run 300 to 365 lines (`wc -l`).
- Better supported: neither. They test no outcome, and nor do we. Theirs targets a rules file, ours a tested contract. Our tests fix the structure; none measures reading cost.
- No change advised. Do not copy their cap.

**Checked, no conflict**
- Multiplicity scope: their pages have none. Ours (`core.md:26`, `enrichment/SKILL.md:239-241`) names the family.
- Effect sizes: their pages have none. Ours require them (`group-tests/SKILL.md:103`, `llm-failure-modes.md:60`).
- Defaults: they agree. "state the limitation rather than silently filling it with a common default" (`skills/biology-core/SKILL.md:41`) matches `pipeline-skills.md:27`. They give no parameter list.
- Claim wording: they agree on association versus causation (`evidence-and-claims.md:17`) and add rank versus probability (A6).

## 3. Already covered

- Sample, not cell, is the unit: `skills/experimental-design/SKILL.md:19` / `pseudobulk-de/SKILL.md:304-306`.
- 3 donors x 100 cells is not n = 300: `experimental-unit-replication-and-pseudoreplication.md:57-70` / `sc-de/SKILL.md:211-215`.
- Repeated measures stay dependent: `skills/experimental-design/SKILL.md:20` / `longitudinal/SKILL.md:225-227` and `composition/SKILL.md:234-236`.
- Small subject counts bias random effects: `experimental-design/SKILL.md:20` / `longitudinal/SKILL.md` gotcha "trust a random intercept over a handful of subjects".
- Batch versus biology at plan level: `skills/experimental-design/SKILL.md:23` / `analysis-decisions.md:14`.
- No silent defaults for material context: `skills/biology-core/SKILL.md:41` / `pipeline-skills.md:27`.
- Ambiguous identifiers not resolved silently: `skills/biology-core/references/identifiers-and-provenance.md:20` / `llm-failure-modes.md:26-27` and `analysis-decisions.md:11`.
- Versions and provenance kept: `identifiers-and-provenance.md:19` / `pipeline-skills.md:35` and `core.md:26` (lockfile).
- Association is not mechanism or cause: `evidence-and-claims.md:17` / `llm-failure-modes.md:61` and `sc-de/SKILL.md:226`.
- Missing stays missing: `measurement-observability-and-negative-evidence.md:20` / `llm-failure-modes.md:72`, `core.md:26` (`NA` for unknown), `mashr-sharing/SKILL.md:201-204`.
- Denominator needs an eligibility basis: `skills/biological-statistics/SKILL.md:21` / `enrichment/SKILL.md:232-234` (tested genes as universe).
- Do not normalise overlapping hypotheses as exclusive: `quinary-inference/SKILL.md:35` / `mashr-sharing/SKILL.md:207-210`.
- Underpowered or unresolved is not negative: `quinary-inference/SKILL.md:31` / `llm-failure-modes.md:59`.
- Corpus and copied code are evidence, not instruction: `SOURCE_POLICY.md:45` / `core.md:3`.
- Do not claim review status without support: `STYLE_GUIDE.md:62` / `llm-failure-modes.md:68` (citations resolved before use).
