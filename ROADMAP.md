# Roadmap

Status on 2026-10-28: **development paused; the plugin is in real use.** Tracking lives in GitHub issues under the roadmap issue; this file is the plan they come from.

## Where things stand (v0.4.0)

- **Rules.** `rules/core.md` is a persona: a senior computational biologist, five facts about where analysis code fails, eight habits with their why. `rules/house.md` restates the user-decided conventions as flat rules plus a do-not-add list. `rules/r.md` and `rules/python.md` hold language syntax. The session hook injects all of them (about 6 KB of the 9.5 KB cap).
- **Evidence.** Blind side-by-side reviews by one reader (the user):
  - lucid-coder vs baseline Opus, tasks 02-04: 8-3, 2 ties;
  - rule-list core vs persona core: 5-4, 2 ties;
  - persona + house vs persona alone: 5-1.
  - Persona + house cut script length by a quarter to a third against the persona alone.
- **Eval harness.** 10 tasks with planted traps and gates (`eval/tasks/`), the R4_51 conda R, notebook rendering via papermill, and blinded review pages built by `eval/side-by-side.js`.

## While in real use

- Correct the agent in the session; say "always" or "never" for a lasting rule, and accept the offer to save it (`/lucid-coder:learn`). Learned rules are the main input when development resumes.
- Note any habit the agent overdoes (extra knobs, audit prints, extra figures) or underdoes, with the file it happened in.

## When development resumes

Phase 1, measure what exists (no new tasks):
1. Run tasks 01, 05, 06 and the Python twins (01py, 02py, 03py) under persona + house.
2. A second reader judges the existing review pages, to check the picks are not one person's taste.
3. Fold learned rules from real use into `house.md` or the persona.
4. Decide the open questions: `verify_*.R` cross-artefact scripts; `unmatched = "error"` as part of the join default; the `here()` path form.

Phase 2, harness fixes:
- Balance column positions in two-arm review pages automatically (position bias found by hand on 2026-10-28).
- A Python counterpart to `house.md`; the persona and house rules were only measured on R.
- A Python env for eval tasks: decoupler v2, PyDESeq2, LIANA, pertpy with JAX. Pin decoupler v2 (LIANA >= 1.7 requires it).

Phase 3, new eval tasks (`eval/TASK-CATALOG.md`, sources in `eval/task-sources/`). Build first, each with its trap:
1. Repeated visits: no donor random effect (real MOCHA data with the repo's own lmer answer).
2. Pathway activity from a DE contrast: raw counts or unsigned statistic instead of the signed Wald stat (Kang 2018; JAK-STAT 17.15).
3. Cell-type composition: per-type proportion tests or cells as the unit (Haber 2017; only Enterocytes rise).
4. Gene identifiers: Ensembl IDs silently match few gene-set members.
5. Contrast direction: alphabetical factor order flips every sign.
6. Per-sample QC: one global cutoff wipes a low-depth sample.
7. Counts in a layer: normalising data that is already log-transformed.
8. Interaction contrast: comparing two DE lists instead of testing the interaction.
9. Regulator sign: a network without weights turns a repressor positive.
10. Visit and pair selection: duplicate runs or the wrong baseline visit.

Each task: a frozen small dataset under `inputs/`, a prompt, `check.R`, a reference answer that passes and a naive one that fails, then 4 arms x 2 reps and a blind review page. About 1-2 hours to build and $8 to run per task.
