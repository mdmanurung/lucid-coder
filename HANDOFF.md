# Handoff — lucid-coder eval: variants run scored, readability review next

**Written:** 2026-11-09 13:29 CET · **Branch:** main @ 45f3465 · **Status:** 07c, 08b, 09b built and run (A and D, 3 reps); all 18 gates pass; arms not separated by gates.

## Check first
Run `git log --oneline 45f3465..HEAD` and `git status --short`. If either shows work this file does not mention, say so before anything else: State, Next action and Read first may be out of date, so re-check each against the code before acting on it.

## Goal
Show whether persona + house (arm D) makes analysis code more readable to a human than baseline Opus (arm A). The user ranks human readability above correctness. Done means the user's blind picks on the newest tasks, unblinded and summarised.

## Next action
Build the blind A vs D page for the variants run (the user has not said yet): `node eval/side-by-side.js eval/results/gated-variants-20261109 eval/results/gated-variants-20261109/review --arms A,D --rep 1`. Ask the user before sending it anywhere.

## State
- Uncommitted: `eval/results/results.tsv` (+18 rows, `gated-variants-20261109/*`); `.claude/` untracked. Committing results.tsv is the user's call.
- Tests run this session and passing: `eval/tests/side-by-side.test.js`, `compare-readers.test.js`, the gates for 07c, 09b, 12, 13 (`gates.test.js`) and 08b (`python.test.js`). The full suite was not rerun.
- Results: 07c A 3/3, D 3/3 (median log10 p shift 0.16, max 0.26, tolerance 0.75). 08b A 3/3, D 3/3 (all six used the donor-paired design: JAK-STAT 14.89). 09b A 3/3, D 3/3 (all six used the Goblet reference). `table fallback: gate pass/fail counts are binary per run`; rows in `results.tsv`.
- Assumed, not run: that the 07c prompt hides the donor from a model beyond what these six runs show; the 08b and 09b gate soft spots below were not stress-tested with other methods.

## Locked decisions
- Human readability matters more than correctness for this plugin (user, 2026-11-01).
- Keep both tiers, persona core + house rules (5-1 blind, 2026-10-28).
- Rule caps: core.md 3327/3350, house.md 1153/1200, r.md 1482/1500; a new rule needs a cut.
- Variant gates (user, 2026-11-08/09): 07c keeps `LOG10_P_TOL <- 0.75` and renames `donor_id` to `pid`; 08b ships the same 2,345 genes as 08 and passes JAK-STAT in 14.3-17.3 plus NFkB positive with padj < 0.05; 09b is 4 vs 4 mice, 6 types, Enterocyte x10, dispersion 0.1, padj 0.05, seed 1; run scope A and D, 3 reps.
- Decisions 2026-11-03: one `verify_<claim>.R` per claim spanning two outputs; R joins state `unmatched = "error"` when every row must match; paths use `here("data/x.tsv")` (see `ROADMAP.md`).

## Dead ends — do not redo
- Gate-only comparison on 07-16 and now 07c/08b/09b hits the ceiling: baseline Opus passes every trap, so gates cannot separate the arms.
- Fold learned rules (#4): one learned rule (r-0001, comments.steps) is not enough yet.
- Pruning `.claude/worktrees/`: all 12 have uncommitted files; check them before removing any.
- 07c with the donor hidden in the sample id: impossible, the real ids (`B024-AP0C2W3`) carry no donor.

## Read first
- `eval/TASK-CATALOG.md` — built table rows 13-15 list each variant's soft spots
- `eval/README.md` — score.js, judge.sh, side-by-side.js (`--seed`, per-task balance), report.js
- `ROADMAP.md` — phases; phases 1-3 mostly done
- `eval/tasks/08b_pathway-from-counts/check.R` — the band and the audit of wrong vectors

## Open items
- #3 second reader: packet at `eval/results/second-reader.zip` (gitignored), compare with `eval/compare-readers.js` — returns when a colleague has answered; lean (not a decision): none.
- #4 fold learned rules — returns after more `/lucid-coder:learn` rules from real use; lean (not a decision): wait.
- 08b band lets four wrong vectors through under the paired fit (|stat| 15.07, -log10 p 15.92, signed -log10 p 16.58, signed -log10 padj 16.59); only the NFkB clause rejects them — returns if a wrong answer ever passes both; lean (not a decision): keep the band.
- 09b passes edgeR with default TMM and its naive crumblr fails on Stem alone (q 0.032 vs 0.05) — returns if a crumblr version flips the verdict; lean (not a decision): keep as is.
- Next real work: readability review of the newest tasks, or the catalog's next tier (multiple-testing family, group-aware CV) — returns on the user's choice; lean (not a decision): readability page first.
