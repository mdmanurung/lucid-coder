# Handoff — lucid-coder readability eval, tasks 07-16

**Updated:** 2026-11-02 · **Branch:** main · **Status:** 07-16 verdict D 12/13; phase 1 verdict D 5/6. Python `main()` guard added to rules/python.md (#7 closed).

## Check first
Run `git log --oneline 80eb63b..HEAD` and `git status --short`. If either shows work this file does not mention, say so before anything else: State, Next action and Read first may be out of date, so re-check each against the code before acting on it.

## Goal
Show whether persona + house (arm D) makes analysis code more readable to a human than baseline Opus (arm A). The user ranks human readability above correctness. Done means the user's blind picks plus the LLM-judge scores for tasks 07-16, unblinded and summarised.

## Next action
#27 (2026-11-02): pixi `dev` env installed (7.6 GB); `pixi run -e dev test` runs node tests then every skill suite. In dev, 20 of 22 suites pass after: with_status concat fix (cec1c1a), demos regenerated under pandas 2.3.3 and stub-Rscript tests for longitudinal/mashr (a35cc21). velocity (2 tests: pseudotime not byte-reproducible, NameError on the GPCCA error path) is being fixed. Suites take 7-56 min each on this file system; run them in parallel per skill.
#23 done 2026-11-02: label args in all 15 pipeline plot calls, Python stat_pwc label_size=7, 21 plots.py copies resynced. Skill pytest suites in the pixi default env fail identically before and after (anndata nullable-string writes, a pandas FutureWarning, velocity deps missing); worth its own issue.
New gated tasks 07b (NAs; wide drop) and 15b (A vs B test; samples are replicates) built 2026-11-02, gates proven in gates.test.js / python.test.js. Run `gated-07b15b-20261102` (A, D x2): 07b all 4 pass (both arms drop NAs per TF). 15b: A 1/2, D 0/2 pass; D ran a cell-level Welch, reasoning that sample is crossed with group and not stated to be a donor. The 15b prompt now names `sample` as the donor; rerun `gated-15b-donor-20261102`: A 2/2, D 2/2 pass (paired test on sample means). So the earlier fails were prompt ambiguity, and 15b is at ceiling once the donor is named.
Py rerun (2026-11-02, `py-main-20261102`, arm D only, 01py-03py x2): the new `main()` guard was taken up 6/6, gates 6/6, step comments in all; one header (02py r2) dropped the `Purpose:` label.
Phase 1 #2 done (2026-11-02): user picked D in 5/6 panels (01, 05, 06, 02py, 03py); A won 01py. Gates 48/48; header A 0/12, D 12/12; checks A 54, D 20; cost A $7.13, D $8.15. Judges at ceiling (Sonnet q2 A 4.67, D 4.92). Python gap from the reasons: the user wants `if __name__ == "__main__":` by default (only A used it), plus D's header and comments. That contradicts `rules/python.md` "flat top-level code"; #7 stays open for it. `review/key-with-picks.tsv`.

Done 2026-11-01: user picked D in 12/13 panels (all 9 whole-script panels and 3 of 4 task-14 sections; A only on 14_fit, "X is verbose but tolerable; W needs comments if extracting >2 results"). `review/picks-unblinded.tsv`. LLM judges at ceiling (q2 A 4.00/4.05 vs D 4.05/4.45). Next: decide on committing results and the open items.

## State
- Uncommitted: `eval/results/results.tsv` (+40 rows for gated-0716); `.claude/` untracked.
- Tests: 69/69 pass (`node --test tests/*.test.js`; npm is not on PATH, node is in `~/.local/bin`).
- Results: gates A 19/20 and D 19/20. Both fails are task 12 (low-depth SRR28679757 kept 0/360 cells). Table only, no figure: `results.tsv` rows `gated-0716-20261101/*`.
- Results: D has the full header 20/20 vs A 0/20, and 28 checks vs A's 65. Cost A $4.17, D $5.09.
- Reviewer read of the code: D reports drops, refuses run-dir overwrite and stops on an unmapped ID. A drops NAs silently (07, 08) and uses cells as replicates in 15 (the gate passed it anyway).
- Review page: https://claude.ai/artifact/QGgGghQP3DQMXPB6v7E1q1 (republished 2026-11-01 under m.d.manurung@lumc.nl) (db capability, `picks` collection empty), 13 panels, rep 1, A vs D. The key is in `review/side-by-side-key.tsv`; do not read it before the user finishes.

## Locked decisions
- Next priority was the 07-16 gated run, preceded by language-neutral house.md (80eb63b). Three-reviewer team discussion (user, 2026-11-01).
- Human readability matters more than correctness for this plugin (user, 2026-11-01).
- Keep both tiers, persona core + house rules (5-1 blind, 2026-10-28).
- Rule caps: core.md 3327/3350, house.md 1153/1200, r.md 1482/1500; a new rule needs a cut.

## Dead ends — do not redo
- Gate-only comparison on 07-16 hits the ceiling: baseline Opus passes the traps, so it cannot separate the arms.
- B/#4 "fold learned rules": first learned rule r-0001 (comments.steps, global) saved 2026-11-01 from the 14_fit pick; one rule is not enough to fold yet.
- Pruning `.claude/worktrees/`: all 12 have uncommitted files (1-35 each); check them before removing any.

## Read first
- `eval/README.md` — score.js, judge.sh, side-by-side.js, report.js usage
- `eval/rubric.md` — readability items the judges score
- `ROADMAP.md` — phases; ceiling result argues for harder or less-hinting prompts

## Open items
- Review page shows "page not found": this session publishes under the `m.d.manurung@lumc.nl` account, and the earlier review pages lived on another account. Fix by switching the browser account, sharing from that account, or opening the HTML file locally (picks go to localStorage). Wait for the user's choice.
- Unblind after the user's picks: join `picks` to `side-by-side-key.tsv`; report picks per arm with reasons.
- Commit `results.tsv`? This is the user's call.
- #7 Python house rules and #8 import time stay open; #2 (tasks 01, 05, 06 and the Python twins) is not run.
- Gates for 07-16 could test what D did differently (cells as replicates in 15, silent drops) — once the readability verdict is in.
