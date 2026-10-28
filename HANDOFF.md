# Handoff — lucid-coder

**Date:** 2026-10-27 · **Branch:** main, not pushed · **Status:** Phases 0-3 committed. On 2026-10-27 the rules were rewritten from the corpus code, R4DS/R Cookbook and Wilke. A review then aligned the skill doc, house-style table, templates, scorer and hooks with them. 80 tests pass. All grid numbers below predate the new rules.

## Goal
Build a ponytail-like Claude Code plugin with three properties:

- It steers agents to readable, reproducible R/Python for bioinformatics and statistics.
- Style stays consistent across repositories.
- It keeps only checks that matter, and it learns from user corrections.

Readability is the goal, not token efficiency. Done means phases 0–3 of `PLAN.md` section 9 pass their checks.

## Next action
1. Both worlds: `rules/core.md` is the persona (judgement + why), `rules/house.md` (new, `93679bb`) restates the user-decided conventions and a do-not-add list as flat rules; readCore() injects both. Rendered prompt 6.0 KB of 9.5. Not yet measured: run D on 01sk (2 reps) and 02-04 from the worktree (so the other session's tree is untouched) and compare blind against the persona-only runs.
2. Done: round 2 head-to-head (results below).
3. Persona core is live (`rules/core.md`, 3,343 bytes): identity, five failure facts, eight habits, authority. The rule-list core it replaced is in git history (`4b14e3b`). Reviewer reports in `skills/lucid-coder/references/persona-panel-20261028/`. Dropped from core in the swap, kept in r.md/python.md or house-style only: no globals, comment toggles, dead code, `_v2`, denominators, distinct seed per replicate, root-relative path wording, "let the package loop", `check_*.R` ban text, separate layers, legend order. Blind test: item 1.
4. The per-file inline packets (`packet-20261027`, `packet-skill-20261027`) are superseded: the user found whole-file rating fatigue-biased.
5. Decisions open for the user: allow `verify_*.R` cross-artefact scripts; `unmatched = "error"` as part of the join default; the `here()` path form; the other session's uncommitted `templates/plots.*`, `bench/`, `skills/*` and house-style rows.

## Rules vs persona v2 on 01sk, blind (2026-10-28)
- Page https://claude.ai/artifact/3doNepJZpSwBf9aGXHqUFi, rep 1, C = rule list, D = persona v2. Unblinded: persona won setup, filter-and-aggregate and the SKILL.md workflow; the rule list won load, fit and the notebook ("more concise in explanations": the persona notebook was 519 lines against 316, with an extra MDS figure); write was a tie.
- Running tally of rule list vs persona, like-for-like panels only: rules 5, persona 4, ties 2. Parity, one reader.

## Persona v2 (2026-10-28)
- The task 02 diff showed the persona adding a `message()` of dims read and a `CONDITION_LEVELS` block the task did not need. Three clauses changed in core (`b46adf3`): report what a step removed, nothing when it removed nothing; a task with no knob has no parameters block; no audit of what went right. Cut to pay for them: the re-plot sentence, "names say what objects hold", "QC and annotation are notebooks too".

## Rules vs persona, blind (2026-10-28)
- Page https://claude.ai/artifact/39ScY8EorPPrspU484YF4y, whole scripts, C = rule-list D (sbs2), D = persona D. Unblinded: rules won 02 and 04, persona won 03. Scripts under the persona were longer (34/137/123 vs 23/97/85).
- Reasons: 02 and 04 "but use suppressPackageStartupMessages" (the rule-list side lacked it; the persona side had it, so no change). 04: four six-line `write.table()` calls lost to the walked list: "if the write table is lengthy like this because so many arguments then I prefer [the list]". Habit 3 now reads: each output gets its own write line, a one-line call or a helper when the call is long; never a walked list. That matches round 1 (four `write_output()` calls beat `iwalk`) and this pick.
- Scope check before the review: 02 and 04 outputs byte-identical across sides; 03 same 40 rows, median line vs mean ± SD; the page carried a scope note per panel (`--notes`).

## Branch note (2026-10-28)
- The other session checked out `docs/clawbio-determinism` in the shared tree while this session was committing; commit `b315a9a` (persona core) landed there and was cherry-picked onto main as `5eb69dc`. Commits from this session now go through a worktree on main. The duplicate on their branch is harmless.

## Round 2 head-to-head results (2026-10-28)
- Unblinded: D won 02 whole, 03 setup ("too many checks" in A), 03 load, 04 load ("concise"), 04 setup; A won 03 figure (plot frame computed before the ggplot call, "offload the wrangling outside the plotting"), 04 filter (6-line `colSums`/`rowSums` matrix filter beat a 45-line long-format pipeline), 04 write (four `write_table()` calls beat `iwalk`, which D used despite the new rule). Totals across both A-vs-D rounds: D 8, A 3, no preference 2.
- Applied: `suppressPackageStartupMessages` (asked twice) to r.md; `r.matrix`, `fig.plot_frame` and `r.startup` rows in house-style; the matrix and plot-frame habits in the persona.

## Tasks 02-04 round 2 (2026-10-28)
- `eval/results/sbs2-0204-20261027`: 4 arms × 1 rep under the header, divider and write rules; gates 12/12. Script lines: 02 A 50, B 28, C 32, D 23; 03 A 120, B 35, C 93, D 97; 04 A 90, B 49, C 91, D 85. side-by-side.js now skips a step that only some versions have (two panels had an empty side).
- ub-mode (github.com/unclebob/ub-mode) read 2026-10-28: a 3,958-byte SKILL.md of ~40 bare imperatives in three phases (make it work, make it right, practice) paired with deterministic gates (crapper, dryer, mutator). Taken: "check the data where it enters, inside trust it" (our check policy), "if the test would still pass when every function returned nothing, throw it out" (checks able to fail), "pin what the code does, change the structure, prove the behaviour did not change" (task 04 gate). Not taken: its shape, a rule list with no voice or why, which is what the user is moving away from.

## Whole-script head-to-head (2026-10-27, night)
- Task 02: D won. Task 04: no preference, "need a mix": A for its descriptive header, full-width dividers and one write call per table; D for UPPER_SNAKE parameters, inline `case_when` without helpers, and comments that explain the next block.
- Like-for-like totals for A vs D: D 3 (02 whole, 03 setup, 03 load), A 0, no preference 2 (03 figure, 04 whole).
- Applied: header lists Inputs and Outputs one file per line with what it holds (core, `comments.header`); every block, parameters too, opens with a divider padded to column 80 (r.md, `code.sections`); one write call per output, never a list walked with `iwalk()` (core ladder 4, `code.writes`). Both templates follow. To fit the caps, core dropped "legend order matches plot order" and the copied-notebooks clause, and r.md shortened the sessionInfo line.

## Head-to-head results (2026-10-27, night)
- First page, unblinded: D won 9 panels, A 1, no preference 1. But arm A wrote tasks 02 and 04 without section banners, so its whole script sat in setup and its other columns were empty: 6 panels compared D against nothing, and 02/04 setup compared a full script against a header. Only the three task-03 panels were like for like: D won setup and load, the figure was a tie. The user flagged the 02 setup panel.
- Fix: side-by-side.js now shows a task as one whole-script panel when any version has no banner; `--tasks` filters tasks. Tasks 02 and 04 went out again as whole scripts.

## Tasks 02-04 for the head-to-head (2026-10-27, night)
- `eval/results/sbs-0204-20261027`: tasks 02-04 × A-D × 1 rep at BUDGET=3, gates 12/12. Script lines: 02 A 57, B 33, C 22, D 28; 03 A 153, B 35, C 71, D 81; 04 A 94, B 58, C 88, D 79. The page compares A and D only, at the user's request for two options at a time.

## Round 2 review results (2026-10-27, night)
- Picks, unblinded: D won setup, load, fit and the notebook; C won filter-and-aggregate and write; A won the SKILL.md Workflow; B none. Plugin arms 6/7 (round 1: 4/7). One rep, one reader.
- Reasons: informative header and parameters up front; key checkpoints, concise; package functions instead of new ones; white space and sectioning with adequate comments; notebook with parameters, settings and libraries first. All match current rules; no rule change.

## 01sk rerun under the step-comment rules (2026-10-27, night)
- `eval/results/skill-nb2-20261027`: gates 8/8 at BUDGET=5, $2.88-5.22 per run (C r2 hit the cap but passed). Script lines A 399/302, B 152/142, C 175/149, D 144/161. UPPER_SNAKE parameters in rep 1: C 4, D 3, A and B 0.

## Side-by-side review results (2026-10-27, night)
- The user answered all 19 panels. Agent panels, unblinded: D won setup, load, fit and the SKILL.md Workflow; B won write and the notebook; A won filter-and-aggregate; C won none. One rep, one reader: a lean, not an effect.
- Applied as rules: top-of-script parameters in UPPER_SNAKE (reverses the R4DS lowercase rule; rubric H6 flipped); a short comment naming each step above its block, labels on long pipe verbs, a blank line between steps (core "Comments", house-style `comments.steps`); SKILL.md workflows keyword first (`docs.workflow`). To fit core, two lines moved out: the `check_*.R`/`anyNA()` examples and the separate-layers line.
- Reinforced, no change: native `|>` (flagged on two human picks), one `stopifnot(identical(...))` on cell IDs vs count columns, no one-use functions, concise notebooks.
- Not applied: "save key parameters in result.json, but toned down" has no clear threshold yet.

## Notebook split and 01sk rerun (2026-10-27, evening)
- Core rule: scripts prep data, fit models and write tables; figures and reports go in a notebook that reads those tables. 01sk now asks for `report.ipynb` (R kernel `ir`), run by papermill and rendered to `report.md` by nbconvert.
- The harness, scorer and gate tests use the R4_51 conda env (R 4.5.1, tidyverse, IRkernel, papermill, pandoc 3.4). Its bin holds Claude Code 2.0.73, so run.sh pins `claude` first. Gate tests 15/15 pass; fixtures regenerate identically.
- `eval/results/skill-nb-20261027`: 4 arms × 2 reps at BUDGET=5 (at $2, two arms ran out). Gates 8/8. Cost $3.26-5.15 per run. Script lines (score.js) fell from 179-445 under the old contract to 106-280.

| Arm | Script lines r1 / r2 |
|---|---|
| A | 228 / 280 |
| B | 161 / 140 |
| C | 155 / 106 |
| D | 188 / 174 |

## Grid 5 under the current rules (2026-10-27, night)
- `eval/results/grid5-20261027/`: tasks 01-04 × A-D × 2 reps, run after `52f5d2b`. All 32 runs succeeded, gates 32/32, reruns 32/32. Cost: A $2.88, B $2.27, C $1.67, D $2.12.

| Arm | air passes | Header | Lines (median) | Code lines | Checks | OUT of budget | Scripts with OUT = 0 |
|---|---|---|---|---|---|---|---|
| A | 5/8 | 0/8 | 98 | 65 | 36 | 23 | 2/8 |
| B | 1/8 | 0/8 | 58.5 | 44.5 | 31 | 18 | 2/8 |
| C | 8/8 | 8/8 | 77.5 (rerun 4: 121.5) | 56 (78) | 19 | 4 | 4/8 |
| D | 8/8 | 8/8 | 74 (rerun 4: 116) | 51 (71.5) | 18 | 2 | 6/8 |

- The plugin arms are now shorter than baseline; under the old rules they were longer. Arm D followed its learned rules (`out_dir` and `options(warn = 1)`) in 8/8 runs. Baseline Opus now uses `|>` in 8/8, so the pipe rule no longer separates the arms. A and B spend most checks on input audits at load (23 and 18 OUT).
- Check classification (Opus agent, strict): C/D 6 OUT of 37 (16.2%), median 0. Five of the six are a second guard on one join, three of them `unmatched = "error"` next to `relationship =`. One is an `identical()` that cannot fail. If `unmatched = "error"` counts as part of the join default, C/D are 3 OUT of 34 (8.8%), as in reruns 3-4. Open question: `r.md` names only `relationship =`; `unmatched = "error"` is what stops silent drops, so it arguably belongs in the default.
- Gaps: 3 C/D scripts (02_C_r1, 04_C_r1, 02_D_r2) have no overwrite guard. C/D have no usage guard (A/B 8/8). 26 of 34 C/D lints are `nonportable_path_linter` from `here::here("data/x.tsv")`; `here::here("data", "x.tsv")` avoids them, but `r.md` has 10 characters free.

## Review and alignment (2026-10-27, late)
- A review found the injected rules (`rules/`) were current but everything around them was stale. Fixes, one commit each:
  - `2ce4e5c`: SKILL.md defers to `rules/` (347 to 186 lines); house-style rows for checks, exclusions, notebooks, automation, constants, width, errors, joins, recode and Python args revised.
  - `015b7ec`: template examples rewritten to the new style, outputs byte-identical.
  - `52f5d2b`: `header_complete` no longer needs Notes; new `join_guard_count` column; agents may run `air`/`ruff` with the template configs; R prompts allow stringr and forcats. Rescoring rerun 4 changed no stored value.
  - `b954475`: `/lucid-coder lang` is per repo (JSON map keyed by repo id); language-scoped learned rules follow the languages in use; on/off/lang replies say they apply from the next session start or `/clear`.
- Host fact: the R 4.5.2 module has no tidyverse meta-package, `qs`, `limma`, `janitor` or `broom`. So no MPL demo script could run here; all were judged by reading.
- Caps: core 3,332/3,350 and `r.md` 1,490/1,500 leave no room; a new rule needs a cut.
- Another session works in this repo concurrently (`bench/`, `skills/pseudobulk-de/`, `templates/plots.*`, a `fig.text` row). Commits from this session stage named files only.

## Before/after demo on a user script (2026-10-27)
- Live A-vs-C refactor of `Senegal_phipseq/.../src/_MPL_differential_analysis.R` (100 lines, limma trend test on MPL-stimulated subcluster frequencies). Same prompt, Opus, `--safe-mode`, no data; C got `arm-prompt.js C r` (4,580 chars). Files in `eval/demo/mpl-refactor/`: `original.R`, `prompt.md`, `A_refactored.R` (144 lines, $0.94, 9 turns), `C_refactored.R` (146 lines, $1.98, 14 turns).
- A kept the script's shape and added explanation: a purpose comment, a `tidy_toptable()` helper, `REFERENCE_SAMPLE`, and a comment warning that `matCounts` columns and `matMeta` rows align by position. It kept `%>%`, `remove = F`, dot names, 4-space indent, all 12 libraries and three `source()` calls.
- C changed what a reader can verify: header with producers (Zenodo id from the README, `[gap]` for the metadata file) and consumers (Fig4-6); `here::here()`; `|>`; snake_case; the three libraries the script uses and one `source()`; constants `N_THREADS`, `REFERENCE_SAMPLE`, `ROBUST_MAXIT` with a reason each; `message()` with the count at each drop; `coef = "group.L"`/`"group.Q"` by name (verified: `weightedcontrasts.R` names the columns `.L`, `.Q`, `.C`); writes grouped at the end; `sessionInfo()`. Exactly one check, the default "IDs identical before a positional copy", as a `stop()`; A only commented on the same risk.
- Caveats: C dropped the unused `groupCat = group3` column from the design table (safe, nothing reads it) and calls LD413 a "reference (control) sample", which is its guess (the other scripts only drop it). Neither version was executed.
- Write-from-scratch twin (`eval/demo/mpl-from-scratch/`): the same analysis given as a five-step spec with the two data files described and `weightedcontrasts.R` present, no source script. A: 100 lines, $0.57, 4 turns; C: 196 lines, $0.85, 8 turns. Both arms named the coefficients (`"group2.L"`), aligned the design to the matrix with `match()` + an identity check, and used `droplevels()`, so those are Opus defaults, not plugin effects. A's three checks are `stopifnot()` (one-batch-per-patient, order, `nlevels == 3`); C's five are `stop()` with a fix (two column checks at load, one-batch-per-patient, order, outputs non-empty) plus `relationship = "many-to-one"`/`"one-to-one"` on every join. C's additions over A: the header (with honest `[gap]`s and a Notes block), `here()`, `|>`, four named constants with a reason each (LD413 is "excluded at the analyst's request; reason not recorded", not guessed this time), a `message()` listing the patients removed by `drop_na()`, writes grouped at the end, `session_info.txt`. Cost of the additions: twice the line count; the Notes block spends four lines on why `clean_names()` was not applied. Neither passes `air format --check` (no over-width lines in C; both differ in call layout). Neither was executed.
- Arm D on the same spec (`D_MPL_trend_test.R`, 171 lines, $0.93, 6 turns): the three learned rules took (`options(warn = 1)` on line 1, `out_dir`, no plot to title). Checks: the same five as C plus a refuse-if-outputs-exist guard (the core's "new run dir" default, applied to fixed paths) and `relationship =` on joins; it dropped C's outputs-non-empty check. Shorter than C because the drop logging is two `message()` calls and the Notes block states facts (FDR family, medians use all cells) rather than a `clean_names()` apology. `LINEAGE_PATTERN <- "_.*$"` as a constant is over-naming. Fails `air format --check` on call layout, no over-width lines.
- Pipe rule (`83d4b80`): `r.md` now asks for a new line after every `|>`. D rerun on the same spec (`D2_pipe-rule_MPL_trend_test.R`, 178 lines, $0.25): no step after `|>` on the same line (old D had several). It added one spare check: a row-count `stop()` after a join that already has `relationship = "many-to-one"`. HTML report of the examples: `eval/demo/mpl-refactor-examples.html`, published at https://claude.ai/artifact/KsHsD2dAij1dRJV7TMZKSR (private).
- air as the layout authority (`5025e08`): `r.md` now tells agents to run `air format <file>` after every edit and to build `message()`/`stop()` text with `sprintf()` (air puts comma-pasted pieces one per line); the "one argument per line" clause went, since air does it. `air` keeps short pipe chains on one line, so the pipe rule stays. From-scratch rerun with `Bash(air:*)` allowed (`eval/demo/mpl-from-scratch-air/`): all three pass `air format --check` (C and D ran air; A matched by chance). lintr: A 34 (27 `%>%`, 5 path strings, 2 names), C 0, D 1. Substance left: A has 4 bare `stopifnot()`, no join guards, silent drops, no constants, no session info; C/D have messages with fixes, `relationship =` on joins, logged drops by ID, named constants. New plugin faults: C strips the lineage prefix from `cluster` (breaks joins on the full name downstream) and again adds a row-count check after a `relationship =` join; D refuses to overwrite fixed output paths. Report republished (version 2) at the same URL.
- Corpus rules (`d2ef5da`): three Opus agents read the code of 17 repos (`cc_learnings/repos/`, `repos_new/`; exemplars immgenT for R, sound-life and aifi for Python) for brevity, project shape and notebook use. Evidence in `skills/lucid-coder/references/corpus-patterns.md`. Core lost the load-time column checks, outputs-non-empty and per-drop message blocks; checks are one-line named `stopifnot()` / `validate=`; added split-map-combine, a function only when a block repeats, package getters, compute steps writing files with one figure script or notebook each, source data from the plotted object. Open conflict: immgenT's `verify_*.R` against core's `check_*.R` ban.
- Book rules (R4DS 2e, R Cookbook 2e, Wilke, data-to-viz; `e970574`, `3d5cc3f`), arm D rerun (`eval/demo/mpl-books-rules/`, one run each). R: the literal R4DS "named arguments one per line" made the agent split every call, so the script was 157 lines and 121 code lines (`D_R_named-args-rule_157-lines.R`). The rule now reads "one argument per line only when a call overflows", and the rerun is 106 lines, 78 code lines, air-clean. It uses `library(tidyverse)`, `join_by()`, `.by =`, `separate_wider_delim()`, `fct()` and lowercase constants. It splits `cluster` into lineage and annotation, as R arm C did twice. Read that as an ambiguous spec ("lineage, cluster") rather than a plugin fault. It also drops `dir.create(out_dir)`. Python: 92 lines, 74 code lines, ruff-clean, `pd.crosstab` for counts, a dict comprehension of fits. It again leaves lineage-absent subclusters as NaN, now stated in the docstring.
- Tightened the demo specs (`eval/demo/mpl-from-scratch/prompt.md`, `eval/demo/mpl-python/prompt.md`): `cluster` is the full subcluster name; a patient with no cells in a whole lineage gets 0, as the original's `values_fill = 0` does; create `outputs/`. Arm D rerun (`eval/demo/mpl-tight-spec/`): R 100 lines (70 code), air-clean; Python 94 lines (77 code), ruff-clean. Both meet the three points. New minor issue: R renames `subclusterPG_merged` to `cluster` on load, so the median table's key column changes name.
- Rerun under the corpus rules, from-scratch MPL spec (`eval/demo/mpl-concise-rules/`, `eval/demo/mpl-python/`, one run per arm). R lines: A 112, C 97 (old rules 183), D 102 (old 189). Python lines: A 162 (8 functions + `main()`), C 116, D 98. R C again strips the lineage prefix from `cluster` and has one 102-column line; Python D leaves lineage-absent subclusters as NaN instead of 0. Python spec uses `eval/demo/mpl-python/weighted_contrasts.py`, a port of `contr.poly.weighted()` checked equal to the R output on a 3/5/4 factor. Report version 4 at the same URL.

## Rerun 4 and the human packet (2026-10-26, night)
- Second grilling: rerun 4 before the packet; name the factor-NA not-check exactly; add a ruff width line to `python.md`; keep the traps as evidence only; R-only packet of 16 (pilot A/B rep 1 + rerun-4 C/D rep 1), judges on the same 16, Python packet only if the R packet shows an effect.
- Rerun 4 (`eval/results/rerun4-cd-20261026/`, C/D × 4 × 3): 24/24 gates, headers 24/24, `|>` 24/24, D rules 12/12. The `stop()` width fix worked: over-width lines 89 → 28, over-width `stop()` 52 → 11, lint median 4 → 1 (11 scripts at 0). `air` still fails 23/24 for a second, mechanical reason: agents put several arguments per line in a split call; air wants one. `r.md` now says so (`5b70501`), unmeasured. Agents cannot run air in the harness, so `air_ok` measures hand-formatting.
- Check classification, rerun 4: 122 checks, 102 default, 8 extra, 12 OUT (9.8%), median 0; rerun 3 was 11 OUT (8.5%), median 0. The literal `anyNA()`-after-`factor()` pattern is gone; the same intent survives as `setdiff(x, LEVELS)` guards (4, from 6). C carries most OUT (10 vs D 2): redundant second `if` on one join (3) and duplicate-row audits (2). So the budget holds at median 0 but C's wording still leaves about one spare check in half its scripts.
- Packet: `eval/results/packet-20261026/packet/P01.R`-`P16.R` + `sheet.tsv`; key at `packet-20261026/packet-key.tsv` (gitignored; do not open before the sheet is filled). Judges (Opus, Sonnet) ran on the same 16 → `packet-20261026/judges.tsv` (committed; keys stay ignored). Proxy only: q1 (one read) is at the ceiling for every arm (A 5.0/4.75, B 4.5/4.75, C 5.0/5.0, D 5.0/5.0, Opus/Sonnet means of 4); q2 (find decisions) edges to C/D (Opus: A 3.75, B 3.75, C 4.25, D 4.25; Sonnet: A 4.5, B 4.25, C 5.0, D 5.0); q3 "would change the checks" is yes for 30 of 32, so it does not separate arms. The human sheet decides.

## Tiers built, rerun 3 and the traps (2026-10-26, evening)
- Built and committed: two-tier injection (`rules/core.md` 3,291 chars, `rules/r.md`, `rules/python.md`; the 2,500 core cap was dropped by the user, since every remaining line is a rule), tasks `05_recode-trap` and `06_grep-trap`, the Python scoring path and twins `01py`-`03py`. 74 tests pass.
- Rerun 3 (`eval/results/rerun3-cd-20261026/`, C/D × 4 tasks × 3 reps): all gates pass, headers 24/24, `|>` 23/23, D rules 12/12. Check classification (Opus agent, strict): 130 checks, 111 default, 8 extra, 11 OUT (8.5%), median OUT 0 per script in both arms. Rerun 1 was 16 OUT of 95 (16.8%), median 1. So the "Not checks" line met its target. The leak that remains is `anyNA()` after `factor(levels = )` (6 of 11 OUT), the exact case the line names.
- Regression in rerun 3: `air` failed on 24/24 scripts (reruns 1-2: 0). The one-line error template in `r.md` was copied as single long `stop()` lines (52 over-width lines; rerun 2 had 1). Fixed in `9e6e632`: the template now says to keep the message short because `air` does not wrap strings. Not yet re-measured.
- Python grid (`eval/results/python-20261026/`, `01py`-`03py` × 4 arms × 2 reps): 24/24 gates pass. Headers C/D 12/12, A/B 0/12. Checks per script: A median 7.5 (02py 10, 10), B 4, C 5, D 3. D's `out_dir` rule 6/6; no-title rule on 03py D 2/2, C 1/2, A 0/2. `ruff format --check` passes only 2/24 (C): the agents write 89-95 column lines that ruff wraps; lint medians C/D 1, A/B 2 (E501 and B905 `zip()` without `strict=`). Python tier text is only 717 chars, so a format-width line there is cheap if wanted.
- Traps (`eval/results/traps-20261026/`, 05 and 06 × 4 arms × 3 reps): 24/24 pass, baseline included. The traps do not discriminate because the agent reads the data before writing (arm A normalised `"BR 1"`/`"br2"` with `toupper(gsub())`; every arm used `!grepl()`). They still show the budget: C/D 6 checks on 05 (all inside the budget) vs A 4-5, B 3-4; on 06 C/D 4 vs A 1-2, B 1.

## Second C/D rerun and the grilled design (2026-10-26)
- `eval/results/rerun2-cd-20261026/`: 16 runs, all gates pass. Check medians: C 6 (rerun 1: 5), D 5 (rerun 1: 6); pooled over both reruns C 6, D 6. So the right figure for the budget rule is 8 → 6, not 8 → 5. Per task, rerun 2, C: 01 8,7; 02 5,4; 03 7,7; 04 5,6. D: 01 5,7; 02 5,4; 03 8,7; 04 7,4.
- Check classification over the 16 rerun-1 scripts (95 checks): 73 default-five, 6 allowed extra, 16 outside the budget, median 1 per script. The overruns are unknown factor level (6), one label per donor (4), duplicate rows after a join (3), file exists (2), NA audit (1). Raw count is the wrong metric: the budget allows 6 for a script with a join, so a median of 6 can be clean.
- Design settled in a grilling session: (1) core shrinks to ≤ 2,500 chars and gains one line naming the three top overruns as not-checks; language tiers `rules/r.md` and `rules/python.md` (≤ 1,500 each) load by file-extension scan under cwd to depth 3, override with `/lucid-coder lang r|python|both|auto`; learned rules unchanged; same text to subagents. (2) Two R trap tasks, 05 recode with a catch-all and 06 `x[-grep()]` on a pattern with no match, run at 4 arms × 3 reps. (3) Python twins of tasks 01-03 with ruff in place of air/lintr and `uv run` for the rerun, 4 arms × 2 reps. (4) Human reading by the user alone on a stratified half (16 scripts, one rep per task × arm), LLM judges as the second rater.

## Per-file audit (2026-10-26, after the C/D rerun)
- Three Opus agents audited the house style against the primary digests, `cc_learnings/repos_summary/` (17 files) and `cc_learnings/per_file/` (198 shards), which no earlier audit had read. A Fable agent spot-checked 12 cited lines and rejected six audit claims, among them "zero `oob = squish`" (MM and immgenT use it) and "`call. = FALSE` in MM only" (immgenT too).
- Applied in `house-style.md`: about 25 evidence corrections and relabels (`r.errors` now credits immgenT's `fn():` prefix and `call. = FALSE`; `r.indent` and `r.assign` are split, not dominant; `fig.theme` names MM, MOCHA and IHA; `comments.header` "no author line" is new, since every corpus header has one; `py.paths` is new; `here()` has 0 of 17 witnesses); 6 new rows: `code.constants`, `code.by_name`, `py.views`, `data.write`, `fig.panel_id`, `fig.once`; 7 clause merges. Dropped as covered or single-lineage: bare-except as a separate row, notebook parameterising, pure helpers, stats outside the plot, scanpy palette form.
- Applied in `SKILL.md`: the same rules plus four before/after examples (constant swallowed by `...`, recode without a catch-all, Python bare `except`, AnnData view mutation).
- `rules/core.md` is unchanged at 3,990 chars: the reviewer judged core should shrink, not grow, because agents still over-check. So the eval arms, which inject only core, have not seen today's additions.
- A second C/D rerun, `eval/results/rerun2-cd-20261026/`, started after these edits. It replicates the first rerun (same core text), taking n to 4 per task and arm.

## C/D rerun and corpus audit (2026-10-26, after the pilot)
- Rerun: `eval/results/rerun-cd-20261026/`, tasks 01-04 x arms C,D x 2 reps, all 16 completed. Report in `eval/results/report.md`; the pilot report is kept as `eval/results/report-pilot-20261026.md`.
- Check budget (D8) took C from a median of 8 checks to 5 in this rerun (6 pooled with rerun 2); D stayed at 6. Baseline is 3. Per task, C: 01 8,6; 02 5,4; 03 6,4; 04 4,7.
- Every new house rule took in 8/8: `|>` only, `theme_bw(base_size = 8)`, zero lints. Arm D followed its three new learned rules (`out_dir`, no titles, `options(warn = 1)`) in 8/8; C in 0/8.
- Arm D's rules were replaced before the rerun: the old `%>%` rule contradicted D1 and the old theme rule had become the default. `results.tsv` gained three columns (`native_pipe`, `rule_out_dir`, `rule_warn1`); pilot rows hold NA there.
- Corpus audit: three Opus agents compared the house style with `cc_learnings/final/` and `repos_new/`; a Fable agent then reviewed the 16 proposed additions. Applied: the seed rule corrected (distinct seed per replicate, never one fixed seed in a repeated helper), exact matching and named joins, recode without a catch-all, logged skips, functions and helper namespacing, no blanket warning suppression, no reported numbers from notebooks, shebang and strict shell, shadowing and full argument names, deprecated args, provenance comments, heatmap/stars/labels, string "NA" and trailing-space labels, conda exports without `prefix:`. Dropped on review: replication unit and paired designs (sc-conventions pack owns them), min sizes and alpha (values split across repos), run receipts (mycelium-extra `verify`), per-group tryCatch logging (one repo).
- Citation fixes in `house-style.md`: source data is IHA, IFN, MOCHA only; `py.paths`/`py.args` now cite cns; `here()` is labelled external-only; the corpus is about 8 lineages, not 9-10.
- `rules/core.md` is at 3,990 chars; the R dispatch detail moved to the house style to make room.

## Phase 3 pilot (2026-10-26)
- Run: `eval/results/pilot-20261026/`, 4 tasks × 4 arms × 2 reps, Opus 5.5, all 32 completed. Usage $11.11. Scores in `eval/results/results.tsv`, report in `eval/results/report.md`.
- Isolation: `claude --safe-mode --no-session-persistence` on the user's subscription. `--bare` needs an API key, and that account had no credit. Plugin hooks don't run in an isolated session, so arms B-D get their text through `--append-system-prompt-file` (`eval/arm-prompt.js`).
- Results, n = 8 per arm:
  - Every gate passed in every arm, both traps included. Opus 5.5 avoids pseudoreplication and the misaligned join without help, so these tasks can't show a safety gain.
  - Complete headers: 0/8 in A and B, 8/8 in C and D.
  - Median checks: 3.5 (A), 3 (B), 8 (C), 6 (D). The plugin over-checked; D8 now caps it.
  - Arm D took up all three learned rules; no other arm did.
- LLM judges hit the ceiling: a median of 5/5 for every arm on one-read clarity, and weak agreement (kappa 0.64, 0.23, 0.00). Treat them as uninformative.
- Human packet: `eval/results/pilot-20261026/packet/` (32 blinded scripts plus `sheet.tsv`; the key is in `packet-key.tsv`). It is gitignored, so it exists only on this machine. One arm B script contains a `ponytail:` comment that reveals its arm; it was left unchanged.
- `run.sh --out NAME` lets one process per arm share a results folder, so the arms run in parallel.

## State
- Phases 0-2 and the Phase 3 harness and pilot are committed on `main`. Nothing is pushed.
- Phase 2 added `templates/` (`air.toml`, `.lintr`, `ruff.toml`, `palettes.{R,py}`, `theme.{R,py}`, `.here`, and `examples/01_fit-contrasts.{R,py}`), `skills/align/`, `skills/review/` and `tests/templates.test.js`. Two parallel Opus agents built it.
- Tests: `node --test tests/*.test.js` runs 36 cases, all passing, none skipped. The template tests run air, ruff, lintr, both example scripts and a scratch-repo align simulation. That simulation is the Phase 2 done-when.
- Tools: air 0.12.0 comes from conda-forge (`pixi global install air`), at `~/.pixi/bin/air` with a symlink at `~/.local/bin/air`. Posit's installer refuses glibc 2.28. R is `module load statistical/R/4.5.2/gcc.8.5.0`, which has lintr, ggplot2 and here but not styler or janitor. Python examples run through `uv run --python 3.12 --with matplotlib,pandas,scipy`.
- Phase 2 live checks, 2026-10-26:
  - `/lucid-coder:review` on a deliberately bad R script gave 10 findings grouped by topic, each with `file:line`.
  - `/lucid-coder:review templates/examples` found real bugs in the examples: two non-empty checks that could never fire, unknown condition labels silently becoming `NA`, and R paths not using `here::here()`. All three are fixed. The Python example then gained scipy p-values and BH q-values; they match the R example to floating-point precision.
  - The `allowed-tools` line on `review` lets `rules.js list` run without a prompt (checked in `-p` mode).
- `templates/.here` marks `templates/` as the root for `here::here()` when the examples run in place.
- Phase 1 added `scripts/rules.js` (the only store writer), `skills/learn/`, `skills/rules/` and `tests/{rules,prompt}.test.js`, and extended `hooks/lib.js` and `hooks/prompt.js`. Two parallel Opus agents built it against one CLI contract.
- Tests: `node --test tests/*.test.js` runs 27 cases, all passing.
- Phase 1 live checks, 2026-10-26:
  - A repo-scope rule made with `rules.js propose` then `approve` was injected in a session in this repo, and absent in a session run from `/tmp`.
  - `/lucid-coder:rules list` ran the CLI through the skill.
  - The nudge hook was run against 9 real transcripts. It fired only on the one whose last turn ran Write.
  - The probe store was deleted afterwards.
- Data dir: hooks under `--plugin-dir` get `CLAUDE_PLUGIN_DATA=~/.claude/plugins/data/lucid-coder-inline`. The Bash tool may lack that variable, so SessionStart writes its path to `~/.claude/.lucid-coder-data` and `rules.js` reads it. `${CLAUDE_PLUGIN_ROOT}` substitutes inside skill text (checked live).
- Phase 0 files: `.claude-plugin/plugin.json`, `hooks/{hooks.json,lib.js,activate.js,subagent.js,prompt.js}`, `rules/core.md` (3,869 chars), `tests/hooks.test.js`, `package.json`.
- Tests: `node --test tests/*.test.js` runs 12 cases, all passing. `npm` is not installed on this host; `package.json` holds the same command. `claude plugin validate .` passes, with directory warnings for a missing README, LICENSE and icon.
- Live check 2026-10-26: `claude --plugin-dir . -p ...` quoted the core heading and ladder rung 1. With `~/.claude/.lucid-coder-active` set to `off`, a new session reported NOT FOUND. `-p "/lucid-coder off"` wrote the flag and `-p "/lucid-coder on"` removed it. An Explore subagent quoted ladder rung 1 from its SubagentStart context. This covers level 2 of `PLAN.md` section 8.
- `mattpocock-skills@claude-plugins-official` is installed at user scope. `AGENTS.md` and `docs/agents/` hold its config: GitHub issues, default triage labels.

## Locked decisions
- Readability-first persona; brevity only breaks ties. User brief, 2026-10-26.
- Checks are kept only if they stop silent data loss, a silent mislabel or a silent overwrite. There are five defaults and no `check_*.R` files. User's standing stopifnot rule (`SYNTHESIS.md:170`).
- One house style across repositories. A repository departs from it only through an approved repo-scope rule, and each rule has a topic key.
- Hooks are written in Node, because this host's `python3` is 3.6.8. Following ponytail's pattern.
- The rule store lives in `${CLAUDE_PLUGIN_DATA}/rules.json` (fallback `~/.claude/lucid-coder/`). A repository rules file loads only while its SHA-256 matches one the user trusted. Reason: prevent prompt injection.
- Do not rebuild what mycelium-extra already does: scilintr gate, data-contract-check, verify, manuscript-audit.
- v0.1 bundle subsystems are deferred: URL/GitHub ingest, learn-method, Acumen, evaluate, compare. See `PLAN.md` section 11.
- D1–D7 (`PLAN.md` section 10): `|>` for new repos, `air`, width 100, tripwires off by default, correction nudge on, no v0.1 Python CLI, file names `NN_verb-noun`. User, 2026-10-26.
- Phase 0 deviations from `PLAN.md` section 6, copied from ponytail:
  - The hook command is a string, `node "${CLAUDE_PLUGIN_ROOT}/hooks/x.js"`, not the `command`+`args` form.
  - SessionStart and UserPromptSubmit print raw text. SubagentStart prints `hookSpecificOutput` JSON, because raw stdout is dropped there.
  - The mode flag is `~/.claude/.lucid-coder-active`. Absent means on; content `off` means off.
- External sources (Bryan, Tang, Openscapes, Tidyverse course, PLOS 2021, R Cookbook) are in `skills/lucid-coder/references/external-sources.md`, labelled `ext` in `house-style.md`. User request, 2026-10-26.

## Dead ends — do not redo
- Injecting the full `SKILL.md` from a hook: it is 10.6 KB, and `additionalContext` is capped at 10,000 chars. Inject `rules/core.md` instead.
- A v0.1 Python hook on this host fails with `SyntaxError` on 3.6.8 (`from __future__ import annotations`).
- `node --test tests/` fails on Node 26: it treats the directory as a module. Use `node --test tests/*.test.js`.
- A `python3` heredoc edit is blocked by the mycelium-extra gate hook when a path holds a shell variable. Use the Edit tool.
- Re-reading the docs corpus or the ten external URLs: the findings are already in `PLAN.md`, `house-style.md` and `external-sources.md`.

## Read first
- `PLAN.md`: the verdict on v0.1 (section 1), layout (section 3), learning store and schema (section 5), hooks facts (section 6), evaluation (section 8), phases (section 9).
- `skills/lucid-coder/SKILL.md`: the full persona.
- `rules/core.md`: the injected form. Keep it under 4,000 chars; a test enforces this.
- `skills/lucid-coder/references/house-style.md`: topic keys and defaults.
- `hooks/lib.js`: store read, repo id, trust check, render, escape, budget, data-dir pointer.
- `scripts/rules.js`: run it with no arguments for the usage text. Exit codes: 0 ok, 1 validation error, 2 conflict.
- `/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans-lucidcoder/lucid-coder/`: the old v0.1 bundle. Its plan docs were revised 2026-10-26 at the user's request: `AGENT_IMPLEMENTATION_PROMPT.md` points the build at this repository, and the other docs are marked superseded. Its code is reference only. The original is archived at `../lucid-coder-v0.1-orig.tar.gz`.

## Open items
- The data-dir pointer is one file. If a dev copy (`--plugin-dir`) and an installed copy both run, the last session to start wins, and `rules.js` writes to that copy's store. Not tested with an installed copy.
- No conflict check runs between a trusted repo file's rules and the user's own repo-scope rules on the same topic. Both are injected.
- `/lucid-coder:learn` has not run live end to end, from correction through `propose`, `approve` and `supersede`. The CLI side is unit-tested.
- `rules.js` `parseArgs` rejects a value that starts with `-`, such as `--statement "-x"`. `--statement=-x` works.
- `claude plugin validate .` warns about a missing README, LICENSE and icon. These matter only before publishing.
- `/lucid-coder:align` has not run end to end, because `-p` mode cannot answer its prompts. The copy step is simulated in tests.
- `.lintr` allows `SNAKE_CASE` alongside `snake_case`, so the house style's UPPER_SNAKE constants pass.
