# Roadmap

Status on 2026-11-10 (v0.5.0): **the plugin is in real use; phases 1 to 3 are mostly done.** Open: #3 (second reader) and #4 (learned rules) wait on people. Tracking lives in GitHub issues under the roadmap issue; this file is the plan they come from.

## Where things stand (v0.5.0)

- **Rules.** `rules/core.md` is a persona: a senior computational biologist, five facts about where analysis code fails, eight habits with their why. `rules/house.md` restates the user-decided conventions as flat rules plus a do-not-add list. `rules/r.md` and `rules/python.md` hold language syntax. The session hook injects all of them (about 7.4 KB of the 9.5 KB cap).
- **Evidence.** Blind side-by-side reviews by one reader (the user):
  - lucid-coder vs baseline Opus, tasks 02-04: 8-3, 2 ties;
  - rule-list core vs persona core: 5-4, 2 ties;
  - persona + house vs persona alone: 5-1.
  - Persona + house cut script length by a quarter to a third against the persona alone.
- **Eval harness.** 25 tasks with planted traps and gates (`eval/tasks/`), the R4_51 conda R, notebook rendering via papermill, and blinded review pages built by `eval/side-by-side.js`.

## While in real use

- Correct the agent in the session; say "always" or "never" for a lasting rule, and accept the offer to save it (`/lucid-coder:learn`). Learned rules are the main input when development resumes.
- Note any habit the agent overdoes (extra knobs, audit prints, extra figures) or underdoes, with the file it happened in.

## Phases

Phase 1, measure what exists. Done except two items:
1. Done (#2, 2026-11-02): tasks 01, 05, 06 and the Python twins under persona + house. D won 5 of 6 panels.
2. Open (#3): a second reader judges the existing review pages. The packet is `eval/results/second-reader.zip`; compare with `eval/compare-readers.js`.
3. Open (#4): fold learned rules from real use into `house.md` or the persona.
4. Done (#5, 2026-11-03): see "Decisions recorded 2026-11-03".
5. Done (#21): tut counts from the kai notebook corpus added to house-style rows.

Phase 2, harness fixes. Done:
- Two-arm review pages balance the arm sides (`a8ed619`, #6).
- Python house rules in `rules/python.md` (#7).
- Python env for eval tasks: the uv lock `eval/envs/python-tasks.txt` pins decoupler 2.3.1, pydeseq2 0.5.4, liana 1.10.0 and `pertpy[jax,de]` 1.4.0.

Phase 3, new eval tasks. Done: all ten catalog tasks are built as 07 to 16, plus 07b and 15b (`eval/TASK-CATALOG.md`). Arms A and D passed the gates of 07 to 16 almost everywhere, so the gates sit at the ceiling. Next: harder variants with fewer hints, then the catalog's next tier.

Phase 4, skill evidence from the kai notebook corpus (#19): mine how 2,845 scverse tutorials call each pipeline skill's method, then propose gotchas and grill-question values for the user to approve (#20).

## Decisions recorded 2026-11-10

From the Fable review of the package and the exactplot comparison:
- Hardening landed: CI with `claude plugin validate --strict` and a skip-count gate; consistency tests for version and doc drift; time-boxed SessionStart scans; `trust` needs `--sha256`; machine paths replaced by env vars; `meta.json` records the Claude Code version, repo SHA, arm-prompt hash and scored-script hash; MIT `LICENSE`.
- Figure helpers: font hierarchy 8/6/5/10 pt (in-plot labels drop from about 7 to 6 pt), page widths 89 and 183 mm, `compose_figure()`, `get_legend()`, fixed 10 mm embedding arrows in R. Python keeps arrows scaled to the data range. LaTeX layout is out: it clashes with `cairo_pdf`.
- The work email stays in old `HANDOFF.md` history; no rewrite.
- CI skip baseline is 0 until the first CI run sets it.
- `claude plugin eval` port dropped (2026-11-10). The eval sandbox needs `socat`, and its Bash tool cannot write to the run's working directory (`mkdir: Read-only file system`), so cases that run `Rscript` and write `results/run1/` cannot pass. Five ported cases (05, 06, 07, 11, 14) scored 0.09 to 0.86, all from that limit. The graders matched correct output when the script ran outside the sandbox. Plugin eval could still check behaviour (a skill fired, a tool ran); the repo's own harness keeps the output gates.
- Not done: BixBench-style harder tasks; `nbstripout` (conflicts with notebooks as the report).

## Decisions recorded 2026-11-03

From #5 (phase 1 open questions):
- One `verify_<claim>.R` is allowed per claim that spans two output files; the `check_*.R` ban stays for everything else.
- R joins state `unmatched = "error"` next to `relationship =` when every row must find a partner.
- Paths use one string: `here("data/x.tsv")`.

From #26 (mycelium-extra integration), written into `skills/lucid-coder/SKILL.md`, "Inside a Mycelium repo":
- Check budget: ask per repo, save the answer as a repo-scoped learned rule.
- Layout: in Mycelium repos follow `outputs/run_<date>_<slug>/` and `source_data/<panel>.tsv`; `results/` and `.csv` elsewhere.
- `.living/conventions.md`: one repo-scoped learned rule per convention, no blanket override.
- Receipts: tested 2026-11-03 on mycelium-extra 0.16.1. The gate ignores scripts outside the repo (no grant, no receipt), filed as mycelium-extra#5. Until it is fixed, run pipeline skills through a repo wrapper named in the plan.
