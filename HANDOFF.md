# Handoff — lucid-coder

**Date:** 2026-10-26 · **Branch:** main, not pushed · **Status:** Phases 0, 1 and 2 committed and verified; Phase 3 not started

## Goal
Build a ponytail-like Claude Code plugin with three properties:

- It steers agents to readable, reproducible R/Python for bioinformatics and statistics.
- Style stays consistent across repositories.
- It keeps only checks that matter, and it learns from user corrections.

Readability is the goal, not token efficiency. Done means phases 0–3 of `PLAN.md` section 9 pass their checks.

## Next action
Plan Phase 3, the evaluation in `PLAN.md` section 8.

## State
- Phases 0, 1 and 2 are committed on `main`. Nothing is pushed.
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
