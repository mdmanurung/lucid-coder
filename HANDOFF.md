# Handoff — lucid-coder

**Date:** 2026-10-26 · **Branch:** main, not pushed · **Status:** Phase 0 done and verified; Phase 1 not started

## Goal
Build a ponytail-like Claude Code plugin with three properties:

- It steers agents to readable, reproducible R/Python for bioinformatics and statistics.
- Style stays consistent across repositories.
- It keeps only checks that matter, and it learns from user corrections.

Readability is the goal, not token efficiency. Done means phases 0–3 of `PLAN.md` section 9 pass their checks.

## Next action
Start Phase 1 (`PLAN.md` section 9): `scripts/rules.js`, the `learn` and `rules` skills, the correction nudge in `hooks/prompt.js`, and trust-by-hash for repository rule files.

## State
- Phase 0 is committed on `main` in five commits. Nothing is pushed.
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
- `hooks/lib.js`: store read, render, escape, budget. Phase 1's `rules.js` should reuse `escapeRule`.
- `/exports/para-lipg-hpc/mdmanurung/mycelium-extra/plans-lucidcoder/lucid-coder/`: the old v0.1 bundle. Its plan docs were revised 2026-10-26 at the user's request: `AGENT_IMPLEMENTATION_PROMPT.md` points the build at this repository, and the other docs are marked superseded. Its code is reference only. The original is archived at `../lucid-coder-v0.1-orig.tar.gz`.

## Open items
- Before Phase 1 fixes the store path, confirm `CLAUDE_PLUGIN_DATA` is set in the hook environment under both `--plugin-dir` and an installed plugin. If only one sets it, the hooks and `rules.js` could read different files. Ponytail never uses this variable, so there is no prior evidence.
- `claude plugin validate .` warns about a missing README, LICENSE and icon. These matter only before publishing.
- Phase 2: `review` should compare outputs before and after a refactor (PLOS 2021 regression-test rule). Add this to `PLAN.md` section 7 when Phase 2 starts.
