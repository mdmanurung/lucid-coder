# lucid-coder: revised plan

Written 2026-10-26. This plan replaces
`mycelium-extra/plans-lucidcoder/lucid-coder/` (the "v0.1 bundle"). That bundle
is unchanged. The drafted persona is in `skills/lucid-coder/SKILL.md`, and the
house style is in `skills/lucid-coder/references/house-style.md`.

**Goal.** Build a ponytail-like Claude Code plugin with three properties:

- Agents write bioinformatics and statistics R/Python that humans can read on one pass.
- Style stays consistent across repositories.
- The persona keeps only checks that matter, and it learns from the user's corrections.

Token efficiency is not a goal. Readability is the goal.

**Status.** Nothing here is built. Every claim about agent behaviour is
**not assessed** until section 8 runs.

## 1. Verdict on the v0.1 bundle

The bundle is careful about scientific claims. But it optimises the wrong
property, spreads across too many subsystems, and its learning store is unsafe
to inject. Findings, most severe first:

1. **The persona is concise-first, not readable-first.**
   - The README headline is "concise first-pass R/Python", and the skill says "simplest complete implementation".
   - `comparison.py:19` scores lines, characters, functions, classes and imports. All of these reward brevity.
   - No setting, rubric or metric measures readability.
2. **The check policy invites checkpoint sprawl.**
   - `skills/lucid-coder/SKILL.md` lists ten guardrails and `references/scilintr-prevention.md` adds twelve rows.
   - Neither says when a check earns its place.
   - The corpus has that rule (`cc_learnings/final/code_smells_and_checks.md:142`), along with the user's standing instruction to keep `stopifnot` "to the utmost essential" (`conventions_synthesis/SYNTHESIS.md:170`).
3. **No house style, so no cross-repository consistency.**
   - The bundle says "honor the repository's local idioms" (`skills/lucid-coder/SKILL.md:35`) and ranks local conventions above global style (`docs/LEARNING_AND_CONFLICTS.md:35-41`).
   - That reproduces the drift the corpus documents: 16 of 17 repositories hard-code paths, and pipe, assignment, indent and palette handling all differ.
4. **The learning store is a prompt-injection path.**
   - `hooks/inject.py:27` injects any profile marked `"status": "approved"` from `cwd/.lucid-coder/`.
   - A cloned repository can ship that file. Nothing ties "approved" to the user.
   - Rule text is not escaped. A newline plus `#` opens a fake heading (reproduced by the code review).
5. **The hook fails on this host and fails closed.**
   - `hooks.json` runs bare `python3`, which is 3.6.8 here. `from __future__ import annotations` then raises `SyntaxError`, so no guidance loads.
   - On 3.11, one malformed profile (`KeyError`, `TypeError`, `AttributeError`) crashes the hook and removes the core skill too (`inject.py:28` catches only `OSError, ValueError`).
6. **The learning state machine is mostly prose.**
   - The docs describe `observed -> ... -> superseded/archived`, but the code implements only propose, approve and activate.
   - There is no reject or supersede. Two contradictory profiles can both be active.
   - `approve` flips every rule in a profile at once, so the user cannot approve a subset.
   - An over-budget profile is dropped whole with no notice (`store.py:181`).
7. **Scope covers five subsystems the brief does not need.** They are:
   - URL and GitHub ingestion
   - `learn-method`
   - the Acumen bridge
   - `evaluate`
   - the HTML `compare` report

   They carry the defect classes the corpus already found: SSRF (`0x7f.1` and `.localhost` pass `sources.py:69`; the redirect check runs after the request), git credential-helper use, symlink following and loose `float()` comparison (`'0123' == '123'`).
8. **Packaging defects.**
   - Tests fail without `PYTHONPATH=src` (`ModuleNotFoundError`).
   - `scripts/smoke-test.sh` and `scripts/lucid-coder` are mode 644 (`Permission denied`).
   - The full skill is both auto-discovered and injected, so it can appear twice.
   - There is no off switch.

What the bundle gets right, and this plan keeps:

- imported text is data, not instructions
- nothing activates without approval
- scientific correctness outranks style
- evidence labels on every claim
- no silent edits to `AGENTS.md`

## 2. What lucid-coder is now

lucid-coder has three parts. A persona sets how scientific code reads. A house
style keeps that the same across repositories. A small rule store lets the user
correct both.

lucid-coder **owns**:

- the readability ladder, script layout and comment policy
- the check policy: which checks earn their place
- the house style for R, Python, files and figures, with a topic key per rule
- learned rules from explicit user corrections, with approval and supersession
- repository config templates (`air.toml`, `.lintr`, `ruff.toml`) that make the house style mechanical

lucid-coder **uses, and does not rebuild**:

| Need | Owner | Evidence |
|---|---|---|
| scilintr rules (about 27 Python, 44 R) | mycelium-extra `card_lint`, `lint_before_run`, `lint_on_write` | `07_audit_round_2.md:122-131` |
| Data contracts, sample alignment, h5ad layers | mycelium-extra `data-contract-check` | `INTEGRATED_SPEC.md:29-35` |
| Run receipts, provenance, plan fidelity | mycelium-extra `verify` | same |
| Manuscript figure map and table validator | mycelium-extra `manuscript-audit` | same |
| Single-cell method guidance | `sc-conventions-revised` pack, 19 guides | `SYNTHESIS.md:233` |
| Project decisions and conventions memory | Mycelium `.living/` | `INTEGRATED_SPEC.md:27` |
| Coding restraint | ponytail, when installed alongside | |

The persona points at these tools when they are installed. It does not restate
their rules. When ponytail is also active, its ladder decides *whether* code
should exist. lucid-coder decides how the code that does exist *reads*.

**Non-goals:**

- a validator, gate, approval or status system; mycelium-extra owns these (`RUNBOOK.md:177`)
- fetching URLs or cloning repositories; the user clones, and lucid-coder reads local files
- method extraction from tutorials
- benchmark optimisation, Acumen, or an HTML diff report
- a repository scaffold; mycelium-extra's planned repo-template skill owns that (`SYNTHESIS.md:283`)

## 3. Plugin layout

The layout copies ponytail's working pattern:

- Node hooks, fail-open
- a mode flag and `/lucid-coder off`
- SubagentStart injection with an agent-type matcher
- skills for the commands

```text
lucid-coder/
├── .claude-plugin/plugin.json   # name, version, description; no "hooks" key (default path loads)
├── hooks/
│   ├── hooks.json               # SessionStart, SubagentStart, UserPromptSubmit
│   ├── activate.js              # SessionStart: core ruleset + active rules
│   ├── subagent.js              # SubagentStart: same text, filtered by LUCID_SUBAGENT_MATCHER
│   ├── prompt.js                # UserPromptSubmit: mode commands + correction nudge
│   └── lib.js                   # store read, render, escape, budget
├── rules/core.md                # injected ruleset, hard cap 4,000 characters
├── skills/
│   ├── lucid-coder/SKILL.md     # full persona, loaded on demand
│   │   └── references/house-style.md
│   ├── learn/SKILL.md           # /lucid-coder:learn: turn a correction into a proposed rule
│   ├── rules/SKILL.md           # /lucid-coder:rules: list, approve, reject, supersede, retire
│   ├── review/SKILL.md          # /lucid-coder:review: read a diff against house style + rules
│   └── align/SKILL.md           # /lucid-coder:align: install formatter/linter configs into a repo
├── scripts/rules.js             # the only writer of the store; called by the skills
├── templates/                   # air.toml, .lintr, ruff.toml, palettes.R, palettes.py, theme.R
├── tests/                       # node --test
└── README.md
```

Why Node:

- Python 3.6.8 is the default `python3` here, and the corpus shows a declared floor drifting from the real interpreter (`07_audit_round_2.md:144`).
- Node avoids the problem, and ponytail proves it works across hosts.
- One runtime serves the hooks and the store writer.

Other hosts (Codex, OpenCode) are **deferred** until the Claude Code version
passes section 8. Ponytail's `writeHookOutput` shows the per-host output shapes
to copy then.

## 4. Persona and house style

### 4.1 Rules for the persona

The full persona is `skills/lucid-coder/SKILL.md` (8.5 KB). The injected
`rules/core.md` is its compact form, because hook context is capped at 10,000
characters per `additionalContext` (`code.claude.com/docs/en/hooks`). Core rules
plus active learned rules must fit with margin.

The persona makes five changes from v0.1:

1. **The readability ladder replaces the brevity ladder.** Its rungs are:
   - inline single-use logic
   - name intermediates after what they hold
   - make decisions visible as named constants with reasons
   - loop over repeated blocks
   - add a helper file only when two or more scripts call it
   - use the established package

   Brevity only breaks ties.
2. **A fixed script layout.** Every script has a header (Purpose, Inputs with producer or `[gap]`, Outputs, Run, Notes), five folding sections in fixed order, and a run directory argument. There is no `setwd()`. Evidence: immgenT is the corpus's best header (`coding_conventions_R_python.md:76`). Seven of 17 repositories have no usable header (`:94`).
3. **A comment policy.** Comments say why: at a decision, a number's origin, a reference level, a non-obvious API need, a known limit.
   - Banned: comments that restate code, results pasted into comments, tutorial voice, and commented-out code.
   - Stale or wrong comments appear in all 17 repositories (`code_smells_and_checks.md:136`). So the persona updates the comment in the same edit as the code.
4. **The check policy.** A check earns its place only if it stops silent data loss, a silent mislabel or a silent overwrite (`code_smells_and_checks.md:142`). The default five are: new run directory, inputs and columns at load, join cardinality, ID identity before a positional copy, outputs non-empty.
   - The persona names what not to write.
   - Python raises exceptions instead of `assert`.
   - Messages name the object and the fix.
   - This follows the corpus: 7 of 17 repositories have zero checks, and the useful checks cluster at joins and loads (`code_smells_and_checks.md:31,61`).
5. **Figures by default.** One named palette file, one theme set once, sizes in mm, `cairo_pdf`, no `dpi` on vectors, source data per panel, file named by panel.

### 4.2 House style

`references/house-style.md` pins one default per topic, with evidence labels.
It counts 17 repository names as about 9-10 independent lineages. A repository
may depart from the house style only through an approved exception (section 5).
"Local idiom" is no longer a default.

The mechanical half ships as configs: `/lucid-coder:align` copies `air.toml`,
`.lintr` and `ruff.toml` into a repository, after showing the diff and asking.
A formatter makes style identical across repositories without spending
persona context on whitespace.

## 5. Learning from feedback

### 5.1 Flow

```text
user correction ──> applied now, in this task
      │ lasting? ("always", "never", "in this lab")
      ▼
offer once: /lucid-coder:learn ──> proposed rule (shown in full)
      ▼
user approves ──> active ──> injected next session and into subagents
      │
      ├─ a later rule on the same topic + scope ──> old rule becomes superseded
      └─ /lucid-coder:rules retire <id> ──> retired (kept, not injected)
```

The `UserPromptSubmit` hook does not store prompts. It matches a short pattern
list ("always", "never", "don't use", "prefer", "instead of") on prompts that
follow a code edit. On a match it adds one line of context: "If this is a
lasting style correction, offer `/lucid-coder:learn` once." This is a nudge,
not a capture. False positives cost one offer, and the user can switch the
nudge off.

### 5.2 Store

The store has three properties:

- **Location:** `${CLAUDE_PLUGIN_DATA}/rules.json`, in the user's home. It survives plugin updates.
- **No project file can inject rules.** This closes finding 4. A repository that wants team rules commits `lucid-coder.rules.json`, and `/lucid-coder:rules trust` records its SHA-256 in the user store. If the file changes, it stops loading until the user trusts it again.
- **One writer:** `scripts/rules.js`. It validates the store, writes atomically to a temporary file and then renames, and keeps a `.bak` copy.

Rule schema:

```json
{
  "id": "r-0007",
  "topic": "r.pipe",
  "scope": {"kind": "repo", "repo": "github.com/lab/atlas"},
  "statement": "Use %>% in this repository; it predates base R 4.1.",
  "why": "Mixed pipes in one codebase slow review.",
  "good": "counts %>% filter(n_cells >= MIN_CELLS)",
  "bad": "counts |> filter(n_cells >= MIN_CELLS)",
  "status": "active",
  "supersedes": null,
  "source": "user-correction",
  "created": "2026-10-26",
  "approved": "2026-10-26"
}
```

Field rules:

- `topic` is required. Each rule uses a house-style key or a new `user.*` key.
- `scope.kind` is one of `global`, `language`, `repo`, `file-glob`. The other scope key is `language` (`r` or `python`), `repo` (the origin remote as `host/owner/repo`, else the repo root path) or `glob`.
- `status` is one of `proposed`, `active`, `rejected`, `superseded`, `retired`.
- `source` is one of `user-correction`, `user-preference`, `repo-review`. A rule from `repo-review` was proposed by `/lucid-coder:review` reading local code, and it still needs approval.

### 5.3 Conflict and supersession, mechanically

- At most one `active` rule per `(topic, scope)`. Activating a second rule asks to supersede the first. The old rule keeps its record with status `superseded`, and the new rule's `supersedes` points back to it.
- Precedence across scopes is `file-glob > repo > language > global > house style`. This applies only within a topic, and it is below scientific correctness and the current request.
- Recency never resolves a conflict on its own. The user chooses.
- Approval is per rule, never per profile.

### 5.4 What is never stored

The store never holds:

- raw prompts or transcripts
- data values
- file contents
- people's names
- paths outside the repository

`rules.js` rejects statements over 300 characters, statements with newlines,
and the common secret patterns. Those patterns include `Bearer`, `AKIA` and
`secret=`, which the v0.1 filter missed. Rendering escapes rule text to one
line with no leading `#`.

### 5.5 Injection budget

The core rules get 4,000 characters. Learned rules get 4,000 characters, which
leaves 2,000 of margin under the 10,000 cap. Over budget, the hook injects the
narrowest-scope rules first. It then adds one line: "N rules not shown; run
`/lucid-coder:rules`." Nothing is dropped silently.

## 6. Hooks

Facts below come from `code.claude.com/docs/en/hooks` and `plugins-reference`,
checked 2026-10-26.

| Event | Matcher | Does | Output |
|---|---|---|---|
| SessionStart | `startup\|resume\|clear\|compact` | inject core + active rules unless mode is `off` | `hookSpecificOutput.additionalContext` |
| SubagentStart | agent type, via `LUCID_SUBAGENT_MATCHER` (unset = all) | same text; SessionStart context does not reach subagents | same |
| UserPromptSubmit | — | handle `/lucid-coder on\|off\|status`; correction nudge | same |

Hook rules:

- Use the exec form (`"command": "node"`, `"args": ["${CLAUDE_PLUGIN_ROOT}/hooks/activate.js"]`).
- Set a 5 s timeout.
- Every hook wraps its body in try/catch and always emits the core ruleset, even when the store is unreadable. One bad rule must not remove the persona (finding 5).
- Write no state to `${CLAUDE_PLUGIN_ROOT}`, which changes on update.

PostToolUse lint feedback is **deferred**. mycelium-extra's `lint_on_write`
already does it where installed, and a second linter hook would double the messages.

## 7. Commands (skills)

| Command | Reads | Writes |
|---|---|---|
| `/lucid-coder on\|off\|status` | mode flag | mode flag |
| `/lucid-coder:learn "<correction>"` | store | one `proposed` rule, after showing it |
| `/lucid-coder:rules [list\|approve\|reject\|supersede\|retire\|trust]` | store | store, one rule at a time |
| `/lucid-coder:review [path\|diff]` | code, house style, rules | nothing. Report: findings by topic, each with file:line and a fix |
| `/lucid-coder:align` | repo configs | `air.toml`, `.lintr`, `ruff.toml`, after showing the diff and asking |

`review` may propose rules ("this repository uses `%>%` in 40 files; add a
repo exception?"). It never activates them. When a change is a refactor, it
asks whether outputs were compared before and after (Stoudt et al. 2021).

## 8. Evaluation: does the persona change what agents write?

Claims are made by level, and passing one level never proves the next:

1. Hooks parse and emit.
2. The text reaches the session.
3. Output changes.
4. Humans find it more readable.

The current evidence is none of these.

**Task set.** Use eight tasks, each as both an R and a Python variant where it fits:

- pseudobulk differential expression
- metadata join and alignment
- scRNA QC table
- composition with a donor denominator
- ggplot/matplotlib figure with palette
- refactor a messy legacy script without changing results
- add a seeded clustering step
- fix a broken join

Two of the tasks are scientific negative controls: cell-level pseudoreplication
bait and a misaligned metadata order.

**Arms.** Run baseline, ponytail alone, lucid-coder alone, and lucid-coder with
three learned rules. Use the same model and budget, with 3 runs per task and arm.

**Measures.** Each is reported separately; they are never averaged into one score.

- **Hard gates:** the code runs, and the scientific negative controls are caught. A failed gate is reported, never offset by readability.
- **Mechanical style:**
  - `air format --check` / `ruff format --check` / `lintr` pass rate under the house configs
  - header present and complete
  - count of `stopifnot`/`assert`/`raise` per script, and the share that guard one of the default five
  - commented-out code lines
- **Learned-rule uptake:** the share of runs that follow each active learned rule.
- **Human readability:** two blinded readers per script answer three questions:
  - "Can you say what this script does after one read?"
  - "Can you find where each decision is made?"
  - "Would you change any check?"

  Report agreement between the readers.

The evaluation harness is a shell loop plus a TSV of results. It is not a
framework.

## 9. Build phases

Each phase ends with a check that fails if the phase is broken.

| Phase | Work | Time | Done when |
|---|---|---|---|
| 0 | Manifest, `hooks.json`, `activate.js`, `subagent.js`, `rules/core.md`, mode flag | half a day | `claude plugin validate .` passes; `node --test` covers: corrupt store still emits core, `off` emits nothing, output under 10,000 chars; live session shows the rules (manual, step 2 of section 8) |
| 1 | `rules.js`, `learn`, `rules` skills, `prompt.js` nudge, trust-by-hash | 1-2 days | tests: second active rule on one topic is refused without supersede; repo file with changed hash does not load; newline/`#`/secret statements rejected; budget notice appears |
| 2 | `templates/` configs and palette/theme files, `align`, `review` | 1 day | `align` on a scratch repo then `air format --check` and `ruff format --check` pass on the templates' own example scripts |
| 3 | Evaluation run (section 8) | 2-3 days including human reading | results TSV plus a one-page report labelled by level |
| 4 | Codex and OpenCode adapters, only if phase 3 shows an effect | 1 day | per-host live check |

## 10. Decisions

All seven decided 2026-10-26: the user accepted every recommendation.

| ID | Decision | Decided |
|---|---|---|
| D1 | Pipe: `\|>` or `%>%` for new repositories | `\|>`: base R, no dependency; existing repositories keep theirs via a repo rule |
| D2 | R formatter: `air` or `styler` | `air`: fast, one config file, already used in immgenT's newer code (`coding_conventions_R_python.md:157`) |
| D3 | R line width: 80 (air default) or 100 | 100, per the corpus recommendation; set in `air.toml` |
| D4 | Hard-coded count tripwires (immgenT `107L`) | off by default, matching your stopifnot rule; allowed when you ask |
| D5 | Correction nudge on by default | on; one offer per correction, never repeated |
| D6 | Keep v0.1's Python CLI | no; reuse only `store.py`'s ideas in `rules.js` |
| D7 | File names: `NN_verb_noun` or `NN_verb-noun` | `NN_verb-noun`: `_` between fields, `-` between words, per Bryan's naming sequel |

## 11. Deferred from v0.1, with the condition to revisit

| Item | Revisit when |
|---|---|
| URL/GitHub ingestion | never inside the plugin; the user clones and `review` reads local files |
| `learn-method`, analysis-skill template | a method skill is needed that the sc-conventions pack does not cover |
| Acumen bridge, `evaluate` | section 8 shows an effect worth optimising |
| HTML `compare` report | a reviewer asks for side-by-side diffs |
| PostToolUse lint hook | mycelium-extra is not installed and users want inline lint |

## 12. Evidence used

Primary sources, per your direction:

- `cc_learnings/final/`: coding, smells, plotting, analysis and the cross-repo comparison
- `cc_learnings/repos_summary/`: all 17 repositories

Secondary sources:

- `conventions_synthesis_2026-10-26/`: SYNTHESIS, 02, 05, 07
- `mycelium_integrated_implementation_bundle/docs/`
- `mycelium_reproducibility_execution_2026-10-26/phase_A/`
- the manuscript-starter
- a live run of the v0.1 tests, smoke test and hook on this host

External sources, added at the user's request on 2026-10-26. Their rules and
quotes are in `skills/lucid-coder/references/external-sources.md`:

- Jenny Bryan, "Code smells and feels": <https://github.com/jennybc/code-smells-and-feels>
- Jenny Bryan, "How to name files: the sequel": <https://speakerdeck.com/jennybc/how-to-name-files-the-sequel>
- Ming Tang, reproducible computing:
  <https://divingintogeneticsandgenomics.com/files/slides/2025-06-25_short_good_enough_practice_reproducible_computing.pdf>,
  <https://divingintogeneticsandgenomics.com/files/slides/2022-07-09_ISCB_reproducible.pdf>,
  <https://divingintogeneticsandgenomics.com/post/repro-comp/>
- Openscapes: <https://openscapes.github.io/series/how-do-i/refactor-modularize-code.html>,
  <https://openscapes.github.io/series/core-lessons/coding-strategies.html#project-oriented-workflows>
- Wright et al., "Tidyverse Skills for Data Science": <https://jhudatascience.org/tidyversecourse/>
- Stoudt et al., PLOS Comput Biol 2021: <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008770>
- Long and Teetor, "R Cookbook", 2nd edition: <https://rc2e.com/>

The corpus is a static audit. Zero published figures were reproduced, and
every corpus rule is advice unless labelled otherwise.
