---
name: review
description: >-
  Review code against the lucid-coder house style and the user's active
  learned rules. Reports findings by topic with file:line and a fix; writes
  nothing. Use when the user says "review this against the house style",
  "lucid review", "check my script's readability", or "does this follow our
  conventions".
argument-hint: "[path|diff]"
allowed-tools: Bash(node "${CLAUDE_PLUGIN_ROOT}/scripts/rules.js" list:*), Bash(git diff:*), Bash(git ls-files:*)
---

# Review code against the house style

Read-only. The report is the only output.

## What to read

1. The code. With no argument or `diff`, review the working tree change:

   ```bash
   git diff HEAD
   git ls-files --others --exclude-standard
   ```

   Read each untracked file the second command lists. With a path, review
   that file or every R and Python file under that directory.
2. The defaults and topic keys in
   `${CLAUDE_PLUGIN_ROOT}/skills/lucid-coder/references/house-style.md`.
3. The user's active rules:

   ```bash
   node "${CLAUDE_PLUGIN_ROOT}/scripts/rules.js" list --status active --json
   ```

   An active rule replaces the default for its topic. The narrowest scope wins.

## Report

Group findings by topic key. Each finding has three parts: `file:line`, the
default or rule id it breaks, and a one-line fix.

Order the groups by consequence:

1. Scientific risk: `checks.*`, `data.*`, `repro.*`, `stats.*`.
2. Layout, comments and control flow: `files.*`, `comments.*`, `code.*`.
3. Style: `r.*`, `py.*`, `fig.*`.

Anything a formatter fixes (whitespace, quotes, line breaks, indentation) is
one line: "run `air format .`" or "run `ruff format .`". It is never a finding.

Report at most the 10 most consequential findings. Say how many more there
are and which topics they fall under.

### Regression line

If the change is a refactor, meaning behaviour is meant to be unchanged, ask
whether outputs were compared before and after. If a `results/<run_dir>/`
exists, name the command:

```bash
diff -r results/<previous_run> results/<new_run>
```

### Repo-wide patterns

A pattern across many files (`=` assignment, 4-space indent, inline palettes
in N files) is one proposed-rule suggestion, not N findings. Give the exact
command for the user to run or decline, with the file count in `--why`. Wrap
values in single quotes; write a single quote inside a value as `'\''`.
`%>%` is the exception: the house style is `|>` always, so report it as one
`r.pipe` finding with the file count and never propose a rule that keeps it.
Example for `=` assignment in 40 files:

```bash
node "${CLAUDE_PLUGIN_ROOT}/scripts/rules.js" propose --topic 'r.assign' --scope repo \
  --statement 'Keep = for assignment in this repository.' \
  --why 'repo-review: = assignment in 40 files; switching to <- would touch all of them' \
  --source repo-review
```

### mycelium-extra

scilintr is a separate install; mycelium's `analyze` and mycelium-extra's
`verify` run it when it is present. If it is installed (`scilintr` on `PATH`
for Python, the R package `scilintr` for R), do not
repeat its findings (bare `except`, `setwd()`, unseeded randomness, global
warning suppression). Write one line instead: "scilintr covers X; run
`scilintr` or `/mycelium-extra:verify`."

## Never

- Never edit a file.
- Never run `propose`, `approve` or `trust`. Show the command and stop.
- Never claim numerical equivalence from a diff you did not run.
- Never report a finding without a `file:line`.
