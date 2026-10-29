---
name: align
description: >-
  Install the lucid-coder house-style configs into this repository, one file
  at a time, after showing each diff and asking. Use when the user says
  "align this repo", "install the house style configs", "set up air/ruff/lintr",
  or "add the palette file".
argument-hint: "[r|python|all]"
---

# Align a repository with the house style

Copy the templates in `${CLAUDE_PLUGIN_ROOT}/templates/` to the repository
root, only where the user says yes. The argument limits the language;
with no argument, do both.

Templates:

- R: `air.toml`, `.lintr`, `palettes.R`, `theme.R`, `plots.R`
- Python: `ruff.toml`, `palettes.py`, `theme.py`, `plots.py`
- Both: `lucid-coder.yaml`, filled in by step 4 rather than copied as is

## Steps

1. Detect languages. A language applies when the argument allows it and
   the repository holds its files:

   ```bash
   find . -path ./.git -prune -o \( -name '*.R' -o -name '*.Rmd' -o -name '*.qmd' \) -print -quit
   find . -path ./.git -prune -o \( -name '*.py' -o -name '*.ipynb' \) -print -quit
   ```

   If neither prints a path, say so and stop.
2. For each applicable template, compare it with the repository root.
   - Absent: `cat "${CLAUDE_PLUGIN_ROOT}/templates/<file>"` and show it.
   - Present: `diff <file> "${CLAUDE_PLUGIN_ROOT}/templates/<file>"` and show
     the diff. An empty diff means the file is already aligned; say so and
     move on.
   - If `pyproject.toml` has a `[tool.ruff]` section, show it and ask
     whether to move it into `ruff.toml` or keep it. Ruff prefers `ruff.toml`
     when both sit in the same directory.
3. Ask once per file: "Copy `<file>`?" Stop and wait for the answer. On yes:

   ```bash
   cp "${CLAUDE_PLUGIN_ROOT}/templates/<file>" ./<file>
   ```

   A `palettes.R` or `palettes.py` that exists and differs from the template
   holds project colours. Never copy over it. Offer to append the
   `palette_for` helper from the template instead, and only on yes.
4. Offer `lucid-coder.yaml`, the env pointer the SessionStart hook reads. It
   names the env only, never packages. If it exists, show it and ask whether
   it is still right. Otherwise ask which env manager (pixi, conda, uv or
   system) and which env name or prefix. Activate that env, then fill `r`
   and `python` from:

   ```bash
   command -v Rscript
   command -v python
   ```

   Each path must sit under the env prefix; if one does not, the env is not
   active, so say so and ask. Leave a field out rather than guess it. Show
   the filled file and ask before writing it to the repository root. A
   malformed file is ignored by the hook with a one-line warning.
5. Count files that use the magrittr pipe:

   ```bash
   grep -rl '%>%' --include='*.R' --include='*.Rmd' --include='*.qmd' . | wc -l
   ```

   If the count is above 0, say how many files. The house style is `|>`
   always, so offer to convert them as a separate change. Do not convert
   here, because `air` does not rewrite pipes.
6. Report what was copied and what was skipped. Then list the commands to run
   next, for the languages that apply:

   ```bash
   air format .
   Rscript -e 'lintr::lint_dir()'
   ruff format .
   ```

   Say that folders and scaffolding belong to mycelium-extra's repo-template
   skill, not to this one.

   Mention analysis contracts too, without writing one: each analysis gets
   `contracts/<analysis>.yaml` (from `templates/contract.yaml`), agreed before
   any code. A value of `ask` is a question the grill still has to settle.

   ```bash
   node "${CLAUDE_PLUGIN_ROOT}/scripts/contract.js" validate contracts/<analysis>.yaml
   node "${CLAUDE_PLUGIN_ROOT}/scripts/contract.js" asks contracts/<analysis>.yaml
   node "${CLAUDE_PLUGIN_ROOT}/scripts/contract.js" approve contracts/<analysis>.yaml --by <name>
   node "${CLAUDE_PLUGIN_ROOT}/scripts/contract.js" status contracts/<analysis>.yaml
   node "${CLAUDE_PLUGIN_ROOT}/scripts/contract.js" gates contracts/<analysis>.yaml
   node "${CLAUDE_PLUGIN_ROOT}/scripts/contract.js" card <skill>
   ```

   A contract's `skill:` links the method card `skills/<skill>/card.yaml`;
   its `engine:` links `skills/<skill>/variants/<engine>.yaml`. Every card
   `requires` entry needs a `columns` entry in the contract.

   In a mycelium-extra repo, `contract.js mycelium <file> --table <path>
   [--sample-col <col>]` prints a `mycelium-extra.data_contract.v1` JSON for
   `data-contract-check` to stdout. Pipe it in; never save it in the repo
   unasked.

## Never

- Never write outside the repository root, and never change anything under
  `${CLAUDE_PLUGIN_ROOT}`.
- Never copy or overwrite a file without a yes for that file.
- Never replace a `.lintr`, `air.toml` or `ruff.toml` that has different
  values without showing the diff first.
- Never run `air format`, `ruff format` or `lintr` on the user's code without
  asking.
