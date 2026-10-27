# lucid-coder evaluation harness

This harness tests whether lucid-coder changes what agents write (PLAN.md section 8).
It is a shell loop plus a TSV, not a framework.
Every measure is reported on its own; nothing is averaged into one score.

## Arms

All arms run `claude --safe-mode --no-session-persistence` with the same model, budget and tools, from a fresh git repo.
The tools are `Read`, `Write`, `Edit` and `Bash` limited to `Rscript`, `uv`, `ls`, `mkdir`, `cat`, `head`, `air` and `ruff`.
The rules ask for `air format` or `ruff format` after every edit, so run.sh copies `templates/air.toml` and `templates/ruff.toml` into every repo before the run.

| arm | what the model gets |
|---|---|
| A | nothing extra (baseline) |
| B | ponytail's `full` instructions |
| C | lucid-coder's core ruleset, rendered by `hooks/lib.js` |
| D | C plus the three learned rules in `arms/rules-d.json` |

Arms B to D go in through `--append-system-prompt-file`.
Step 0 (2026-10-26) showed that plugin hooks do not run in an isolated session, `--plugin-dir` ones included.
`node eval/arm-prompt.js <A|B|C|D>` prints the exact text; run.sh saves it as `<run_dir>/arm-prompt.md`.

## Auth and isolation

Runs use your Claude Code login, so they count against your subscription usage.
`--safe-mode` turns off your `CLAUDE.md`, installed plugins, hooks, skills and MCP servers.
Step 0 confirmed this: a safe-mode session saw none of ponytail, i-have-adhd or the global `CLAUDE.md`.
run.sh and judge.sh drop any exported `ANTHROPIC_API_KEY`, which would otherwise override the login.
An earlier version used `--bare` with an API key; it was dropped because `--bare` accepts no login.

## Tasks

`eval/tasks/0N_*` are the R tasks. `01py_pseudobulk-de`, `02py_metadata-join` and `03py_figure-palette` are Python twins of 01 to 03: the same inputs, truth and `check.R`, and a prompt that asks for one Python script run through `uv run --python 3.12 --with pandas,scipy,matplotlib python`.
Run them with `--tasks 01py,02py,03py`. The agent gets `Bash(uv:*)` for this; each task's `reference/` holds a house-style answer and `reference/naive/` a wrong one.

## Check the commands first

```bash
eval/run.sh --tasks 01,02,03,04 --arms A,B,C,D --reps 2 --dry-run
```

This prints one line per run (32 here) and creates nothing.

## Step 0 smoke test

Four runs, task 03, one per arm, about $8 at most:

```bash
eval/run.sh --tasks 03 --arms A,B,C,D --reps 1
```

## Pilot

```bash
eval/run.sh --tasks 01,02,03,04 --arms A,B,C,D --reps 2
```

Runs go one at a time into `eval/results/<timestamp>/<task>_<arm>_r<rep>/`.
Each holds `repo/` (the agent's working copy), `arm-prompt.md`, `claude.json`, `claude.stderr` and `meta.json`.
`MODEL` (default `claude-opus-5-5`) and `BUDGET` (default 2 USD per run) override the defaults and land in meta.json.

**Cost cap:** 32 runs × $2 = $64 for the pilot, plus the judges below.

## Score

```bash
node eval/score.js eval/results/<timestamp>
```

This appends one row per run to `eval/results/results.tsv` and skips runs already there.
The script is the first `scripts/*.R` or `scripts/*.py` (never `00_legacy.*`); the language follows from the task id, so there is no language column.
It reruns each script on a clean copy of the repo and calls the task's `check.R` for the gate: `Rscript` with a 120 s limit for `.R`, `uv run --python 3.12 --with pandas,scipy,matplotlib python` with a 600 s limit for `.py` (uv may fetch packages on its first call).
For `.R` it then copies `templates/air.toml` and `templates/.lintr` into the repo and runs `air format --check` and lintr.
For `.py` it copies `templates/ruff.toml` and runs `ruff format --check` (column `air_ok`, reported as "format check pass"), `ruff check` (`lint_count`) and `ruff check --select ERA001` (`commented_code_lines`).
`header_complete` needs Purpose, Inputs, Outputs and Run in the first 20 lines; Notes is optional, as in `rules/core.md`.
`join_guard_count` counts `relationship =` (R) or `validate=` (Python) outside comment lines.
It is the last column, so rows scored before it existed are one cell short; report.js reads the missing cell as NA, and score.js rewrites an old header when it appends.
The R-only columns `uptake_pipe`, `uptake_theme`, `native_pipe` and `rule_warn1` are NA for `.py`; `check_count` counts `raise`, `sys.exit(` and `assert` there.

## LLM judges (proxy)

```bash
eval/judge.sh eval/results/<timestamp>
```

Two judges (`claude-opus-5-5`, `claude-sonnet-5-5`) score each script once, under neutral ids, with tools off.
The cost is 2 calls per script, each capped by `JUDGE_BUDGET` (default 0.5 USD): at most $32 for 32 scripts.
Outputs: `judges.tsv`, the private `judge-key.tsv`, and raw replies in `judge-raw/`.

## Human packet

```bash
node eval/make-packet.js eval/results/<timestamp>
```

This writes `packet/` (blinded `P01.R` or `P01.py`, ... keeping each script's extension, and an empty `sheet.tsv`) and the private `packet-key.tsv` beside it.
Send each reader the `packet/` folder only.
Each reader fills a copy of `sheet.tsv`: q1 and q2 are 1 to 5, q3 is yes or no.

## Report

```bash
node eval/report.js eval/results/<timestamp> --human reader1.tsv --human reader2.tsv [--step0 notes.txt]
```

This writes `eval/results/report.md` with sections L1 to L4 and a separate LLM judge section.
Without `--human` sheets, L4 says pending.
Agreement is weighted Cohen's kappa: quadratic weights for q1 and q2, unweighted for q3.

## Tests

```bash
node --test eval/tests/*.test.js
```

`harness.test.js` needs no API key; it uses stub `claude` binaries and a mock task.
`python.test.js` runs the Python twins' reference and naive scripts through uv and `check.R`, then scores a reference run; it skips when uv or the R module is missing.

## Rubric

`rubric.md` defines how to judge a script objectively. Every item is a count or a yes/no with a decision rule, evidence by line number, two blinded readers and no total score. Part 1 is quality; Part 2 is lucid-coder conformance and is reported separately. Each item says whether it comes from `score.js`, a regex or a reader.

## ClawBio-style task (01sk)

`tasks/01sk_pseudobulk-de/` asks for the task 01 analysis packaged as a ClawBio skill. The agent writes `skills/pseudobulk-de/SKILL.md`, with frontmatter and the 11 required sections, an R script for data prep and model fitting, and `report.ipynb` (R kernel `ir`) for figures and the report. The script runs as `Rscript skills/pseudobulk-de/pseudobulk_de.R --input data --output <dir>` and writes `de.tsv`, `result.json` and `reproducibility/commands.sh`. Papermill runs the notebook with `-p out_dir <dir>`, and `jupyter nbconvert --to markdown --no-input` turns it into `report.md`. The gate checks the science, these files, and a `parameters` cell in the notebook. `score.js` reruns only the script.

All R runs use the R4_51 conda env (R 4.5.1 with tidyverse, limma, edgeR, IRkernel, papermill and pandoc); set `R_ROOT` to override. That env ships an old Claude Code, so `run.sh` resolves `claude` before it changes PATH.

## Side-by-side review

```bash
node eval/side-by-side.js eval/results/<run-dir> eval/results/<run-dir>/review --rep 1 --human <panels.json>
```

This writes one page with a panel per analysis step: setup, load, filter and aggregate, fit, write, the report notebook and the SKILL.md Workflow section. Each panel shows the four arms as columns, shuffled per panel, so the reader compares one step at a time instead of scoring whole files in sequence. The reader picks a column, or no preference, and gives a reason. `side-by-side-key.tsv` maps each panel column to its run and stays out of git. Optional human panels show one step from four corpus repos; the repo names appear after a pick is saved.
