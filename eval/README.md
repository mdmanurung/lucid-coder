# lucid-coder evaluation harness

This harness tests whether lucid-coder changes what agents write (PLAN.md section 8).
It is a shell loop plus a TSV, not a framework.
Every measure is reported on its own; nothing is averaged into one score.

## Arms

All arms run `claude --safe-mode --no-session-persistence` with the same model, budget and tools, from a fresh git repo.

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
It reruns each script on a clean copy of the repo (120 s limit) and calls the task's `check.R` for the gate.
It then copies `templates/air.toml` and `templates/.lintr` into the repo and runs `air format --check` and lintr.

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

This writes `packet/` (blinded `P01.R`, ... and an empty `sheet.tsv`) and the private `packet-key.tsv` beside it.
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
