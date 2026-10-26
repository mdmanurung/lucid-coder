# lucid-coder eval: rerun-cd-20261026

Model claude-opus-5-5, budget $2 per run. Each measure is reported on its own; there is no combined score, and a failed gate is never offset by style.

## L1: the harness ran

16 runs found, 16 scored. claude exited non-zero in 0. No usable claude.json (cost NA) in 0. No script written in 0. Total cost $6.59.

## L2: the text reached the session

Arms delivered by --append-system-prompt-file in --safe-mode sessions, because plugin hooks do not run in an isolated session (step 0, 2026-10-26).
Channel recorded in meta.json: append-system-prompt-file.

## L3: the output changed

n per task × arm cell = 2; n per arm = C: 8, D: 8. Rate cells show k/n with NA excluded.

| measure | C lucid-coder | D lucid-coder + 3 rules |
|---|---|---|
| n runs | 8 | 8 |
| gate pass | 8/8 | 8/8 |
| script found | 8/8 | 8/8 |
| rerun ok | 8/8 | 8/8 |
| air format pass | 2/8 | 2/8 |
| lint count, median | 0 | 0.5 |
| header complete | 8/8 | 8/8 |
| check count, median | 5.5 | 6.5 |
| commented-code lines, median | 0 | 0 |
| uptake: %>% only | 0/8 | 0/8 |
| uptake: no plot title | 2/2 | 2/2 |
| uptake: theme_bw(base_size = 8) | 2/2 | 2/2 |
| house: |> only (pipe used) | 8/8 | 8/8 |
| D rule: out_dir variable | 0/8 | 8/8 |
| D rule: options(warn = 1) | 0/8 | 8/8 |

Gate pass by task (n per task × arm cell = 2):

| task | C lucid-coder | D lucid-coder + 3 rules |
|---|---|---|
| 01_pseudobulk-de | 2/2 | 2/2 |
| 02_metadata-join | 2/2 | 2/2 |
| 03_figure-palette | 2/2 | 2/2 |
| 04_refactor-legacy | 2/2 | 2/2 |

Every gate passed.

## L4: humans find it more readable

Pending: no filled reader sheets yet.

## LLM judges (proxy, not human evidence)

Not run.
