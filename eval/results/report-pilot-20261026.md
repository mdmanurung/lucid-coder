# lucid-coder eval: pilot-20261026

Model claude-opus-5-5, budget $2 per run. Each measure is reported on its own; there is no combined score, and a failed gate is never offset by style.

## L1: the harness ran

32 runs found, 32 scored. claude exited non-zero in 0. No usable claude.json (cost NA) in 0. No script written in 0. Total cost $11.11.

## L2: the text reached the session

Arms delivered by --append-system-prompt-file in --safe-mode sessions, because plugin hooks do not run in an isolated session (step 0, 2026-10-26).
Channel recorded in meta.json: none, append-system-prompt-file.

## L3: the output changed

n per task × arm cell = 2; n per arm = A: 8, B: 8, C: 8, D: 8. Rate cells show k/n with NA excluded.

| measure | A baseline | B ponytail | C lucid-coder | D lucid-coder + 3 rules |
|---|---|---|---|---|
| n runs | 8 | 8 | 8 | 8 |
| gate pass | 8/8 | 8/8 | 8/8 | 8/8 |
| script found | 8/8 | 8/8 | 8/8 | 8/8 |
| rerun ok | 8/8 | 8/8 | 8/8 | 8/8 |
| air format pass | 0/8 | 0/8 | 1/8 | 2/8 |
| lint count, median | 0 | 0.5 | 2 | 1 |
| header complete | 0/8 | 0/8 | 8/8 | 8/8 |
| check count, median | 3.5 | 3 | 8 | 6 |
| commented-code lines, median | 0 | 0 | 0 | 0 |
| uptake: %>% only | 8/8 | 7/8 | 0/8 | 8/8 |
| uptake: no plot title | 1/2 | 2/2 | 2/2 | 2/2 |
| uptake: theme_bw(base_size = 8) | 0/2 | 0/2 | 0/2 | 2/2 |

Gate pass by task (n per task × arm cell = 2):

| task | A baseline | B ponytail | C lucid-coder | D lucid-coder + 3 rules |
|---|---|---|---|---|
| 01_pseudobulk-de | 2/2 | 2/2 | 2/2 | 2/2 |
| 02_metadata-join | 2/2 | 2/2 | 2/2 | 2/2 |
| 03_figure-palette | 2/2 | 2/2 | 2/2 | 2/2 |
| 04_refactor-legacy | 2/2 | 2/2 | 2/2 | 2/2 |

Every gate passed.

## L4: humans find it more readable

Pending: no filled reader sheets yet.

## LLM judges (proxy, not human evidence)

Two blinded model judges; n per task × arm cell = 2; n per arm = A: 8, B: 8, C: 8, D: 8.

**claude-opus-5-5**

| measure | A baseline | B ponytail | C lucid-coder | D lucid-coder + 3 rules |
|---|---|---|---|---|
| n scripts | 8 | 8 | 8 | 8 |
| q1 one read, median (1-5) | 5 | 5 | 5 | 5 |
| q2 find decisions, median (1-5) | 4 | 4 | 4 | 4 |
| q3 would change a check (yes) | 8/8 | 8/8 | 8/8 | 8/8 |

**claude-sonnet-5-5**

| measure | A baseline | B ponytail | C lucid-coder | D lucid-coder + 3 rules |
|---|---|---|---|---|
| n scripts | 8 | 8 | 8 | 8 |
| q1 one read, median (1-5) | 5 | 5 | 5 | 5 |
| q2 find decisions, median (1-5) | 4 | 4 | 5 | 5 |
| q3 would change a check (yes) | 7/8 | 7/8 | 6/8 | 7/8 |

Agreement, claude-opus-5-5 vs claude-sonnet-5-5 (weighted Cohen's kappa; quadratic weights for q1 and q2, unweighted for q3): q1 0.64 (n = 32), q2 0.23 (n = 32), q3 0.00 (n = 32).
