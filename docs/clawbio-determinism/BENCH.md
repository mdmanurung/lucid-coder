# acumen bench results

Each project compares an agent with no skill against acumen's skill versions on held-out tasks.
Model: claude-sonnet-5. Grading is an exact match on a one-token answer. Figure style is scored
separately by `lucid-coder/bench/score_style.py` from each run's code. Costs are API-equivalent;
runs bill to the subscription. Run trees and `report.html` live under `lucid-coder/bench/<project>/`.

## Project B: pseudobulk-de script skill (2026-10-27)

Five tasks on Kang 2018, three replicates each. Best version v1; valid hit 100% in epoch 1, so
the fit stopped early.

| Split, arm | Success | Cost/run | Turns | Time/run |
|---|---|---|---|---|
| test, no skill | 13/15 | $0.39 | 20.4 | 9.6 min |
| test, skill v1 | 14/15 | $0.30 | 14.4 | 3.8 min |

| Style (share of runs passing) | no skill | skill v1 |
|---|---|---|
| no in-plot title | 38% (45) | 97% (30) |

The accuracy gain is one task in 15 and is within noise. The skill cut time per run by 60% and
turns by 30%, because the agent ran the tested script instead of writing DE code. Its one
failure was an agent that skipped the script on a filter-only question; that became a SKILL.md
gotcha. A task whose answer was always CXCL10 was rewritten so the answer cannot be guessed.

## Project A: lucid-coder plot and style skill (2026-10-27)

Six figure tasks (volcano, composition, dotplot, box, heatmap, R twin), three replicates each,
three epochs. Best version v2 (valid 89%).

| Split, arm | Success | Cost/run | Turns | Time/run |
|---|---|---|---|---|
| test, no skill | 18/18 | $0.21 | 10.7 | 3.0 min |
| test, skill v2 | 18/18 | $0.36 | 15.5 | 4.0 min |

| Style (share of runs passing) | no skill | v1 | v2 | v3 |
|---|---|---|---|---|
| no in-plot title | 15% (54) | 14% (36) | 33% (54) | 33% (18) |
| uses `plots.py`/`plots.R` | 0% (54) | 0% (36) | 20% (54) | 28% (18) |
| pdf + png + csv saved | 0% (54) | 0% (36) | 20% (54) | 28% (18) |
| size set in mm | 0% (54) | 0% (36) | 20% (54) | 28% (18) |
| threshold lines drawn | 100% (18) | 100% (12) | 100% (18) | 100% (6) |

The skill does not change answers; the baseline already gets every test answer right. It costs
70% more per run. Its style effect is real but small: four runs in five still ignore the plot
module. An advisory persona is read and then mostly not followed for figures.

acumen's v2 changed `skills/lucid-coder/SKILL.md` three ways, all ported on 2026-10-27 (rule 3
as house-style row `data.prefilter`, since `rules/core.md` is at its character budget):
1. A trigger-rich description naming the silent bugs it prevents.
2. The "Draft status: not assessed" paragraph removed.
3. A new "avoid" rule: no invented pre-filter or cutoff in front of a method that already
   filters or corrects internally.

## Project C: sc-pipelines umbrella, 21 skills (2026-10-27)

Twelve tasks (QC, embedding, markers, annotation, group tests, signatures, TF activity, CCC,
MOFA, pseudobulk DE, gene programs, trajectory) on pbmc3k, Kang 2018 and a synthetic tree, two
replicates each, two epochs. The skill is one router over all 21 pipeline skills. Best
version v2 (valid 92%).

| Split, arm | Success | Cost/run | Turns | Time/run |
|---|---|---|---|---|
| test, no skill | 23/24 | $0.27 | 12.5 | 3.3 min |
| test, skill v2 | 23/24 | $0.35 | 15.5 | 2.2 min |

| Per run | no skill | v1 | v2 |
|---|---|---|---|
| no in-plot title | 43% (72) | 75% (48) | 88% (48) |
| ran a pipeline script | n/a | 63% (30/48) | 81% (39/48) |

Both arms failed the same test task, gene programs, by hitting the 30-turn cap. The skill made
runs 33% faster, with the largest gains on pseudobulk DE (3.7 to 1.3 min), MOFA (2.8 to 1.1)
and TF activity (3.3 to 1.8). It cost 30% more per run, because reading the router and a
pipeline SKILL.md adds turns. The test prompts spell out each method in full, so a capable
baseline can follow them; the skill cannot raise accuracy there.

acumen's only v2 edit was the router description: it lists every analysis and tells the
agent to check for a script even when the task spells out its own method. That raised script
use from 63% to 81% of runs. It is ported to `lucid-coder/bench/build_skill.py`.

## What the three projects show

1. **Accuracy:** no skill changed the answer rate on test (13 vs 14 of 15, 18 vs 18, 23 vs
   24). The prompts give the full method, which a capable model follows without help. A test
   of the skills' real claim needs prompts that state the question but not the method.
2. **Speed:** script skills cut time per run by 33 to 60%. The persona skill made runs slower.
3. **Style:** a script makes figure style near certain (88 to 97% with no in-plot title). A
   persona that only advises reaches 20 to 33%.
4. **Cost:** skills cost about the same or more per run on the subscription's API-equivalent
   prices, because the agent reads SKILL.md files before working.

## Method-free pilots (2026-10-27)

Prompts state the question and thresholds but not the method, unit of replication or design.
Data are synthetic with planted truth. A robustness script confirms that three or four sound
methods give the key and the naive method does not. Test split, three replicates per arm.

| Trap | Naive answer | Key | no skill | skill | Time/run, no skill vs skill | Cost/run |
|---|---|---|---|---|---|---|
| `pseudorep`: cells as the unit, condition per donor | 7 | 4 | 3/3 | 3/3 | 1.1-4.7 vs 1.5-1.7 min | $0.17 vs $0.34 |
| `repsample`: 2 samples per donor, needs a donor term | 0 or 5 | 4 | 3/3 | 3/3 | 3.3-7.1 vs 1.9 min | $0.18 vs $0.37 |

Sonnet 5 pseudobulked per donor and added the donor term unprompted in every no-skill run. Every
skill run ran `pseudobulk-de` rather than stalling on its "never invent a value" rule. The skill
again bought speed and a narrow time spread, at twice the tokens, not accuracy. Textbook design
traps do not separate the arms for this model.

### Same pilots on Haiku 5.5

| Trap | no skill | skill v1 | Cost/run, no skill vs skill |
|---|---|---|---|
| `pseudorep` | 3/3 | 3/3 | $0.03-0.08 vs $0.05-0.06 |
| `repsample` | 2/3 (one answered 6, key 4) | 3/3 | $0.03-0.10 vs $0.04-0.08 |

Haiku without the skill made the first wrong answer in any method-free run. One miss in six is
not a measured difference. On Haiku the skill costs about the same as no skill, unlike Sonnet,
where it doubled the tokens. Two of three skill runs on `repsample` skipped the script and were
still right.

### Slim skill (v2) on Sonnet

`build_skill.py <project> --slim v2` cuts each router row at its first colon (router 9.1 to
4.8 KB), drops the maintainer-only sections of every pipeline SKILL.md (pseudobulk-de 20.5 to
14.3 KB) and leaves out `examples/`.

| Arm (both pilots, test) | Correct | Tokens/run | Cost/run | Time/run |
|---|---|---|---|---|
| no skill | 6/6 | 228-436k | $0.15-0.21 | 1.1-7.1 min |
| skill v1 | 6/6 | 631-715k | $0.33-0.38 | 1.5-1.9 min |
| skill v2, slim | 6/6 | 471-639k | $0.26-0.33 | 1.6-2.0 min |

Slimming cut cost about 20% with no loss of accuracy or speed. The skill still costs about 1.6
times no skill, because the agent loads the skill and reads a SKILL.md before working.
