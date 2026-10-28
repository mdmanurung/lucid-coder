# lucid-coder

A Claude Code plugin for readable, reproducible R and Python in computational
biology. It gives Claude one house style for analysis code and figures, 21
tested pipeline skills for single-cell and statistics work, and rules learned
from your own corrections.

**Docs:** <https://mdmanurung.github.io/lucid-coder/>

## Install

In Claude Code:

```
/plugin marketplace add mdmanurung/lucid-coder
/plugin install lucid-coder@lucid-coder
```

Start a new session. The rules load at session start and into every subagent.

## What is inside

- **Persona and rules** (`rules/`, `skills/lucid-coder/`): headers, knobs at
  the top, checks where data enters, seeded steps, donors as the unit.
- **Learned rules** (`/lucid-coder:learn`, `/lucid-coder:rules`): a lasting
  correction becomes a rule after you approve it.
- **Figure style** (`skills/figure-style/`, `templates/`): one look in ggplot2,
  ComplexHeatmap and plotnine with plotnine-extra.
- **Pipeline skills** (`skills/<name>/`): sc-qc, sc-normalize, sc-embed,
  sc-annotate, sc-markers, pseudobulk-de, sc-de, variance-partition,
  mashr-sharing, composition, milo-da, enrichment, signature-scoring,
  tf-activity, gene-programs, ccc, trajectory, velocity, mofa, group-tests and
  longitudinal. Each is a tested script with no analysis defaults, fixed
  errors and a fixed output tree.

## Commands

| Command | Effect |
|---|---|
| `/lucid-coder on`, `off`, `status` | Turn the persona on or off |
| `/lucid-coder lang r` (`python`, `both`, `auto`) | Choose the language rules |
| `/lucid-coder:review` | Review code against the house style |
| `/lucid-coder:align` | Install the house configs into a repo |

## Develop

```bash
npm test                          # hooks, rules, templates
quarto preview site               # docs site
```

Each pipeline skill's test command is in its `references/maintainer.md`. See
[Develop](https://mdmanurung.github.io/lucid-coder/develop.html) for the full
guide.
