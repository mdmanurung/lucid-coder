# lucid-coder plot bench

An [acumen](https://github.com/scverse/acumen) project for the lucid-coder plot/style skill. `plots/tasks.yaml` holds six corpus figure procedures on Kang 2018 PBMCs (volcano, composition bar, marker dotplot, box with jitter, capped heatmap, and an R volcano). Each task has train, valid and test splits. The source notebook sits in a comment above each task.

## What it measures

acumen grades `answer.md` by exact match after `strip()`, so each task ends in one token: a gene count, a cell type, a gene or a label. That score shows whether the agent ran the specified procedure. It says nothing about the figure. `box` is binary (`stim`/`ctrl`), so a guess scores 50%.

`score_style.py` gives a separate style score from the `fig.*` rows of `house-style.md`: no title, subtitle or caption; cut-off lines on volcanoes; a seed wherever jitter appears; a pdf+png+csv triplet; the plot templates imported or sourced; an mm figure size. A run keeps only `answer.md`, `script.py` and its transcript, so the scorer reads `script.py` plus every Write, Edit and Bash call in `transcript.jsonl`. That is how it sees R code. Figures are not kept; the file check runs only on a run dir that holds them.

## Build the skill

`python build_skill.py plots` copies `skills/lucid-coder/` and `templates/{plots,theme,palettes}.{py,R}` into `plots/skills/v1/`. It also appends a note on importing the templates. acumen reads `skills/v1/SKILL.md` directly, so the files sit in `v1/`, not in `v1/lucid-coder/`. The script refuses to overwrite v1; delete it to rebuild. `fit` keeps an existing v1, so v1 is our skill and epoch 2 writes v2.

## Run

```bash
cd plots
acumen prepare                  # data cache in plots/cache/ (read-only)
acumen check                    # reproducers, then a cheap review agent
acumen fit --max-epochs 3       # costs money: baseline, v1, then improved versions
acumen report                   # report.html
python ../score_style.py        # style table per arm; per-run CSV in runs/style.csv
python ../score_style.py --selftest
```

**Cost caps:** `config.yaml` sets one model (`claude-sonnet-5`), 3 replicates, `max_usd: 1.0` and `max_turns: 30` per run, and 4 runs at once. With `--max-epochs 3`, a worst-case fit is about $15–30.

## R and data

The target `../target` is a stub package. It pins anndata, decoupler, matplotlib, numpy, pandas, pydeseq2, scanpy and scipy, and puts an `Rscript` shim on the agent PATH. The shim runs host R 4.5.2 with DESeq2, ggplot2 and scattermore. Its paths are baked in for this host; override them with `LUCID_R_HOME` and `LUCID_R_LIBS_USER`. Reproducers call plain `Rscript`, because `acumen check` puts the venv `bin` first on PATH.

`common/prepare/` holds the data scripts, and `plots/prepare` links to it. Two entries go into `plots/cache/`:
- `kang2018.h5ad`: 24,673 cells, raw counts. obs `label` (`ctrl`, `stim`), `replicate` (8 patients, e.g. `patient_101`) and `cell_type` (`B cells`, `CD14+ Monocytes`, `CD4 T cells`, `CD8 T cells`, `Dendritic cells`, `FCGR3A+ Monocytes`, `Megakaryocytes`, `NK cells`). decoupler 2.x has no Kang loader and figshare serves a bot challenge, so the script copies a local pertpy `kang_2018.h5ad` (`LUCID_KANG_H5AD`). It falls back to the figshare URL only when that file is missing.
- `pbmc3k.h5ad`: `scanpy.datasets.pbmc3k()` raw counts. No task uses it yet.
