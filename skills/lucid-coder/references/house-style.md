# House style: the default every repository starts from

Each rule has a **topic key**. A learned rule or a repository exception
names the topic it replaces, so two active rules can never disagree on
one topic (see `PLAN.md` section 5).

Evidence labels: **obs** means seen in the 17-repo corpus. **rec** means the
corpus recommends it. **new** means no repo did it and we are deciding.
**ext** means an external source recommends it; see `external-sources.md`.
The corpus is 17 repository names but about 8 independent lineages.
Four are AIFI, two are REFLEX versions of one codebase, three are Haniffa
lab (skin_fib, spatial, cns), and dreamlet and variancePartition share an
author. Counts below are repository names, not independent witnesses.

Sources are abbreviated:

- `C` = `cc_learnings/final/coding_conventions_R_python.md`
- `S` = `cc_learnings/final/code_smells_and_checks.md`
- `P` = `cc_learnings/final/plotting_conventions.md`
- `RS/<repo>` = `cc_learnings/repos_summary/<repo>.md`
- `A` = `cc_learnings/final/analysis_conventions.md`
- `X` = `cc_learnings/final/cross_repo_comparison.md`
- `02` = `cc_learnings/02_coding_R_python.md`
- `SYNTHESIS.md` = `conventions_synthesis_2026-10-26/SYNTHESIS.md`
- `PF/<repo>/shard_NN` = `cc_learnings/per_file/<repo>/shard_NN.md`

## R

| Topic | Default | Evidence |
|---|---|---|
| `r.assign` | `<-`, never `=` at top level | obs clean in immgenT, aifi, sound-life; mixed in dreamlet, MOCHA, IHA, MM, aging (`RS/dreamlet_analysis:53`, `RS/MOCHA_Manuscript:69`, `RS/IHA-Figure:130`, `RS/multiple-myeloma:67`, `RS/single_cell_aging:66`) |
| `r.indent` | 2 spaces | obs split: 2 in immgenT, MM, aging (`C:155`, `RS/multiple-myeloma:67`, `RS/single_cell_aging:66`); 4 in IFN, reflex and the function bodies of aifi and sound-life (`RS/IFN_manuscript:86`, `PF/reflex_manuscript/shard_01:99`, `RS/aifi-healthy-pbmc-reference:125`, `RS/sound-life-scrna-analysis:84`); mixed in ALTRA, MOCHA |
| `r.pipe` | `\|>` always, never `%>%`, in new and existing repositories; a new line after every `\|>`, one step per line | user rule 2026-10-26 (D1 revised in `PLAN.md`). `%>%` appears in 16 of 17 (`C:154`); `\|>` is base R 4.1+ and needs no dependency |
| `r.tidyverse` | tidyverse first for data frames: dplyr, tidyr, purrr, readr, stringr, forcats; never `data.table`, base loops, `apply()` or `merge()` on data frames; base R only for matrices and modelling APIs (`limma`, `colSums()`) | user rule 2026-10-27; corpus mixes data.table (MOCHA `.N`), reshape2 (MM) and base `apply`/`expand.grid` (immgenT, dreamlet), so corpus base-R idioms are translated: `expand_grid()`, `set_names()`, `group_split()`/`map()`, `count() \|> filter()` |
| `r.r4ds` | R4DS 2e idioms: keep the input name, save results under a new name; `count(x, sort = TRUE)`; `slice_max(x, n = 1, with_ties = FALSE)`; `pivot_longer(names_to =, values_to =, values_drop_na = TRUE)`, `".value"` with `names_sep`; `join_by()`, never `by = c("a" = "x")`; `anti_join()` to find unmatched keys; `unmatched = "error"` on inner joins; `fct(x, levels = )` (errors on an unknown level, `factor()` makes it NA); `fct_infreq() \|> fct_rev()`, `fct_recode(new = "old")`, `fct_lump_n()`; `if_else()` over `ifelse()`; name logical conditions in `mutate()` before `filter()`; `near()` for doubles, `is.na()` never `== NA`; `separate_wider_delim()` over `strsplit()`, `too_few = "debug"`; `str_c()` keeps NA, `str_glue()` writes "NA" (labels and paths only); raw strings `r"{...}"`; `complete()`, `count(.drop = FALSE)`, `na_if()` at read; verb function names, `{{ var }}` and `pick()`; `\(x)` lambdas; `set_names(basename) \|> map() \|> list_rbind(names_to =)`; `walk2()` for writes; matrix code `drop = FALSE`, `vapply()`; models `lm(y ~ x, data = df)` then `coef()`, `confint()`, `broom::augment()`. Not adopted: `possibly(..., NULL)` (drops failed files silently) | ext R4DS (`data-transform`, `data-tidy`, `joins`, `factors`, `logicals`, `strings`, `regexps`, `missing-values`, `functions`, `iteration`, `base-R` chapters at r4ds.hadley.nz); R Cookbook 2e recipes 6.7, 11.1, 11.3, 11.15 (rc2e.com). Where they differ, R4DS wins: RC2E uses `%>%`, `TRUE ~`, base strings and `by = c()` |
| `r.quotes` | double quotes | new: repos mix them (`RS/reflex_manuscript:62`) |
| `r.bool` | `TRUE`/`FALSE`, never `T`/`F` | rec (`C:209`) |
| `r.width` | soft limit 100 columns, set as `line-width = 100` in `air.toml`; one argument per line only when a call overflows | decided 2026-10-26 (D3); rec (`C:209`); no repo had a rule (`C:156`); revised 2026-10-27 to match rules/r.md (corpus rules) |
| `r.formatter` | `air` | decided 2026-10-26 (D2). No repo had a formatter config (`cc_learnings/README.md:48`) |
| `r.library` | `library()` at the top, only for what the script uses; `pkg::fn()` for one-off calls; `suppressPackageStartupMessages()` is fine (startup messages are not warnings); no `quiet_library()` wrapper | rec (`02:99`); the wrapper ban is new: the wrapper appears in 5 repos (`C:142`) |
| `r.names` | snake_case objects; parameters at the top in UPPER_SNAKE (user pick 2026-10-27, reversing R4DS lowercase constants), related names share a prefix, plot objects by panel (`p_fig2c`); never `c`, `T`, `F`, `df` or `data`; argument names spelled in full | rec (`02:84`); obs six repos snake_case (`02:77`, `C:65`, `C:968`); shadowing in IFN (`02:81`) and immgenT (`C:1024`); partial matching in 5 repos (`S:122`) |
| `r.helpers` | helper files use `pkg::fn()` only, never `library()` or `require()` inside a function; a `source()` line ends with a comment naming what it provides | obs immgenT (`PF/immgenT-GP-analysis/shard_07:95`, 8 of 12 helpers `C:1032`); anti-patterns GTestimate `require()` in helpers (`C:698`), MOCHA `require()` mid-file (`RS/MOCHA_Manuscript:73`) |
| `r.errors` | scripts check with a one-line named `stopifnot("<what holds>" = cond)`; `stop(sprintf(...), call. = FALSE)` only inside helpers, message prefixed `fn():` | new, a synthesis of MM `call. = FALSE` plus fix-in-message (`RS/multiple-myeloma:79`, `:181`, `S:40`), immgenT `call. = FALSE` (`PF/immgenT-GP-analysis/shard_06:124`) and function prefix (`RS/immgenT-GP-analysis:74`); revised 2026-10-27 to match rules/r.md (corpus rules) |
| `r.dispatch` | `switch()` with `match.arg()` to dispatch on a string; `dplyr::case_when()` over nested `ifelse()`; `%\|\|%` for fallbacks; `is.numeric(x)`, never `class(x)[[1]] == "numeric"`; S3 methods when branching on class | ext: Bryan, code smells |

## Any language

| Topic | Default | Evidence |
|---|---|---|
| `code.conditionals` | guard clauses first: `stop()` or `return()` early, no `else` after an exit; a fussy condition becomes a named predicate (`is_paired_sample()`); a few small functions over one long one; less nesting | ext: Bryan, code smells |
| `code.toggles` | never comment code in and out to change behaviour; use a named argument or constant | ext: Bryan, code smells |
| `code.functions` | a named function only when called 3+ times in a script, or shared by 2+ scripts through a helper file; used once or twice, the code stays inline, and an inline `\(x)` goes in `map()`; inputs as arguments, one value out, no globals, defined before first use (user rule 2026-10-27) | rec (`02:70`, `C:208`); globals in ALTRA, MOCHA, IHA (`C:145`); use before definition in MOCHA (`C:146`) |
| `code.matching` | exact matches: `%in%`, `fixed = TRUE`, anchored patterns; never `grep(paste(genes, collapse = "\|"))`; filter with `!grepl()`, never `x[-grep()]` or `x[-which()]` | obs 7 repos (`S:125`); empty negative index in MM, variancePartition, MOCHA (`S:125`, `RS/variancePartition_analysis:139`, `PF/multiple-myeloma/shard_01:88`) |
| `code.join_by` | `join_by()` in every join plus `relationship = "many-to-one"` (or the true cardinality) instead of a check block; pandas `validate="m:1"` on every merge | `by =`: obs natural joins in 9+ repos (`S:68`), ALTRA `merge` drops (`S:288`); pandas `validate=` is new (`S:565` records its absence); revised 2026-10-27 to match rules/core.md (corpus rules) |
| `code.warnings` | never `suppressWarnings()`, `options(warn = -1)` or `filterwarnings("ignore")`; silence one category with a one-line reason; in Rmd no `opts_chunk$set(warning = FALSE, message = FALSE)` | rec (`02:114`); obs aifi, IHA, sound-life, spatial (`C:378`, `C:936`, `C:1524`, `C:1823`); ALTRA `options(warn = -1)` in 12 notebooks (`RS/ALTRA-manuscript:102`) |
| `code.constants` | values a reader may change (cutoffs, levels, formula, seed) sit at the top; one-use paths and patterns stay inline; UPPER_SNAKE in R and Python (user pick in the 2026-10-27 side-by-side review: the case marks what a reader may change); prose names the value, never repeats it; shared across scripts means one sourced `R/constants.R` | obs drift in 5 lineages (`RS/reflex_manuscript_v2:203`, `RS/IFN_manuscript:194`, `RS/ALTRA-manuscript:185`, `PF/multiple-myeloma/shard_04:213`, `PF/single_cell_aging/shard_04:50`); model immgenT (`RS/immgenT-GP-analysis:76`); revised 2026-10-27 to match rules/core.md (corpus rules) |
| `code.by_name` | select rows, columns and coefficients by name, never by position; assign by key, never positionally after a join | obs IHA (`PF/IHA-Figure/shard_12:9`), MOCHA (`RS/MOCHA_Manuscript:63`), aging (`RS/single_cell_aging:179`), reflex (`RS/reflex_manuscript:72`, `RS/reflex_manuscript_v2:215`) |
| `code.deprecated` | `linewidth` not `size` for lines; no `aes_string()`, `spread()`, `is_categorical_dtype()` | obs 10 repos (`S:130`) |
| `code.shebang` | `#!/usr/bin/env Rscript` or `#!/usr/bin/env python3` on runnable scripts; shell and SLURM wrappers use `set -euo pipefail` and `${1:?usage}` | obs immgenT (`C:981`, `RS/immgenT-GP-analysis:38`); cns minority, 9 of 36 shell files (`RS/cns-haematopoiesis:89`), `${1:?}` once (`:110`); rec (`C:216`) |

## Python

| Topic | Default | Evidence |
|---|---|---|
| `py.formatter` | `ruff format` (black-compatible, 88 columns) | rec: black at 88 (`C:209`); ruff is one tool for format and lint; obs MM pins ruff/black (`PF/multiple-myeloma/shard_03:14`) |
| `py.strings` | f-strings | new: repos mix `.format()` (`RS/sound-life-scrna-analysis:81`) |
| `py.paths` | `pathlib.Path`, `ROOT = Path(__file__).resolve().parents[1]` | new (`C:247-250` shows argparse `Path`, not a root constant) |
| `py.args` | flat top-level code with paths and parameters as a few variables at the top; `argparse` only for a value that varies per job; no `main()` scaffolding | obs cns (`RS/cns-haematopoiesis:47`, `C:109`); revised 2026-10-27 to match rules/python.md (corpus rules) |
| `py.errors` | raise `ValueError`/`KeyError` with the fix in the message; no `assert` outside tests | rec (`SYNTHESIS.md:172`); obs cns, spatial (`RS/cns-haematopoiesis:189`, `RS/spatial_skin_atlas:815`, `:821`) |
| `py.views` | `.copy()` before assigning into a subset's `.obs` or `.var`; never mutate an AnnData view | obs 38 Haniffa shard files (`RS/skin_fibroblast_atlas:175`, `RS/spatial_skin_atlas:1431`); one lineage, but the bug is universal to AnnData |
| `py.imports` | one import block at the top; never re-import mid-notebook | obs: cns re-imports scanpy 3-6 times (`C:160`) |
| `py.seed` | `random_state` as a function argument | rec (`SYNTHESIS.md:146`) |

## Files and layout

| Topic | Default | Evidence |
|---|---|---|
| `files.names` | `NN_verb-noun.R`: `_` between fields, `-` between words; lowercase, no spaces or other punctuation; zero-padded prefix, with gaps such as `90_` for odd-one-out scripts; ISO 8601 `YYYY-MM-DD` first when order is by date | decided 2026-10-26 (D7, `PLAN.md`); ext: Bryan, naming sequel; `SYNTHESIS.md:109` recommends the underscore form `NN_verb_noun`; all nine digested repositories fail either (`SYNTHESIS.md:111`) |
| `files.layout` | `data/` raw and read-only (`chmod u-w -R data/`), `scripts/`, `results/<run_dir>/`; one project per repository; a README at the root that lists any manual step with its exact command | ext: Tang, Openscapes, Tidyverse course; chmod: 0 corpus witnesses |
| `files.paths` | `here::here()` in R; Python: see `py.paths`; every path relative to the project root; save source, never the workspace; each script runs in a fresh process | ext: Openscapes, Tidyverse course, Tang. 0 of 17 repos use `here()` (`C:668`, `C:871`, `RS/reflex_manuscript:38`, `RS/IHA-Figure:70`, `RS/multiple-myeloma:44`); the corpus recommends running from the repo root with a root guard (`C:205`), which immgenT does (`RS/immgenT-GP-analysis:54`) |
| `files.versions` | version through git and run folders; never `_v2`, `_final`, `-Copy1` | obs failures (`SYNTHESIS.md:111`) |
| `files.outputs` | `<run_dir>/` passed as an argument; files named by content or panel | rec (`C:211`); 8 repos had names that drift from panels (`S:120`) |
| `files.helpers` | a helper file only when 2 or more scripts call it; helpers live in `R/<topic>.R` (immgenT, MM) or one function per file | rec (`cc_learnings/02_coding_R_python.md:72`); obs GTestimate (`C:136`) |
| `files.notebooks` | scripts prep data, fit models and write tables; figures and reports go in a notebook (`.ipynb`, R kernel `ir`, or `.qmd`) that reads those tables, never report text pasted together with `writeLines()`/`sprintf()` in a script; QC review and annotation are notebooks too, with markdown headings instead of banner comments; a notebook takes its run dir from a `parameters` cell (papermill); runs top to bottom from a fresh kernel; strip `.ipynb` outputs before commit; no `.ipynb_checkpoints` committed | rec (`phase_A/CONVENTIONS.md`); ext: Tang; revised 2026-10-27 to match rules/core.md (corpus rules) |

## Reproducibility and data

| Topic | Default | Evidence |
|---|---|---|
| `repro.env` | commit a lockfile: `renv.lock`, `uv.lock` or `environment.yml`; export conda environments without the machine `prefix:` line; a container only when the lockfile is not enough; Python session info is `sc.logging.print_versions()` or `session_info.show()` | ext: Tang, PLOS 2021; obs conda exports with `prefix:` and missing packages (`A:74`, `A:82`, `A:263`); session info: Haniffa 15 files (`RS/spatial_skin_atlas:10`) |
| `repro.automation` | a script for every reported number; batch steps are scripts, repeated routines in a shell script or a Snakemake or `{targets}` rule; judgement steps may be notebooks; no number is typed by hand | ext: Tang; rec (`C:203`, `X:145`); hand-typed numbers in GTestimate, immgenT, spatial, MM (`A:222`); revised 2026-10-27 to match rules/core.md (corpus rules) |
| `data.recode` | one named vector drives the recode and `fct(x, levels = )` (never `factor()`); no catch-all arm (`TRUE ~ "Other"`, `.fillna("Other")`), so an unmapped label surfaces as `NA` | obs cns (`S:342`), IHA (`S:601`), MM trailing-space keys (`S:761`), spatial (`S:1154`); revised 2026-10-27 to match rules/r.md (corpus rules) |
| `data.tidy` | long format before analysis; one value per cell; `NA` for unknown, never blank, an undeclared `-999` or the string `"NA"` (`S:169`, `S:601`); `trimws()` label values at load (`S:761`); `janitor::clean_names()` on import; one case-consistent code per category; plain-text data with ISO 8601 dates; the raw-count layer name is one project constant (e.g. `layers["counts"]`), never two spellings; inputs are never overwritten in place | ext: Tidyverse course, Tang; `trimws`: obs MM (`PF/multiple-myeloma/shard_07:187`); `clean_names()` is ext only; layers: obs Haniffa `raw_counts` 31 vs `counts` 19, inputs overwritten in 10 files (`RS/cns-haematopoiesis:62`) |
| `data.prefilter` | no invented pre-filter or cutoff in front of a method that already filters or corrects internally (DESeq2 independent filtering, edgeR `filterByExpr`, a tool's own multiple-testing step): the two compound and shift the result; apply the method once, with the user's threshold, or let it decide | acumen bench 2026-10-27, plots project: skill v2 added it from failed runs (`bench/plots/skills/v2/SKILL.md`); not in rules/core.md, which is at its 3350-character budget |
| `data.write` | `readr::write_tsv()` or `write.csv(row.names = FALSE)`, rownames moved to a column first; never write a file and read it back within one script | obs IHA ~60 `write.csv` then `select(-X)` (`RS/IHA-Figure:88`), MOCHA `read.csv()[,-1]` (`RS/MOCHA_Manuscript:20`, `PF/MOCHA_Manuscript/shard_05:58`), MM (`RS/multiple-myeloma:50`), ALTRA rds handoffs (`RS/ALTRA-manuscript:176`) |
| `data.exclusions` | report removals in one line per filter step or one counts table, with the reason; a loop that skips a group says which and why | obs sound-life logs every QC step and label edit to a table (`PF/sound-life-scrna-analysis/shard_04:41`, `:106`), aifi prints the count removed per filter (`PF/aifi-healthy-pbmc-reference/shard_01:98`); ext: PLOS 2021; silent skips in IFN, aging, spatial (`S:470`, `S:25`, `S:74`); model MM `warning()` + `next` (`S:774`); revised 2026-10-27 to match rules/core.md (corpus rules) |
| `repro.seed` | seed each stochastic step in the script, a distinct seed per replicate; never one fixed seed inside a repeated helper; `RNGkind("L'Ecuyer-CMRG")` or seed in the worker for parallel code; `scvi.settings.seed`; pass `seed =` to `geom_text_repel()` and `position_jitter()` | rec (`A:257`); obs IHA, IFN per-iteration seeds (`A:54-55`); failures: MOCHA fixed seed in a helper (`A:67`, `A:259`), variancePartition `mclapply` (`A:65`), unseeded scVI in cns and spatial (`A:51`, `A:64`); RNGkind: 0 corpus witnesses; `seed =` in geoms: immgenT (`RS/immgenT-GP-analysis:182`) |

## Checks and comments

These topics let the user tune the two things they care about most.
One example correction is "stop checking column presence in this repo".

| Topic | Default | Evidence |
|---|---|---|
| `checks.default_five` | only checks that stop silent data loss, mislabel or overwrite: ID identity before a positional copy; join cardinality via `relationship =` / `validate=`; a new run dir is not overwritten; each one line | rec (`S:142`); ext: PLOS 2021 "data gut checks" after merges; run dir refusal: 0 corpus witnesses; MM `relationship = "many-to-one"` on 6 joins (`RS/multiple-myeloma:84`); revised 2026-10-27 to match rules/core.md (corpus rules) |
| `checks.extra` | only for a known scientific assumption that would fail silently | rec (`S:142`); user rule (`SYNTHESIS.md:170`) |
| `checks.budget` | the defaults plus at most one extra per script; skip what fails loudly anyway (a missing file or column); no input audits, type checks, expected counts or re-derived results; in a notebook a printed `dim()` or crosstab is fine | eval pilot 2026-10-26: without a budget, lucid-coder wrote a median of 8 checks per script against 3 without it; user rule; revised 2026-10-27 to match rules/core.md (corpus rules) |
| `checks.tripwires` | off: no hard-coded expected counts unless asked | decided 2026-10-26 (D4); immgenT uses them deliberately (`RS/immgenT-GP-analysis:110`) |
| `checks.style` | R: one-line named `stopifnot("meta rows match counts" = identical(...))`; `stop(sprintf(...), call. = FALSE)` only in helpers; Python: one statement `if not ok: raise ValueError(f"...")`; `assert` only in tests | new, a synthesis of MM `call. = FALSE` plus fix-in-message (`RS/multiple-myeloma:79`, `:181`, `S:40`), immgenT `call. = FALSE` (`PF/immgenT-GP-analysis/shard_06:124`) and function prefix (`RS/immgenT-GP-analysis:74`); revised 2026-10-27 to match rules/core.md (corpus rules) |
| `checks.files` | no separate `check_*.R` files | user rule (`SYNTHESIS.md:170`) |
| `code.sections` | every block, parameters included, opens with `# Name ` padded with `-` to column 80 | user pick 2026-10-27 (arm A dividers, task 04) |
| `code.writes` | one explicit write call per output file (a helper when 3+ files), never a named list walked with `iwalk()`/`walk2()` | user pick 2026-10-27: "repeated write_output to me is much more intuitive" (task 04) |
| `r.startup` | `suppressPackageStartupMessages({ library(...) })` around the library block | user pick 2026-10-28 (round 2 head-to-head, 02 and 04 setup: "please use suppressPackageStartupMessages") |
| `r.matrix` | a count matrix stays a matrix: filter with `rowSums()`, `colSums()` and `[`; pivot to long form only for a table or a plot | user pick 2026-10-28: a 6-line matrix filter beat a 45-line long-format pipeline ("use available, widely used functions") |
| `fig.plot_frame` | build the plot data frame and any caption text before the `ggplot()` call; the ggplot chain only maps and styles | user pick 2026-10-28: "offload the wrangling outside the plotting ... instead of very lengthy ggplot2 scripts" |
| `comments.steps` | a short comment above each block names the step (`# merge annotations with sample metadata`); a long pipe may label its verbs; a blank line between steps; dense code with no blank lines fails | user picks 2026-10-27: dreamlet `Kfoury_CancerCell_2021_get_data.Rmd:105-124`, sound-life `15-R_assemble_type_frequencies.ipynb` cell 23, ALTRA `RA_06-...CON1_vs_ARI.Rmd:200-214`; rejected a dense ggplot block with no blank lines |
| `docs.workflow` | a SKILL.md or README workflow lists steps keyword first, detail second (`1. **Aggregate**: sum raw counts per donor`) | user pick 2026-10-27 (01sk Workflow panel) |
| `comments.why` | comment decisions, number origins, reference levels, API needs, limits, and the source of a marker list or copied code; a platform UUID or file-id literal names its file; never restate code | obs dreamlet, MM, immgenT (`C:600`, `C:1159`, `C:1032`); provenance in cns, aifi, MOCHA (`S:365`, `C:377`, `S:735`); UUID literals: ALTRA (`RS/ALTRA-manuscript:63`) |
| `comments.header` | Purpose and Run, then Inputs [producer] and Outputs one file per line with what it holds, then Notes; no blank fields, no author line (user pick 2026-10-27: the descriptive file list of arm A, task 04) | fields rec (`C:206`); no author line is new: every corpus header carries one (`RS/ALTRA-manuscript:39`, `RS/immgenT-GP-analysis:37`, `RS/MOCHA_Manuscript:26`) |
| `comments.dead_code` | delete commented-out code; no parked code behind `knit_exit()` or a bare `STOP` cell | rec (`C:210`); 0 corpus witnesses; parked code in dreamlet (`PF/dreamlet_analysis/shard_01:32`), skin_fib (`RS/skin_fibroblast_atlas:38`) |

## Figures

| Topic | Default | Evidence |
|---|---|---|
| `fig.palette` | one `palettes.R` or `palettes.py` of named maps keyed by label; unknown label fails | rec (`P:449-464`); 12 repos define colours inline (`P:344`); immgenT stops on unknown label in `level2_group_palette.R` (`RS/immgenT-GP-analysis:102`) but falls back to grey elsewhere (`:116`); MM returns silent NA (`RS/multiple-myeloma:58`) |
| `fig.theme` | `theme_set(theme_bw(base_size = 8))` once per script | user rule 2026-10-26 (D9); obs: `theme_classic` in 12 of 13 ggplot repos (`P:62`); named themes in MM, MOCHA, IHA (`P:62`, `RS/MOCHA_Manuscript:103`, `PF/IHA-Figure/shard_11:64`) |
| `fig.size` | width and height in mm; 89 mm single column, 183 mm double | rec (`P:449-464`) |
| `fig.export` | `ggsave(..., units = "mm", device = cairo_pdf)` plus a PNG for review; no `dpi` on vector files; `pdf.fonttype = 42` and `ps.fonttype = 42`; Python: `fig.savefig(run_dir / name, bbox_inches = "tight")` then `plt.close(fig)`; one `sc.settings.figdir = run_dir` | obs 6 repos set fonttype (`P:348`); IHA `dpi = 1000` on PDF did nothing (`cc_learnings/03_plotting.md:75`); MM `cairo_pdf` 34 uses (`RS/multiple-myeloma:116`) |
| `fig.raster` | rasterise dense point layers only | rec (`P:456`) |
| `fig.source_data` | one source-data table per panel from the frame passed to the plot | obs IHA, IFN, MOCHA (`P:429-435`) |
| `fig.text` | no subtitle and no caption inside the plot; a title only when the panel stands alone; the caption sits in the report under the figure, one or two sentences of at most 60 words, opening with what is shown ("Box plot of ..."); the caption also says how the view was made: reduction (PCs, neighbours, `min_dist`, seed), downsampling, transforms, and the genome build for genomic panels | ext immgenT figure site (aguerozz.github.io/immgenT-GP-analysis): no `subtitle` or `caption` in any `labs()`, `labs(title = "")` in Figure 2E, page captions of 15-27 words in Figure 2 and up to 61 in Figure 1a; ext biology-skills v0.6.0 (github.com/biologyskills/biology-skills) `hadley-viz` embeddings `:20`, `:54`, aggregation `:30`; adopted 2026-10-27, their in-plot subtitle rejected |
| `fig.wilke` | each donor a point over any summary, never bar plus error bar (box with `outlier.shape = NA` under jittered points; `sns.boxplot(showfliers=False)` + `stripplot`); no violin for few points; name SD, SE or CI in the caption; bars start at 0, log-scale bars at 1, otherwise dots; size maps to area (`scale_size_area()`); axis and legend titles carry the unit; unordered categories ordered by value (`fct_reorder()`), natural orders kept, palette order kept for groups; legend order matches the plot, direct labels past ~8 categories (`legend_loc="on data"`); colour only when it carries meaning, alpha on group fills, colour plus shape, no rainbow for continuous, no red-green; small multiples share axes or the caption says not, panels in a meaningful order, lowercase tags (`tag_levels = "a"`); overplotting: smaller points, then alpha, then `geom_hex()`/`hexbin`; filled shapes over outlines, alpha on overlapping densities, try several bin widths; no 3D; check at print width, PDF for print, PNG for review, never JPEG; paired donors get an x = y line or a slopegraph. Chart choice: forest plot for effects with intervals (reference at 0, or 1 on log); stacked bars per sample for many proportions, faceted box with points for inference; row-scaled heatmap, stated, with a light-centred diverging palette; few points means under 10 per group (`FEW_POINTS` in `plots.*`): points plus a median, no box or violin; many subjects over time may add a mean line with a named interval (SD, SE or CI) beside the per-subject lines; never draw a line across a missing time point, and note dropped subjects | ext Wilke, Fundamentals of Data Visualization (clauswilke.com/dataviz: boxplots-violins, visualizing-uncertainty, proportional-ink, figure-titles-captions, visualizing-amounts, redundant-coding, color-pitfalls, color-basics, multi-panel-figures, overlapping-points, avoid-line-drawings, no-3d, small-axis-labels, image-file-formats, visualizing-associations, visualizing-proportions, histograms-density-plots); data-to-viz.com caveats (error_bar, boxplot, radius_or_area, order_data, color_com_nothing, overplotting, stacking, multi_distribution, heatmap); r-graph-gallery.com and python-graph-gallery.com boxplot-with-jitter code; ext biology-skills v0.6.0 (github.com/biologyskills/biology-skills) `hadley-viz` distributions `:16`, `:38`, longitudinal `:18`, `:31`, missingness `:17`; user decision 2026-10-27 (n < 10) |
| `fig.heatmap` | diverging palette centred at 0; symmetric limits held in a named constant; `oob = scales::squish`; the caption names the order basis: clustered (metric and linkage) or fixed by label | obs centred at 0 in six repos (`P:529-937`); named limit in ALTRA, immgenT, MM (`PF/ALTRA-manuscript/shard_12:137`, `PF/immgenT-GP-analysis/shard_06:113`, `P:859`); `squish` in MM, immgenT; IFN and reflex_v2 use literal or asymmetric limits (`RS/IFN_manuscript:149`, `RS/reflex_manuscript_v2:130`); ext biology-skills v0.6.0 (github.com/biologyskills/biology-skills) `hadley-viz` heatmaps `:17`, `:37` |
| `fig.stars` | `***` < 0.001, `**` < 0.01, `*` < 0.05 on adjusted p; omit "ns" | obs MM drops ns (`P:869`); cns keeps `NS` (`RS/cns-haematopoiesis:158`); spatial stars raw p (`RS/spatial_skin_atlas:1544`); ALTRA and IFN differ (`A:223`, `RS/IFN_manuscript:187`, `:213`) |
| `fig.labels` | group sizes in tick labels `(n = 12)`; `bquote()` for maths; never a raw column or variable name | rec (`P:460`); obs spatial, variancePartition (`P:410-414`) |
| `fig.colourblind` | check the palette with a colour-vision-deficiency simulator before release | **new**: no repo documents one (`cc_learnings/03_plotting.md:41`); 0 corpus witnesses |
| `fig.panel_id` | one named vector, `panels <- c(fig2c = "umap-donor")`, drives the plot object name, the file name and the source-data name | fails in 5 lineages (`RS/reflex_manuscript:12`, `RS/IFN_manuscript:22`, `RS/ALTRA-manuscript:24`, `RS/IHA-Figure:294`, `PF/MOCHA_Manuscript/shard_02:183`); model immgenT (`PF/immgenT-GP-analysis/shard_08:126`) |
| `fig.once` | build each plot once and pass it explicitly, `ggsave(plot = p_fig2c, ...)`; never a second copy of the plot code for export, never `last_plot()`; display and save size from the same constants | obs reflex_v2 (`RS/reflex_manuscript_v2:112`), IFN (`RS/IFN_manuscript:139-140`), ALTRA `last_plot()` in 7 notebooks (`RS/ALTRA-manuscript:115`), IHA size mismatch in 14 files (`RS/IHA-Figure:177`); same-name overwrites in cns, skin_fib (`RS/cns-haematopoiesis:140`, `RS/skin_fibroblast_atlas:118`) |

## Gaps the corpus does not settle

Choose these per project and record them in the palette file header:

- dotplot defaults
- volcano thresholds; one pair is seen (`Q < 0.1`, `|log2FC| > 0.5`)
- heatmap anchors; one example is seen (`-4/0/4`)
- a font family
