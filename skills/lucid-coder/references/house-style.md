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

## R

| Topic | Default | Evidence |
|---|---|---|
| `r.assign` | `<-`, never `=` at top level | obs: dominant in most repos, `=` in dreamlet and MOCHA (`RS/dreamlet:53`, `RS/MOCHA:69`) |
| `r.indent` | 2 spaces | obs: most R repos; 4 in reflex (`C:155`) |
| `r.pipe` | `\|>` always, never `%>%`, in new and existing repositories | user rule 2026-10-26 (D1 revised in `PLAN.md`). `%>%` appears in 16 of 17 (`C:154`); `\|>` is base R 4.1+ and needs no dependency |
| `r.quotes` | double quotes | new: repos mix them (`RS/reflex_manuscript:62`) |
| `r.bool` | `TRUE`/`FALSE`, never `T`/`F` | rec (`C:209`) |
| `r.width` | soft limit 100 columns, set as `line-width = 100` in `air.toml`; one argument per line for long calls | decided 2026-10-26 (D3); rec (`C:209`); no repo had a rule (`C:156`) |
| `r.formatter` | `air` | decided 2026-10-26 (D2). No repo had a formatter config (`cc_learnings/README.md:48`) |
| `r.library` | `library()` at the top, only for what the script uses; `pkg::fn()` for one-off calls; `suppressPackageStartupMessages()` is fine (startup messages are not warnings); no `quiet_library()` wrapper | rec (`02:99`); the wrapper ban is new: the wrapper appears in 5 repos (`C:142`) |
| `r.names` | snake_case objects, UPPER_SNAKE constants, plot objects by panel (`p_fig2c`); never `c`, `T`, `F`, `df` or `data`; argument names spelled in full | rec (`02:84`); obs six repos snake_case (`02:77`, `C:65`, `C:968`); shadowing in IFN (`02:81`) and immgenT (`C:1024`); partial matching in 5 repos (`S:122`) |
| `r.helpers` | helper files use `pkg::fn()` only, never `library()` or `require()` inside a function; a `source()` line ends with a comment naming what it provides | obs immgenT (`C:1032`, `C:1035`), MOCHA (`02:66`); anti-pattern GTestimate `require()` in 11 of 16 helpers (`C:698`) |
| `r.errors` | `stop("<object>: <problem>; <fix>.", call. = FALSE)` | obs MM, IFN (`RS/multiple-myeloma:79`, `RS/IFN_manuscript:98`) |
| `r.dispatch` | `switch()` with `match.arg()` to dispatch on a string; `dplyr::case_when()` over nested `ifelse()`; `%\|\|%` for fallbacks; `is.numeric(x)`, never `class(x)[[1]] == "numeric"`; S3 methods when branching on class | ext: Bryan, code smells |

## Any language

| Topic | Default | Evidence |
|---|---|---|
| `code.conditionals` | guard clauses first: `stop()` or `return()` early, no `else` after an exit; a fussy condition becomes a named predicate (`is_paired_sample()`); a few small functions over one long one; less nesting | ext: Bryan, code smells |
| `code.toggles` | never comment code in and out to change behaviour; use a named argument or constant | ext: Bryan, code smells |
| `code.functions` | inputs as arguments, one value out, no globals, defined before first use | rec (`02:70`, `C:208`); globals in ALTRA, MOCHA, IHA (`C:145`); use before definition in MOCHA (`C:146`) |
| `code.matching` | exact matches: `%in%`, `fixed = TRUE`, anchored patterns; never `grep(paste(genes, collapse = "\|"))`; filter with `!grepl()`, never `x[-grep()]` or `x[-which()]` | obs 7 repos (`S:125`); empty-result bugs in MM, cns, spatial (`S:40`, `S:45`, `S:1172`) |
| `code.join_by` | name `by =` in every join; pandas `validate=` | obs natural joins in 9+ repos (`S:68`); ALTRA `merge` drops (`S:288`) |
| `code.warnings` | never `suppressWarnings()`, `options(warn = -1)` or `filterwarnings("ignore")`; silence one category with a one-line reason | rec (`02:114`); obs aifi, IHA, sound-life, spatial (`C:378`, `C:936`, `C:1524`, `C:1823`) |
| `code.deprecated` | `linewidth` not `size` for lines; no `aes_string()`, `spread()`, `is_categorical_dtype()` | obs 10 repos (`S:130`) |
| `code.shebang` | `#!/usr/bin/env Rscript` or `#!/usr/bin/env python3` on runnable scripts; shell and SLURM wrappers use `set -euo pipefail` and `${1:?usage}` | obs cns (`C:86`, `C:109`), immgenT (`C:981`); rec template (`C:216`) |

## Python

| Topic | Default | Evidence |
|---|---|---|
| `py.formatter` | `ruff format` (black-compatible, 88 columns) | rec: black at 88 (`C:209`); ruff is one tool for format and lint |
| `py.strings` | f-strings | new: repos mix `.format()` (`RS/sound-life-scrna-analysis:81`) |
| `py.paths` | `pathlib.Path`, `ROOT = Path(__file__).resolve().parents[1]` | rec template (`C:247-250`) |
| `py.args` | `argparse` with typed `Path` arguments; `main()` behind an `if __name__ == "__main__"` guard | obs cns: argparse in 32 of 46 scripts (`C:109`, `C:204`, `repos_new/cns-haematopoiesis`) |
| `py.errors` | raise `ValueError`/`KeyError` with the fix in the message; no `assert` outside tests | rec (`SYNTHESIS.md:172`); obs cns, spatial (`RS/cns-haematopoiesis:189`, `RS/spatial_skin_atlas:819`) |
| `py.imports` | one import block at the top; never re-import mid-notebook | obs: cns re-imports scanpy 3-6 times (`C:160`) |
| `py.seed` | `random_state` as a function argument | rec (`SYNTHESIS.md:146`) |

## Files and layout

| Topic | Default | Evidence |
|---|---|---|
| `files.names` | `NN_verb-noun.R`: `_` between fields, `-` between words; lowercase, no spaces or other punctuation; zero-padded prefix, with gaps such as `90_` for odd-one-out scripts; helpers `helperNN_verb-noun.R`; ISO 8601 `YYYY-MM-DD` first when order is by date | decided 2026-10-26 (D7); ext: Bryan, naming sequel; rec (`SYNTHESIS.md:109`); all nine digested repositories fail it (`SYNTHESIS.md:111`) |
| `files.layout` | `data/` raw and read-only (`chmod u-w -R data/`), `scripts/`, `results/<run_dir>/`; one project per repository; a README at the root that lists any manual step with its exact command | ext: Tang, Openscapes, Tidyverse course |
| `files.paths` | `here::here()` in R; Python: see `py.paths`; every path relative to the project root; save source, never the workspace; each script runs in a fresh process | ext: Openscapes, Tidyverse course, Tang. No corpus repo uses `here()` (`C:668`, `C:871`); the corpus recommends running from the repo root with a root guard (`C:205`) |
| `files.versions` | version through git and run folders; never `_v2`, `_final`, `-Copy1` | obs failures (`SYNTHESIS.md:111`) |
| `files.outputs` | `<run_dir>/` passed as an argument; files named by content or panel | rec (`C:211`); 8 repos had names that drift from panels (`S:120`) |
| `files.helpers` | a helper file only when 2 or more scripts call it; one function per file in `R/`, or one file per topic (`R/<topic>.R`) | rec (`cc_learnings/02_coding_R_python.md:72`); obs GTestimate (`C:136`) |
| `files.notebooks` | runs top to bottom from a fresh kernel; prefer `.qmd`/`.Rmd` or marimo, which diff in git; strip `.ipynb` outputs before commit; no `.ipynb_checkpoints` committed | rec (`phase_A/CONVENTIONS.md`); ext: Tang |

## Reproducibility and data

| Topic | Default | Evidence |
|---|---|---|
| `repro.env` | commit a lockfile: `renv.lock`, `uv.lock` or `environment.yml`; export conda environments without the machine `prefix:` line; a container only when the lockfile is not enough | ext: Tang, PLOS 2021; obs conda exports with `prefix:` and missing packages (`A:74`, `A:82`, `A:263`) |
| `repro.automation` | a script for every step; repeated routines in a shell script or a Snakemake or `{targets}` rule; a notebook only explores: no reported number comes from one or is typed by hand | ext: Tang; rec (`C:203`, `X:145`); hand-typed numbers in GTestimate, immgenT, spatial, MM (`A:222`) |
| `data.recode` | one named level vector drives both the recode and `factor()`; no catch-all arm (`TRUE ~ "Other"`, `.fillna("Other")`), so an unmapped label surfaces as `NA` | obs cns (`S:342`), IHA (`S:601`), MM trailing-space keys (`S:761`), spatial (`S:1154`) |
| `data.tidy` | long format before analysis; one value per cell; `NA` for unknown, never blank, an undeclared `-999` or the string `"NA"` (`S:169`, `S:601`); `trimws()` label values at load (`S:761`); `janitor::clean_names()` on import; one case-consistent code per category; plain-text data with ISO 8601 dates | ext: Tidyverse course, Tang |
| `data.exclusions` | log every removed sample, cell or gene with the reason and the count removed, and every group a loop skips, with why | ext: PLOS 2021; silent skips in IFN, aging, spatial (`S:470`, `S:25`, `S:74`); model MM `warning()` + `next` (`S:774`) |
| `repro.seed` | seed each stochastic step in the script, a distinct seed per replicate; never one fixed seed inside a repeated helper; `RNGkind("L'Ecuyer-CMRG")` or seed in the worker for parallel code; `scvi.settings.seed` | rec (`A:257`); obs IHA, IFN per-iteration seeds (`A:54-55`); failures: MOCHA fixed seed in a helper (`A:67`, `A:259`), variancePartition `mclapply` (`A:65`), unseeded scVI in cns and spatial (`A:51`, `A:64`) |

## Checks and comments

These topics let the user tune the two things they care about most.
One example correction is "stop checking column presence in this repo".

| Topic | Default | Evidence |
|---|---|---|
| `checks.default_five` | new run dir; inputs and columns at load; join cardinality; ID identity before positional copy; outputs non-empty | rec (`S:142`); ext: PLOS 2021 "data gut checks" after merges |
| `checks.extra` | only for a known scientific assumption that would fail silently | rec (`S:142`); user rule (`SYNTHESIS.md:170`) |
| `checks.budget` | the default five plus at most one extra per script; one `if` per input at load and per join; skip a check when the next step already fails loudly; no audits of input content unless the input is untrusted, except raw integer counts before a count model | eval pilot 2026-10-26: without a budget, lucid-coder wrote a median of 8 checks per script against 3 without it; user rule |
| `checks.tripwires` | off: no hard-coded expected counts unless asked | decided 2026-10-26 (D4); immgenT uses them deliberately (`RS/immgenT-GP-analysis:110`) |
| `checks.style` | R `if (...) stop(..., call. = FALSE)`; Python `raise`; message names object and fix | obs (`RS/multiple-myeloma:79`, `RS/cns-haematopoiesis:189`) |
| `checks.files` | no separate `check_*.R` files | user rule (`SYNTHESIS.md:170`) |
| `comments.why` | comment decisions, number origins, reference levels, API needs, limits, and the source of a marker list or copied code; never restate code | obs dreamlet, MM, immgenT (`C:600`, `C:1159`, `C:1032`); provenance in cns, aifi, MOCHA (`S:365`, `C:377`, `S:735`) |
| `comments.header` | Purpose, Inputs [producer], Outputs, Run, Notes; no blank fields, no author line | rec (`C:206`) |
| `comments.dead_code` | delete commented-out code | rec (`C:210`) |

## Figures

| Topic | Default | Evidence |
|---|---|---|
| `fig.palette` | one `palettes.R` or `palettes.py` of named maps keyed by label; unknown label fails | rec (`P:449-464`); 12 repos define colours inline (`P:344`) |
| `fig.theme` | `theme_set(theme_bw(base_size = 8))` once per script | user rule 2026-10-26 (D9); obs: `theme_classic` in 12 of 13 ggplot repos (`P:62`); only MM defines a named theme (`P:62`) |
| `fig.size` | width and height in mm; 89 mm single column, 183 mm double | rec (`P:449-464`) |
| `fig.export` | `ggsave(..., units = "mm", device = cairo_pdf)` plus a PNG for review; no `dpi` on vector files; `pdf.fonttype = 42` and `ps.fonttype = 42` | obs 6 repos set fonttype (`P:348`); IHA `dpi = 1000` on PDF did nothing (`cc_learnings/03_plotting.md:75`) |
| `fig.raster` | rasterise dense point layers only | rec (`P:456`) |
| `fig.source_data` | one source-data table per panel from the frame passed to the plot | obs IHA, IFN, MOCHA (`P:429-435`) |
| `fig.heatmap` | diverging palette centred at 0; symmetric limits held in a named constant; `oob = scales::squish` | obs ALTRA, dreamlet, IFN, immgenT, MM, reflex_v2 (`P:529-937`); the limit itself varies (±2 in 4 of 10) |
| `fig.stars` | `***` < 0.001, `**` < 0.01, `*` < 0.05 on adjusted p; omit "ns" | obs MM, cns, spatial (`A:223`, `P:869`, `P:413`); ALTRA differs |
| `fig.labels` | group sizes in tick labels `(n = 12)`; `bquote()` for maths; never a raw column or variable name | rec (`P:460`); obs spatial, variancePartition (`P:410-414`) |
| `fig.colourblind` | check the palette with a colour-vision-deficiency simulator before release | **new**: no repo documents one (`cc_learnings/03_plotting.md:41`) |

## Gaps the corpus does not settle

Choose these per project and record them in the palette file header:

- dotplot defaults
- volcano thresholds; one pair is seen (`Q < 0.1`, `|log2FC| > 0.5`)
- heatmap anchors; one example is seen (`-4/0/4`)
- a font family
