# Code rubric

A rubric for judging an analysis script by what can be counted or checked, not by taste. Each item has a decision rule that two readers apply the same way. They record evidence as line numbers, and they compute agreement before results are reported.

## Principles

- **Counts and yes/no only.** There are no 1-5 scales. The LLM judges put every arm at the top of their 1-5 readability scale (`judges.tsv`, q1), so a scale does not separate arms.
- **Every mark carries evidence.** A "no" or a count lists the line numbers. A mark without them does not count.
- **Report items separately.** Never sum them into one score. A short script with a silent mislabel is not "mostly good". Results are reported per item and per arm, as `report.js` already does for L3.
- **Score blind.** Readers see neutral ids (`make-packet.js`), never the arm.
- **Quality and house style are separate.** Part 1 measures properties any good analysis script has. Part 2 measures conformance to lucid-coder's own choices (`|>`, `theme_bw(8)`, learned rules). Only Part 1 says whether code is better. Part 2 only says whether the plugin was followed.
- **"n/a" is not a pass.** An item that does not apply, such as no joins in the script, is marked n/a and is left out of that item's denominator.

Measure codes: **A** = a `score.js` column (regex items come from `eval/measures.js`); **G** = a regex not yet a column; **M** = manual, by a blinded reader using the rule below.

## Part 1: quality

### Correctness

| Id | Item | Decision rule | Measure |
|---|---|---|---|
| Q1 | Task gate passes | `check.R` returns `pass` | A `gate` |
| Q2 | Clean rerun | the script runs from a fresh run dir with exit 0 | A `rerun_ok` |
| Q3 | Inference unit is right | the test's rows are the replication unit the design names, such as one value per donor; pseudo-replication is a "no" | M |
| Q4 | Contrasts and coefficients by name | every coefficient, column or list element that drives a result is selected by name; one `coef = 2`, `[, 3]` or `$d` partial match is a "no" | M, partly G |

### Silent-failure guards

| Id | Item | Decision rule | Measure |
|---|---|---|---|
| S1 | Join keys named | every join names its keys (`by =`, `join_by()`, `on=`); count joins without them | A `joins_unkeyed` |
| S2 | Join cardinality guarded | every join whose expected cardinality is many-to-one or one-to-one has `relationship =` / `validate=`, or one check right after it; count unguarded joins | A `join_guard_count`, plus M to count joins |
| S3 | Positional copy guarded | every positional pairing is preceded by an identity check or done by name. Examples: `cbind` of two tables, `colnames<-`, a design matrix built from a table and matched to a count matrix. Count unguarded pairings | M |
| S4 | Removals reported | every step that drops samples, cells or genes prints how many it dropped (a `message()`, `print()` or one counts table); count silent drops | M |
| S5 | No silent NA or warning suppression | no `na.rm = TRUE` on a result that should not have NA, no fill of NA without a comment, and no blanket `suppressWarnings()`, `options(warn = -1)` or `filterwarnings("ignore")`; count | A `warning_suppress` for suppression, M for NA |
| S6 | Stochastic steps seeded | every random step (jitter, repel, sampling, UMAP, clustering) has a seed; n/a if none; count unseeded | M, partly G |
| S7 | Outputs cannot overwrite | the script refuses an existing run dir, or writes only into a new one; yes/no | M |

### Check economy

The budget is `rules/core.md` "## Checks", applied as in the grid 5 classification (`HANDOFF.md`).

| Id | Item | Decision rule | Measure |
|---|---|---|---|
| E1 | Checks out of budget | count checks classified OUT: input audits, column checks at load, type checks, expected counts, re-derived results, a second check on one join, file-exists, outputs-non-empty | M |
| E2 | Checks that cannot fail | count checks whose condition is true by construction, such as `identical()` on two objects built from the same names | M |
| E3 | Defaults present | of the defaults that apply (ID identity before a positional copy, join cardinality, new run dir), count the missing ones | M |

### Readability, measured

| Id | Item | Decision rule | Measure |
|---|---|---|---|
| R1 | Code lines | non-blank, non-comment lines; compare within a task, never across tasks | A `code_lines` |
| R2 | Longest pipe | the most steps in one pipe chain; flag above 10 (R4DS) | A `longest_pipe` |
| R3 | Deepest nesting | the most levels of `{`, or of Python indentation, inside a function or loop; flag above 3 | A `max_nesting` |
| R4 | Uninformative names | count assigned names in `tmp`, `temp`, `df`, `df2`, `data`, `x`, `res2`, `foo`, or a single letter outside a lambda or math index | A `vague_names` |
| R5 | Repeated literals | count numeric or string literals used 2 or more times that drive a decision (cutoffs, levels, ids) and are not named | M, partly G |
| R6 | Pasted blocks | count blocks of 4 or more lines repeated with only a name or literal changed, which a loop over a named vector would replace | M |
| R7 | Dead code | count commented-out code lines | A `commented_code_lines` |
| R8 | Header | Purpose, Inputs, Outputs and Run are all present in the first 20 lines | A `header_complete` |
| R9 | Decisions explained | of the named cutoffs, reference levels and formulas, count those with no reason (a comment or a header line) | M |
| R10 | Wrappers used once | count functions defined and called exactly once, outside an iterator or lambda (a function passed to `map()` or `group_modify(~ f(.x))` is reused) | A `wrappers_once` |

### Reproducibility

| Id | Item | Decision rule | Measure |
|---|---|---|---|
| P1 | Portable paths | no absolute paths, `setwd()` or `~`; count | A via `lint_count` (absolute and nonportable path linters), G |
| P2 | Versions recorded | session info is written to the run dir (`sessionInfo()`, `print_versions()`); yes/no | A `session_info` |
| P3 | Writes stay in the run dir | every write targets the run dir, except outputs the task names; count stray writes | M |
| P4 | Loads at the top | every `library()`/`require()` or `import` sits above the first other code line (setup such as `options()` may precede it); count late loads. A package used without loading it fails Q2 | A `late_loads` |

### Format

| Id | Item | Decision rule | Measure |
|---|---|---|---|
| F1 | Formatter clean | `air format --check` (R) or `ruff format --check` (Python) passes with the template config | A `air_ok` |
| F2 | Lints | lintr / ruff findings with the template config | A `lint_count` |

## Part 2: house-style conformance (lucid-coder only)

These items say whether the plugin was followed. They do not say whether the code is better.

| Id | Item | Decision rule | Measure |
|---|---|---|---|
| H1 | Native pipe | `|>` used, `%>%` absent | A `native_pipe` |
| H2 | House theme | `theme_bw(base_size = 8)` set; n/a without a plot | A `uptake_theme` |
| H3 | No plot titles | no `ggtitle()` or `labs(title =)`; n/a without a plot | A `uptake_no_title` |
| H4 | Learned rules (arm D) | `out_dir` named; `options(warn = 1)` first | A `rule_out_dir`, `rule_warn1` |
| H5 | Tidyverse first | no `data.table`, `apply()` on a data frame, `merge()` or base loops over rows; count | A `base_df_calls` |
| H6 | Lowercase R constants | no UPPER_SNAKE assignment in R; count | A `upper_constants` |

## Protocol

1. Run `node eval/score.js <results dir>`; this fills the A items.
2. The regex items are `score.js` columns too (`eval/measures.js`, tests in `eval/tests/measures.test.js`). They blank strings and comments first. They are regexes, not a parser, so check a sample against a reader before trusting a surprising count.
3. Make a blinded packet with `node eval/make-packet.js`. Two readers each mark every M item with line numbers, working separately.
4. Agreement: Cohen's kappa for yes/no items and weighted kappa for counts, binned 0 / 1 / 2+. `report.js` already computes weighted kappa. Report any item below 0.6 as unreliable, and reword its rule before the next round.
5. Report per arm: for each count the median and the number of scripts with 0; for each yes/no the rate. Put n per cell next to every table. There is no total.

An LLM judge may mark the M items as a second rater. It must quote the line for every mark. It replaces neither human reader.
