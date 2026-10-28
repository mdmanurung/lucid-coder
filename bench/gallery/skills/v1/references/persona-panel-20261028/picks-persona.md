# The coder this user wants to read: picks, contradictions, persona

Source: three blinded rounds (user-picks.md) read against the actual panel code. One line per pick on what the code differed in; "(inferred)" where the user gave no reason.

## 1. What differed, pick by pick

| Round, panel | Picked vs rejected: the concrete code difference |
|---|---|
| R1 setup | D: 4-line arg parse, `stop()`; A: `parse_args()` while loop + `die()` + `sh()` helpers (70 lines). |
| R1 load | D: `inner_join(..., relationship = "one-to-one", unmatched = "error")` in 14 lines; A: 30 `die()` lines. The praised `stopifnot(identical(cells[["cell_id"]], colnames(counts_mat)))` is in C, not A (user-picks.md misattributes it). |
| R1 filter/agg | A: 6 lines, `t(rowsum(t(m), cells$donor_id))`; D and C: 100-line summaries and `result` list. A was short because its gene filter sat in the fit panel. |
| R1 fit | D: `left_join` back onto all genes so filtered genes keep NA, three `# ... ----` blocks; C: tested genes only; A: adds `lmFit` sensitivity test. |
| R1 write | B: writes + `result` list + commands.sh, no helpers; D: adds `sh_arg <- \(x) ...`; A: 190 lines. |
| R1 notebook | B: `md()`, `tbl()`, `rd()` helpers each used many times, 90 lines; D, A, C: 300+ lines of `display_markdown(sprintf(paste(...)))`. |
| R1 workflow | D: `1. **Read and match.** ...`; others plain sentences. |
| Corpus R setup | immgenT: `source("code/R/plot_utils.R") # lineage_colors()`; wants ALTRA's `###### LIBRARIES ######` added. |
| Corpus R load | dreamlet: base `merge`, `match`, one-line comment per block. |
| Corpus R wrangle | sound-life: one pipe, a comment above each verb (uses `%>%`). |
| Corpus R model | ALTRA: `group_by(labels) |> group_modify(~freq_glm(.x, ...))`, commented; wants `|>`. |
| Corpus R plot | none; immgenT's one-line 300-char `theme(...)` called "too dense". |
| Corpus R write | ALTRA: `#### export results ####` banner, `dir.create` guard, `write_csv`; banner "too crowded". |
| Corpus R report | GTestimate: plot, theme tweak, `ggsave`, no dead code (inferred; IFN has commented-out `pdf()`/`dev.off()`). |
| Corpus Py setup | cns: imports, then `# ---- Defaults ----` block of `N_LATENT = 25`. |
| Corpus Py load | cns: numbered `# 2. ... # 6.` comments, reports missing rows (inferred). |
| Corpus Py model | aifi: plain `for res_num, res_name in tasks:`; rejected `ThreadPoolExecutor` and commented-out code (inferred). |
| Corpus Py plot | skin atlas: groupby to proportions, then one stacked bar with `n=` labels (inferred). |
| R2 setup | D: Purpose/Inputs/Outputs/Run header, `MIN_CELLS_PER_DONOR <- 10 # why`, `stopifnot("usage: ..." = ...)`; A: 90-line parser. |
| R2 load | D: visible `stopifnot(identical(sort(colnames(counts_mat)), sort(cells$cell_id)))`; C: same check hidden in `unmatched = "error"`. |
| R2 filter/agg | C: `sumTechReps(count_mat, ID = cells$donor_id)` + `filterByExpr` + `normLibSizes` in one place; A: `add_filter` closure; B: `sapply(... rowSums)`. |
| R2 fit | D: three sectioned blocks, `glmQLFit(robust = TRUE)`, all-gene `left_join`; A: 70 lines with sensitivity `lmFit` and NA-filled `detail` frame. |
| R2 write | C: 27 lines of pure `write_tsv`/`write_json`, `result` built earlier; D: 30-line `result <- list(...)` inline (inferred). |
| R2 notebook | D: parameters cell, then libraries, `options`, `theme_set`, `N_TOP_GENES`, `CONDITION_COLOURS` before any read; C: `withCallingHandlers` muffle block first. |
| R2 workflow | A: eight one-line verb-first steps; D: bold keywords, wrapped lines, trailing paragraph (inferred). |
| R3 02 whole | D: 43 lines, one pipe `pivot_longer |> summarise(.by) |> inner_join(unmatched = "error")`; A: 83 lines, five `if/stop` blocks (inferred). |
| R3 03 setup | D: `CELL_TYPE_LABELS <- c(t_cell = "T cell", ...)` upfront; A: six validation `if/stop` blocks in setup (inferred). |
| R3 03 load | D: `pull(hex, name = label)` with `# fct() errors on any condition ... not listed above`; A: explicit `cols(...)` and `setNames`. |
| R3 03 figure | none; D bare `theme_bw`, A over-styled with `pretty_cell_type()` and 10-line `theme()`. |
| R3 04 whole | none: likes A's per-file header and `# Parameters ------` dividers, A's `write_output()` called 4x over D's `iwalk(\(table, name) ...)`; likes D's `MIN_LIB_SIZE`, D's `filter(sum(count) > MIN_LIB_SIZE, .by = sample_id)`, D's comment-then-block. |

## 2. Contradictions and how to resolve them

1. **Checks.** "I am done with the checks" and "to the point", yet twice (R1, R2 load) he wanted a named `stopifnot` even when `unmatched = "error"` already enforced it. Resolve: one visible named assertion at each trust boundary; no validation walls; a check hidden in a function argument does not count as visible.
2. **Aggregation.** R1 picked base `t(rowsum(t(m), ...))` "to the point"; R2 picked `sumTechReps` as "existing function" over `rowsum`. Resolve: length wins first; at equal length the domain package's named function wins.
3. **Workflow format.** R1 picked bold keywords; R2 picked plain verb-first lines over bold keywords. Resolve: verb first, one line per step; bold is optional.
4. **Headers.** Corpus write banner "too crowded"; R3 preferred A's longer header and full-width dividers over D's compact four lines. Resolve: no `####` walls; yes to `# Section ------` to the margin and one line per input/output file.
5. **Functions.** "No helper unless used three times", yet he liked B's notebook helpers (each used 5+ times, consistent) and A's `write_output()` (4 calls) over an `iwalk` lambda. Resolve: a named helper called repeatedly beats an anonymous lambda in an iterator, even when the lambda is shorter.
6. **Tidyverse-first** vs picks of base `merge`/`match` (dreamlet) and `rowsum` (A). Resolve: tidyverse for the data-frame layer and pipes; matrix steps take the one-line base or package call.
7. **Write section.** R1 "saving key parameters, tone down"; R2 picked the arm whose write block held no `result` list. Resolve: result.json yes, built before the write section, which holds writes only.

## 3. Which kind of preference drove picks

1. **Readability of text** (headers, sectioning, whitespace, comment placement, UPPER_SNAKE parameters upfront): about 14 reasoned picks.
2. **Shape of code** (no one-off helpers, existing package function, plain loop over iterator, pipe over if-chain): about 9.
3. **Scientific procedure**: about 3, but each is a veto. "Lengthy but necessary to keep the whole gene lists" beat brevity (R1 fit); the alignment check beat concision (R1, R2 load).

## 4. Never commented on

`col_types` strictness, error-message wording, `robust = TRUE`, the BH family (all input genes vs tested genes varied across arms), A's `lmFit` sensitivity check (never praised, A lost), reproducibility files (commands.sh, sessionInfo, md5), figure styling (mm, cairo, palettes, point shapes), reference lists in notebooks, CLI `--help`, run-time. He named exactly one check across three rounds.

## 5. Persona

A computational immunologist who reads analysis scripts top to bottom the way a reviewer reads a methods section. He wants to know in the first screen what the script does, what goes in, what comes out, and which numbers are tunable; then he wants each step to arrive as a short, labelled block he can skim by its comment line. He does not trust cleverness, iterators, or wrappers; he trusts a named edgeR or dplyr call he already knows, repeated plainly if need be. He does not want defensive code, but he wants one sentence-length assertion at every join that would otherwise fail silently, and he will accept extra lines only when they protect the dataset's completeness (every gene stays in the table) or the replication unit. He optimises for the time it takes a colleague to understand and rerun the script, not for line count, elegance, or exhaustiveness.

## 6. Habits

1. **Open with Purpose / Inputs / Outputs / Run, one line per file.** He called this "informative" (R2 setup) and preferred A's per-file descriptions over a crammed block (R3 04).
2. **Put every tunable at the top as `MIN_CELLS_PER_DONOR <- 10 # why`.** UPPER_SNAKE with a trailing reason "makes it clear" (R1 setup, R2 notebook, corpus Py setup `N_LATENT = 25`).
3. **Divide the script with `# Load data ----` style lines and a blank line between blocks.** "Sectioned well by linebreak smartly" (R1 fit); "white space good" (R2 fit); dense one-liners lose (corpus R plot).
4. **Write one comment line above each block saying what it does, not how.** "Comments informative and well positioned" (R3 03 load: `# fct() errors on any condition ... not listed above`); sound-life's comment-per-verb pipe won.
5. **Parse arguments in four lines, no `parse_args()`, no `die()`.** `input_dir <- args[match("--input", args) + 1]` plus a named `stopifnot` beat every 70-line parser (R1, R2 setup).
6. **Call the package's own function instead of rebuilding it.** `sumTechReps(count_mat, ID = cells$donor_id)` over `sapply(... rowSums)` (R2 prep); `glmQLFit`/`topTags` over hand-rolled frames; "I would use limma/edgeR".
7. **Keep one visible named assertion after every join.** `stopifnot("every count column has a metadata row" = identical(...))` was wanted even when `unmatched = "error"` covered it (R1, R2 load).
8. **Repeat a plain call rather than iterate with a lambda.** `write_output(group_means, run_dir, "group_means.tsv")` four times is "directly understandable"; `list(...) |> iwalk(\(table, name) ...)` is not (R3 04).
9. **Keep the whole gene table.** `tibble(gene = rownames(count_mat)) |> left_join(de_tested, ...)` so filtered genes carry NA was "lengthy but necessary" (R1 fit).
10. **Finish the script with writes only; put plots in a notebook.** The 27-line pure-write block won (R2 write); the notebook that opened with parameters, libraries and settings won (R2 notebook); end-of-script report text is "very distracting".
