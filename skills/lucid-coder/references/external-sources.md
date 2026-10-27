# External sources

These sources were added at the user's request on 2026-10-26. `house-style.md`
labels the rules they support **ext**. Quotes are the authors' own words.
Rules the persona already had (no `setwd()`, `library()` at the top, delete
commented-out code, seeds, session info) are not repeated here.

## Jenny Bryan, "Code smells and feels" (useR! 2018)

<https://github.com/jennybc/code-smells-and-feels>

- "Do not comment and uncomment sections of code to alter behaviour."
- "A few little functions >> a monster function."
- "Extract fussy condition logic into well-named function", for example `is_parental()`.
- Guard clauses: "Hoist quick `stop()`s and `return()`s to the top." "An `if` block that ends with `stop()` or `return()` does not require `else`."
- "less indentation >> lots of indentation".
- "`switch()` is ideal if you need to dispatch different logic, based on a string." Pair it with `match.arg()`.
- "`dplyr::case_when()` is ideal if you need to dispatch different data, based on data", instead of nested `ifelse()`.
- Use `is.numeric(x)`, not `class(x)[[1]] == "numeric"`. When conditions branch on class, write S3 methods.
- `%||%` for fallback chains.

## Jenny Bryan, "How to name files: the sequel" (NormConf 2022)

<https://speakerdeck.com/jennybc/how-to-name-files-the-sequel>

- "Pick a filename convention. Any convention. Just pick one."
- Names are "machine readable, human readable, sorted in a useful way".
- "avoid spaces, punctuation, accented characters, case sensitivity".
- "`_` Underscore to delimit fields"; "`-` Hyphen to delimit words within fields".
- "embrace the slug": `01_marshal-data.R`, not `01.R`.
- Prefix helpers: `helper01_load-counts.R`, not one `helpers.R`.
- "left pad numbers with zeros"; leave gaps such as `90_` for odd-one-out scripts.
- "use the ISO 8601 standard for dates": `YYYY-MM-DD`.

## Ming Tang, reproducible computing (2022 ISCB slides, 2025 slides, blog post)

- <https://divingintogeneticsandgenomics.com/files/slides/2022-07-09_ISCB_reproducible.pdf>
- <https://divingintogeneticsandgenomics.com/files/slides/2025-06-25_short_good_enough_practice_reproducible_computing.pdf>
- <https://divingintogeneticsandgenomics.com/post/repro-comp/>

- Folders "/data (read-only), /scripts, /results"; "chmod u-w -R data/".
- "Always generate the output/intermediate files/figures in the results folder using the scripts in the scripts folder."
- "Use relative path or better use here() to construct file path"; "pyhere" in Python.
- "Pin package versions!": renv, uv, conda; containers when needed.
- "Write scripts for everything unless it is not possible. (manual editing, document, document, document!)"
- "The best documentation is automation"; repeated tasks go in a shell script or Snakemake, Nextflow or `{targets}`.
- "Jupyter notebook is not git friendly"; "Quarto is Git-friendly"; marimo for Python.
- "Use the right Null values"; `janitor::clean_names()`; "Be consistent" with codes and identifiers.

## Openscapes, "Refactor and modularize code" and "Coding strategies for future us"

- <https://openscapes.github.io/series/how-do-i/refactor-modularize-code.html>
- <https://openscapes.github.io/series/core-lessons/coding-strategies.html#project-oriented-workflows>

- Refactoring changes "how it works, not what it does".
- Remove hard-coded file paths; move "away from long for-loops" to small functions and `purrr`.
- "Save the source, not the workspace"; restart R often; every script runs in a fresh process.
- Avoid `rm(list = ls())`; use `here()`; one project = one repository; a project "should just work" on another computer.

## Wright, Ellis, Hicks and Peng, "Tidyverse Skills for Data Science" (2021)

<https://jhudatascience.org/tidyversecourse/>

- "Each variable you measure should be in one column"; each observation in its own row; one table per type of data.
- One code per value, consistent case; `NA` for unknown, no empty cells; one thing per cell (`weight_lbs`).
- Plain-text data files; ISO 8601 dates; `here()` for paths; never `setwd()`.
- Convert to long format before analysis; one pipeline step per line.
- Check duplicates after tidying, for example `janitor::get_dupes()`.

## Stoudt, Vásquez and Martinez, "Principles for data analysis workflows" (PLOS Comput Biol 2021)

<https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008770>

- "Data gut checks": dimensions after merges; nulls, zeros and negatives; how missing values are coded.
- Document removed data points and why.
- Record methodological choices and the calculations behind parameter values.
- Separate mechanistic code (functions) from the analysis narrative; shared helpers live in one file.
- Regression tests compare outputs before and after code changes.
- Record software versions; use a package manager.

## Long and Teetor, "R Cookbook", 2nd edition (2019)

<https://rc2e.com/>

Relevant recipes: 3.2 (one RStudio project per nontrivial piece of work),
3.3 and 3.16 (do not rely on the saved workspace), 3.14 (`Rscript` batch
runs), 8.4 (`set.seed`), 12.20 (code sections), 15.10 (reusable functions in
a sourced file; a package for larger work).

## Wickham, Çetinkaya-Rundel and Grolemund, "R for Data Science", 2nd edition (2023)

<https://r4ds.hadley.nz/workflow-style.html> is the R code-style authority (user rule, 2026-10-27).
Taken into `rules/r.md`: lowercase snake_case names, constants included (its "Avoid" example is `SHORTFLIGHTS`); variants share a prefix, not a suffix; a space before `|>` and `|>` at the end of the line; named arguments one per line; pipelines over about 10-15 lines split into named sub-tasks, especially after reshaping or summarising; ggplot2 `+` follows the pipe rules; sectioning comments `# Load data ----`. Spacing rules are left to `air`.
Also taken (R4DS idiom chapters and R Cookbook 2e, 2026-10-27): the `r.r4ds` row in `house-style.md`. R Cookbook predates `|>`, `join_by()` and `.default`; where the two books differ, R4DS wins. R Cookbook still supplies the modelling advice: fit with `data =`, then read results with accessors.

## Wilke, "Fundamentals of Data Visualization"; data-to-viz caveats; R and Python graph galleries

<https://clauswilke.com/dataviz/>, <https://www.data-to-viz.com/caveats.html>, <https://r-graph-gallery.com/>, <https://python-graph-gallery.com/>. These are the figure authority (user rule, 2026-10-27). Core `## Figures` now has points over summaries, named error bars, bars from 0, units in axis titles, shared panel axes and legend order. The full list and a chart-choice table are in the `fig.wilke` row of `house-style.md`. The 8 pt base size is the house's choice, not Wilke's. Wilke only says axis labels are usually too small.
