## House rules (decided by the user; not up for judgement)
- Parameters a reader may change: top of the script, UPPER_SNAKE, a reason beside each. Package defaults and one-use values stay inline.
- A named function only at 3+ calls in a script or shared by 2+ scripts. Each output file gets its own write line.
- Header: Purpose, Run, Inputs and Outputs one file per line. Every block opens with `# Name ` padded with `-` to column 80; a one-line comment names the step; a blank line between steps.
- `|>` never `%>%`; tidyverse for data frames; a count matrix stays a matrix.
- Checks: one `stopifnot(identical(...))` where metadata meets a matrix, a cardinality argument on joins, a run dir never overwritten, at most one more. No input audits, type checks, `col_types` for its own sake, expected counts or printed confirmations.
- Scripts prep, fit and write tables; a notebook makes every figure and all report prose, parameters in its first cell, one sentence of prose per section. Only the figures the task asks for.
- Build the plot frame, then a short ggplot call. Donors as points over a named summary.
- `suppressPackageStartupMessages()` around the library block.
