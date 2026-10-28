## House rules (the user decided these)
- Parameters a reader may change: right after the header, UPPER_SNAKE, a reason each. Package defaults and one-use values stay inline.
- A named function only at 3+ calls in a script or shared by 2+ scripts. Each output file gets its own write line.
- Header: Purpose, Run, Inputs and Outputs one file per line. Every block opens with `# Name ` padded with `-` to column 80; a one-line comment names the step; a blank line between steps.
- `|>` never `%>%`; tidyverse for data frames; a count matrix stays a matrix.
- Checks: one `stopifnot(identical(...))` where metadata meets a matrix, a cardinality argument on joins, a run dir never overwritten, at most one more. No input audits, type checks, expected counts or printed confirmations.
- Scripts prep, fit and write tables; a notebook makes every figure and all report prose, parameters in its first cell, one sentence of prose per section. Draw the figures a reviewer expects for the method (DE: donor MDS, volcano, top genes per donor).
- Build the plot frame, then a short ggplot call. Donors as points over a named summary.
- `suppressPackageStartupMessages()` around the library block.
