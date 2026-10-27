const test = require('node:test');
const assert = require('node:assert');
const { codeMeasures, strip } = require('../measures');

const R = `# Purpose: toy; this comment mentions merge( and |> and must not count
options(warn = 1)
library(dplyr)
CUTOFF <- 0.05
df <- read.csv("a|>b.csv")
x <- 1
fit_one <- function(d) {
  if (nrow(d) > 0) {
    for (i in seq_len(nrow(d))) {
      if (d$p[i] < CUTOFF) {
        message("hit")
      }
    }
  }
}
mapped <- function(d) d
res <- df |>
  filter(p < CUTOFF) |>
  mutate(q = p) |>
  left_join(meta, by = "id")
res2 <- merge(res, meta)
other <- left_join(res, meta)
fit_one(res)
per_group <- function(d) d
grouped <- res |> group_modify(~ per_group(.x))
out <- lapply(split(res, res$g), mapped)
summary <- res |>
  summarise(
    n = n(),
    by = 1
  )
suppressWarnings(as.numeric("1"))
library(tidyr)
sessionInfo()
`;

test('strip blanks comments and string contents but keeps lines', () => {
  const s = strip('a <- "x # y" # note\nb', false);
  assert.strictEqual(s, 'a <- "" \nb');
  assert.strictEqual(strip('"""doc\n# not code"""\nx = 1', true), '""""""\n\nx = 1');
});

test('R measures', () => {
  const m = codeMeasures(R, false);
  assert.strictEqual(m.joins_unkeyed, 2, 'merge() and left_join() without by');
  assert.strictEqual(m.warning_suppress, 1);
  assert.strictEqual(m.longest_pipe, 4, 'res chain: three |> join four steps');
  assert.strictEqual(m.max_nesting, 4, 'function > if > for > if');
  assert.strictEqual(m.vague_names, 3, 'df, x, res2; a named argument on its own line is not counted');
  assert.strictEqual(m.wrappers_once, 1, 'fit_one is called once; mapped is passed to lapply');
  assert.strictEqual(m.session_info, true);
  assert.strictEqual(m.late_loads, 1, 'library(tidyr) below code; options() above library() is setup');
  assert.strictEqual(m.base_df_calls, 3, 'merge(), lapply() and a for loop over nrow()');
  assert.strictEqual(m.upper_constants, 1);
  assert.ok(m.code_lines > 20);
});

const PY = `"""Purpose: toy."""
import warnings
import pandas as pd

warnings.filterwarnings("ignore")
warnings.simplefilter("ignore", category=FutureWarning)
CUTOFF = 0.05
tmp = pd.read_csv("x.csv")
merged = tmp.merge(meta, how="left")
keyed = tmp.merge(meta, on="id", validate="m:1")
summary = (
    tmp.groupby("g")
    .agg(n=("x", "size"))
    .reset_index()
)


def helper(d):
    return d


helper(tmp)
import numpy as np
`;

test('Python measures', () => {
  const m = codeMeasures(PY, true);
  assert.strictEqual(m.joins_unkeyed, 1);
  assert.strictEqual(m.warning_suppress, 1, 'a category filter is allowed');
  assert.strictEqual(m.longest_pipe, 3, 'two .method lines after the first');
  assert.strictEqual(m.vague_names, 1, 'tmp');
  assert.strictEqual(m.wrappers_once, 1);
  assert.strictEqual(m.session_info, false);
  assert.strictEqual(m.late_loads, 1, 'import numpy below code');
  assert.strictEqual(m.base_df_calls, null);
  assert.strictEqual(m.upper_constants, null);
});
