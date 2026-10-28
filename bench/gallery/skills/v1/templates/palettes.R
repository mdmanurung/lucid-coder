# Project palettes and figure choices the house style leaves open.
# Fill each field once per project; every figure script reads this file.
#   dotplot defaults:    [fill in: size range, colour scale, dot.min]
#   volcano thresholds:  [fill in; one seen pair is Q < 0.1 and |log2FC| > 0.5]
#   heatmap anchors:     [fill in; one seen example is -4 / 0 / 4]
#   font family:         [fill in: e.g. Helvetica]
#
# Palettes are named vectors keyed by biological label, never by position.
# The two below are placeholders; the project replaces them with its own labels.
pal_condition <- c(healthy = "#4477AA", disease = "#EE6677")
pal_cell_type <- c(t_cell = "#228833", b_cell = "#CCBB44", monocyte = "#AA3377")

# Colours for a vector of labels. An unknown label stops the script, so a
# mislabelled or renamed group can never fall through as an unplanned colour.
palette_for <- function(labels, pal) {
  # Character first: a factor would index pal by level number, not by name.
  labels <- as.character(labels)
  missing <- setdiff(unique(labels), names(pal))
  if (length(missing) > 0) {
    stop(
      "palette_for: labels not in the palette: ",
      toString(missing),
      "; add them to palettes.R or fix the labels.",
      call. = FALSE
    )
  }
  pal[labels]
}

if (sys.nframe() == 0) {
  stopifnot(identical(unname(palette_for("disease", pal_condition)), "#EE6677"))
  stopifnot(inherits(try(palette_for("unknown", pal_condition), silent = TRUE), "try-error"))
}
