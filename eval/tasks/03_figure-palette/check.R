# Purpose: gate for task 03. The figure and its source data exist; the data has one row per
#          donor x cell type, with the input condition and proportion for each; the PNG shows
#          the lab colour of both conditions and none of the cell-type colours.
# Inputs:  <run_repo>/results/run1/proportions_by-condition.{pdf,png,csv},
#          <run_repo>/data/{proportions,palette}.tsv, <task_dir>/truth/donor_cell_types.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).

# ---- Setup ----
PANEL <- "proportions_by-condition"
CONDITIONS <- c("healthy", "disease")
# Per-channel tolerance (0-255) when matching a pixel to a palette colour; anti-aliasing blends edges.
PIXEL_TOLERANCE <- 6
# Fewest pixels of a colour that count as drawn.
MIN_PIXELS <- 20
# Horizontal gap, as a share of image width, that starts a new column of points.
GROUP_GAP <- 0.08
# Fewest pixels of each condition colour in a column for it to count (drops legend keys).
MIN_GROUP_PIXELS <- 60

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  message("usage: Rscript check.R <run_repo> <task_dir>")
  quit(status = 1)
}
run_repo <- args[[1]]
task_dir <- args[[2]]

gate <- function(pass, detail) {
  verdict <- list(gate = if (pass) "pass" else "fail", detail = detail)
  cat(jsonlite::toJSON(verdict, auto_unbox = TRUE), "\n", sep = "")
  quit(status = 0)
}

# ---- Load truth ----
truth <- read.delim(file.path(task_dir, "truth", "donor_cell_types.tsv"))

# ---- Check the panel files ----
panel_files <- paste0(PANEL, c(".pdf", ".png", ".csv"))
panel_paths <- file.path(run_repo, "results", "run1", panel_files)
# file.size is NA for a missing file.
is_absent <- is.na(file.size(panel_paths)) | file.size(panel_paths) == 0
if (any(is_absent)) {
  gate(FALSE, paste("missing or empty in results/run1:", toString(panel_files[is_absent])))
}
source_data <- tryCatch(read.csv(panel_paths[[3]]), error = function(e) NULL)
if (is.null(source_data)) {
  gate(FALSE, paste(panel_files[[3]], "is not a readable CSV"))
}
if (nrow(source_data) != nrow(truth)) {
  gate(
    FALSE,
    sprintf(
      "%s has %d rows, expected %d (one per donor x cell type)",
      panel_files[[3]],
      nrow(source_data),
      nrow(truth)
    )
  )
}

# ---- Check the data against the inputs ----
inputs <- read.delim(file.path(run_repo, "data", "proportions.tsv"))
needed <- c("donor_id", "condition", "cell_type", "proportion")
if (!all(needed %in% names(source_data))) {
  gate(FALSE, paste(panel_files[[3]], "lacks columns", toString(setdiff(needed, names(source_data)))))
}
matched <- merge(source_data[needed], inputs, by = c("donor_id", "cell_type"), suffixes = c("", ".input"))
if (nrow(matched) != nrow(truth)) {
  gate(FALSE, "CSV donor x cell type rows do not match the input rows")
}
wrong_condition <- matched$condition != matched$condition.input
if (any(wrong_condition)) {
  gate(FALSE, paste("wrong condition for:", toString(unique(matched$donor_id[wrong_condition]))))
}
if (any(abs(matched$proportion - matched$proportion.input) > 1e-6)) {
  gate(FALSE, "CSV proportions differ from data/proportions.tsv")
}

# ---- Check the palette in the PNG ----
palette <- read.delim(file.path(run_repo, "data", "palette.tsv"))
hex <- setNames(palette$hex, palette$label)
image <- png::readPNG(panel_paths[[2]])
pixels <- round(matrix(image[, , 1:3], ncol = 3) * 255)
pixels_of <- function(colour) {
  target <- as.integer(col2rgb(colour))
  sum(abs(pixels[, 1] - target[1]) <= PIXEL_TOLERANCE &
    abs(pixels[, 2] - target[2]) <= PIXEL_TOLERANCE &
    abs(pixels[, 3] - target[3]) <= PIXEL_TOLERANCE)
}
drawn <- vapply(hex, pixels_of, numeric(1))
missing_condition <- CONDITIONS[drawn[CONDITIONS] < MIN_PIXELS]
if (length(missing_condition) > 0) {
  gate(FALSE, paste("PNG lacks the lab colour of:", toString(missing_condition)))
}
# Grey (dendritic_cell) is skipped: axes and text are grey too.
cell_type_hex <- hex[!(names(hex) %in% CONDITIONS) & palette$hex != "#BBBBBB"]
wrong_colour <- names(cell_type_hex)[drawn[names(cell_type_hex)] >= MIN_PIXELS]
if (length(wrong_colour) > 0) {
  gate(FALSE, paste("PNG uses cell-type colours, expected condition colours:", toString(wrong_colour)))
}

# ---- Check which condition sits where ----
# Swapped colours keep both hexes in the PNG, so compare heights: group the coloured pixels into
# columns by horizontal gaps, then ask in each group which condition's points sit higher.
pixel_mask <- function(colour) {
  target <- as.integer(col2rgb(colour))
  abs(pixels[, 1] - target[1]) <= PIXEL_TOLERANCE &
    abs(pixels[, 2] - target[2]) <= PIXEL_TOLERANCE &
    abs(pixels[, 3] - target[3]) <= PIXEL_TOLERANCE
}
pixel_x <- rep(seq_len(ncol(image)), each = nrow(image))
pixel_y <- rep(seq_len(nrow(image)), times = ncol(image))
healthy_mask <- pixel_mask(hex[["healthy"]])
disease_mask <- pixel_mask(hex[["disease"]])
xs <- sort(unique(pixel_x[healthy_mask | disease_mask]))
group_id <- cumsum(c(1, diff(xs) > GROUP_GAP * ncol(image)))
group_of_x <- setNames(group_id, xs)
group_of_pixel <- group_of_x[as.character(pixel_x)]
mean_height <- function(mask, g) mean(pixel_y[mask & group_of_pixel %in% g])
# Both conditions need enough pixels in a group; a lone legend key falls out here.
groups <- unique(group_id)
usable <- vapply(
  groups,
  \(g) sum(healthy_mask & group_of_pixel %in% g) >= MIN_GROUP_PIXELS &&
    sum(disease_mask & group_of_pixel %in% g) >= MIN_GROUP_PIXELS,
  logical(1)
)
# Smaller row index is higher on the page.
healthy_higher <- vapply(
  groups[usable],
  \(g) mean_height(healthy_mask, g) < mean_height(disease_mask, g),
  logical(1)
)
type_means <- tapply(inputs$proportion, list(inputs$cell_type, inputs$condition), mean)
expected_share <- mean(type_means[, "healthy"] > type_means[, "disease"])
if (length(healthy_higher) < 3) {
  gate(FALSE, sprintf("PNG shows %d point columns with both condition colours, expected one per cell type", length(healthy_higher)))
}
if ((mean(healthy_higher) > 0.5) != (expected_share > 0.5)) {
  gate(FALSE, "condition colours look swapped: healthy and disease sit the wrong way round against the data")
}
gate(TRUE, sprintf("PDF, PNG and CSV present; CSV has %d rows matching the inputs; PNG uses the condition colours in the right places", nrow(truth)))
