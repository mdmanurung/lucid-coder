# Figure theme and panel export, shared by every figure script.
library(ggplot2)

# base_size 7: journal panels are read at print size, where 7 pt is the usual minimum.
theme_set(theme_classic(base_size = 7))

# Write one panel as <run_dir>/<name>.pdf (vector, for the figure), .png (for
# review) and .csv (the source data, the same frame that was passed to ggplot()).
save_panel <- function(p, data, run_dir, name, width_mm, height_mm) {
  if (!dir.exists(run_dir)) {
    stop(
      "save_panel: run_dir ",
      run_dir,
      " does not exist; create it before saving.",
      call. = FALSE
    )
  }
  stem <- file.path(run_dir, name)
  # No dpi on the PDF: vector output ignores it.
  ggsave(
    paste0(stem, ".pdf"),
    p,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = cairo_pdf
  )
  ggsave(
    paste0(stem, ".png"),
    p,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 300
  )
  write.csv(data, paste0(stem, ".csv"), row.names = FALSE)
  invisible(stem)
}
