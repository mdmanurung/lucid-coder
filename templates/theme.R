# Figure theme and panel export, shared by every figure script.
library(ggplot2)

# One font hierarchy in pt, read by the theme, text geoms and panel letters
# (user-approved 2026-11-10, after exactplot's xp_init): main 8 for axis and
# legend text, small 6 for in-plot labels (the print floor), tiny 5 only for
# dense secondary text such as a number-at-risk table, large 10 for panel letters.
FONT_PT <- c(main = 8, small = 6, tiny = 5, large = 10)
# Page widths in mm for compose_figure() and save_panel(): one journal column or
# the full double-column page (user-approved 2026-11-10).
PAGE_MM <- c(single = 89, double = 183)
# Embedding corner arrows (exactplot's xp_umap_axis default): 10 mm long, so every UMAP
# in a figure carries the same arrow; the pad leaves room for the axis labels.
ARROW_MM <- 10
ARROW_PAD_MM <- 4

# theme_pubr at the main size: the lab default (user pick 2026-10-28); panel
# letters bold at the large size.
theme_set(
  ggpubr::theme_pubr(base_size = FONT_PT[["main"]]) +
    theme(
      plot.tag = element_text(size = FONT_PT[["large"]], face = "bold"),
      # The letter gets its own margin, so it never overlaps a panel title.
      plot.tag.location = "margin"
    )
)
# geom_text() and geom_label() default to 3.88 mm (about 11 pt); use the small size.
update_geom_defaults("text", list(size = FONT_PT[["small"]] / .pt))
update_geom_defaults("label", list(size = FONT_PT[["small"]] / .pt))

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
