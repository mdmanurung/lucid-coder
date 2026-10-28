# The seven common plot forms, drawn the same way by every figure script.
# Source theme.R and palettes.R first; each function returns a ggplot object and
# never saves it, so the caller passes it to save_panel().
# No title, subtitle or caption is drawn; the caption goes in the report (fig.text).
# Columns are passed as strings. No statistics beyond the lm fit in
# plot_scatter_fit(): stars and p-values are the caller's job.
# Colours, theme and dashed threshold lines are fixed house style. Every tunable
# value is a required argument; a grill asks the user for it per project:
#   padj_max:     which adjusted-p cutoff marks a gene significant?
#   lfc_min:      which absolute log2 fold change marks a gene significant?
#   limit:        at which absolute value should the heatmap colour scale saturate?
#   jitter_width: how wide should the jitter spread points within a box?
#   seed:         which seed fixes the jitter so reruns draw the same panel?
#   max_size:     what is the largest dot size in the dotplot, in mm?

# fig.wilke: under this many points per group a box misleads; draw points + median.
FEW_POINTS <- 10

# Stop with the function name when a required column is absent.
check_cols <- function(data, cols, fn) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop(sprintf("%s: missing columns: %s", fn, toString(missing)), call. = FALSE)
  }
  invisible(data)
}

# Factor in palette order; palette_for() stops first on any label the palette lacks.
factor_by_pal <- function(x, pal) {
  palette_for(x, pal)
  factor(x, levels = names(pal))
}

# Significant = padj < padj_max and |lfc| > lfc_min; both cut-offs are required.
plot_volcano <- function(data, lfc, padj, padj_max, lfc_min, label = NULL) {
  check_cols(data, c(lfc, padj, label), "plot_volcano")
  data$significant <- data[[padj]] < padj_max & abs(data[[lfc]]) > lfc_min
  p <- ggplot(data, aes(.data[[lfc]], -log10(.data[[padj]]), colour = significant)) +
    geom_point(size = 0.5) +
    geom_hline(yintercept = -log10(padj_max), linetype = "dashed") +
    geom_vline(xintercept = c(-lfc_min, lfc_min), linetype = "dashed") +
    scale_colour_manual(values = c(`TRUE` = "darkred", `FALSE` = "grey60"), guide = "none") +
    labs(x = "log2 fold change", y = "-log10 adjusted p")
  if (is.null(label)) {
    return(p)
  }
  p +
    geom_text(
      aes(label = .data[[label]]),
      data = data[data$significant, ],
      vjust = -0.5,
      size = 2.5
    )
}

# Box per group with the rows on top; pair joins the same id across groups.
plot_box <- function(data, group, value, pal, jitter_width, seed, value_label, pair = NULL) {
  check_cols(data, c(group, value, pair), "plot_box")
  data[[group]] <- factor_by_pal(data[[group]], pal)
  n <- table(data[[group]])
  p <- ggplot(data, aes(.data[[group]], .data[[value]])) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_colour_manual(values = pal, guide = "none") +
    scale_x_discrete(labels = \(g) sprintf("%s (n = %d)", g, as.integer(n[g]))) +
    labs(x = NULL, y = value_label)
  # fig.wilke: under FEW_POINTS per group a box misleads; draw points and a median.
  p <- if (min(n) < FEW_POINTS) {
    p +
      stat_summary(
        aes(colour = .data[[group]]),
        fun = median,
        fun.min = median,
        fun.max = median,
        geom = "crossbar",
        width = 0.5
      )
  } else {
    p + geom_boxplot(aes(fill = .data[[group]]), outlier.shape = NA)
  }
  if (is.null(pair)) {
    return(
      p + geom_point(position = position_jitter(width = jitter_width, height = 0, seed = seed))
    )
  }
  # No jitter with pairs, so each line meets its points.
  p +
    geom_line(aes(group = .data[[pair]]), colour = "grey60") +
    geom_point()
}

# UMAP or similar; scattermore rasterises dense point clouds when installed.
# `axis_name` names the embedding ("UMAP"); axes read "UMAP 1" and "UMAP 2".
plot_embedding <- function(data, x, y, colour, pal, axis_name) {
  check_cols(data, c(x, y, colour), "plot_embedding")
  data[[colour]] <- factor_by_pal(data[[colour]], pal)
  points <- if (requireNamespace("scattermore", quietly = TRUE)) {
    scattermore::geom_scattermore(pointsize = 2)
  } else {
    geom_point(size = 0.3)
  }
  ggplot(data, aes(.data[[x]], .data[[y]], colour = .data[[colour]])) +
    points +
    scale_colour_manual(values = pal) +
    coord_equal() +
    labs(x = paste(axis_name, 1), y = paste(axis_name, 2), colour = colour) +
    theme(axis.ticks = element_blank(), axis.text = element_blank())
}

# Diverging fill centred at 0 with symmetric limits; values beyond limit are squished.
# Row and column order are the caller's factor levels: no clustering here.
plot_heatmap <- function(data, row, col, value, limit) {
  check_cols(data, c(row, col, value), "plot_heatmap")
  ggplot(data, aes(.data[[col]], .data[[row]], fill = .data[[value]])) +
    geom_tile() +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-limit, limit),
      oob = scales::squish
    ) +
    labs(x = col, y = row, fill = value) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}

# Dot size = fraction expressing (0-1), colour = mean expression.
plot_dotplot <- function(data, feature, group, mean, frac, max_size) {
  check_cols(data, c(feature, group, mean, frac), "plot_dotplot")
  ggplot(data, aes(.data[[feature]], .data[[group]])) +
    geom_point(aes(size = .data[[frac]], colour = .data[[mean]])) +
    scale_size_area(max_size = max_size) +
    scale_colour_gradient(low = "grey90", high = "#B2182B") +
    labs(x = feature, y = group, size = frac, colour = mean) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}

# Stacked bar of each label's share of the count per sample.
plot_composition <- function(data, sample, label, count, pal) {
  check_cols(data, c(sample, label, count), "plot_composition")
  data[[label]] <- factor_by_pal(data[[label]], pal)
  ggplot(data, aes(.data[[sample]], .data[[count]], fill = .data[[label]])) +
    geom_col(position = "fill") +
    scale_fill_manual(values = pal) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = sample, y = sprintf("%% of %s", count), fill = label)
}

# Points with a linear fit and its 95% confidence band.
plot_scatter_fit <- function(data, x, y, label = NULL) {
  check_cols(data, c(x, y, label), "plot_scatter_fit")
  p <- ggplot(data, aes(.data[[x]], .data[[y]])) +
    geom_point(size = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, level = 0.95) +
    labs(x = x, y = y)
  if (is.null(label)) {
    return(p)
  }
  p + geom_text(aes(label = .data[[label]]), vjust = -0.5, size = 2.5)
}

if (sys.nframe() == 0) {
  source("theme.R")
  source("palettes.R")
  # A null device, so no stray Rplots.pdf lands in templates/.
  pdf(NULL)
  donors <- data.frame(
    donor_id = rep(sprintf("D%02d", 1:4), times = 2),
    condition = rep(names(pal_condition), each = 4),
    log2_expr = c(5.1, 4.8, 5.4, 5.0, 8.2, 7.9, 8.5, 8.1)
  )
  de <- data.frame(
    gene = c("IFI27", "ACTB", "CD3E", "MX1"),
    log2_fc = c(3.1, 0.1, -1.2, 2.0),
    q_value = c(0.001, 0.9, 0.04, 0.2)
  )
  cells <- data.frame(
    umap_1 = c(0.1, 0.5, 1.2, 2.0, 2.4, 3.1),
    umap_2 = c(1.0, 0.2, 2.2, 1.5, 0.4, 2.8),
    cell_type = rep(names(pal_cell_type), times = 2)
  )
  grid <- tidyr::expand_grid(
    gene = c("IFI27", "ACTB"),
    cell_type = names(pal_cell_type)
  )
  grid$z <- c(-5, -1, 0, 1, 2, 5)
  grid$frac <- c(0.1, 0.5, 0.9, 0.3, 0.7, 1.0)
  counts <- data.frame(
    sample = rep(c("S1", "S2"), each = 3),
    cell_type = rep(names(pal_cell_type), times = 2),
    n_cells = c(10, 20, 30, 5, 5, 40)
  )
  panels <- list(
    volcano = list(plot_volcano(de, "log2_fc", "q_value", 0.05, 1, label = "gene"), de),
    box = list(
      plot_box(donors, "condition", "log2_expr", pal_condition, 0.2, 1, "log2 expression"),
      donors
    ),
    box_paired = list(
      plot_box(
        donors,
        "condition",
        "log2_expr",
        pal_condition,
        0.2,
        1,
        "log2 expression",
        pair = "donor_id"
      ),
      donors
    ),
    embedding = list(
      plot_embedding(cells, "umap_1", "umap_2", "cell_type", pal_cell_type, "UMAP"),
      cells
    ),
    heatmap = list(plot_heatmap(grid, "cell_type", "gene", "z", 4), grid),
    dotplot = list(plot_dotplot(grid, "gene", "cell_type", "z", "frac", 4), grid),
    composition = list(
      plot_composition(counts, "sample", "cell_type", "n_cells", pal_cell_type),
      counts
    ),
    scatter_fit = list(plot_scatter_fit(de, "log2_fc", "q_value", label = "gene"), de)
  )
  run_dir <- tempfile()
  dir.create(run_dir)
  for (name in names(panels)) {
    save_panel(panels[[name]][[1]], panels[[name]][[2]], run_dir, name, 89, 60)
    stopifnot(file.exists(file.path(run_dir, paste0(name, ".pdf"))))
  }
  stopifnot(inherits(try(plot_volcano(de, "lfc", "q_value", 0.05, 1), silent = TRUE), "try-error"))
  donors$condition[[1]] <- "unknown"
  stopifnot(inherits(
    try(
      plot_box(donors, "condition", "log2_expr", pal_condition, 0.2, 1, "log2 expression"),
      silent = TRUE
    ),
    "try-error"
  ))
  invisible(dev.off())
}
