#' Scatter two samples with their binding cutoffs drawn as a box
#'
#' Plots one sample against another, overlays the rectangle bounded by each
#' sample's `cutoff_50` and upper-tail cutoff (the "bound, not tail" window),
#' highlights a gene set, labels named genes, and adds marginal densities.
#'
#' @param data Data frame of signal with one column per sample, plus
#'   `diff_gene` and `mgi_symbol` columns.
#' @param x_col,y_col Column names to plot.
#' @param x_condition,y_condition Replicate names used to look the cutoffs up
#'   in `annotation`.
#' @param annotation Cutoff table with `replicate_names`, `cutoff_50` and
#'   `CrossingSmoothedValue`, e.g. from [chip_binarize_tail_removed()].
#' @param label_genes Character vector of `mgi_symbol` values to label. Was
#'   previously read from the global environment, so the function only worked
#'   if you happened to have a `label_genes` object lying around.
#' @param mark_genes Genes drawn with an x marker and always labelled.
#' @param mark_category Only genes whose `max_cat` equals this get the hollow
#'   circle marker. Pass `NULL` to skip that layer, which also lets `data`
#'   omit a `max_cat` column.
#' @param winsorize_max Upper quantile for [DescTools::Winsorize()].
#' @param highlight_colour,other_colour Point colours for the two
#'   `diff_gene` levels.
#'
#' @return A `ggExtra` object (a `gtable`), as returned by
#'   [ggExtra::ggMarginal()]. Print it to draw.
#'
#' @examples
#' \dontrun{
#' chip_scatter_cutoffs(
#'   signal, "H3K27me3_P7_rep1", "H3K27me3_P56_rep1",
#'   "H3K27me3_P7_rep1", "H3K27me3_P56_rep1",
#'   annotation = breaks_table,
#'   label_genes = c("Gli1", "Ccnd1", "Mycn")
#' )
#' }
#' @import ggplot2
#' @importFrom rlang .data
#' @export
chip_scatter_cutoffs <- function(data, x_col, y_col, x_condition, y_condition,
                                 annotation,
                                 label_genes = character(0),
                                 mark_genes = character(0),
                                 mark_category = NULL,
                                 winsorize_max = 0.9,
                                 highlight_colour = "red",
                                 other_colour = "grey") {

  require_pkg("DescTools", "ggExtra", "ggrepel")

  cutoff_window <- function(condition) {
    row <- annotation[annotation$replicate_names == condition, , drop = FALSE]
    if (nrow(row) == 0L) {
      stop("no annotation row for replicate ", sQuote(condition), call. = FALSE)
    }
    c(min = as.numeric(row$cutoff_50[1]),
      max = as.numeric(row$CrossingSmoothedValue[1]))
  }

  x_cutoffs <- cutoff_window(x_condition)
  y_cutoffs <- cutoff_window(y_condition)

  wins <- function(v) DescTools::Winsorize(v, maxval = winsorize_max)

  plot_data <- data.frame(
    x_value    = wins(data[[x_col]]),
    y_value    = wins(data[[y_col]]),
    diff_gene  = as.factor(data$diff_gene),
    mgi_symbol = data$mgi_symbol,
    stringsAsFactors = FALSE
  )
  if (!is.null(mark_category)) {
    plot_data$max_cat <- data$max_cat
  }

  axis_label <- function(nm) gsub("_", " ", gsub("rep.*", "", nm))

  colour_values <- c(highlight_colour, other_colour)
  names(colour_values) <- levels(plot_data$diff_gene)[seq_along(colour_values)]

  p <- ggplot(plot_data, aes(x = .data$x_value, y = .data$y_value,
                             fill = .data$diff_gene, colour = .data$diff_gene)) +
    annotate("rect",
             xmin = x_cutoffs[["min"]], xmax = x_cutoffs[["max"]],
             ymin = y_cutoffs[["min"]], ymax = y_cutoffs[["max"]],
             fill = NA, colour = "black") +
    geom_point(aes(alpha = .data$diff_gene, size = .data$diff_gene)) +
    scale_colour_manual(values = colour_values) +
    scale_fill_manual(values = colour_values) +
    scale_alpha_manual(values = stats::setNames(c(1, 0.3)[seq_along(colour_values)],
                                                names(colour_values))) +
    scale_size_manual(values = stats::setNames(c(1, 0.5)[seq_along(colour_values)],
                                               names(colour_values))) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = axis_label(x_col), y = axis_label(y_col))

  if (length(label_genes)) {
    p <- p + ggrepel::geom_text_repel(
      data = plot_data[plot_data$mgi_symbol %in% label_genes, , drop = FALSE],
      aes(label = .data$mgi_symbol),
      colour = "black", point.padding = 0.2, box.padding = 0.2,
      min.segment.length = 0, size = 4, seed = 42, alpha = 1, force = 10
    )
  }

  if (!is.null(mark_category) && length(label_genes)) {
    circled <- plot_data[plot_data$mgi_symbol %in% label_genes &
                           plot_data$max_cat == mark_category, , drop = FALSE]
    if (nrow(circled)) {
      p <- p + geom_point(data = circled, colour = "black", shape = 1,
                          size = 2, alpha = 1, inherit.aes = TRUE)
    }
  }

  if (length(mark_genes)) {
    marked <- plot_data[plot_data$mgi_symbol %in% mark_genes, , drop = FALSE]
    if (nrow(marked)) {
      p <- p +
        geom_point(data = marked, colour = "black", shape = 4, size = 2, alpha = 1) +
        ggrepel::geom_text_repel(
          data = marked, aes(label = .data$mgi_symbol),
          colour = "black", point.padding = 0.2, box.padding = 0.2,
          min.segment.length = 0, size = 4, seed = 42, alpha = 1, force = 10
        )
    }
  }

  ggExtra::ggMarginal(p, type = "density", groupFill = TRUE)
}
