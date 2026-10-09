#' Row-normalise a count table for heatmap plotting
#'
#' Divides each gene by its own row mean and takes `log2`, giving a
#' fold-change-from-average matrix that is symmetric about zero — which is
#' what makes a diverging heatmap palette read correctly. Genes whose maximum
#' falls below `cutoff` are dropped first.
#'
#' Zeros become `-Inf` under the log, so they are clamped: either to the
#' negative of the matrix maximum (the default, which keeps the colour scale
#' symmetric) or to the smallest finite non-zero value present.
#'
#' @param data Genes x samples data frame of expression or signal.
#' @param cutoff Genes whose row maximum is below this are dropped.
#' @param normbymax If `TRUE` (default) clamp `-Inf` and anything below it to
#'   `-max(log2 values)`, keeping the palette symmetric. If `FALSE`, clamp to
#'   the smallest finite non-zero value instead.
#' @param verbose If `TRUE` print the clamp value. The original always printed.
#'
#' @return A genes x samples matrix of row-normalised `log2` values, ready to
#'   hand to `pheatmap()`.
#'
#' @section Fixes:
#' `max(norm_log2)` was taken over a matrix that still held `-Inf`, and the
#' row maximum was computed without `na.rm`, so a single `NA` dropped the gene
#' via the cutoff rather than via the `complete.cases()` step. The nested
#' `min_finite()` also shadowed the package function of the same name.
#'
#' @seealso [wins_norm_histones()] for the winsorise-and-rescale alternative,
#'   [heatmap_table()]
#'
#' @examples
#' \dontrun{
#' pheatmap::pheatmap(plot_norm_heatmap(tpm, cutoff = 5))
#' }
#' @export
plot_norm_heatmap <- function(data, cutoff = 1, normbymax = TRUE, verbose = FALSE) {

  # Local helper, kept rather than calling the package's min_finite(): this
  # one also excludes exact zeros, which matters because a zero here would
  # become the new floor for every -Inf cell.
  min_finite_nonzero <- function(m) {
    vals <- m[is.finite(m) & m != 0]
    if (length(vals) == 0L) {
      stop("no finite non-zero values after normalisation", call. = FALSE)
    }
    min(vals)
  }

  # Work on a matrix throughout. The original stayed on a data frame, where
  # max() happens to work but is.finite() does not, so any attempt to screen
  # out -Inf before taking the maximum failed with "default method not
  # implemented for type 'list'".
  mat <- as.matrix(data)

  # drop genes whose maximum is below the cutoff, then any incomplete rows
  row_max <- apply(mat, 1, max, na.rm = TRUE)
  mat <- mat[row_max > cutoff, , drop = FALSE]
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  if (nrow(mat) == 0L) {
    stop("no gene has a maximum above cutoff = ", cutoff, call. = FALSE)
  }

  norm_log2 <- log2(sweep(mat, 1, rowMeans(mat), FUN = "/"))

  if (normbymax) {
    # Clamp symmetrically about zero so a diverging palette reads correctly.
    finite_vals <- norm_log2[is.finite(norm_log2)]
    if (length(finite_vals) == 0L) {
      stop("no finite values after normalisation", call. = FALSE)
    }
    max_val <- -1 * max(finite_vals)
    if (verbose) print(max_val)
    norm_log2[!is.finite(norm_log2) & norm_log2 < 0] <- max_val
    norm_log2[norm_log2 < max_val] <- max_val
  } else {
    min_val <- min_finite_nonzero(norm_log2)
    norm_log2[!is.finite(norm_log2) & norm_log2 < 0] <- min_val
  }

  norm_log2
}
