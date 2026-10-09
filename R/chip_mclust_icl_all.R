#' Cluster every column of a signal matrix with mclust, choosing k by ICL
#'
#' For each column, fits Gaussian mixtures over `G_range` components, selects
#' the model with the best integrated complete-data likelihood (ICL), and
#' returns the per-observation classification. Values above a per-column tail
#' cutoff are held out of the fit and assigned their own top class, so that a
#' long right tail of very highly bound genes does not dominate the mixture.
#'
#' @param data Data frame or matrix of signal; one column per sample.
#' @param crossing_points_df Tail cutoffs as returned by [trim_upper_tail()],
#'   with columns `Column` and `CrossingSmoothedValue`. Pass `NULL` to cluster
#'   every value.
#' @param G_range Component counts to compare.
#'
#' @return `data` with one extra `<column>_class` integer column per input
#'   column.
#'
#' @seealso [trim_upper_tail()] to produce `crossing_points_df`,
#'   [chip_segment()] for the two steps in one call.
#'
#' @examples
#' \dontrun{
#' cutoffs <- trim_upper_tail(histone_rlog)
#' classes <- chip_mclust_icl_all(histone_rlog, cutoffs)
#' }
#' @export
chip_mclust_icl_all <- function(data, crossing_points_df = NULL, G_range = 2:8) {

  require_pkg("mclust")

  # mclust::Mclust() and mclust::mclustICL() build a call to mclustBIC and
  # resolve it in the CALLER's frame, so they fail with
  # `could not find function "mclustBIC"` unless mclust is attached. Binding it
  # here makes them work with mclust in Suggests rather than Imports.
  mclustBIC <- mclust::mclustBIC

  results <- list()

  for (column_name in names(data)) {
    column_data <- data[[column_name]]

    tail_cutoff <- NULL
    if (!is.null(crossing_points_df)) {
      hit <- crossing_points_df[crossing_points_df$Column == column_name,
                                "CrossingSmoothedValue"]
      if (length(hit) > 0 && is.finite(as.numeric(hit[1]))) {
        tail_cutoff <- as.numeric(hit[1])
      }
    }

    above <- if (is.null(tail_cutoff)) {
      integer(0)
    } else {
      which(column_data > tail_cutoff)
    }

    # `x[-integer(0)]` returns an EMPTY vector, not x. The original relied on
    # negative indexing here, so any column whose values all fell below the
    # tail cutoff was clustered on nothing and came back entirely NA.
    below <- setdiff(seq_along(column_data), above)
    data_to_process <- column_data[below]
    data_to_process <- data_to_process[is.finite(data_to_process)]

    if (length(data_to_process) < max(G_range)) {
      warning("column ", column_name, ": only ", length(data_to_process),
              " usable values, skipping", call. = FALSE)
      results[[column_name]] <- rep(NA_integer_, length(column_data))
      next
    }

    mclust_data <- data.frame(column_data = data_to_process)
    icl_results <- mclust::mclustICL(mclust_data, G = G_range)

    best <- which(icl_results == max(icl_results, na.rm = TRUE), arr.ind = TRUE)[1, ]
    best_g <- as.integer(dimnames(icl_results)[[1]][best[["row"]]])
    best_model <- dimnames(icl_results)[[2]][best[["col"]]]

    # `model =` is not an Mclust argument; it was silently swallowed by `...`
    # and the model family was re-selected by BIC instead of being fixed to
    # the one ICL chose.
    fitted <- mclust::Mclust(mclust_data, G = best_g,
                             modelNames = best_model, verbose = FALSE)

    final <- rep(NA_integer_, length(column_data))
    final[below[is.finite(column_data[below])]] <- fitted$classification
    if (length(above) > 0) {
      final[above] <- max(final, na.rm = TRUE) + 1L
    }

    results[[column_name]] <- final
  }

  results_df <- as.data.frame(results)
  colnames(results_df) <- paste(names(results), "class", sep = "_")

  cbind(data, results_df)
}
