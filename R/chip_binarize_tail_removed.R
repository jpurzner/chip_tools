#' Binarize each column after dropping its upper tail
#'
#' Runs [binarize_counts()] per column on the values *below* that column's
#' upper-tail cutoff, and returns the cutoff table joined to the resulting
#' mixture summaries. Fitting the mixture without the tail keeps a handful of
#' very highly bound features from pulling the signal component upward.
#'
#' @param data_frame Data frame of signal; one column per sample.
#' @param tail_cutoff_table Cutoffs from [trim_upper_tail()], with columns
#'   `Column` and `CrossingSmoothedValue`.
#' @param plot If `TRUE` (default) draw the per-column mixture fit.
#'
#' @return `tail_cutoff_table` with `cutoff_50`, `cutoff_75`, `mean_upper`,
#'   `mean_lower`, `loglik`, `aic` and `bic` joined on `Column`. Columns whose
#'   fit failed are dropped, with a warning naming them.
#'
#' @seealso [trim_upper_tail()], [binarize_counts()]
#'
#' @examples
#' \dontrun{
#' cutoffs <- trim_upper_tail(histone_rlog)
#' breaks  <- chip_binarize_tail_removed(histone_rlog, cutoffs)
#' }
#' @export
chip_binarize_tail_removed <- function(data_frame, tail_cutoff_table, plot = TRUE) {

  results <- list()
  failed <- character(0)

  for (column_name in names(data_frame)) {
    cutoff_value <- tail_cutoff_table[tail_cutoff_table$Column == column_name,
                                      "CrossingSmoothedValue"]

    if (length(cutoff_value) == 0L || !is.finite(as.numeric(cutoff_value[1]))) {
      failed <- c(failed, column_name)
      next
    }

    trimmed <- data_frame[[column_name]]
    trimmed <- trimmed[is.finite(trimmed) & trimmed < as.numeric(cutoff_value[1])]

    if (length(trimmed) < 10L) {
      failed <- c(failed, column_name)
      next
    }

    max_val <- max(trimmed)

    # A failed mixture fit for one sample should not abort the whole table.
    fit <- tryCatch(
      binarize_counts(trimmed,
                      return_what = "summary",
                      exclude_extreme = TRUE,
                      label_name = column_name,
                      plot = plot,
                      mu = c(0.5, max_val - 0.2)),
      error = function(e) {
        warning("column ", column_name, ": ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )

    if (is.null(fit)) {
      failed <- c(failed, column_name)
      next
    }
    results[[column_name]] <- fit
  }

  if (length(failed)) {
    warning("no mixture summary for: ", paste(failed, collapse = ", "),
            call. = FALSE)
  }
  if (length(results) == 0L) {
    stop("no column could be binarized", call. = FALSE)
  }

  summaries <- as.data.frame(do.call(rbind, results))
  summaries$Column <- names(results)
  rownames(summaries) <- NULL

  tail_cutoff_table$Column <- as.character(tail_cutoff_table$Column)
  merge(tail_cutoff_table, summaries, by = "Column")
}
