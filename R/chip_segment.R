#' Segment ChIP signal into classes, end to end
#'
#' Convenience wrapper that runs the two-step segmentation used throughout the
#' GNP/MB ChIP-seq analysis: find the per-column upper-tail cutoff with
#' [trim_upper_tail()], then cluster the values below that cutoff with
#' [chip_mclust_icl_all()], giving the tail its own top class.
#'
#' @param data_frame Data frame of signal; one column per sample, one row per
#'   feature.
#' @param deriv_cutoff Derivative threshold handed to [trim_upper_tail()].
#' @param G_range Component counts handed to [chip_mclust_icl_all()].
#' @param plot If `TRUE` (default) draw the rank and histogram diagnostics
#'   from [trim_upper_tail()].
#'
#' @return A list with `cutoffs` (the `trim_upper_tail()` table) and
#'   `classes` (`data_frame` plus one `<column>_class` column per sample).
#'
#' @examples
#' \dontrun{
#' seg <- chip_segment(histone_rlog)
#' seg$cutoffs
#' head(seg$classes)
#' }
#' @export
chip_segment <- function(data_frame,
                         deriv_cutoff = 3e-4,
                         G_range = 2:8,
                         plot = TRUE) {

  # The original body ignored both of its arguments and read a global called
  # `histone_rlog_prune_mean_only`, referenced `combined_data` fifty lines
  # before creating it, and carried a second, stale copy of
  # chip_mclust_icl_all() nested inside itself. It could not run as written.
  # This version is the pipeline it was reaching for.

  require_pkg("mclust")
  stopifnot(is.data.frame(data_frame))

  cutoffs <- trim_upper_tail(data_frame,
                             derivative_threshold = deriv_cutoff,
                             generate_plots = plot)

  classes <- chip_mclust_icl_all(data_frame,
                                 crossing_points_df = cutoffs,
                                 G_range = G_range)

  list(cutoffs = cutoffs, classes = classes)
}
