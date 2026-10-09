#' Winsorise and 0-1 normalise histone-mark intensities
#'
#' Clips each column at a fixed *number of points* from each end rather than a
#' quantile, then rescales to 0-1. Clipping by point count keeps the scaling
#' stable when one sample has a handful of extreme outliers, which quantile
#' winsorising does not.
#'
#' @param histone_rlog Features x samples data frame of rlog or VST values.
#' @param median_subtract If `TRUE`, subtract each column's median before
#'   winsorising, to align distributions that sit at different baselines.
#' @param low,high Number of points to clip from the bottom and top of each
#'   column. Set either to `0` to leave that end alone.
#' @param normalize If `TRUE` (default) rescale each column to 0-1 after
#'   winsorising. Set `FALSE` to get the winsorised values on their original
#'   scale, which is the only way `median_subtract` changes the result — see
#'   below.
#'
#' @return A data frame the same shape as the input: 0-1 scaled by default, or
#'   winsorised on the original scale if `normalize = FALSE`.
#'
#' @section Fixes:
#' `median_subtract` was applied to the wrong matrix: the median-subtracted
#' values were computed into `histone_rlog_medsubt` and then the *original*
#' `histone_rlog` was winsorised, so the argument was ignored. The
#' normalisation also ignored `NA`.
#'
#' @section On median_subtract:
#' Even applied correctly, `median_subtract` cannot change the default output.
#' Subtracting a column median is a constant shift, and the 0-1 min-max
#' rescaling that follows is shift-invariant, so the two cancel exactly. The
#' argument only has an effect with `normalize = FALSE`. It is kept, correctly
#' wired, for that case and for call compatibility.
#'
#' @seealso [trim_upper_tail()], [vsd_norm()]
#'
#' @examples
#' \dontrun{
#' normed <- wins_norm_histones(histone_rlog, median_subtract = TRUE)
#' }
#' @export
wins_norm_histones <- function(histone_rlog, median_subtract = FALSE,
                               low = 500, high = 20, normalize = TRUE) {
  require_pkg("pheatmap", "DescTools")


  

  
  
  
# from http://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode



estimate_mode <- function(x) {
  d <- stats::density(x)
  d$x[which.max(d$y)]
}

col_metrics <- function(x) {
  range_df <- data.frame(mins = apply(x, 2, min), 
                         maxes = apply(x, 2, max),
                         median = apply(x, 2, stats::median),
                         mode = apply(x, 2, estimate_mode))
  return(range_df)
}                        

# winsorize the data 

winsorize_by_points <- function(df, low_wins, high_wins) {
  # Handle edge case: if low_wins <= 0, don't clip lower end
  if (low_wins > 0) {
    low_cut <- apply(df, 2, function(x) sort(x)[low_wins])
  } else {
    low_cut <- apply(df, 2, min, na.rm = TRUE)
  }

  # Handle edge case: if high_wins <= 0, don't clip upper end
  if (high_wins > 0) {
    high_cut <- apply(df, 2, function(x) sort(x, decreasing = TRUE)[high_wins])
  } else {
    high_cut <- apply(df, 2, max, na.rm = TRUE)
  }

  # Use pmax/pmin for manual clipping (avoids DescTools API changes)
  wins_df <- sapply(1:ncol(df), function(i) {
    x <- df[,i]
    pmin(pmax(x, low_cut[i]), high_cut[i])
  })
  colnames(wins_df) <- colnames(df)
  row.names(wins_df) <- row.names(df)

  return(wins_df)

}

if (median_subtract) {
  # median subtract the data to better align the transformed values 
  histone_rlog_medsubt <- sweep(histone_rlog, 2, apply(histone_rlog, 2, stats::median), `-`)
} else {
  histone_rlog_medsubt <- histone_rlog
}

# The original winsorized `histone_rlog`, not `histone_rlog_medsubt`, so the
# median_subtract argument had no effect on the returned values at all.
histone_rlog_medsubt_wins <- winsorize_by_points(histone_rlog_medsubt, low, high)
if (normalize) {
  unit_scale <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / diff(rng)
  }
  out <- apply(histone_rlog_medsubt_wins, 2, unit_scale)
} else {
  out <- histone_rlog_medsubt_wins
}

return(as.data.frame(out))

} # end of function 
