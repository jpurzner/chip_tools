#' Map t50 from the interpolated axis back onto sampled time points
#'
#' A t50 computed on the dense interpolated time axis is in units of
#' interpolated index, which is awkward to plot against a handful of real
#' sampling ages. This builds the mapping from the full axis to the compressed
#' axis of actually-sampled points and applies it.
#'
#' @param t50 Numeric vector of t50 values on the interpolated axis.
#' @param time_ref The `time_ref` logical frame passed to [time_expr_interp()].
#'   Columns one and two are the two data sets' sampling masks.
#'
#' @return `t50` expressed on the compressed axis.
#'
#' @seealso [expr_half_max()], [time_expr_interp()]
#'
#' @examples
#' \dontrun{
#' t50_compressed <- t50_to_compressed(t50, time_ref[, c("scott", "hatten")])
#' }
#' @export
t50_to_compressed <- function (t50, time_ref) {
  
  t_ref <- time_ref
  t_ref$i <- 1:nrow(time_ref)
  t_ref$short_i <- cumsum(time_ref[,1] | time_ref[,2])
  t_ref$short_i[!(time_ref[,1] | time_ref[,2])] <- NA
  t_ref$short_i <- zoo::na.approx(t_ref$short_i)
  af <- stats::approxfun(t_ref$i, t_ref$short_i)
  return(af(t50))
} 