#' Smallest finite value in a data frame
#'
#' Row-wise minimum over finite values only, then the minimum of those.
#' Convenient for picking a plot floor from log-transformed data where zeros
#' have become `-Inf`.
#'
#' @param df Data frame or matrix of numbers.
#'
#' @return A single number, or `Inf` if no value is finite.
#'
#' @examples
#' min_finite(data.frame(a = c(1, -Inf), b = c(3, 2)))
#' @export
min_finite <- function(df) {
  vals <- as.matrix(df)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) {
    warning("no finite values in df", call. = FALSE)
    return(Inf)
  }
  min(vals)
}
