#' Thin out the zero spike in a count distribution
#'
#' Count data over features often has a zero bin tall enough to swamp every
#' other bin, which both hides the shape of the signal and dominates a mixture
#' fit. This replaces a random subset of the zeros with `NA` so the zero bin is
#' at most `scale_max` times the tallest non-zero bin.
#'
#' @param counts Numeric vector of counts.
#' @param scale_max Target height of the zero bin as a multiple of the tallest
#'   non-zero bin.
#' @param breaks Passed to [graphics::hist()].
#' @param plot If `TRUE` draw the histogram used to size the bins. The original
#'   always drew it, with no way to turn it off.
#'
#' @return `counts` with some zeros replaced by `NA`. If the zero bin is
#'   already below the target, `counts` is returned unchanged.
#'
#' @examples
#' set.seed(1)
#' x <- c(rep(0, 500), rpois(200, 5))
#' sum(is.na(soften_zero(x, plot = FALSE)))
#' @export
soften_zero <- function(counts, scale_max = 1.5, breaks = 100, plot = FALSE) {

  histogram <- graphics::hist(counts, plot = plot, breaks = breaks)

  if (length(histogram$counts) < 2L) {
    warning("histogram has fewer than two bins; returning counts unchanged",
            call. = FALSE)
    return(counts)
  }

  max_value <- max(histogram$counts[-1])
  zero_indices <- which(counts == 0)
  n_remove <- length(zero_indices) - (max_value * scale_max)

  # `sample(x, n)` errors on a negative n, so the original blew up whenever
  # the zero bin was already shorter than scale_max * the tallest other bin.
  if (n_remove <= 0) {
    return(counts)
  }

  counts[sample(zero_indices, floor(n_remove))] <- NA
  counts
}
