#' Smooth a closed polygon with a periodic spline
#'
#' Splines the x and y coordinates of a polygon boundary separately, wrapping
#' `k` vertices around each end so the result closes smoothly. Used to tidy
#' the jagged cluster outlines that come out of [outline_tsne()].
#'
#' @param xy An n x 2 matrix of boundary coordinates in order, with `n >= k`.
#' @param vertices Number of spline points to generate (some are clipped from
#'   the ends, so the result has fewer rows than this).
#' @param k Number of vertices wrapped around each end.
#' @param ... Passed to [stats::spline()].
#'
#' @return A two-column matrix of smoothed boundary coordinates.
#'
#' @source Adapted from whuber's answer at
#'   <https://gis.stackexchange.com/questions/24827>.
#'
#' @seealso [outline_tsne()], [remove_intersect_all()]
#'
#' @examples
#' square <- cbind(c(0, 1, 1, 0), c(0, 0, 1, 1))
#' head(spline.poly(square, vertices = 50))
#' @export

spline.poly <- function(xy, vertices, k=3, ...) {
  #
  # Splining a polygon.
  #
  #   The rows of 'xy' give coordinates of the boundary vertices, in order.
  #   'vertices' is the number of spline vertices to create.
  #              (Not all are used: some are clipped from the ends.)
  #   'k' is the number of points to wrap around the ends to obtain
  #       a smooth periodic spline.
  #
  #   Returns an array of points. 
  # 
  
  # Assert: xy is an n by 2 matrix with n >= k.
  
  # Wrap k vertices around each end.
  n <- dim(xy)[1]
  if (k >= 1) {
    data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
  } else {
    data <- xy
  }
  
  # Spline the x and y coordinates.
  data.spline <- stats::spline(1:(n+2*k), data[,1], n=vertices, ...)
  x <- data.spline$x
  x1 <- data.spline$y
  x2 <- stats::spline(1:(n+2*k), data[,2], n=vertices, ...)$y
  
  # Retain only the middle part.
  cbind(x1, x2)[k < x & x <= n+k, ]
}  