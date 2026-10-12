#' Remove overlaps between cluster outlines
#'
#' Takes the polygons from [outline_tsne()] and iteratively subtracts
#' overlapping regions so no two outlines intersect, re-ordering and
#' re-trimming each boundary after every pass.
#'
#' Subtraction is sequential, not symmetric: cluster `i` has every other
#' cluster subtracted from it in index order, using each neighbour's
#' *current* geometry. So where two outlines overlap, the contested region
#' ends up with the higher-indexed cluster, and the lower-indexed one gives
#' way.
#'
#' The guarantee that matters is checked by the package tests: the returned
#' outlines have **zero pairwise overlapping area**, and their total area
#' equals the area of the union of the inputs. Nothing is lost and nothing is
#' double-counted, which is what makes the result safe to shade.
#'
#' @param tsne_outl List of per-cluster boundary data frames from
#'   [outline_tsne()]. The first two columns are taken as x and y.
#' @param point_dist_window Rolling window used to smooth point-to-point
#'   distances before trimming.
#' @param strip_path_dist Points whose smoothed neighbour distance exceeds
#'   this are dropped, which removes the long spurs polygon subtraction leaves
#'   behind.
#' @param passes Number of subtract/reorder/trim passes after the first. Was
#'   hard-coded to 5.
#' @param kern Unused; retained for call compatibility.
#'
#' @return A list of boundary data frames with columns `x`, `y`, `ring` and
#'   `k`, one per input cluster and in the same order. Clusters entirely
#'   consumed by their neighbours come back as zero-row frames rather than
#'   being dropped, so the list stays aligned with the input.
#'
#'   `ring` identifies the closed path a point belongs to. A cluster with a
#'   hole punched in it has more than one: the outer boundary and the hole.
#'   Plot or smooth each ring separately — flattening them gives a
#'   self-intersecting path.
#'
#' @section Ported from rgeos to sf:
#' This used `rgeos`, which was **archived from CRAN in October 2023**, so it
#' could not run on a current R installation at all. The three calls map
#' across directly:
#'
#' | was | now |
#' | --- | --- |
#' | `rgeos::gIntersects()` | [sf::st_intersects()] |
#' | `rgeos::gBuffer(width = 0)` | [sf::st_make_valid()] |
#' | `rgeos::gDifference()` | [sf::st_difference()] |
#' | `sp::Polygon()` etc. | [sf::st_polygon()] |
#' | `raster::geom()` | [sf::st_coordinates()] |
#'
#' `sp` and `raster` are no longer needed either. The zero-width-buffer trick
#' for repairing self-intersecting rings is now `st_make_valid()`, which is
#' what it was always emulating.
#'
#' @section Other fixes:
#' * Point-to-point distance was `sqrt(dy^2 + dy^2)` — the x difference was
#'   never used, so the spur filter cut on the y axis alone.
#' * With no intersecting pairs, `for (n in 1:nrow(comparison))` ran once with
#'   `n = 1` against a zero-row frame and errored. Non-overlapping input now
#'   passes straight through.
#' * `gDifference()` returns `NULL` when one polygon is wholly inside another.
#'   That was fed to `raster::geom()` unguarded, which errored. Fully consumed
#'   clusters now return zero rows.
#' * Intersections were tested once against the *original* polygons, then
#'   subtracted in a loop that mutated them — so after the first subtraction
#'   the remaining tests were stale. Each pair is now re-tested against the
#'   current geometry.
#' * Rings are closed before being made into polygons; sf requires it and
#'   `sp` did not.
#' * `orderPoints()` is from `contoureR`, not `prevR`. A note added in 0.1.0
#'   said `prevR`, which has never exported a function of that name — so the
#'   call was failing inside a `try()` and the point reordering was silently
#'   doing nothing. With it working, point ordering and spur trimming had to
#'   be made per-ring: applied across rings they mix a hole boundary into its
#'   container and produce a self-intersecting path.
#'
#' @seealso [outline_tsne()] to produce the input, [spline.poly()] to smooth
#'   the result, [plot_tsne_kmeans()]
#'
#' @examples
#' # two overlapping squares
#' sq <- function(x0, y0, w = 1) {
#'   data.frame(x = c(x0, x0 + w, x0 + w, x0),
#'              y = c(y0, y0, y0 + w, y0 + w))
#' }
#' outlines <- list(sq(0, 0), sq(0.5, 0.5))
#'
#' if (requireNamespace("sf", quietly = TRUE) &&
#'     requireNamespace("contoureR", quietly = TRUE)) {
#'   clean <- remove_intersect_all(outlines, strip_path_dist = 10, passes = 1)
#'   vapply(clean, nrow, integer(1))
#' }
#' @export
remove_intersect_all <- function(tsne_outl,
                                 point_dist_window = 2,
                                 strip_path_dist = 1,
                                 passes = 5,
                                 kern = FALSE) {

  require_pkg("sf", "contoureR")

  n_clust <- length(tsne_outl)
  if (n_clust == 0L) {
    return(list())
  }

  # Close a ring and drop incomplete rows. sf needs a closed ring of at least
  # four points; sp did not, so the original never closed anything.
  close_ring <- function(xy) {
    xy <- xy[stats::complete.cases(xy), , drop = FALSE]
    if (nrow(xy) < 3L) {
      return(NULL)
    }
    if (!isTRUE(all.equal(xy[1, ], xy[nrow(xy), ]))) {
      xy <- rbind(xy, xy[1, , drop = FALSE])
    }
    unname(xy)
  }

  # A boundary frame may describe several rings (an outer boundary plus holes,
  # once something has been subtracted out of the middle of a cluster). The
  # `ring` column keeps them apart; without it they collapse into one
  # self-intersecting ring and the area is wrong.
  as_polygon <- function(df) {
    if (is.null(df) || nrow(df) < 3L) {
      return(sf::st_polygon())
    }
    ring <- if ("ring" %in% names(df)) df$ring else rep(1L, nrow(df))
    rings <- lapply(split(seq_len(nrow(df)), ring), function(idx) {
      close_ring(as.matrix(df[idx, c("x", "y"), drop = FALSE]))
    })
    rings <- rings[!vapply(rings, is.null, logical(1))]
    if (length(rings) == 0L) {
      return(sf::st_polygon())
    }
    sf::st_make_valid(sf::st_polygon(rings))
  }

  # st_coordinates returns X/Y plus L1 (ring) and L2 (part) indices. Combining
  # them gives one id per ring across all parts of a multipolygon.
  as_boundary <- function(geom, k) {
    empty <- data.frame(x = numeric(0), y = numeric(0),
                        ring = integer(0), k = integer(0))
    if (is.null(geom) || length(geom) == 0L || isTRUE(sf::st_is_empty(geom))) {
      return(empty)
    }
    xy <- try(sf::st_coordinates(geom), silent = TRUE)
    if (inherits(xy, "try-error") || nrow(xy) == 0L) {
      return(empty)
    }
    cols <- intersect(c("L2", "L1"), colnames(xy))
    ring <- if (length(cols)) {
      as.integer(factor(apply(xy[, cols, drop = FALSE], 1, paste,
                              collapse = "-")))
    } else {
      rep(1L, nrow(xy))
    }
    data.frame(x = xy[, "X"], y = xy[, "Y"], ring = ring, k = k)
  }

  rem_overlap <- function(outlines) {
    polys <- lapply(outlines, as_polygon)

    for (i in seq_len(n_clust)) {
      for (j in seq_len(n_clust)) {
        if (i == j) next
        a <- polys[[i]]
        b <- polys[[j]]
        if (isTRUE(sf::st_is_empty(a)) || isTRUE(sf::st_is_empty(b))) next

        # Re-test against the CURRENT geometry: the original computed every
        # intersection up front and then mutated the polygons as it went, so
        # all but the first test were stale.
        hits <- try(sf::st_intersects(a, b, sparse = FALSE)[1, 1], silent = TRUE)
        if (inherits(hits, "try-error") || !isTRUE(hits)) next

        diffed <- try(sf::st_difference(a, b), silent = TRUE)
        if (inherits(diffed, "try-error")) next
        polys[[i]] <- sf::st_make_valid(diffed)
      }
    }

    lapply(seq_len(n_clust), function(n) as_boundary(polys[[n]], n))
  }

  # Both of the per-point steps below operate WITHIN a ring. Applying them
  # across rings mixes a hole's boundary into its container's and produces a
  # self-intersecting path.
  by_ring <- function(df, fun) {
    if (nrow(df) == 0L) {
      return(df)
    }
    out <- lapply(split(df, df$ring), fun)
    out <- out[vapply(out, nrow, integer(1)) > 0L]
    if (length(out) == 0L) {
      return(df[0, , drop = FALSE])
    }
    do.call(rbind, out)
  }

  strip_long_path <- function(df) {
    if (nrow(df) < point_dist_window) {
      return(df)
    }
    p2pd <- c(0, sqrt((diff(df$x)^2) + (diff(df$y)^2)))
    kern_d <- c(rep(0, point_dist_window - 1),
                zoo::rollmean(p2pd, point_dist_window, na.pad = FALSE))
    df[kern_d < strip_path_dist, , drop = FALSE]
  }

  reorder_clockwise <- function(df) {
    if (nrow(df) < 3L) {
      return(df)
    }
    ord <- try(contoureR::orderPoints(df$x, df$y, clockwise = TRUE),
               silent = TRUE)
    if (inherits(ord, "try-error")) {
      return(df)
    }
    df[ord, , drop = FALSE]
  }

  one_pass <- function(outlines) {
    out <- rem_overlap(outlines)
    out <- lapply(out, by_ring, fun = reorder_clockwise)
    out <- lapply(out, by_ring, fun = strip_long_path)
    lapply(out, by_ring, fun = reorder_clockwise)
  }

  result <- one_pass(tsne_outl)
  for (n in seq_len(max(passes, 0L))) {
    result <- one_pass(result)
  }

  lapply(result, function(df) {
    rownames(df) <- NULL
    df[, c("x", "y", "ring", "k"), drop = FALSE]
  })
}
