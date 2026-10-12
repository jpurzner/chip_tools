# remove_intersect_all() was unrunnable before 0.4.0 -- it needed rgeos,
# archived from CRAN in October 2023. These tests pin the contract of the sf
# port: zero pairwise overlapping area, and total area equal to the union of
# the inputs.

skip_unless_spatial <- function() {
  skip_if_not_installed("sf")
  skip_if_not_installed("contoureR")
}

sq <- function(x0, y0, w = 1) {
  data.frame(x = c(x0, x0 + w, x0 + w, x0),
             y = c(y0, y0, y0 + w, y0 + w))
}

# Turn a returned boundary frame back into a polygon so we can measure it.
# Must respect the `ring` column: a cluster with a hole has two rings, and
# flattening them into one gives a self-intersecting path and the wrong area.
as_poly <- function(df) {
  if (nrow(df) < 3L) return(sf::st_polygon())
  rings <- lapply(split(seq_len(nrow(df)), df$ring), function(idx) {
    m <- as.matrix(df[idx, c("x", "y"), drop = FALSE])
    if (nrow(m) < 3L) return(NULL)
    if (!isTRUE(all.equal(m[1, ], m[nrow(m), ]))) m <- rbind(m, m[1, ])
    unname(m)
  })
  rings <- rings[!vapply(rings, is.null, logical(1))]
  if (!length(rings)) return(sf::st_polygon())
  sf::st_make_valid(sf::st_polygon(rings))
}

poly_area <- function(df) {
  p <- as_poly(df)
  if (sf::st_is_empty(p)) 0 else as.numeric(sf::st_area(p))
}

total_overlap_area <- function(result) {
  ps <- lapply(result, as_poly)
  total <- 0
  for (i in seq_along(ps)) {
    for (j in seq_along(ps)) {
      if (i >= j) next
      if (sf::st_is_empty(ps[[i]]) || sf::st_is_empty(ps[[j]])) next
      inter <- suppressWarnings(sf::st_intersection(ps[[i]], ps[[j]]))
      if (length(inter) && !sf::st_is_empty(inter)) {
        total <- total + as.numeric(sf::st_area(inter))
      }
    }
  }
  total
}

test_that("overlapping outlines come back with zero overlapping area", {
  skip_unless_spatial()
  res <- remove_intersect_all(list(sq(0, 0), sq(0.5, 0.5)),
                              strip_path_dist = 10, passes = 1)

  expect_length(res, 2L)
  expect_lt(total_overlap_area(res), 1e-9)

  # total area is the union of two unit squares overlapping on a 0.5 x 0.5 patch
  expect_equal(sum(vapply(res, poly_area, numeric(1))), 1.75, tolerance = 1e-6)
})

test_that("area is conserved across a chain of three outlines", {
  skip_unless_spatial()
  res <- remove_intersect_all(list(sq(0, 0), sq(0.5, 0), sq(1, 0)),
                              strip_path_dist = 10, passes = 1)

  expect_lt(total_overlap_area(res), 1e-9)
  # union of unit squares at x = 0, 0.5, 1 spans x in [0, 2], y in [0, 1]
  expect_equal(sum(vapply(res, poly_area, numeric(1))), 2, tolerance = 1e-6)
})

test_that("non-overlapping input passes straight through", {
  skip_unless_spatial()
  # `for (n in 1:nrow(comparison))` used to run once with n = 1 against a
  # zero-row frame and error out here.
  res <- remove_intersect_all(list(sq(0, 0), sq(5, 5)),
                              strip_path_dist = 10, passes = 1)
  expect_length(res, 2L)
  expect_equal(vapply(res, poly_area, numeric(1)), c(1, 1), tolerance = 1e-6)
  expect_lt(total_overlap_area(res), 1e-9)
})

test_that("a nested outline does not error and conserves area", {
  skip_unless_spatial()
  # gDifference() returned NULL for a wholly-contained polygon, and that NULL
  # went unguarded into raster::geom(). This is also the case that needs
  # per-ring handling: the outer cluster ends up as a ring with a hole.
  res <- remove_intersect_all(list(sq(0, 0, 4), sq(1, 1, 1)),
                              strip_path_dist = 20, passes = 1)
  expect_length(res, 2L)
  expect_lt(total_overlap_area(res), 1e-9)
  expect_equal(sum(vapply(res, poly_area, numeric(1))), 16, tolerance = 1e-6)
})

test_that("degenerate and empty input are handled", {
  skip_unless_spatial()
  expect_length(remove_intersect_all(list()), 0L)

  # fewer than three points is not a polygon
  res <- remove_intersect_all(list(data.frame(x = c(0, 1), y = c(0, 1)), sq(0, 0)),
                              strip_path_dist = 10, passes = 1)
  expect_equal(nrow(res[[1]]), 0L)
  expect_gt(nrow(res[[2]]), 0L)
})

test_that("the result keeps input order and carries a cluster index", {
  skip_unless_spatial()
  res <- remove_intersect_all(list(sq(0, 0), sq(0.5, 0.5), sq(9, 9)),
                              strip_path_dist = 10, passes = 1)
  expect_length(res, 3L)
  for (n in seq_along(res)) {
    expect_setequal(names(res[[n]]), c("x", "y", "ring", "k"))
    if (nrow(res[[n]])) expect_equal(unique(res[[n]]$k), n)
  }
})

test_that("passes = 0 still does one subtraction pass", {
  skip_unless_spatial()
  res <- remove_intersect_all(list(sq(0, 0), sq(0.5, 0.5)),
                              strip_path_dist = 10, passes = 0)
  expect_lt(total_overlap_area(res), 1e-9)
})

test_that("remove_intersect_all no longer needs rgeos", {
  skip_unless_spatial()
  # The whole point of the port: this used to be impossible on current R.
  expect_false("rgeos" %in% loadedNamespaces())
  expect_no_error(
    remove_intersect_all(list(sq(0, 0), sq(0.5, 0.5)),
                         strip_path_dist = 10, passes = 1)
  )
})

test_that("spline.poly smooths an outline returned by remove_intersect_all", {
  skip_unless_spatial()
  res <- remove_intersect_all(list(sq(0, 0), sq(0.5, 0.5)),
                              strip_path_dist = 10, passes = 1)
  keep <- res[[which.max(vapply(res, nrow, integer(1)))]]
  keep <- keep[keep$ring == keep$ring[1], ]   # outer ring only
  sm <- spline.poly(as.matrix(keep[, c("x", "y")]), vertices = 50)
  expect_equal(ncol(sm), 2L)
  expect_true(all(is.finite(sm)))
})

test_that("rings are kept separate when a cluster ends up with a hole", {
  skip_unless_spatial()
  # Subtracting a small square from the middle of a big one leaves a ring.
  # Ordering or trimming points ACROSS the two rings mixes the hole boundary
  # into the outer boundary and produces a self-intersecting path -- which is
  # exactly what happened once orderPoints() started working.
  res <- remove_intersect_all(list(sq(0, 0, 4), sq(1, 1, 1)),
                              strip_path_dist = 20, passes = 1)
  outer <- res[[1]]

  expect_gt(length(unique(outer$ring)), 1L)
  # the outer ring spans the full 4x4 extent, the hole only 1x1
  spans <- tapply(outer$x, outer$ring, function(v) diff(range(v)))
  expect_equal(sort(as.vector(round(spans, 6))), c(1, 4))

  # and the reconstructed polygon is valid, i.e. not self-intersecting
  expect_true(sf::st_is_valid(as_poly(outer)))
})
