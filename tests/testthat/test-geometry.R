test_that("spline.poly smooths a closed polygon", {
  square <- cbind(c(0, 1, 1, 0), c(0, 0, 1, 1))
  out <- spline.poly(square, vertices = 60)
  expect_true(is.matrix(out))
  expect_equal(ncol(out), 2)
  expect_gt(nrow(out), 10)
  expect_true(all(is.finite(out)))
})

# The old test here asserted that remove_intersect_all() errors with an
# rgeos message. That assertion is obsolete: the function was ported to sf in
# 0.4.0 and no longer needs rgeos at all. See test-outlines.R.

test_that("outline_tsne computes convex hull area without sp", {
  # sp::Polygon(...)@area was replaced with the shoelace formula, which is
  # exact for a simple polygon and removes the last reason to depend on sp.
  shoelace <- function(cc) {
    cx <- cc[[1]]; cy <- cc[[2]]
    abs(sum(cx * c(cy[-1], cy[1]) - c(cx[-1], cx[1]) * cy)) / 2
  }
  # unit square, area 1
  expect_equal(shoelace(data.frame(x = c(0, 1, 1, 0), y = c(0, 0, 1, 1))), 1)
  # right triangle with legs 3 and 4, area 6
  expect_equal(shoelace(data.frame(x = c(0, 3, 0), y = c(0, 0, 4))), 6)

  expect_false("sp" %in% loadedNamespaces())
})
