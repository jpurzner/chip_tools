test_that("spline.poly smooths a closed polygon", {
  square <- cbind(c(0, 1, 1, 0), c(0, 0, 1, 1))
  out <- spline.poly(square, vertices = 60)
  expect_true(is.matrix(out))
  expect_equal(ncol(out), 2)
  expect_gt(nrow(out), 10)
  expect_true(all(is.finite(out)))
})

test_that("remove_intersect_all reports its missing dependency clearly", {
  # rgeos was archived from CRAN in October 2023. Previously this failed deep
  # inside the function with "could not find function gIntersects".
  skip_if(requireNamespace("rgeos", quietly = TRUE),
          "rgeos is installed, so the guard does not fire")
  expect_error(remove_intersect_all(list()), "rgeos")
  expect_error(remove_intersect_all(list()), "archived from CRAN")
})
