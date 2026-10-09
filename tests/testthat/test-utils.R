# Pins the fixes to the data-wrangling and normalisation helpers.

test_that("grouped_col_mean keeps single-replicate groups", {
  # The old version split counts into single- and multi-replicate halves, then
  # indexed the multi-replicate half with positions from the UNSPLIT grouping
  # vector, and dropped the single-replicate groups from the output entirely.
  counts <- data.frame(a1 = c(1, 2), a2 = c(3, 4), b1 = c(10, 20))
  meta <- data.frame(condition = c("a", "a", "b"))

  out <- grouped_col_mean(counts, meta, "condition")

  expect_equal(colnames(out), c("a", "b"))
  expect_equal(out$a, c(2, 3))      # mean of a1 and a2
  expect_equal(out$b, c(10, 20))    # single replicate carried through
})

test_that("grouped_col_mean rejects mismatched metadata", {
  counts <- data.frame(a1 = 1:2, a2 = 3:4)
  expect_error(grouped_col_mean(counts, data.frame(condition = "a"), "condition"),
               "correspond one to one")
  expect_error(grouped_col_mean(counts, data.frame(condition = c("a", "a")), "nope"),
               "not a column of metadata")
})

test_that("df_rowfilt returns the aligned frame, not its input", {
  # The old version built `df_new` and then returned `df_filt`, so it was a
  # no-op that printed three dimensions.
  filt <- data.frame(score = c(5, 6), row.names = c("g1", "g2"))
  main <- data.frame(expr = c(1, 2, 3), row.names = c("g1", "g2", "g3"))

  out <- df_rowfilt(filt, main)

  expect_equal(nrow(out), 3)
  expect_equal(colnames(out), "score")
  expect_true(is.na(out["g3", "score"]))
  expect_equal(out["g1", "score"], 5)
})

test_that("df_rowmatch zero-fills by default", {
  filt <- data.frame(score = 5, row.names = "g1")
  main <- data.frame(expr = c(1, 2), row.names = c("g1", "g2"))

  expect_equal(df_rowmatch(filt, main)["g2", "score"], 0)
  expect_true(is.na(df_rowmatch(filt, main, na2zero = FALSE)["g2", "score"]))
})

test_that("min_finite ignores infinities and does not print", {
  expect_silent(out <- min_finite(data.frame(a = c(1, -Inf), b = c(3, 2))))
  expect_equal(out, 1)
  expect_warning(expect_equal(min_finite(data.frame(a = Inf)), Inf),
                 "no finite values")
})

test_that("soften_zero returns input unchanged rather than erroring", {
  # sample(x, n) errors on a negative n, so the old version blew up whenever
  # the zero bin was already shorter than scale_max * the tallest other bin.
  set.seed(1)
  few_zeros <- c(rep(0, 3), stats::rpois(500, 20))
  expect_no_error(out <- soften_zero(few_zeros, plot = FALSE))
  expect_equal(out, few_zeros)
})

test_that("soften_zero thins a dominant zero spike", {
  set.seed(1)
  x <- c(rep(0, 2000), stats::rpois(300, 20))
  out <- soften_zero(x, scale_max = 1.5, plot = FALSE)
  expect_gt(sum(is.na(out)), 0)
  expect_lt(sum(out == 0, na.rm = TRUE), 2000)
  # only zeros are ever removed
  expect_true(all(x[is.na(out)] == 0))
})

test_that("wins_norm_histones scales each column to 0-1", {
  skip_if_not_installed("DescTools")
  skip_if_not_installed("pheatmap")
  set.seed(1)
  m <- data.frame(s1 = stats::rnorm(200, 0), s2 = stats::rnorm(200, 50))

  out <- wins_norm_histones(m, low = 5, high = 5)

  expect_equal(dim(out), dim(m))
  expect_true(all(unlist(out) >= 0 & unlist(out) <= 1))
  expect_equal(unname(apply(out, 2, min)), c(0, 0))
  expect_equal(unname(apply(out, 2, max)), c(1, 1))
})

test_that("median_subtract reaches the winsorised values but cancels under 0-1 scaling", {
  skip_if_not_installed("DescTools")
  skip_if_not_installed("pheatmap")
  # The argument used to be ignored outright -- the median-subtracted matrix
  # was computed and then the ORIGINAL matrix was winsorised. It is now wired
  # up, but min-max rescaling is shift-invariant, so it can only be observed
  # with normalize = FALSE. See ?wins_norm_histones.
  set.seed(1)
  m <- data.frame(s1 = stats::rnorm(200, 0), s2 = stats::rnorm(200, 50))

  raw_plain <- wins_norm_histones(m, median_subtract = FALSE, low = 5, high = 5,
                                  normalize = FALSE)
  raw_subtr <- wins_norm_histones(m, median_subtract = TRUE, low = 5, high = 5,
                                  normalize = FALSE)
  expect_false(isTRUE(all.equal(raw_plain, raw_subtr)))
  expect_lt(abs(stats::median(raw_subtr$s2)), abs(stats::median(raw_plain$s2)))

  # ... and cancels once rescaled
  expect_equal(
    wins_norm_histones(m, median_subtract = FALSE, low = 5, high = 5),
    wins_norm_histones(m, median_subtract = TRUE,  low = 5, high = 5)
  )
})

test_that("opti_map pairs the labellings up", {
  # x group 1 overlaps y group 2, x 2 -> y 3, x 3 -> y 1
  out <- opti_map(c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                  c(2, 2, 2, 3, 3, 3, 1, 1, 1))
  expect_equal(nrow(out), 3)
  expect_setequal(out$x, 1:3)
  expect_setequal(out$y, 1:3)
})

test_that("counts2tpm columns sum to a million", {
  counts <- data.frame(s1 = c(100, 200, 300), s2 = c(50, 400, 150),
                       row.names = paste0("g", 1:3))
  lens <- data.frame(len = c(1000, 2000, 500))
  tpm <- counts2tpm(counts, lens, read_length = 100)
  expect_equal(unname(colSums(tpm)), c(1e6, 1e6), tolerance = 1e-6)
  expect_equal(rownames(tpm), rownames(counts))
})

test_that("chip_mclust_icl_all survives a column with nothing above the cutoff", {
  skip_if_not_installed("mclust")
  # `x[-integer(0)]` returns an EMPTY vector, so any column whose values all
  # fell below the tail cutoff was clustered on nothing and came back all NA.
  set.seed(1)
  dat <- data.frame(sample_a = c(stats::rnorm(150, 1, 0.3),
                                 stats::rnorm(150, 4, 0.3)))
  cutoffs <- data.frame(Column = "sample_a", CrossingSmoothedValue = 1e6)

  out <- chip_mclust_icl_all(dat, cutoffs, G_range = 2:3)

  expect_true("sample_a_class" %in% names(out))
  expect_false(all(is.na(out$sample_a_class)))
  expect_equal(nrow(out), nrow(dat))
})
