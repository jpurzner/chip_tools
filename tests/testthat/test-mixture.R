# These tests pin the mixture-model fixes. Each one fails against the
# pre-restructure code; see REVIEW.md for the matching finding.

skip_if_no_mixtools <- function() {
  testthat::skip_if_not_installed("mixtools")
}

bimodal <- function(n = 600, seed = 1) {
  set.seed(seed)
  c(stats::rnorm(n * 0.6, 0.30, 0.08),
    stats::rnorm(n * 0.4, 0.70, 0.08))
}

test_that("mixture_fit_stats counts parameters and uses the real sample size", {
  # BIC used to be computed with n = length(data) -- `data` being the base R
  # function, so length() was 1 and log(n) was 0, making BIC == AIC always.
  model <- list(mu = c(0.3, 0.7), sigma = c(0.1, 0.1),
                lambda = c(0.6, 0.4), loglik = -100)
  fit <- chiptools:::mixture_fit_stats(model, n = 600)

  # 2 means + 1 shared variance + 1 free mixing proportion
  expect_equal(fit[["n_par"]], 4)
  expect_equal(fit[["aic"]], -2 * -100 + 2 * 4)
  expect_equal(fit[["bic"]], -2 * -100 + log(600) * 4)
  expect_gt(fit[["bic"]], fit[["aic"]])
})

test_that("binarize_counts returns calls for every input value", {
  skip_if_no_mixtools()
  x <- bimodal()
  calls <- binarize_counts(x, plot = FALSE)

  expect_s3_class(calls, "data.frame")
  expect_equal(nrow(calls), length(x))
  expect_setequal(names(calls), c("counts", "bin"))
  expect_true(all(calls$bin %in% c(0L, 1L)))
  # the two components should not both be called the same way
  expect_gt(sum(calls$bin), 0)
  expect_lt(sum(calls$bin), length(x))
})

test_that("binarize_counts works when the input is an expression, not a variable", {
  skip_if_no_mixtools()
  # The old code indexed a column literally named "counts", so passing
  # anything other than a variable called `counts` produced an error.
  expect_no_error(binarize_counts(bimodal() * 1, plot = FALSE))
})

test_that("binarize_counts applies exclude_extreme to the transformed values", {
  skip_if_no_mixtools()
  # With log_data and exclude_extreme both set, the old code overwrote the
  # log-transformed x with a subset of the RAW counts, silently discarding the
  # transform. On raw counts in the thousands the 0.1 < x < 1 filter leaves
  # nothing, so this now errors clearly instead of fitting the wrong scale.
  raw <- bimodal() * 1000
  expect_error(
    binarize_counts(raw, log_data = TRUE, exclude_extreme = TRUE, plot = FALSE),
    "fewer than 10 usable values"
  )
})

test_that("binarize_counts summary has the names the cutoff tables expect", {
  skip_if_no_mixtools()
  s <- binarize_counts(bimodal(), return_what = "summary", plot = FALSE)
  expect_setequal(names(s), c("cutoff_50", "cutoff_75", "mean_upper",
                              "mean_lower", "loglik", "aic", "bic"))
  expect_lte(s[["mean_lower"]], s[["mean_upper"]])
})

test_that("split_counts compares k and returns the BIC winner", {
  skip_if_no_mixtools()
  # The ~100 lines after the old `return(best_model)` were unreachable, and
  # BIC was penalised by k rather than the parameter count.
  cmp <- split_counts(bimodal(), min_k = 2, max_k = 3,
                      plot = FALSE, return_what = "comparison")
  expect_true(all(c("k", "loglik", "n_par", "aic", "bic") %in% names(cmp)))
  expect_equal(cmp$k, 2:3)
  expect_true(all(cmp$n_par > cmp$k))   # more parameters than components

  model <- split_counts(bimodal(), min_k = 2, max_k = 3, plot = FALSE)
  expect_true(!is.null(model$mu))
})

test_that("trinarize_counts gives three states and ordered cutoffs", {
  skip_if_no_mixtools()
  # The old version hit `return(model)` before any of its documented return
  # options, so it could only ever hand back the raw mixEM object.
  set.seed(42)
  x <- c(stats::rnorm(300, 0.2, 0.05),
         stats::rnorm(300, 0.5, 0.05),
         stats::rnorm(300, 0.8, 0.05))

  cutoffs <- trinarize_counts(x, plot = FALSE, return_what = "cutoffs")
  expect_setequal(names(cutoffs), c("low", "high"))
  expect_lt(cutoffs[["low"]], cutoffs[["high"]])

  calls <- trinarize_counts(x, plot = FALSE)
  expect_equal(nrow(calls), length(x))
  expect_true(all(calls$state %in% 0:2))
  expect_equal(length(unique(calls$state)), 3)
})

test_that("normal_mixture_cutoff brackets the root between the component means", {
  model <- list(mu = c(0.3, 0.7), sigma = c(0.08, 0.08), lambda = c(0.5, 0.5))
  cutoff <- chiptools:::normal_mixture_cutoff(model, proba = 0.5, i = 1)
  expect_gt(cutoff, 0.3)
  expect_lt(cutoff, 0.7)
  # equal weights and equal variances put the 0.5 crossing at the midpoint
  expect_equal(cutoff, 0.5, tolerance = 1e-6)
})

test_that("require_pkg names the missing package and the caller", {
  expect_error(
    chiptools:::require_pkg("definitelyNotAnInstalledPackage"),
    "definitelyNotAnInstalledPackage"
  )
  expect_silent(chiptools:::require_pkg("stats"))
})
