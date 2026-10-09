#' Internal helpers
#'
#' Shared machinery for the mixture-model and plotting functions. Not exported.
#'
#' @name chiptools-internal
#' @keywords internal
NULL

#' Require an optional package
#'
#' Most of the heavy dependencies here (mixtools, mclust, DESeq2, topGO, ...)
#' are only needed by a handful of functions, so they live in `Suggests` rather
#' than `Imports`. This gives a single actionable error when one is missing,
#' instead of the `could not find function` error you would otherwise get
#' several frames deep.
#'
#' @param ... Package names.
#' @param reason Optional explanation appended to the error.
#' @noRd
require_pkg <- function(..., reason = NULL) {
  pkgs <- c(...)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0L) {
    return(invisible(TRUE))
  }
  caller <- deparse(sys.call(-1)[[1]])
  stop(sprintf(
    "%s needs %s, which %s not installed.%s\n  install.packages(c(%s))",
    caller,
    paste(sQuote(missing), collapse = ", "),
    if (length(missing) == 1L) "is" else "are",
    if (is.null(reason)) "" else paste0(" ", reason),
    paste(sprintf('"%s"', missing), collapse = ", ")
  ), call. = FALSE)
}

#' Standard mixture-fit diagnostics
#'
#' AIC and BIC for a fitted `mixtools` normal mixture.
#'
#' A 2-component univariate normal mixture with `arbvar = FALSE` (one shared
#' variance) has `k` means + 1 variance + `k - 1` free mixing proportions
#' parameters. The original scripts computed this as
#' `length(lambda) + 2 * length(mu)`, double counted the variance, and passed
#' `length(data)` — the length of the base R `data()` function, i.e. 1 — as the
#' sample size, so every BIC was wrong. Both are fixed here.
#'
#' @param model A `mixEM` object from [mixtools::normalmixEM()].
#' @param n Number of observations the model was fitted to.
#' @return Named numeric vector of `loglik`, `n_par`, `aic`, `bic`.
#' @noRd
mixture_fit_stats <- function(model, n) {
  k <- length(model$mu)
  n_sigma <- length(unique(model$sigma))
  n_par <- k + n_sigma + (k - 1L)
  c(loglik = model$loglik,
    n_par  = n_par,
    aic    = -2 * model$loglik + 2 * n_par,
    bic    = -2 * model$loglik + log(n) * n_par)
}

#' House style for the mixture histogram + component density plot
#'
#' `binarize_counts()`, `split_counts()` and `trinarize_counts()` all drew the
#' same figure with slightly different code. This is that figure, once.
#'
#' @param df Data frame with column `x` plus one column per component density.
#' @param component_cols Density column names, in plotting order.
#' @param component_colours Line colours, recycled to `component_cols`.
#' @param cutoffs Numeric vector of vertical-line positions.
#' @param x_label,title Axis label and optional title.
#' @param bins Histogram bin count.
#' @return A `ggplot` object.
#' @import ggplot2
#' @noRd
mixture_density_plot <- function(df,
                                 component_cols,
                                 component_colours = c("dodgerblue1", "seagreen4", "red"),
                                 cutoffs = numeric(0),
                                 x_label = "counts",
                                 title = NULL,
                                 bins = 100) {
  component_colours <- rep_len(component_colours, length(component_cols))

  p <- ggplot(df, aes(x = .data$x)) +
    geom_histogram(aes(y = after_stat(density)), bins = bins,
                   colour = "black", fill = "grey", alpha = 0.3)

  for (i in seq_along(component_cols)) {
    p <- p + geom_line(aes(y = .data[[component_cols[i]]]),
                       colour = component_colours[i],
                       linewidth = 1, alpha = 0.8)
  }

  if (length(cutoffs)) {
    p <- p + geom_vline(xintercept = cutoffs, linetype = "dotdash", alpha = 0.5)
  }

  p <- p +
    labs(x = x_label, y = "density", title = title) +
    theme(axis.title.y = element_text(size = 16, angle = 90),
          axis.text.y  = element_text(size = 16),
          axis.text.x  = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border     = element_blank(),
          panel.background = element_blank())
  p
}

#' Posterior-probability cutoff between two mixture components
#'
#' Solves for the value of `x` at which the posterior probability of belonging
#' to component `i` equals `proba`, given a two-component normal mixture.
#'
#' @param model A `mixEM` object.
#' @param proba Target posterior probability.
#' @param i Index of the component the probability refers to.
#' @return A single numeric cutoff.
#' @noRd
normal_mixture_cutoff <- function(model, proba = 0.5, i = which.min(model$mu)) {
  f <- function(v) {
    dens <- model$lambda * stats::dnorm(v, mean = model$mu, sd = model$sigma)
    proba - (dens[i] / sum(dens))
  }

  # Bracket the root between the two component means rather than between an
  # arbitrary data point and max(x); the original code took `low` from
  # `x[which.min(f(x))]`, which is not guaranteed to bracket a sign change and
  # was the usual cause of the "failed to find root" fallback firing.
  lower <- min(model$mu)
  upper <- max(model$mu)

  root <- tryCatch(stats::uniroot(f, lower = lower, upper = upper)$root,
                   error = function(e) NULL)
  if (is.null(root)) {
    warning("uniroot failed to find a cutoff; using the midpoint of the ",
            "component means", call. = FALSE)
    return(mean(range(model$mu)))
  }
  root
}

#' Null-coalescing operator
#'
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x
