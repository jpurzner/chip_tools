#' Split counts into three states with a 3-component Gaussian mixture
#'
#' Fits a three-component normal mixture to a vector of counts and derives the
#' two cutoffs that separate the components, giving a low / mid / high call for
#' every observation. This is the three-state counterpart to
#' [binarize_counts()].
#'
#' Cutoffs are placed where the posterior probability of belonging to the lower
#' of an adjacent pair of components equals `proba` (0.5 by default), solved
#' with [stats::uniroot()]. If the root solver fails the midpoint between the
#' two component means is used instead and a warning is signalled.
#'
#' @param counts Numeric vector of counts (or any continuous signal).
#' @param log_data If `TRUE`, fit on `log10(counts + 1)` and back-transform the
#'   cutoffs onto the original scale before classifying.
#' @param exclude_extreme If `TRUE`, fit on the subset `0.1 < x < 1` only. The
#'   returned classification still covers all of `counts`.
#' @param mu,lambda,sigma Starting values for the three component means,
#'   mixing proportions and standard deviations passed to
#'   [mixtools::normalmixEM()].
#' @param proba Posterior probability defining each cutoff.
#' @param label_name Optional plot title.
#' @param plot If `TRUE` (default) draw the fit.
#' @param return_what One of `"calls"` (default) for a data frame of counts and
#'   their 0/1/2 state, `"cutoffs"` for the two cutoffs, `"means"` for the
#'   component means, or `"all"` for a list holding the calls, cutoffs, the
#'   fitted model, fit statistics and the plot.
#'
#' @return Depends on `return_what`; see above.
#'
#' @seealso [binarize_counts()], [split_counts()]
#'
#' @examples
#' \dontrun{
#' states <- trinarize_counts(rlog_values, log_data = TRUE)
#' table(states$state)
#' }
#' @import ggplot2
#' @export
trinarize_counts <- function(counts,
                             log_data = FALSE,
                             exclude_extreme = FALSE,
                             mu = c(0.3, 0.5, 0.7),
                             lambda = c(1/3, 1/3, 1/3),
                             sigma = c(0.3, 0.2, 0.25),
                             proba = 0.5,
                             label_name = NULL,
                             plot = TRUE,
                             return_what = c("calls", "cutoffs", "means", "all")) {

  require_pkg("mixtools")
  return_what <- match.arg(return_what)

  x_all <- if (log_data) log10(counts + 1) else counts

  # Subset used for fitting only; `counts` is still what gets classified.
  x_fit <- x_all[is.finite(x_all)]
  if (exclude_extreme) {
    x_fit <- x_fit[x_fit > 0.1 & x_fit < 1]
  }
  if (length(x_fit) < 10) {
    stop("fewer than 10 usable values to fit a 3-component mixture", call. = FALSE)
  }

  model <- mixtools::normalmixEM(x = x_fit, k = 3, mu = mu, lambda = lambda,
                                 sigma = sigma, arbvar = FALSE, maxit = 10000,
                                 maxrestarts = 30, epsilon = 1e-20)

  ord <- order(model$mu)

  # Cutoff between an adjacent pair of components (lo, hi are indices into
  # the unordered model parameters), at the point where the posterior for the
  # lower component equals `proba`.
  pair_cutoff <- function(lo, hi) {
    f <- function(v) {
      dens <- model$lambda * stats::dnorm(v, mean = model$mu, sd = model$sigma)
      proba - (dens[lo] / (dens[lo] + dens[hi]))
    }
    lower <- model$mu[lo]
    upper <- model$mu[hi]
    root <- tryCatch(stats::uniroot(f, lower = lower, upper = upper)$root,
                     error = function(e) NULL)
    if (is.null(root)) {
      warning("uniroot failed between components ", lo, " and ", hi,
              "; using the midpoint of their means", call. = FALSE)
      return(mean(c(model$mu[lo], model$mu[hi])))
    }
    root
  }

  cutoffs <- c(low = pair_cutoff(ord[1], ord[2]),
               high = pair_cutoff(ord[2], ord[3]))

  if (plot) {
    plot_df <- data.frame(x = x_fit)
    for (i in seq_len(3)) {
      plot_df[[paste0("nd", i)]] <-
        model$lambda[ord[i]] * stats::dnorm(plot_df$x,
                                            model$mu[ord[i]],
                                            model$sigma[ord[i]])
    }
    p <- mixture_density_plot(
      plot_df,
      component_cols = c("nd1", "nd2", "nd3"),
      component_colours = c("dodgerblue1", "seagreen4", "red"),
      cutoffs = cutoffs,
      x_label = if (log_data) "log10(counts + 1)" else "counts",
      title = label_name
    )
    print(p)
  } else {
    p <- NULL
  }

  # Back-transform the cutoffs so classification happens on the input scale.
  cutoffs_out <- if (log_data) (10^cutoffs) - 1 else cutoffs

  state <- cut(counts,
               breaks = c(-Inf, cutoffs_out[["low"]], cutoffs_out[["high"]], Inf),
               labels = c(0L, 1L, 2L), right = FALSE)
  calls <- data.frame(counts = counts, state = as.integer(as.character(state)))

  switch(return_what,
         calls   = calls,
         cutoffs = cutoffs_out,
         means   = stats::setNames(model$mu[ord], c("mean_low", "mean_mid", "mean_high")),
         all     = list(calls   = calls,
                        cutoffs = cutoffs_out,
                        model   = model,
                        fit     = mixture_fit_stats(model, length(x_fit)),
                        plot    = p))
}
