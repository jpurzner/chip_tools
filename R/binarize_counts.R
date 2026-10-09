#' Separate signal from background with a 2-component Gaussian mixture
#'
#' Fits a two-component normal mixture to a vector of counts and uses the
#' posterior-probability crossing point as the threshold between background and
#' signal. Written for ChIP-seq count data over gene bodies or promoters, where
#' the distribution is reliably bimodal.
#'
#' The approach follows the answer at
#' <https://stats.stackexchange.com/questions/57993>, adapted for
#' next-generation sequencing counts.
#'
#' @param counts Numeric vector of counts (or any continuous signal).
#' @param log_data If `TRUE`, fit on `log10(counts + 1)` and back-transform the
#'   cutoffs onto the original scale before classifying.
#' @param exclude_extreme If `TRUE`, restrict the *fit* to `0.1 < x < 1`. The
#'   returned classification still covers every element of `counts`.
#' @param mu Starting values for the two component means.
#' @param label_name Optional plot title.
#' @param plot If `TRUE` (default) draw the fitted mixture.
#' @param return_what One of `"calls"` (default) for a data frame of counts and
#'   their 0/1 call, `"cutoffs"` for the two cutoffs, `"mean_upper"` for the
#'   mean of the signal component, `"summary"` for the named vector the
#'   downstream cutoff tables expect, or `"all"` for everything plus the model
#'   and plot.
#'
#' @return Depends on `return_what`. `"summary"` gives a named numeric vector
#'   of `cutoff_50`, `cutoff_75`, `mean_upper`, `mean_lower`, `loglik`, `aic`
#'   and `bic` — this is the form [chip_binarize_tail_removed()] consumes.
#'
#' @seealso [trinarize_counts()] for three states, [split_counts()] to let BIC
#'   choose the component count.
#'
#' @examples
#' \dontrun{
#' calls <- binarize_counts(counts_vector, log_data = TRUE)
#' table(calls$bin)
#'
#' # the row that goes into a cutoff table
#' binarize_counts(counts_vector, log_data = TRUE, return_what = "summary")
#' }
#' @import ggplot2
#' @export
binarize_counts <- function(counts,
                            log_data = FALSE,
                            exclude_extreme = FALSE,
                            mu = c(0.3, 0.65),
                            label_name = NULL,
                            plot = TRUE,
                            return_what = c("calls", "cutoffs", "mean_upper",
                                            "summary", "all")) {

  require_pkg("mixtools")
  return_what <- match.arg(return_what)

  x_all <- if (log_data) log10(counts + 1) else counts

  # Fit on the (optionally trimmed) transformed values. The original code
  # subset `counts` here rather than `x`, which silently threw away the log
  # transform whenever both `log_data` and `exclude_extreme` were set.
  x_fit <- x_all[is.finite(x_all)]
  if (exclude_extreme) {
    x_fit <- x_fit[x_fit > 0.1 & x_fit < 1]
  }
  if (length(x_fit) < 10) {
    stop("fewer than 10 usable values to fit a 2-component mixture", call. = FALSE)
  }

  model <- mixtools::normalmixEM(x = x_fit, k = 2, mu = mu,
                                 lambda = c(0.6, 0.4), sigma = c(0.3, 0.2),
                                 arbvar = FALSE, maxit = 10000,
                                 maxrestarts = 200, epsilon = 1e-20)

  index_lower <- which.min(model$mu)
  index_upper <- which.max(model$mu)

  cutoffs <- c(cutoff_50 = normal_mixture_cutoff(model, 0.50, index_lower),
               cutoff_75 = normal_mixture_cutoff(model, 0.75, index_lower))

  if (plot) {
    plot_df <- data.frame(
      x   = x_fit,
      nd1 = model$lambda[1] * stats::dnorm(x_fit, model$mu[1], model$sigma[1]),
      nd2 = model$lambda[2] * stats::dnorm(x_fit, model$mu[2], model$sigma[2])
    )
    p <- mixture_density_plot(
      plot_df,
      component_cols = c("nd1", "nd2"),
      cutoffs = cutoffs[["cutoff_50"]],
      x_label = if (log_data) "log10(counts + 1)" else "counts",
      title = label_name
    )
    print(p)
  } else {
    p <- NULL
  }

  # Back-transform so the call is made on the scale the user passed in.
  cutoffs_out <- if (log_data) (10^cutoffs) - 1 else cutoffs

  # The original indexed the result by a column literally named "counts",
  # which broke whenever the caller passed an expression rather than a
  # variable called `counts`. Build the column explicitly instead.
  calls <- data.frame(counts = counts,
                      bin = as.integer(counts >= cutoffs_out[["cutoff_75"]]))

  fit <- mixture_fit_stats(model, length(x_fit))

  summary_vec <- c(cutoffs_out,
                   mean_upper = model$mu[index_upper],
                   mean_lower = model$mu[index_lower],
                   loglik = fit[["loglik"]],
                   aic    = fit[["aic"]],
                   bic    = fit[["bic"]])

  switch(return_what,
         calls      = calls,
         cutoffs    = cutoffs_out,
         mean_upper = model$mu[index_upper],
         summary    = summary_vec,
         all        = list(calls   = calls,
                           cutoffs = cutoffs_out,
                           summary = summary_vec,
                           model   = model,
                           fit     = fit,
                           plot    = p))
}
