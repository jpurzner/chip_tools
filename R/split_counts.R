#' Choose the number of mixture components by BIC
#'
#' Like [binarize_counts()], but instead of assuming two components it fits
#' `k = min_k:max_k` normal mixtures and keeps the one with the lowest BIC.
#' Use it when you do not know in advance how many states the signal has.
#'
#' @param counts Numeric vector of counts (or any continuous signal).
#' @param log_data If `TRUE`, fit on `log10(counts + 1)`.
#' @param exclude_extreme If `TRUE`, restrict the fit to `0.1 < x < 1`.
#' @param min_k,max_k Range of component counts to try.
#' @param label_name Optional plot title.
#' @param plot If `TRUE` (default) draw the winning fit.
#' @param return_what `"model"` (default) for the selected `mixEM` object,
#'   `"comparison"` for the per-`k` fit table, or `"all"` for both plus the
#'   plot.
#'
#' @return Depends on `return_what`.
#'
#' @seealso [binarize_counts()], [trinarize_counts()]
#'
#' @examples
#' \dontrun{
#' fit <- split_counts(counts_vector, log_data = TRUE, max_k = 5)
#' fit$mu
#' split_counts(counts_vector, log_data = TRUE, return_what = "comparison")
#' }
#' @import ggplot2
#' @export
split_counts <- function(counts,
                         log_data = FALSE,
                         exclude_extreme = FALSE,
                         min_k = 2L,
                         max_k = 4L,
                         label_name = NULL,
                         plot = TRUE,
                         return_what = c("model", "comparison", "all")) {

  require_pkg("mixtools")
  return_what <- match.arg(return_what)
  stopifnot(min_k >= 1L, max_k >= min_k)

  x_all <- if (log_data) log10(counts + 1) else counts
  x_fit <- x_all[is.finite(x_all)]
  if (exclude_extreme) {
    x_fit <- x_fit[x_fit > 0.1 & x_fit < 1]
  }
  if (length(x_fit) < 10) {
    stop("fewer than 10 usable values to fit a mixture", call. = FALSE)
  }

  ks <- seq.int(min_k, max_k)
  fits <- vector("list", length(ks))
  stats_rows <- vector("list", length(ks))

  for (j in seq_along(ks)) {
    model <- tryCatch(
      mixtools::normalmixEM(x = x_fit, k = ks[j], arbvar = FALSE,
                            maxit = 10000, maxrestarts = 30, epsilon = 1e-20),
      error = function(e) NULL
    )
    if (is.null(model)) next
    fits[[j]] <- model
    # Note: both AIC and BIC are penalised by the actual parameter count here.
    # The original penalised BIC by `k`, the number of components, which made
    # the two criteria incomparable and biased selection toward larger k.
    stats_rows[[j]] <- c(k = ks[j], mixture_fit_stats(model, length(x_fit)))
  }

  keep <- !vapply(fits, is.null, logical(1))
  if (!any(keep)) {
    stop("no mixture converged for k in ", min_k, ":", max_k, call. = FALSE)
  }

  comparison <- as.data.frame(do.call(rbind, stats_rows[keep]))
  best_model <- fits[keep][[which.min(comparison$bic)]]
  best_k <- comparison$k[which.min(comparison$bic)]

  if (plot) {
    ord <- order(best_model$mu)
    plot_df <- data.frame(x = x_fit)
    comp_cols <- paste0("nd", seq_len(best_k))
    for (i in seq_len(best_k)) {
      plot_df[[comp_cols[i]]] <-
        best_model$lambda[ord[i]] * stats::dnorm(x_fit,
                                                 best_model$mu[ord[i]],
                                                 best_model$sigma[ord[i]])
    }
    p <- mixture_density_plot(
      plot_df,
      component_cols = comp_cols,
      component_colours = grDevices::hcl.colors(best_k, "Dark 3"),
      x_label = if (log_data) "log10(counts + 1)" else "counts",
      title = label_name %||% sprintf("Best fit by BIC: k = %d", best_k)
    )
    print(p)
  } else {
    p <- NULL
  }

  switch(return_what,
         model      = best_model,
         comparison = comparison,
         all        = list(model = best_model, k = best_k,
                           comparison = comparison, plot = p))
}
