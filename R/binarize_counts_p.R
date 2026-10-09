#' Separate signal from background with a 2-component gamma mixture
#'
#' The gamma counterpart to [binarize_counts()]. Counts bounded below at zero
#' and right-skewed are often fitted better by a gamma mixture than a normal
#' one, particularly before any log transform.
#'
#' Starting values matter a great deal for [mixtools::gammamixEM()], so by
#' default they are estimated from the data quantiles (method of moments:
#' `beta = var/mean`, `alpha = mean/beta`) and the fit is attempted
#' `maxrestarts` times with jittered starts, keeping the highest likelihood.
#'
#' @param counts Numeric vector of counts.
#' @param log_data If `TRUE`, fit on `log10(counts + 1)`.
#' @param exclude_extreme If `TRUE`, restrict the fit to `0 < x < 1`.
#' @param alpha,beta,lambda Optional starting shape, scale and mixing
#'   proportions. Left `NULL` (the default) they are estimated from the data.
#' @param maxit,epsilon Passed to [mixtools::gammamixEM()].
#' @param maxrestarts Number of jittered restarts to try.
#' @param auto_init If `TRUE` (default) always derive starting values from the
#'   data, ignoring `alpha`/`beta`/`lambda`.
#' @param label_name Optional plot title.
#' @param plot If `TRUE` (default) draw the fitted mixture.
#' @param verbose Print the starting values, per-restart likelihoods and
#'   fitted parameters.
#' @param return_what `"calls"` (default) for a data frame of counts and their
#'   0/1 call, `"cutoffs"` for the two cutoffs, or `"all"` for everything plus
#'   the model and plot.
#'
#' @return Depends on `return_what`.
#'
#' @section Zeros:
#' The gamma density is zero at zero, so exact zeros cannot be fitted. They
#' are replaced by one tenth of the smallest non-zero value for the *fit*
#' only; classification still uses the values as supplied.
#'
#' @seealso [binarize_counts()] for the normal-mixture version.
#'
#' @examples
#' \dontrun{
#' calls <- binarize_counts_p(max_expr, verbose = TRUE)
#' table(calls$bin)
#' }
#' @import ggplot2
#' @export
binarize_counts_p <- function(counts,
                              log_data = FALSE,
                              exclude_extreme = FALSE,
                              alpha = NULL,
                              beta = NULL,
                              lambda = NULL,
                              maxit = 10000,
                              maxrestarts = 20,
                              epsilon = 1e-8,
                              auto_init = TRUE,
                              label_name = NULL,
                              plot = TRUE,
                              verbose = FALSE,
                              return_what = c("calls", "cutoffs", "all")) {

  require_pkg("mixtools")
  return_what <- match.arg(return_what)

  x_all <- if (log_data) log10(counts + 1) else counts

  # Fit on the transformed values. The original subset `counts` here, which
  # discarded the log transform when both flags were set.
  x_fit <- x_all[is.finite(x_all) & x_all >= 0]
  if (exclude_extreme) {
    x_fit <- x_fit[x_fit > 0 & x_fit < 1]
  }
  if (length(x_fit) < 10) {
    stop("fewer than 10 usable non-negative values to fit a gamma mixture",
         call. = FALSE)
  }

  # Gamma support is (0, Inf): nudge exact zeros off the boundary.
  nonzero <- x_fit[x_fit > 0]
  if (length(nonzero) == 0L) {
    stop("all values are zero; nothing to fit", call. = FALSE)
  }
  x_fit[x_fit == 0] <- min(nonzero) / 10

  if (auto_init || is.null(alpha) || is.null(beta) || is.null(lambda)) {
    q25 <- stats::quantile(x_fit, 0.25, names = FALSE)
    q75 <- stats::quantile(x_fit, 0.75, names = FALSE)
    data_var <- stats::var(x_fit)

    # Method of moments for each component: mean = alpha*beta, var = alpha*beta^2
    moment_start <- function(mean_hat, var_hat) {
      if (!is.finite(var_hat) || var_hat <= 0) var_hat <- 0.01
      if (!is.finite(mean_hat) || mean_hat <= 0) mean_hat <- 0.01
      b <- max(var_hat / mean_hat, 1e-3)
      c(alpha = max(mean_hat / b, 0.1), beta = b)
    }

    lower_start <- moment_start(q25, (q25^2) / 2)
    upper_start <- moment_start(q75, data_var)

    if (is.null(alpha))  alpha  <- c(lower_start[["alpha"]], upper_start[["alpha"]])
    if (is.null(beta))   beta   <- c(lower_start[["beta"]],  upper_start[["beta"]])
    if (is.null(lambda)) lambda <- c(0.5, 0.5)

    if (verbose) {
      cat("Auto-estimated starting parameters:\n")
      cat("  alpha: ", paste(signif(alpha, 4), collapse = ", "), "\n", sep = "")
      cat("  beta:  ", paste(signif(beta, 4), collapse = ", "), "\n", sep = "")
      cat("  lambda:", paste(signif(lambda, 4), collapse = ", "), "\n", sep = "")
    }
  }

  best_model <- NULL
  best_loglik <- -Inf

  for (restart in seq_len(maxrestarts)) {
    if (restart == 1L) {
      alpha_try <- alpha; beta_try <- beta; lambda_try <- lambda
    } else {
      jitter <- stats::runif(2, 0.5, 1.5)
      alpha_try <- alpha * jitter
      beta_try  <- beta * jitter
      w <- stats::runif(1, 0.3, 0.7)
      lambda_try <- c(w, 1 - w)
    }

    model <- tryCatch(
      suppressWarnings(
        mixtools::gammamixEM(x = x_fit, k = 2, alpha = alpha_try,
                             beta = beta_try, lambda = lambda_try,
                             maxit = maxit, maxrestarts = 3,
                             epsilon = epsilon, verb = FALSE)
      ),
      error = function(e) {
        if (verbose) cat("Restart", restart, "failed:", conditionMessage(e), "\n")
        NULL
      }
    )

    if (!is.null(model) && is.finite(model$loglik %||% NA) &&
        model$loglik > best_loglik) {
      best_loglik <- model$loglik
      best_model <- model
      if (verbose) cat("Restart", restart, "- loglik:", model$loglik, "\n")
    }
  }

  if (is.null(best_model)) {
    stop("failed to fit a gamma mixture in ", maxrestarts, " attempts; try ",
         "supplying alpha/beta/lambda with auto_init = FALSE", call. = FALSE)
  }
  model <- best_model

  # mixtools stores gamma.pars as a 2 x k matrix of (alpha, beta); mean = alpha * beta
  comp_mean <- model$gamma.pars[1, ] * model$gamma.pars[2, ]
  index_lower <- which.min(comp_mean)

  if (verbose) {
    cat("\nFitted model:\n")
    for (i in seq_len(2)) {
      cat(sprintf("  Component %d - alpha: %.4g beta: %.4g lambda: %.4g mean: %.4g\n",
                  i, model$gamma.pars[1, i], model$gamma.pars[2, i],
                  model$lambda[i], comp_mean[i]))
    }
    cat("  Lower component:", index_lower, "\n")
  }

  gamma_posterior_cutoff <- function(proba, i) {
    f <- function(v) {
      dens <- model$lambda * stats::dgamma(v, shape = model$gamma.pars[1, ],
                                           rate = 1 / model$gamma.pars[2, ])
      proba - (dens[i] / sum(dens))
    }
    # Scan for a sign change rather than guessing a bracket; gamma mixtures
    # can be flat over long stretches and uniroot() needs f(lower)*f(upper) < 0.
    grid <- seq(min(x_fit), max(x_fit), length.out = 1000)
    vals <- vapply(grid, f, numeric(1))
    crossings <- which(diff(sign(vals)) != 0)
    if (length(crossings) == 0L) {
      warning("no posterior crossing found; using the midpoint of the ",
              "component means", call. = FALSE)
      return(mean(comp_mean))
    }
    tryCatch(
      stats::uniroot(f, lower = grid[crossings[1]],
                     upper = grid[crossings[1] + 1])$root,
      error = function(e) mean(comp_mean)
    )
  }

  cutoffs <- c(cutoff_50 = gamma_posterior_cutoff(0.50, index_lower),
               cutoff_75 = gamma_posterior_cutoff(0.75, index_lower))

  if (plot) {
    plot_df <- data.frame(x = x_fit)
    for (i in seq_len(2)) {
      plot_df[[paste0("nd", i)]] <-
        model$lambda[i] * stats::dgamma(x_fit,
                                        shape = model$gamma.pars[1, i],
                                        rate = 1 / model$gamma.pars[2, i])
    }
    p <- mixture_density_plot(
      plot_df,
      component_cols = c("nd1", "nd2"),
      cutoffs = cutoffs,
      x_label = if (log_data) "log10(counts + 1)" else "counts",
      title = label_name %||%
        sprintf("Gamma mixture (loglik: %s)", round(model$loglik, 2))
    )
    print(p)
  } else {
    p <- NULL
  }

  cutoffs_out <- if (log_data) (10^cutoffs) - 1 else cutoffs

  calls <- data.frame(counts = counts,
                      bin = as.integer(counts >= cutoffs_out[["cutoff_75"]]))

  switch(return_what,
         calls   = calls,
         cutoffs = cutoffs_out,
         all     = list(calls = calls, cutoffs = cutoffs_out,
                        model = model, plot = p))
}
