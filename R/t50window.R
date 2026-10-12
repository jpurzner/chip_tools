#' Bin genes into equal-sized windows along their t50
#'
#' Splits a table of genes and their half-max times (t50) by group, orders each
#' group by t50, and chops it into `window_number` equal-sized windows. Each
#' window becomes a column named after its mean t50, so the result can be fed
#' straight to a per-window GO enrichment.
#'
#' @param df Data frame with gene name in column 1, t50 in column 2, and the
#'   grouping factor in column 3.
#' @param overlap Extra rows added to each window, giving overlapping windows.
#' @param window_number Number of windows per group.
#'
#' @return A list with one element per group: a character matrix of gene names,
#'   `window_number` columns wide, column names being each window's mean t50.
#'
#' @section Note:
#' `floor(n / window_number)` rows go in each window, so up to
#' `window_number - 1` of the highest-t50 genes are dropped. With
#' `overlap > 0` the final window indexes past the end of the vector and those
#' cells come back `NA`.
#'
#' Empty and undersized groups are dropped with a warning. An unused factor
#' level previously produced a zero-row matrix and then failed several steps
#' later with "length of dimnames 2 not equal to array extent" -- so
#' passing a filtered column whose factor still carried its original levels
#' was enough to break the call.
#'
#' @seealso [expr_half_max()] to compute t50, [topGO_timeseries()]
#'
#' @examples
#' \dontrun{
#' windows <- t50window(t50_table, window_number = 10)
#' }
#' @export
t50window <- function(df, overlap = 0, window_number = 10) {

  stopifnot(ncol(df) >= 3, window_number >= 1)

  # droplevels() first: an unused factor level otherwise becomes a zero-row
  # group that fails downstream rather than here.
  group <- droplevels(as.factor(df[[3]]))
  max_list <- split(df, f = group)
  max_list <- lapply(max_list, function(x) x[order(x[[2]]), ])

  # A group needs at least one row per window to be binned at all.
  too_small <- vapply(max_list, nrow, integer(1)) < window_number
  if (any(too_small)) {
    warning("dropping group(s) with fewer than ", window_number, " rows: ",
            paste(names(max_list)[too_small], collapse = ", "), call. = FALSE)
    max_list <- max_list[!too_small]
  }
  if (length(max_list) == 0L) {
    stop("no group has at least ", window_number, " rows", call. = FALSE)
  }

  # Uniform window size, so gene names fit a (rows x window_number) matrix.
  vect2mat <- function(v, overlap = 0, window_number = 10) {
    r <- floor(length(v) / window_number) + overlap
    m <- matrix(NA_character_, r, window_number)
    for (i in seq_len(window_number)) {
      idx <- seq.int(r * (i - 1) + 1, length.out = r)
      idx[idx > length(v)] <- NA_integer_
      m[, i] <- as.character(v[idx])
    }
    m
  }

  gene_list <- lapply(max_list, function(x)
    vect2mat(x[[1]], overlap = overlap, window_number = window_number))
  t50_list <- lapply(max_list, function(x)
    vect2mat(x[[2]], overlap = overlap, window_number = window_number))

  t50_summary <- lapply(t50_list, function(y) {
    nums <- matrix(as.numeric(y), nrow = nrow(y))
    colMeans(nums, na.rm = TRUE)
  })

  out <- lapply(seq_along(gene_list), function(n) {
    m <- gene_list[[n]]
    colnames(m) <- t50_summary[[n]]
    m
  })
  names(out) <- names(gene_list)
  out
}
