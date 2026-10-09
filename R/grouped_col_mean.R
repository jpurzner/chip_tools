#' Average count columns within groups
#'
#' Collapses a features x samples count matrix to features x groups by taking
#' the row mean of the columns belonging to each group. Groups with a single
#' sample are carried through unchanged rather than dropped.
#'
#' @param counts Data frame or matrix of counts; one column per sample.
#' @param metadata Data frame with one row per column of `counts`, in the same
#'   order, holding the grouping factor. See [load_metadata()].
#' @param meta_factor Name of the column in `metadata` to group by.
#'
#' @return A data frame with one column per group, in the order the groups
#'   first appear, and the row names of `counts`.
#'
#' @section Fixes:
#' The original split `counts` into single-replicate and multi-replicate
#' halves, then indexed the multi-replicate half with positions taken from the
#' *unsplit* grouping vector, so columns were misaligned whenever any group had
#' exactly one sample. Single-replicate groups were also computed into an
#' unused variable and silently dropped from the result, and the output column
#' names were built from a vector of the wrong length. This version groups in
#' one pass.
#'
#' @examples
#' \dontrun{
#' meta <- load_metadata("metadata.txt")
#' grouped_col_mean(counts, meta, "condition")
#' }
#' @export
grouped_col_mean <- function(counts, metadata, meta_factor) {

  if (!meta_factor %in% colnames(metadata)) {
    stop(sQuote(meta_factor), " is not a column of metadata", call. = FALSE)
  }
  if (nrow(metadata) != ncol(counts)) {
    stop("metadata has ", nrow(metadata), " rows but counts has ",
         ncol(counts), " columns; they must correspond one to one",
         call. = FALSE)
  }

  group_index <- as.character(metadata[[meta_factor]])
  groups <- unique(group_index)

  counts_avg <- vapply(
    groups,
    function(g) {
      cols <- which(group_index == g)
      if (length(cols) == 1L) {
        as.numeric(counts[[cols]])
      } else {
        rowMeans(as.matrix(counts[, cols, drop = FALSE]), na.rm = TRUE)
      }
    },
    numeric(nrow(counts))
  )

  counts_avg <- as.data.frame(counts_avg, stringsAsFactors = FALSE)
  colnames(counts_avg) <- groups
  rownames(counts_avg) <- rownames(counts)
  counts_avg
}
