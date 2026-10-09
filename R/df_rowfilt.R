#' Restrict one data frame to the row names of another
#'
#' Reorders and pads `df_filt` so its rows line up with `df_main`'s, filling
#' missing rows with `NA`. Useful for getting two tables into the same row
#' order before `cbind`.
#'
#' @param df_filt Data frame whose columns you want to keep.
#' @param df_main Data frame whose row names define the output rows.
#'
#' @return `df_filt`'s columns, with one row per row of `df_main`.
#'
#' @section Fixes:
#' The original built the aligned frame as `df_new` and then returned
#' `df_filt`, so the function was a no-op that just printed three dimensions.
#'
#' @seealso [df_rowmatch()], which also zero-fills by default.
#'
#' @examples
#' \dontrun{
#' aligned <- df_rowfilt(peak_scores, expression_table)
#' }
#' @export
df_rowfilt <- function(df_filt, df_main) {
  tmp_merge <- merge(x = df_filt, y = df_main, by.x = 0, by.y = 0, all.y = TRUE)
  df_new <- tmp_merge[, seq_len(ncol(df_filt) + 1L), drop = FALSE]
  rownames(df_new) <- df_new$Row.names
  df_new$Row.names <- NULL
  colnames(df_new) <- colnames(df_filt)
  df_new
}
