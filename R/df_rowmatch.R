#' Align one data frame onto another's row names, zero-filling gaps
#'
#' Reorders `df_filt` so its rows correspond to `df_main`'s, inserting rows for
#' any row name present in `df_main` but not `df_filt`. Unlike [df_rowfilt()],
#' the inserted values default to zero rather than `NA`, which is usually what
#' you want for counts.
#'
#' @param df_filt Data frame (or object coercible to one) whose columns are kept.
#' @param df_main Data frame whose row names define the output rows.
#' @param na2zero If `TRUE` (default) replace `NA` with `0` after the join.
#'
#' @return `df_filt`'s columns with one row per row of `df_main`.
#'
#' @seealso [df_rowfilt()]
#'
#' @examples
#' \dontrun{
#' peaks_aligned <- df_rowmatch(peak_counts, expression_table)
#' }
#' @export
df_rowmatch <- function(df_filt, df_main, na2zero = T) {

  # check data is df and convert if it isn't 
  if (!is.data.frame(df_filt)) {
    df_filt <- as.data.frame(df_filt)  
  }
  if (!is.data.frame(df_main)) {
    df_main <- as.data.frame(df_main)
  }
  #print(head(df_filt))
  #print(dim(df_filt))
  #print(dim(df_main))

  tmp_merge <- merge(x = df_filt, y = df_main, by.x = 0, by.y = 0, all.y = T)

  #print(dim(df_filt))
  df_new <- tmp_merge[,1:(dim(df_filt)[2]+1)]
  row.names(df_new) <- df_new$Row.names
  df_new$Row.names <- NULL
  colnames(df_new) <- colnames(df_filt)
  if (na2zero) {
    df_new[is.na(df_new)] <- 0  
  } 
  #print(dim(df_new))
  return(df_new)
  
}