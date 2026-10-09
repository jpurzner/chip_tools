#' Bin genes into equal-sized windows along their t50
#'
#' Splits a table of genes and their half-max times (t50) by group, orders each
#' group by t50, and chops it into `window_number` equal-sized windows. Each
#' window becomes a column named after its mean t50, so the result can be fed
#' straight to a per-window GO enrichment.
#'
#' @param df Data frame with gene name in column 1, t50 in column 2 and the
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
#' @seealso [expr_half_max()] to compute t50, [topGO_timeseries()]
#'
#' @examples
#' \dontrun{
#' windows <- t50window(t50_table, window_number = 10)
#' }
#' @export
t50window <-  function(df, overlap = 0, window_number= 10) {

  # split list of data-frame by groups 
  max_list <- split( df , f = df[,3] )
  # order the data frames by t50 
  max_list <- lapply(max_list, function(x) x[order(x[,2]),] )

  
  # since uniform size can store gene names as a data-frame with window_size rows and window_number columns
  vect2mat <- function(v, overlap = 0 , window_number = 10 ) {
    r <- floor(length(v)/window_number)
    r <- r + overlap
    m <- matrix(0, r, window_number)
    for (i in c(1:window_number)) {
      begin <- (r*(i-1)+1)
      stop <- (r*(i))
      m[,i] <- as.character(v[c(begin:stop)])  
    } # end of for loop
  return(m)  
  } # end of vect2mat
  
  gene_list <- lapply(max_list, function(x) vect2mat(x[,1], overlap = overlap, window_number = window_number ))
  t50_list <- lapply(max_list, function(x) vect2mat(x[,2], overlap = overlap, window_number = window_number ))
  t50_list <- lapply(t50_list, function(y) data.frame(apply(y, 2, function(x) as.numeric(as.character(x)))))
  t50_summary <-  lapply(t50_list, function(x) as.vector(colMeans(x)))
  gene_list <- lapply(c(1:length(gene_list)), function(n) {
    df = gene_list[[n]]
    colnames(df) <- t50_summary[[n]]
    return(df)
    })
  
  return(gene_list)
    
  
}
