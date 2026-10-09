#' Average aligned time courses across data sets
#'
#' Element-wise mean across a list of genes x timepoints matrices that share
#' row and column layout, ignoring `NA`. Used to pool interpolated time courses
#' from different studies into one trajectory per gene.
#'
#' @param expr_align_list List of genes x timepoints matrices, all the same
#'   shape, as returned by [time_expr_interp()].
#'
#' @return A genes x timepoints matrix of means.
#'
#' @seealso [time_expr_interp()], [time_expr_align()]
#'
#' @examples
#' \dontrun{
#' pooled <- time_expr_avg(time_expr_interp(expr_list, time_ref))
#' }
#' @export
time_expr_avg <- function(expr_align_list) {

gene_num <-  dim(expr_align_list[[1]])[1]
expr_avg <- t(sapply(c(1:gene_num), 
                     function(j) apply( sapply(1:length(expr_align_list), 
                                               function(k) expr_align_list[[k]][j,]), 1, 
                                        function (i) mean(i, na.rm = TRUE ))))

row.names(expr_avg) <- row.names(expr_align_list[[1]]) 
return(expr_avg)


}