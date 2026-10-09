#' Group GO terms into clusters by semantic similarity
#'
#' Walks a table of GO terms in order and assigns each term the index of the
#' first term it is semantically similar to, so that redundant terms (the
#' parent/child pairs GO enrichment returns in bulk) collapse into one cluster.
#'
#' @param df Data frame of GO results, one row per term.
#' @param GO_sim_mat Square semantic-similarity matrix with GO IDs as both row
#'   and column names, e.g. from `GOSemSim::mgoSim()`.
#' @param GO_name Name of the column in `df` holding the GO ID.
#' @param cutoff Similarity above which two terms are called the same cluster.
#'
#' @return `df` with an integer `GO_cluster` column.
#'
#' @section Fixes:
#' `cutoff` was accepted and then ignored; the comparison used a hard-coded
#' `0.5`, so the argument did nothing.
#'
#' @section Note:
#' This is O(n^2) in the number of terms and assigns greedily in row order, so
#' the clustering depends on how `df` is sorted. Sort by p-value or gene count
#' first if you want the most significant term to name each cluster.
#'
#' @seealso [go2sym()], [topGO_wrap()]
#'
#' @examples
#' \dontrun{
#' sim <- outer(res$GO, res$GO, Vectorize(function(x, y)
#'   GOSemSim::goSim(x, y, semData = mmGO, measure = "Wang")))
#' dimnames(sim) <- list(res$GO, res$GO)
#' clustered <- go_cluster(res, sim, cutoff = 0.7)
#' }
#' @export
go_cluster <- function(df,  GO_sim_mat, GO_name = "GO", cutoff = 0.5) {

df$GO_cluster <- 0

for (i in 1:nrow(df)) {
  curr_go <- df[i,GO_name]
  # `cutoff` was accepted and then ignored - 0.5 was hard-coded here, so
  # passing any other value had no effect.
  assoc_go <- row.names(GO_sim_mat)[GO_sim_mat[match(curr_go, colnames(GO_sim_mat)), ] > cutoff]
  #print(assoc_go)
  assoc_idx <- which(df[,GO_name] %in% assoc_go)
  #print(assoc_idx)
  for (j in assoc_idx) {
    df$GO_cluster[j] <-  ifelse (df$GO_cluster[j] == 0, 
                                 i, 
                                 df$GO_cluster[j])
    
  }
} 
return(df)
} # end of go_cluster