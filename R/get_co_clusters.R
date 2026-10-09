#' Index genes by their combination of cluster assignments
#'
#' Given fuzzy-clustering membership matrices from several data sets, assigns
#' each gene its hard cluster in each data set, pastes those into a single
#' co-cluster ID, and returns the ID table ranked by how many genes share it.
#' Genes with the same co-cluster ID behave the same way across every data set.
#'
#' @param memb_list Named list of genes x clusters membership matrices, one per
#'   data set, all with the same genes in the same order.
#'
#' @return A list of `index` (co-clusters ranked by gene count) and `bygene`
#'   (each gene's co-cluster ID and per-data-set cluster).
#'
#' @section Fixes:
#' Read an undefined global `fuzz_list` for the data set names, the loop length
#' and the row names, so it only ran if such an object happened to be bound in
#' the caller. Everything now comes from `memb_list`.
#'
#' @seealso [graph_cluster_plot()], [plot_cluster_line()]
#'
#' @examples
#' \dontrun{
#' co <- get_co_clusters(list(scott = m1$membership, hatten = m2$membership))
#' head(co$index)
#' }
#' @export
get_co_clusters <- function(memb_list) {

  # The original read a global called `fuzz_list` for the names, the length
  # and the row names, so it only ran if you happened to have one bound in
  # the calling environment. Everything now comes from `memb_list`.
  ds_names <- names(memb_list)
  if (is.null(ds_names) || anyNA(ds_names) || any(ds_names == "")) {
    stop("memb_list must be a named list, one element per data set",
         call. = FALSE)
  }

  
  # return the co cluster index dataframe ranked by the number of genes in each cluster
  
  idmaker = function(vec){
    return(paste(vec, collapse="-"))
  }
  
  max_ind <- lapply(seq_along(memb_list), function(i) apply(memb_list[[i]], 1, which.max)) 
  memb_ind <-  as.data.frame(matrix(unlist(max_ind), nrow=length(unlist(max_ind[1]))))
  colnames(memb_ind) <- ds_names
  row.names(memb_ind) <- names(max_ind[[1]])

  # determine the co-id for each gene
  co_id <- apply(as.matrix(memb_ind[, ds_names, drop = FALSE]), 1, idmaker)
  memb_ind <- cbind(co_id, memb_ind)
  co_cluster_index <- plyr::ddply(memb_ind, ds_names, "nrow")
  co_cluster_index <- co_cluster_index[order(co_cluster_index$nrow, decreasing = T), ]
  co_id_label <- apply(as.matrix(co_cluster_index[, ds_names, drop = FALSE]), 1, idmaker) 
  row.names(co_cluster_index) <- co_id_label
  co_cluster_index <- co_cluster_index[,c(1:(dim(co_cluster_index)[2])-1)]
  
  return(list("index" = co_cluster_index, "bygene"  = memb_ind))
  
} # end get_co_clusters