#' Membership product across every cross-data-set cluster combination
#'
#' Given fuzzy-clustering membership matrices from several data sets, forms
#' every combination of one cluster per data set and scores each gene by the
#' product of its memberships across that combination. A gene scores highly
#' only if it belongs to the corresponding cluster in *every* data set, which
#' is what makes the product the right operator here.
#'
#' @param all_memb_matched Named list of genes x clusters membership matrices,
#'   one per data set, all with the same genes in the same row order and with
#'   cluster labels already matched across data sets (see [opti_map()]).
#'
#' @return A genes x combinations matrix. Column names are the per-data-set
#'   cluster indices joined with `-`, e.g. `"2-1-3"`.
#'
#' @section Note:
#' The number of columns is the product of the per-data-set cluster counts, so
#' this grows fast: five data sets at eight clusters each is 32,768 columns.
#'
#' @seealso [memb2graph_edgeeat()] to turn this into a graph,
#'   [get_co_clusters()] for the hard-assignment equivalent, [opti_map()]
#'
#' @examples
#' \dontrun{
#' memb_prod <- multi_memb_prod(matched_membership_list)
#' }
#' @export
multi_memb_prod <- function(all_memb_matched) {
  require_pkg("pbapply")

  # product of number of clusters for each dataset 
  co_cluster_number <- sapply(1:length(all_memb_matched), function (d) dim(all_memb_matched[[d]])[2])
  # generate the index matrix 
  co_cluster_index <- expand.grid(lapply(co_cluster_number, function(x) seq(1:x)))
  # create hyphenated names 
  idmaker = function(vec){
    return(paste(vec, collapse="-"))
  }
  co_id_label <- apply(as.matrix(co_cluster_index), 1, idmaker)
  
  memb_prod <- function(i) {
    ival <- sapply(1:(dim(co_cluster_index)[2]), 
                function (ds) all_memb_matched[[ds]][,co_cluster_index[i,ds]])
    return(apply(ival, 1, prod))
  }
  
  memb_prod <- pbapply::pbsapply(1:dim(co_cluster_index)[1], function(i) memb_prod(i))
  colnames(memb_prod) <- co_id_label
  row.names(memb_prod) <- row.names(all_memb_matched[[1]])
  return(memb_prod)
  
} # end of multi_overlap