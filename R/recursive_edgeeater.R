#' Collapse a membership graph down to a target vertex count
#'
#' Repeatedly deletes the lowest-ranked vertex with [edgeeat()], transferring
#' its gene membership to its neighbours, until the graph is down to `end_num`
#' vertices. Rank is in-strength plus the vertex's own weight, so the vertex
#' carrying the least membership goes first.
#'
#' @param g An igraph graph from [memb2graph_edgeeat()], with `weight`,
#'   `genes` and `gene_membership` on both vertices and edges.
#' @param end_num Stop once the graph has this many vertices.
#' @param memb_cutoff Gene membership below which a gene is dropped rather
#'   than transferred. Now forwarded to [edgeeat()] (see Fixes).
#' @param verbose If `TRUE` print each deleted vertex index.
#'
#' @return A list of `graph` (the collapsed graph) and `metric_df` (total
#'   membership and gene count after each deletion, for choosing `end_num`).
#'
#' @section Fixes:
#' The returned list element was spelled `metrif_df`, so any caller reaching
#' for `$metric_df` silently got `NULL`. `memb_cutoff` was also accepted and
#' then hard-coded to `0` on the [edgeeat()] call, so it never took effect.
#'
#' @seealso [edgeeat()], [memb2graph_edgeeat()], [multi_memb_prod()]
#'
#' @examples
#' \dontrun{
#' res <- recursive_edgeeater(g, end_num = 12)
#' plot(res$metric_df$gene_num_byV)
#' }
#' @export
recursive_edgeeater <- function(g, end_num, memb_cutoff = 0, verbose = FALSE) {
  # iterates through the 
  
  v_num <- length(igraph::V(g))
  memb_byV <- numeric(length = v_num-end_num)
  gene_num_byV <- numeric(length = v_num-end_num)
  count = 1 
  while (v_num >= end_num) { 
    g_rank <- igraph::strength(g, mode ="in") + igraph::V(g)$weight
    vdel <- which.min(g_rank)
    if (verbose) print(vdel)
    # `memb_cutoff` was accepted and then hard-coded to 0 on this call
    g <- edgeeat(g, vdel, memb_cutoff = memb_cutoff) 
    memb_byV[count] <- sum(igraph::V(g)$weight)
    gene_num_byV[count] <- length(unlist(igraph::V(g)$genes))
    # test internal consistency of membership values
    #if (test_graph_membership(g)) {
    #  warning('vertex weight not consistent with gene membership')  
    #}
    v_num <- length(igraph::V(g))
    count = count + 1; 
  }
  metric_df <- data.frame(memb_byV = memb_byV, gene_num_byV = gene_num_byV)
  # was `metrif_df`, so callers reaching for `$metric_df` got NULL
  return(list(graph = g, metric_df = metric_df))
}