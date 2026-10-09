#' Plot cluster expression lines for a community graph, with metadata
#'
#' As [graph_cluster_plot()], but passes sample metadata and a pseudotime axis
#' through to [plot_cluster_line_EBseq()].
#'
#' @param g An igraph graph whose vertices carry a `genes` attribute and a `name`.
#' @param expr_list List of genes x timepoints expression matrices.
#' @param meta_list List of per-data-set sample metadata frames.
#' @param p_time Optional pseudotime axis.
#'
#' @return A `ggplot` object.
#'
#' @seealso [plot_cluster_line_EBseq()], [get_co_clusters()]
#'
#' @examples
#' \dontrun{
#' graph_cluster_plot_EBseq(community_graph, expr_list, meta_list)
#' }
#' @export
graph_cluster_plot_EBseq <- function (g, expr_list, meta_list, p_time = NULL) {
  require_pkg("igraph")

  # comm_obj: community from igraph 
  # g: graph object from igraph 
  # comm_n: number of communities (set to NULL to not cut)
  # co_id: either the results of get_co_clusters or a vector with co_clusters as value and genes as names
  # expr_list: list of expression objects  
  
  v_genes <- V(g)$genes
  v_names <-  V(g)$name
  names(v_genes) <- v_names
  
  flat_v_names <- unlist(lapply(1:length(v_genes), 
                                function(n) rep(names(v_genes[n]), length(v_genes[[n]]))))
  
  flat_v_names <- unlist(flat_v_names)

  gene2cluster <- data.frame(gene = unlist(v_genes), cluster_id =  flat_v_names)

  #print(head(gene2cluster))
  #print(table(gene2cluster$cluster_id))
  p <- plot_cluster_line_EBseq(expr_list, 
                         meta_list, 
                         gene2cluster,
                         p_time) 
  return(p)
} # end community cut   