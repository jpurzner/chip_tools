#' Plot cluster expression lines for an igraph community graph
#'
#' Flattens the per-vertex gene lists of a graph whose vertices are co-clusters
#' into a gene-to-cluster mapping, then hands it to [plot_cluster_line()].
#'
#' @param g An igraph graph whose vertices carry a `genes` attribute (a list of
#'   gene IDs) and a `name`.
#' @param expr_list List of genes x timepoints expression matrices.
#'
#' @return A `ggplot` object.
#'
#' @seealso [graph_cluster_plot_EBseq()] for the version that also takes
#'   sample metadata and a pseudotime axis; [get_co_clusters()].
#'
#' @examples
#' \dontrun{
#' graph_cluster_plot(community_graph, expr_list)
#' }
#' @export
graph_cluster_plot <- function (g, expr_list) {
  
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

  print(head(gene2cluster))
  print(table(gene2cluster$cluster_id))
  p <- plot_cluster_line(expr_list, 
                         c("Scott_GNP", "Scott_MB", "Frank_Cerebellum"), 
                         gene2cluster, 
                         c("E15_GNP", "P1_GNP", "P7_GNP", "P7", "P14_GNP", "P14", "P28_GNP", "P60", "MB_Ptch_het")) 
  return(p)
} # end community cut   