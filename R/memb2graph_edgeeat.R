#' Build a directed co-cluster graph from membership products
#'
#' Turns the genes x combinations matrix from [multi_memb_prod()] into a
#' directed graph. Each gene's top-ranked co-cluster becomes a vertex and its
#' 2nd..`steps` ranked co-clusters become a chain of edges, so an edge carries
#' the membership that a gene would contribute to its next-best home. Vertex
#' and edge weights are those memberships summed over genes.
#'
#' The point of the construction is that [recursive_edgeeater()] can then
#' collapse the graph: the chain structure is what lets membership flow to the
#' next-best cluster when a vertex is deleted.
#'
#' @param memb_prod Genes x combinations membership matrix from
#'   [multi_memb_prod()].
#' @param memb_prod_rank Same shape as `memb_prod`, holding each gene's rank
#'   of each combination (1 = best).
#' @param steps How many ranks to chain per gene.
#' @param plot If `TRUE` draw the membership-by-rank violin diagnostic, which
#'   shows how fast membership falls off with rank and so how large `steps`
#'   needs to be. The original always drew it.
#'
#' @return An igraph graph with `weight`, `genes` and `gene_membership` on
#'   both vertices and edges. Edges to co-clusters that are no gene's top rank
#'   are pruned.
#'
#' @seealso [multi_memb_prod()] for the input, [recursive_edgeeater()] and
#'   [edgeeat()] to collapse the result, [graph_cluster_plot()] to plot it
#'
#' @examples
#' \dontrun{
#' memb_prod <- multi_memb_prod(matched_membership_list)
#' rank_mat <- t(apply(-memb_prod, 1, rank, ties.method = "first"))
#' g <- memb2graph_edgeeat(memb_prod, rank_mat, steps = 5, plot = TRUE)
#' }
#' @import ggplot2
#' @importFrom rlang .data
#' @export
memb2graph_edgeeat <- function(memb_prod, memb_prod_rank, steps = 5, plot = FALSE) {
  require_pkg("igraph")

  
  
  # trying to find out why NA values showing up in graph
  #print(memb_prod)
  
  # generate a df of membership by rank for each gene 
  memb_prod_byrank <- lapply(1:(dim(memb_prod_rank)[1]), function(n) memb_prod[n,match(1:steps,memb_prod_rank[n,])])
  memb_prod_byrank <- data.frame(matrix(unlist(memb_prod_byrank), nrow=dim(memb_prod_rank)[1], byrow=T))
  row.names(memb_prod_byrank) <- row.names(memb_prod)
  colnames(memb_prod_byrank) <- seq(1:steps)
  byrank_molten <- reshape2::melt(memb_prod_byrank)
  # The rank-membership violin is a diagnostic, not the function's output;
  # it used to be drawn unconditionally.
  if (plot) {
    print(ggplot(byrank_molten, aes(x = .data$variable, y = .data$value)) +
            geom_violin() +
            labs(x = "co-cluster rank", y = "membership product"))
  }
 
  
  # generate a df of coclust by rank for each gene
  coclust_index <- lapply(1:(dim(memb_prod_rank)[1]), function(n) colnames(memb_prod_rank)[match(1:steps,memb_prod_rank[n,])])
  coclust_index <- data.frame(matrix(unlist(coclust_index), nrow=dim(memb_prod_rank)[1], byrow=T))
  colnames(coclust_index) <- colnames(memb_prod_byrank)
  row.names(coclust_index) <- row.names(memb_prod_byrank)
  #print(head(coclust_index))
  
  #Trouble shooting in case the row.names give error
  #print(length(row.names(memb_prod_byrank)))
  #print(length(c((1+((2-1)*dim(memb_prod_byrank)[1])):(dim(memb_prod_byrank)[1]*2))))
  
  # generate a flat df with node1 node2 and membership transfer for each gene 
  
  # in this approach we are going represent each edge as a contribution of membership to another node 
  # to do this each gene will contribute the value of the membership for rank 2,3,4,5 
  # these will then be summed as before to determine the total membership contributin of the edge 
  
  edge_df <- lapply(c(2:steps), 
                    function(n) data.frame(
                      gene = row.names(memb_prod_byrank), 
                      from = coclust_index[,n-1],
                      to = coclust_index[,n],
                      membership = memb_prod_byrank[,n],  
                      row.names = c((1+((n-1)*dim(memb_prod_byrank)[1])):(dim(memb_prod_byrank)[1]*n))))
  
  edge_df <- do.call("rbind", edge_df)
  
  # collapse by sum along all unige N1 and N2 
  edge_df_collapse <- plyr::ddply(edge_df, c("from", "to"), summarise, 
                            weight = sum(membership),
                            genes = paste(gene, collapse =  ","),
                            gene_membership = paste(membership, collapse =  ","))
  
  # set vertices as the gene number * membership 
  vertex_df <-  data.frame(
                      gene = row.names(memb_prod_byrank), 
                      name = coclust_index[,1],
                      membership = memb_prod_byrank[,1], 
                      row.names = c(1:(dim(memb_prod_byrank)[1])))
  
  
  
  vertex_df_collapse <- plyr::ddply(vertex_df, "name", summarise, 
                              weight = sum(membership),
                              genes = paste(gene, collapse =  ","),
                              gene_membership = paste(membership, collapse =  ","))
  #return(vertex_df_collapse)
  
  vertex_genes <- strsplit(vertex_df_collapse[,3], ",")
  vertex_gene_membership <- strsplit(vertex_df_collapse[,4], ",")
  vertex_gene_membership <- lapply(vertex_gene_membership, as.numeric)
  
  vertex_df_collapse$genes <- NULL
  vertex_df_collapse$gene_membership <- NULL
  
  #print(head(vertex_df_collapse))
    
  # prune all edges that are not connected represented in the rank1 only vertexes. 
  edge_df_collapse <- subset(edge_df_collapse, (edge_df_collapse$from %in% vertex_df_collapse$name) & 
                               (edge_df_collapse$to %in% vertex_df_collapse$name))
  
  edge_genes <- strsplit(edge_df_collapse[,4], ",")
  edge_gene_membership <- strsplit(edge_df_collapse[,5], ",")
  edge_gene_membership <- lapply(edge_gene_membership, as.numeric)
  edge_df_collapse$genes <- NULL
  edge_df_collapse$gene_membership <- NULL
  
  g <- igraph::graph_from_data_frame(edge_df_collapse, directed=TRUE, vertices=vertex_df_collapse)
  g <- igraph::set_vertex_attr(g, "genes", value = vertex_genes)
  g <- igraph::set_vertex_attr(g, "gene_membership", value = vertex_gene_membership)
  g <- igraph::set_edge_attr(g, "genes", value = edge_genes)
  g <- igraph::set_edge_attr(g, "gene_membership", value = edge_gene_membership)
  
  
  return(g)
  
}