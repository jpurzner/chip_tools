#' topGO_wrap: Perform Gene Ontology (GO) Enrichment Analysis
#'
#' This function performs GO enrichment analysis using the `topGO` package on a specified subset of genes.
#' It organizes results based on semantic similarity of GO terms and can output data in either wide or long format.
#' Additionally, it allows the removal of GO terms with p-values above a specified threshold.
#'
#' @param data A data frame containing gene information, with columns for grouping criteria and gene identifiers.
#' @param set1 Column name for the primary grouping criteria.
#' @param name Column name containing gene identifiers.
#' @param set2 Optional second column for additional grouping criteria (default is NULL).
#' @param min_genes Minimum number of genes in a group to be included in analysis (default is 10).
#' @param max_genes Maximum number of genes in a group to be included in analysis (default is 2000).
#' @param ontology GO ontology to use ("BP", "MF", or "CC") (default is "MF").
#' @param topNodes Number of top GO nodes to retrieve (default is 70).
#' @param test If TRUE, processes only the first two groups for testing purposes (default is FALSE).
#' @param order_by_similarity If TRUE, orders GO terms by semantic similarity (default is TRUE).
#' @param genes_as_row If TRUE, formats results with genes as rows and GO terms as columns (default is FALSE).
#' @param max_pval Maximum p-value threshold for filtering GO terms (default is 0.01).
#' @param pval_column Column to use for p-value filtering; options are "classicFisher", "elimFisher", "topgoFisher" (default is "elimFisher").
#'
#' @return A data frame containing GO enrichment results, with options for ordering by semantic similarity and filtering by p-value.
#'
#' @examples
#' \dontrun{
#' # needs org.Mm.eg.db and real mouse gene symbols
#' go <- topGO_wrap(gene_table, set1 = "cluster_id", name = "mgi_symbol")
#' }
#'
#' @export
topGO_wrap <- function(data, set1, name, 
                       set2 = NULL,  
                       min_genes = 10, 
                       max_genes = 2000,
                       ontology = "MF", 
                       topNodes = 70, 
                       test = FALSE,
                       order_by_similarity = TRUE,
                       genes_as_row = FALSE,
                       max_pval = 0.01,                 # Maximum p-value threshold for filtering GO terms
                       pval_column = "elimFisher") {    # Column to use for p-value filtering, default is elimFisher
  
  # Description:
  # topGO_wrap performs Gene Ontology (GO) enrichment analysis using the `topGO` package on a specified subset of genes.
  # This function organizes results based on semantic similarity of GO terms and can output data in either wide or long format.
  # Additionally, it allows the removal of GO terms with p-values above a specified threshold.
  
  # Verify that the 'name' column exists in 'data'
  require_pkg("org.Mm.eg.db", "GO.db", "reshape", "stringr", "parallel",
              "topGO", "GOSemSim", "cluster")

  # Verify that the 'name' column exists in 'data'
  if (!name %in% colnames(data)) {
    stop("the 'name' column ", sQuote(name), " does not exist in 'data'",
         call. = FALSE)
  }
  
  # Load required packages
  
  # Set up parallel processing
  no_cores <- parallel::detectCores() - 2
  cl <- parallel::makeCluster(no_cores)
  # Workers need topGO and the annotation db, not GOstats - gostarter()
  # below calls topGO. Loading the wrong package here made every worker fail.
  parallel::clusterEvalQ(cl, {
    library(topGO)
    library(org.Mm.eg.db)
  })
  
  # Create a contingency table based on the provided columns (set1 and optional set2)
  if (!is.null(set2)) {
    trans_mat <- table(data[, set1], data[, set2])
    trans_mat <- reshape::melt(trans_mat)
    gene_list <- lapply(seq_len(nrow(trans_mat)), function(x) {
      data[(data[[set1]] == trans_mat[x, 1]) & (data[[set2]] == trans_mat[x, 2]), name]
    })
    names(gene_list) <- paste(trans_mat[, 1], trans_mat[, 2], sep = "_")
    trans_mat$pair_id <- paste(trans_mat[, 1], trans_mat[, 2], sep = "_")
    colnames(trans_mat)[1] <- set1
    colnames(trans_mat)[2] <- set2
    colnames(trans_mat)[3] <- "gene_num"
  } else {
    trans_mat <- table(data[, set1])
    trans_mat <- data.frame(set1 = names(trans_mat), gene_num = as.vector(trans_mat))
    gene_list <- lapply(trans_mat$set1, function(x) data[data[[set1]] == x, name])
    names(gene_list) <- trans_mat$set1
    trans_mat$pair_id <- trans_mat$set1
  }
  
  # Filter gene bins based on min_genes and max_genes
  bin_size <- sapply(names(gene_list), function(x) length(gene_list[[x]]))
  gene_list <- gene_list[bin_size >= min_genes & bin_size <= max_genes]
  if (length(gene_list) == 0) {
    stop("No gene groups meet the min_genes and max_genes criteria.")
  }
  
  universe <- data[[name]]  # Define the universe of genes
  
  # Helper function for topGO enrichment analysis
  gostarter <- function(group, universe) {
    geneList <- factor(as.integer(universe %in% group))
    names(geneList) <- universe
    GOdata <- new("topGOdata", ontology = ontology,
                  allGenes = geneList,
                  description = "test",
                  annot = annFUN.org,
                  mapping = "org.Mm.eg.db",
                  ID = "Symbol")
    resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    resultElim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
    resultTopgo <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
    allRes <- GenTable(GOdata, classicFisher = resultClassic, elimFisher = resultElim,
                       topgoFisher = resultTopgo, orderBy = "elimFisher",
                       ranksOf = "elimFisher", topNodes = topNodes)
    allRes$genes <- sapply(allRes$GO.ID, function(go_term) {
      genes_in_term <- genesInTerm(GOdata, go_term)[[1]]
      test_genes <- genes_in_term[genes_in_term %in% group]
      paste(test_genes, collapse = ", ")
    })
    return(allRes)
  }
  
  # Run gostarter function on each gene group in parallel
  parallel::clusterExport(cl, varlist = c("gene_list", "universe", "gostarter"), envir = environment())
  all_go <- if (test) parallel::parLapply(cl, names(gene_list[1:2]), function(x) gostarter(gene_list[[x]], universe))
  else parallel::parLapply(cl, names(gene_list), function(x) gostarter(gene_list[[x]], universe))
  parallel::stopCluster(cl)
  
  # Combine and name the GO results
  names(all_go) <- names(gene_list)[1:length(all_go)]
  go_df <- all_go[sapply(1:length(all_go), function(n) dim(all_go[[n]])[1]) > 0]
  go_fdf <- do.call("rbind", go_df)
  flat_vtx_out <- unlist(lapply(1:length(go_df), function(n) rep(names(go_df[n]), dim(go_df[[n]])[1])))
  go_fdf <- cbind(flat_vtx_out, go_fdf)
  colnames(go_fdf)[1] <- "pair_id"
  go_fdf <- merge(x = go_fdf, y = trans_mat, by.x = 1, by.y = "pair_id", all.x = TRUE)
  
  # Filter based on max p-value threshold in the specified p-value column
  go_fdf <- go_fdf[go_fdf[[pval_column]] <= max_pval, ]
  
  # Order by semantic similarity if specified
  if (order_by_similarity) {
    go_terms <- as.character(go_fdf$GO.ID)
    valid_go_terms <- go_terms[sapply(go_terms, function(go) {
      go_info <- GOTERM[[go]]
      !is.null(go_info) && Ontology(go_info) == ontology
    })]
    sem_data <- godata("org.Mm.eg.db", ont = ontology, keytype = "GO", annoDb = "org.Mm.eg.db")
    sim_matrix <- mgeneSim(valid_go_terms, semData = sem_data, measure = "Wang")
    dissimilarity <- as.dist(1 - sim_matrix)
    hc <- hclust(dissimilarity, method = "average")
    ordered_go_terms <- valid_go_terms[hc$order]
    go_fdf <- go_fdf[match(ordered_go_terms, go_fdf$GO.ID), ]
  }
  
  # Format results with genes as rows and GO terms as columns if genes_as_row is TRUE
  if (genes_as_row) {
    long_data <- go_fdf %>%
      tidyr::separate_rows(genes, sep = ",\\s*") %>%
      dplyr::mutate(term_combined = paste(GO.ID, Term, sep = ": ")) %>%
      dplyr::select(genes, term_combined) %>%
      tidyr::pivot_wider(names_from = term_combined, values_from = term_combined, values_fill = NA)
    all_genes <- unique(unlist(gene_list))
    long_data <- dplyr::full_join(data.frame(genes = all_genes), long_data, by = "genes")
    go_fdf <- organize_genes_and_terms(long_data, ordered_go_terms)
    go_fdf <- collapse_GO_columns(go_fdf)
  }
  
  return(go_fdf)
} # End of topGO_wrap
