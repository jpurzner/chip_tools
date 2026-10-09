#' Order a gene-by-term table by semantic term order and gene similarity
#'
#' Reorders the columns of a gene-by-term table to follow a supplied GO term
#' ordering (e.g. from GOSemSim semantic clustering), then reorders the rows by
#' hierarchical clustering of the genes on binary (Jaccard) distance over their
#' term membership. Genes with no remaining terms are appended at the end.
#'
#' @param long_data Data frame whose first column is `genes` and whose
#'   remaining column names look like `"GO:0007399: nervous system development"`.
#' @param ordered_go_terms Character vector of GO IDs in the desired order.
#'
#' @return `long_data` with reordered rows and columns.
#'
#' @seealso [go_to_gene_table()], [topGO_wrap()]
#'
#' @examples
#' \dontrun{
#' organize_genes_and_terms(gene_term_table, semantic_order)
#' }
#' @export

organize_genes_and_terms <- function(long_data, ordered_go_terms) {
  # Step 1: Extract GO.ID from go_columns
  go_columns <- colnames(long_data)[-1]  # Exclude "genes" column
  go_ids_in_columns <- sub(": .*", "", go_columns)   # Extract GO.ID portion
  
  # Step 2: Reorder columns based on ordered GO terms
  # Match extracted GO.IDs with ordered_go_terms to maintain semantic ordering
  matched_indices <- match(ordered_go_terms, go_ids_in_columns, nomatch = 0)
  reordered_columns <- c("genes", go_columns[matched_indices[matched_indices > 0]])
  
  # Subset and reorder long_data by reordered_columns
  ordered_long_data <- long_data[, reordered_columns, drop = FALSE]
  
  # Step 3: Identify and retain rows with at least one non-NA value in GO term columns
  non_na_rows <- rowSums(!is.na(ordered_long_data[, -1, drop = FALSE])) > 0
  if (!any(non_na_rows)) {
    warning("All GO term columns are empty after filtering. Returning original data structure with NA.")
    #  return(ordered_long_data)
  }
  
  # Retain only rows with at least one non-NA value in the GO term columns
  filtered_data <- ordered_long_data[non_na_rows, , drop = FALSE]
  
  # Step 4: Convert GO terms to a binary presence/absence matrix for clustering genes
  go_term_matrix <- filtered_data[, -1]  # Remove genes column for binary conversion
  binary_matrix <- ifelse(!is.na(go_term_matrix), 1, 0)
  rownames(binary_matrix) <- filtered_data$genes  # Set gene names as rownames
  
  # Step 5: Calculate a binary distance matrix (e.g., Jaccard distance)
  gene_distance <- stats::dist(binary_matrix, method = "binary")
  
  # Step 6: Perform hierarchical clustering on genes
  gene_clustering <- stats::hclust(gene_distance, method = "average")
  
  # Step 7: Reorder rows based on clustering
  ordered_genes <- rownames(binary_matrix)[gene_clustering$order]
  organized_matrix <- filtered_data[match(ordered_genes, filtered_data$genes), ]
  
  # Step 8: Re-add rows with only NA values, if needed, at the end
  if (any(!non_na_rows)) {
    na_rows <- ordered_long_data[!non_na_rows, ]
    organized_matrix <- rbind(organized_matrix, na_rows)
  }
  
  return(organized_matrix)
}