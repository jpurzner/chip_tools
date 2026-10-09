#' Pack sparse GO term columns into as few columns as possible
#'
#' A gene-by-term table is mostly empty: each gene carries a handful of terms
#' out of hundreds of columns. This packs the term columns into the smallest
#' number of columns that have no two terms competing for the same gene, by
#' interval-graph style greedy placement with rarer terms placed first.
#'
#' @param sparse_matrix Data frame whose first column is gene names and whose
#'   remaining columns each hold a term label or `NA`.
#' @param prefix Column-name prefix for the packed columns.
#'
#' @return A data frame of `genes` plus `<prefix>1..n`.
#'
#' @seealso [go_to_gene_table()], which calls this.
#'
#' @examples
#' \dontrun{
#' collapse_GO_columns(sparse_go, prefix = "BP_term_")
#' }
#' @export
collapse_GO_columns <- function(sparse_matrix, prefix = "GO_term_") {
  # Separate the gene names column from the rest of the matrix
  gene_names <- sparse_matrix[, 1]
  go_term_matrix <- sparse_matrix[, -1]
  
  # Step 1: Sort GO term columns by increasing frequency (rarer terms first)
  term_frequencies <- colSums(!is.na(go_term_matrix))
  sorted_indices <- order(term_frequencies, decreasing = FALSE)
  go_term_matrix <- go_term_matrix[, sorted_indices]
  
  # Initialize an empty list to store final columns
  final_columns <- list()
  
  # Step 2: Place each column in the first available non-overlapping position
  for (col_idx in 1:ncol(go_term_matrix)) {
    term_data <- go_term_matrix[, col_idx]
    placed <- FALSE
    
    # Try to place term_data into an existing final column if no overlap
    for (final_idx in seq_along(final_columns)) {
      if (all(is.na(final_columns[[final_idx]]) | is.na(term_data))) {
        # No overlap, so place term_data in this final column
        final_columns[[final_idx]] <- ifelse(is.na(final_columns[[final_idx]]), term_data, final_columns[[final_idx]])
        placed <- TRUE
        break
      }
    }
    
    # If no suitable column was found, create a new column
    if (!placed) {
      final_columns[[length(final_columns) + 1]] <- term_data
    }
  }
  
  # Step 3: Convert the list of final columns to a data frame
  collapsed_go_term_matrix <- as.data.frame(final_columns)
  
  # Step 4: Rename columns with the specified prefix, e.g., "BP_Term_1", "KEGG_Term_1"
  colnames(collapsed_go_term_matrix) <- paste0(prefix, seq_len(ncol(collapsed_go_term_matrix)))
  
  # Step 5: Combine gene names and the collapsed GO term matrix
  collapsed_matrix <- cbind(genes = gene_names, collapsed_go_term_matrix)
  return(collapsed_matrix)
}


