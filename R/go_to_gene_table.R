#' Turn a topGO result table into a gene-by-term table
#'
#' Explodes the comma-separated `genes` column of a GO enrichment result into
#' one row per gene, then collapses the resulting sparse term columns into as
#' few columns as possible with [collapse_GO_columns()]. The output is meant
#' for printing beside a gene list, where you want each gene's terms without
#' one column per term.
#'
#' @param go_fdf A GO result data frame with `GO.ID`, `Term` and a `genes`
#'   column holding comma-separated gene symbols.
#' @param prefix Column-name prefix for the collapsed term columns.
#'
#' @return A data frame with a `genes` column and `<prefix>1..n` term columns.
#'
#' @seealso [collapse_GO_columns()], [topGO_wrap()], [organize_genes_and_terms()]
#'
#' @examples
#' \dontrun{
#' go_to_gene_table(bp_results, prefix = "BP_term_")
#' }
#' @export
go_to_gene_table <- function(go_fdf, prefix) {
  # Step 1: Separate multiple genes into individual rows
  go_fdf_long <- go_fdf %>%
    tidyr::separate_rows(genes, sep = ",\\s*") %>%  # Split by comma and optional space
    dplyr::mutate(term_combined = paste(GO.ID, Term, sep = ": ")) %>%
    dplyr::select(genes, term_combined)
  
  # Step 2: Pivot the data to create the sparse matrix structure
  go_sparse <- go_fdf_long %>%
    tidyr::pivot_wider(names_from = term_combined, values_from = term_combined, values_fill = NA)
  
  # Step 3: Create a unique list of all genes to ensure all are represented
  all_genes <- unique(go_fdf_long$genes)
  
  # Step 4: Full join with all genes to fill in any missing rows
  sparse_matrix <- dplyr::full_join(data.frame(genes = all_genes), go_sparse, by = "genes")
  
  # Step 5: Collapse the GO term columns with the specified prefix
  collapsed_matrix <- collapse_GO_columns(sparse_matrix, prefix)
  
  return(collapsed_matrix)
}
