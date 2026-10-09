#' Gene symbols annotated to a GO term and its descendants
#'
#' Looks up every mouse gene symbol annotated to one or more GO terms,
#' optionally walking down the ontology to include descendant terms as well —
#' which is usually what you want, since genes are annotated to the most
#' specific term that applies rather than to the parent you searched for.
#'
#' @param go_name One or more GO IDs. IDs with no annotated mouse gene are
#'   dropped silently.
#' @param ontology `"BP"`, `"MF"` or `"CC"`.
#' @param relationship Ontology relation to expand along, pasted onto
#'   `"GO<ontology>"` to name a GO.db map — `"OFFSPRING"` (default) for all
#'   descendants, `"CHILDREN"` for one level, or `NULL` for the term alone.
#' @param return_df If `TRUE` return a GO-to-symbol data frame with the
#'   *queried* term in the `GO` column (not the descendant the gene was
#'   actually annotated to). If `FALSE` return a character vector of symbols.
#'
#' @return A character vector of symbols, or a data frame of `GO` and
#'   `SYMBOL`. With several `go_name` values and `return_df = FALSE` you get a
#'   named list.
#'
#' @section Note:
#' Mouse only — `org.Mm.eg.db` is hard-wired. The whole ENTREZID-to-GO table
#' is pulled on every call, which is slow; hoist the call out of a loop.
#'
#' @seealso [go_cluster()], [topGO_wrap()], [go_to_gene_table()]
#'
#' @examples
#' \dontrun{
#' # genes under "nervous system development" and everything below it
#' go2sym("GO:0007399", ontology = "BP")
#' }
#' @export
go2sym <- function(go_name, ontology = "BP", relationship = "OFFSPRING", return_df = FALSE) {
  require_pkg("GO.db", "org.Mm.eg.db", "AnnotationDbi")


  
cols <- c("GO", "SYMBOL")
db <- getExportedValue("org.Mm.eg.db", "org.Mm.eg.db")
all_keys <- AnnotationDbi::keys(db)
go_df <- AnnotationDbi::select(db, all_keys, cols, keytype = "ENTREZID")



go_name <- go_name[go_name %in% unique(go_df$GO)]

go2sym1 <- function(go_name1, ontology = ontology, relationship = relationship) {
  all_go <- as.character(go_name1)

  if (!is.null(relationship)) {
    go_str <- paste0("GO", ontology, relationship) 
    #print(go_str)
    all_go <-  c(all_go, get(as.character(go_name1), eval(parse(text = go_str))))
  } # end of if 

  #print(all_go)
  if (return_df) {
    all_sym <- go_df[(go_df$GO  %in% all_go) & (go_df$ONTOLOGY == ontology),]
    all_sym$EVIDENCE <- NULL
    all_sym$ENTREZID <- NULL
    all_sym$GO <- go_name1
    all_sym <- all_sym[!duplicated(all_sym),]
    # overwrite the GO term to the parent term 
  } else {
  # return just the genes names
    # column 5 was positional; name it so a schema change in org.Mm.eg.db
  # cannot silently return the wrong column
    all_sym <- unique(go_df[(go_df$GO %in% all_go) & (go_df$ONTOLOGY == ontology), "SYMBOL"])
  }
} # end go2sym1


if (length(go_name)  == 1) {
  all_sym <- go2sym1(go_name, ontology = ontology, relationship = relationship)
} else {
  all_sym <- lapply(go_name, function (x)  go2sym1(x, ontology = ontology, relationship = relationship) )
  names(all_sym) <- go_name
  all_sym <- dplyr::bind_rows(all_sym) 
}


return(all_sym)

} # end of go2sym 