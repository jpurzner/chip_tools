#' Collapse transcripts whose ends are close together
#'
#' Within one gene's transcripts, drops transcripts whose end coordinates are
#' within `dist_cutoff` of a retained one, keeping the higher-signal member of
#' each pair and recording what was merged into it. Used to stop one gene's
#' many annotated TSSs from being counted repeatedly.
#'
#' @param gene Data frame of one gene's transcripts, with `end`, `max_val` and
#'   `ensembl_name` columns.
#' @param dist_cutoff Distance in bases below which two ends are "the same".
#'
#' @return `gene` with redundant rows removed and a `contained_transcript`
#'   column naming the discarded transcripts.
#'
#' @section Note:
#' The original collection also had a `prune_gene.R` defining a second, older
#' copy of `prune_close()` — there was never a `prune_gene()` function. Since
#' both files were sourced, whichever loaded last won. Only this version is
#' carried over.
#'
#' @examples
#' \dontrun{
#' pruned <- do.call(rbind, lapply(split(tx, tx$gene_id), prune_close))
#' }
#' @export
prune_close <- function(gene, dist_cutoff = 250) {
  
  # handle case where too few entries to reduce 
  if (nrow(gene) > 2) { 
    # sort by max value 
    gene <- gene[order(gene$max_val),]
    
    #print(gene)
    
    td <- dist(gene$end, upper = TRUE) 
    td <- as.matrix(td)
    # binarize the distance matrix
    td <- td < dist_cutoff
    # simplify the distance matrix 
    td <- td[rowSums(td) > 1,rowSums(td) > 1]
    # extract index of pairs 
    td <- reshape2::melt(td)
    td <- td[td$value == TRUE,]
    td <- td[!(td$Var1 == td$Var2),]
    
    
    td <- data.frame(t(apply(td[,c(1,2)],1,sort)))
    td <- unique(td)
    # handles the case where there is no overlap  
    if (nrow(td) > 0) {  
      # remove rows that contain information redundant to rows before
      redundant <- c("FALSE", sapply(2:nrow(td), function(n) (td[n,1] %in% unlist(td[1:(n-1),])) & (td[n,2] %in% unlist(td[1:(n-1),]))))
      
      # periodically getting NA's in second row, unsure why but filtering out 
      td <- td[!as.logical(redundant),]
      td <- td[complete.cases(td),]
      colnames(td) <- c("keep", "discard")
      #print(td)
      if (nrow(td) > 0) {
        # record the contained transcripts that are being eliminated
        gene$contained_transcript <- NA
        gene$contained_transcript[td$keep] <- paste(gene$ensembl_name[td$discard], collapse=", ") 
        new_gene <-  gene[-td$discard, ]
      } else {
        new_gene <- gene 
      }
    } else 
      new_gene <- gene 
  } else {
    new_gene <- gene 
  } # end if 
  return(new_gene) 
}