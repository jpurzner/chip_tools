prune_close <- function(gene, dist_cutoff = 250) {
  
  # handle case where too few entries to reduce 
  if (nrow(gene) > 2) { 
    # sort by max value 
    gene <- gene[order(gene$max_val),]
    
    td <- dist(gene$end, upper = TRUE) 
    td <- as.matrix(td)
    # binarize the distance matrix
    td <- td < dist_cutoff
    # simplify the distance matrix 
    td <- td[rowSums(td) > 1,rowSums(td) > 1]
    # extract index of pairs 
    td <- melt(td)
    td <- td[td$value == TRUE,]
    td <- td[!(td$Var1 == td$Var2),]
    
    
    td <- data.frame(t(apply(td[,c(1,2)],1,sort)))
    td <- unique(td)
    # handles the case where there is no overlap  
    if (nrow(td) > 0) {  
      # remove lines that contain information redundant to lines before
      redundant <- c("FALSE", sapply(2:nrow(td), function(n) (td[n,1] %in% unlist(td[1:(n-1),])) & (td[n,2] %in% unlist(td[1:(n-1),]))))
      td <- td[!as.logical(redundant),]
      colnames(td) <- c("keep", "discard")
      
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