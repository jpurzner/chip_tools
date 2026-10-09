#' GO enrichment for each level of a grouping column
#'
#' Splits a gene table by one column and runs a topGO enrichment per group,
#' one group per core.
#'
#' @param t50_df Data frame with the gene identifier in column 1 and the
#'   grouping column named by `split_cat`.
#' @param all_genes The gene universe to test against.
#' @param split_cat Name or index of the column to split on.
#' @param ont GO ontology: `"BP"`, `"MF"` or `"CC"`.
#' @param test If `TRUE` process only the first two groups.
#'
#' @return A list of topGO result tables, one per group.
#'
#' @section Fixes:
#' The split was driven by an undefined global (`gnp_multi_split[, split_cat]`)
#' rather than the function's own `t50_df`, so it only ran if such an object
#' happened to be bound in the caller. The parallel workers also loaded
#' `GOstats` while the worker function calls topGO, so every worker failed.
#'
#' @seealso [topGO_wrap()] for the maintained version with semantic ordering
#'   and p-value filtering; [topGO_timeseries()]
#'
#' @examples
#' \dontrun{
#' go <- topGO_split(gene_table, universe, split_cat = "cluster_id")
#' }
#' @export
topGO_split <- function(t50_df, all_genes, split_cat,  ont = "BP",  test = FALSE) {
  require_pkg("org.Mm.eg.db", "GO.db", "reshape", "stringr", "parallel", "topGO")


  # set1 and set2 are the  column names or index that define the selected   
  # name is the column name or index that stores the gene identifier 
  
  # examples 
  # tt <- transGO(gnp, "P7_GNP_MB_norep_nok36_chromHMM_GNP_P7_10_gencode.vM1.tss.1000", 
  #                      "P7_GNP_MB_norep_nok36_chromHMM_MB_10_gencode.vM1.tss.1000", 4)
  
  
  
  # setup cores for parallel 
  no_cores <- parallel::detectCores() - 2
  # Initiate cluster
  cl <- parallel::makeCluster(no_cores)  
  

  # load packages to cluster 
  # Workers need topGO and the annotation db, not GOstats - gostarter()
  # below calls topGO. Loading the wrong package here made every worker fail.
  parallel::clusterEvalQ(cl, {
    library(topGO)
    library(org.Mm.eg.db)
  })
  
  # Was `split(t50_df, f = gnp_multi_split[, split_cat])`: the grouping came
  # from an undefined global, so this only ran if you happened to have a
  # `gnp_multi_split` bound in the calling environment.
  gene_list <- split(t50_df, f = t50_df[, split_cat])
  
  first_col <- function(x) {
    return(x[,1])
  }
  
  gene_list <- lapply( gene_list, first_col)
  
  
  universe <- all_genes
  
  ## troubleshooting ---------
  #print(head(universe))
  #print(head(gene_list))
  #group <- gene_list[[1]]
  #geneList <- factor(as.integer(universe %in% group))
  #names(geneList) <- universe
  #print(table(geneList))
  #axon_gene <- go_object['GO:0007411']
  #axon_gene <- unique(unlist(axon_gene, use.names=F))
  
  gostarter <- function(group, universe, ont = ont) {
    geneList <- factor(as.integer(universe %in% group))
    names(geneList) <- universe
    message(head(geneList))
       
    GOdata <- new("topGOdata",ontology = ont,
                  allGenes = geneList,
                  description ="test",
                  annot=annFUN.org, 
                  mapping="org.Mm.eg.db", 
                  ID="Symbol")
    
    # run the Fisher's exact tests
    resultClassic <- runTest(GOdata, algorithm="classic", statistic="fisher")
    resultElim <- runTest(GOdata, algorithm="elim", statistic="fisher")
    resultTopgo <- runTest(GOdata, algorithm="weight01", statistic="fisher")
    
    allRes <- GenTable(GOdata, 
                       classicFisher = resultClassic, 
                       elimFisher = resultElim, 
                       topgoFisher = resultTopgo, 
                       orderBy = "elimFisher", 
                       ranksOf = "elimFisher", topNodes = 400)
    
    print(".")
    return(allRes)
  } # end of gostater
  

  # load variables to cluster 
  clusterExport(cl, varlist = c("gene_list", "universe", "ont", "gostarter"), envir=environment())
  
  if (test == T) {
    all_go <- parallel::parLapply(cl ,gene_list[1:2], function(x) gostarter(x, universe, ont = ont))
  } else {
    all_go <- parallel::parLapply(cl, gene_list, function(x) gostarter(x, universe, ont = ont))  
  }
  
  # release cores
  parallel::stopCluster(cl)
  names(all_go) <- names(gene_list)
  
  go_df <- bind_rows(all_go, .id = "time_window")
  return(go_df)
  
} # end of topGO_wrap