#' GO enrichment for each window of a time-course gene matrix
#'
#' Runs a topGO enrichment on each row of a matrix whose rows are time windows
#' and whose cells are gene names, as produced by [t50window()], one window per
#' core.
#'
#' @param t50_df Matrix or data frame whose rows are windows of gene names.
#' @param all_genes The gene universe to test against.
#' @param ont GO ontology: `"BP"`, `"MF"` or `"CC"`.
#' @param test If `TRUE` process only the first two windows.
#'
#' @return A list of topGO result tables, one per window.
#'
#' @section Fixes:
#' The parallel workers loaded `GOstats` while the worker function calls
#' topGO, so every worker failed with a missing-function error.
#'
#' @seealso [t50window()] to build the input, [topGO_wrap()], [topGO_split()]
#'
#' @examples
#' \dontrun{
#' go <- topGO_timeseries(t50window(t50_table)[[1]], universe)
#' }
#' @export
topGO_timeseries <- function(t50_df, all_genes, ont = "BP",  test = FALSE) {
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
  
  gene_list <- split(t(t50_df), seq(nrow(t(t50_df))))
  
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
  names(all_go) <- colnames(t50_df)
  
  go_df <- bind_rows(all_go, .id = "time_window")
  return(go_df)
  
} # end of topGO_wrap