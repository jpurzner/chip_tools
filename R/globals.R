# Column names used with non-standard evaluation (dplyr verbs, ggplot2 aes)
# and names resolved at call time from Suggests packages. R CMD check cannot
# see these, so they are declared here rather than left as 60-odd NOTEs.
#
# Rewriting the plotting functions to use `.data[[ ]]` and `{{ }}` would make
# most of this list unnecessary; see "Worth doing next" in REVIEW.md.
utils::globalVariables(c(
  # dplyr / tidyr verbs and the pipe, used bare in the older functions
  "%>%", ".", "across", "all_of", "arrange", "as_tibble", "bind_rows",
  "case_when", "desc", "distinct", "ends_with", "everything", "filter",
  "group_by", "if_else", "left_join", "mutate", "pivot_longer", "pull",
  "rowwise", "slice", "summarise", "ungroup",

  # plyr / reshape2
  "ddply", "dcast", "melt", "mapvalues", "membership",

  # column names referenced inside dplyr and ggplot2 expressions
  "CrossingRank", "CrossingSmoothedValue", "Dataset", "GO.ID", "Ontology",
  "Rank", "Smoothed", "SmoothedDeriv1", "Term", "Var.1", "Var.2", "avg_t50",
  "cell_line", "chr", "cluster_id", "condition", "cutoff_50", "data_colors",
  "data_k", "diff_gene", "dist2center", "ds", "ds1", "gene", "gene_id",
  "genes", "group1", "group2", "initial_k", "key", "mean_avg_t50",
  "mean_lower", "mean_upper", "median_avg_t50_cor", "mgi_symbol", "p.adj",
  "p adj", "p_value", "pair", "position", "ptime", "quant",
  "replicate_cat", "replicate_names", "t50_cat", "term_combined", "tsne_k",
  "name", "tukey", "value", "variable", "x", "x_r", "y", "y.position",
  "y_r",

  # resolved at call time from Suggests packages
  "GOTERM", "GenTable", "Polygon", "V", "annFUN.org", "clusterEvalQ",
  "clusterExport", "col2hex", "densityMclust", "detectCores", "exprs",
  "extractDBSCAN", "genesInTerm", "geom_text_repel", "godata", "gpar",
  "gtable_add_cols", "gtable_add_grob",
  "gtable_filter", "gtable_matrix", "gtable_remove_grobs", "kNN", "locpoly",
  "makeCluster", "mgeneSim", "optics", "parLapply", "pheatmap",
  "plot_layout", "rollapply", "rowMax", "rowMin", "runTest",
  "specClust", "stat_pvalue_manual", "stopCluster", "str_c",
  "DESeq", "DESeqDataSetFromMatrix", "counts", "estimateSizeFactors",
  "results", "brewer.pal", "brewer.pal.info"
))
