#' Variance-stabilising transform, averaged by condition
#'
#' Applies the DESeq2 variance-stabilising transformation to a count table and
#' returns the per-condition means of the transformed values.
#'
#' @param data Features x samples raw count data frame. The last five rows are
#'   dropped, which is the HTSeq convention for its `__no_feature` and related
#'   summary lines.
#' @param meta Data frame with a `name` column matching the column names of
#'   `data` and one further column giving the condition.
#'
#' @return A data frame of VST values averaged within condition.
#'
#' @section Note:
#' The five trailing rows are removed unconditionally. If your count table has
#' no HTSeq summary rows you will lose five real features.
#'
#' @seealso [grouped_col_mean()], [wins_norm_histones()]
#'
#' @examples
#' \dontrun{
#' vsd <- vsd_norm(htseq_counts, meta)
#' }
#' @export
vsd_norm  <- function(data, meta) {
  require_pkg("DESeq2", "SummarizedExperiment")



  row.names(meta) <- meta$name
  meta$name <- NULL
  colnames(meta) <- "condition"
  meta$condition <- factor(meta$condition)
  #print(meta)
  #removes the last 5 lines for HTSeq counting
  data <- data[c(1:(dim(data)[1]-5)),]
  #print(head(data))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ condition)
  dds <- dds[ apply(counts(dds), 1, max) > 25, ]
  
  rld <- DESeq2::rlog(dds)
  vsd <- DESeq2::varianceStabilizingTransformation(dds)
  varstbl <- SummarizedExperiment::assay(vsd)
  group_varstbl <- grouped_col_mean(varstbl, meta, "condition")
    
  group_varstbl$GENE.ID  <- row.names(group_varstbl) 
  # move last to first column
  group_varstbl <- group_varstbl[,c((dim(group_varstbl)[2]),1:(dim(group_varstbl)[2]-1))]
  return(group_varstbl)
  
}