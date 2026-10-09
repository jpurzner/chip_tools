#' Subset a count table and its metadata to two groups
#'
#' Pulls the samples belonging to `labels` out of a count table and builds the
#' matching `colData` frame, with the factor levels in the order given. Doing
#' the comparison on a trimmed table rather than with DESeq2 contrasts keeps
#' the size-factor estimation from being influenced by libraries that are not
#' part of the comparison.
#'
#' @param labels Length-two character vector of group names, in the order you
#'   want them as factor levels (the first is the reference).
#' @param metadata Data frame with one row per sample, row names matching the
#'   column names of `counts`.
#' @param counts Features x samples count data frame.
#' @param labels_col Index of the grouping column in `metadata`.
#'
#' @return A list of `metadata` (one factor column, ready as DESeq2 `colData`)
#'   and `counts` (the matching columns).
#'
#' @seealso [load_metadata()], [chip_de_rpm()]
#'
#' @examples
#' \dontrun{
#' trimmed <- count_table_trim(c("GNP", "MB"), meta, counts)
#' chip_de_rpm(trimmed$counts, trimmed$metadata$condition, "GNP_vs_MB")
#' }
#' @export
# splits off metadata and counts from two df and creates a 
# list of 2 split dataframes, which are formated for DEseq 
# can also do this with contrasts or using factors but with 
# contrasts I think that the normalization is effected by other 
# libraries, which may not be ideal.

# labels  c("group1", "group2")
# metadata df 
# counts df 


count_table_trim <- function (labels, metadata, counts, labels_col = 3) {

	meta <- metadata[metadata[,labels_col] %in% labels,]
	meta_trim <- as.data.frame(factor(meta[,labels_col], levels = labels))
	row.names(meta_trim) <- row.names(meta)
	colnames(meta_trim)[1] <- colnames(meta)[labels_col] 
	counts_trim <- counts[,colnames(counts) %in% row.names(meta_trim)]
	out_list <- list(metadata = meta_trim, counts = counts_trim)
	return(out_list)
} # end of furntion seq_list