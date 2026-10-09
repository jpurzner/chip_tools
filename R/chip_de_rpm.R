#' Differential binding between two conditions, with RPM columns
#'
#' Runs DESeq2 on a ChIP-seq count table over peaks or genes, calls regions up
#' or down at two fold-change thresholds, and returns the raw counts, per-sample
#' RPM, per-condition mean RPM and the DESeq2 results in one wide table.
#'
#' @param count_df Features x samples raw count data frame; row names are
#'   peaks, genes or regions.
#' @param meta Character vector of condition labels, one per column of
#'   `count_df`, e.g. `c("GNP_H3K4me3", "GNP_H3K4me3", "MB_H3K4me3", "MB_H3K4me3")`.
#'   The first level is the DESeq2 reference.
#' @param de_prefix String prepended to the DESeq2 result column names, so
#'   several marks can be joined into one table.
#' @param verbose If `TRUE` print the `colData` and the call tables. The
#'   original always printed these.
#'
#' @return A data frame of `count_df`, its RPM columns, the per-condition mean
#'   RPM columns, and the prefixed DESeq2 results including `subset_2F` and
#'   `subset_1.5F` calls at 2-fold and 1.5-fold.
#'
#' @section Note:
#' `chip_de_tpm()` in the original script collection was a byte-for-byte copy
#' of this function under a different name, computing RPM despite the name. It
#' is not carried over; use this function.
#'
#' @seealso [RNA_de_tpm_3()] for the RNA-seq equivalent, [count_table_trim()]
#'   to build the inputs.
#'
#' @examples
#' \dontrun{
#' res <- chip_de_rpm(peak_counts,
#'                    c("GNP", "GNP", "MB", "MB"),
#'                    de_prefix = "H3K4me3")
#' }
#' @export
# process chip-seq data to indentify differentially bound regions 
# count_df (row.names = genes, peaks, etc, col.names = datasets)
# meta (names of categories example: c("GNP_H3K4me3", "GNP_H3K4me3", "MB_H3K4me3", "MB_H3K4me3")
# de_prefix (the text to prepend the results of DEseq)

chip_de_rpm <- function(count_df, meta, de_prefix, verbose = FALSE) {
  require_pkg("DESeq2", "stringr")


# setup metadata data.frame
meta <- factor(meta)
colData <- colnames(count_df)
colData <- as.data.frame(colData)
row.names(colData) <- colnames(count_df)
colData[,1] <- meta
colnames(colData) <- "condition"

# calculate RPM and average 
count_sums <- colSums(count_df/1000000)
count_RPM <- sweep(count_df, MARGIN = 2, count_sums, "/")
colnames(count_RPM) <- paste0(colnames(count_RPM), "_RPM")
count_RPM_avg <- grouped_col_mean(count_RPM, colData, "condition")
colnames(count_RPM_avg) <- paste0(colnames(count_RPM_avg), "_RPM_avg")

# run DEseq 
if (verbose) print(colData)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_df, colData = colData, design = ~ condition)
dds <- DESeq2::DESeq(dds)
DE_results <- as.data.frame(DESeq2::results(dds))

# determine subset for CHIP abundance change (up / down)
DE_results$subset_2F <- "unchanged"
DE_results[ !is.na(DE_results$padj) & DE_results$padj < 0.05 & DE_results$log2FoldChange < -1, 7] <- paste0(levels(meta)[1], "_high")
DE_results[ !is.na(DE_results$padj) & DE_results$padj < 0.05 & DE_results$log2FoldChange > 1, 7] <- paste0(levels(meta)[2], "_high")
DE_results$subset_1.5F <- "unchanged"
DE_results[ !is.na(DE_results$padj) & DE_results$padj < 0.05 & DE_results$log2FoldChange < -log2(1.5), 8] <- paste0(levels(meta)[1], "_high")
DE_results[ !is.na(DE_results$padj) & DE_results$padj < 0.05 & DE_results$log2FoldChange > log2(1.5), 8] <- paste0(levels(meta)[2], "_high")
if (verbose) print(table(DE_results$subset_2F))
colnames(DE_results) <- paste0(de_prefix, "_", colnames(DE_results))


# create data_frame with all data
all_df <- merge(x = count_df, y = count_RPM, by.x = 0 , by.y = 0, all = F)
all_df <- merge(x = all_df, y = count_RPM_avg, by.x = 1 , by.y = 0, all = F)
all_df <- merge(x = all_df, y = DE_results, by.x = 1 , by.y = 0, all = F)
row.names(all_df) <- all_df$Row.names
all_df$Row.names <- NULL



return(all_df)
}