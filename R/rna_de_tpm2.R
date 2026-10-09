#' Differential expression with TPM columns (version 2)
#'
#' Runs DESeq2 on an RNA-seq count table and returns the counts, RPM, RPKM,
#' TPM, per-condition means and DESeq2 results in one wide table.
#'
#' @param count_df Genes x samples raw count data frame.
#' @param meta Condition labels, one per column of `count_df`.
#' @param de_prefix String prepended to the DESeq2 result column names.
#' @param exon_length Data frame with an `exon_length` column, row-matched to
#'   `count_df`.
#' @param read_cutoff Genes whose maximum count is below this are dropped.
#'
#' @return A wide data frame of counts, normalised values and DESeq2 results.
#'
#' @section Note:
#' Superseded by [RNA_de_tpm_3()], which computes TPM with [counts_to_tpm()]
#' (proper effective-length handling) instead of the inline formula here.
#' Kept because published figures were made with it.
#'
#' @section Fixes:
#' `count_TPM` was computed with `mapply()` over the TPM formula and then
#' immediately overwritten by a `sweep()` of the *raw* counts, so the first
#' computation was dead code and the returned TPM came from the sweep.
#'
#' @seealso [RNA_de_tpm_3()], [counts_to_tpm()], [chip_de_rpm()]
#'
#' @examples
#' \dontrun{
#' res <- RNA_de_tpm_2(counts, meta$condition, "GNP_MB", exon_lengths)
#' }
#' @export
RNA_de_tpm_2 <- function(count_df, meta, de_prefix = NULL, exon_length = NULL, read_cutoff = 5) {
  require_pkg("stringr", "DESeq2", "pheatmap")


# preparing the count_df
# 1. if using HTSeq then make sure you remove the last 5 lines
#    count_df <- count_df[c(1:(dim(count_df)[1]-5)),]

# trouble shooting: 
# 1. Error in colSums(test) : 'x' must be numeric
# 	make sure that the count_df has the identifier in row names and the remainder of 
#   of the columns are numeric 
  
  
  
  #print("troubleshooting")
  #print(meta)
  # setup metadata data.frame
  meta <- factor(meta)
  colData <- colnames(count_df)
  colData <- as.data.frame(colData)
  row.names(colData) <- colnames(count_df) 
  colData[,1] <- meta
  colnames(colData) <- "condition"
  #print(colData)
  
  # ensure exon_length matches count_df
  ## Need to handle the case where exon_length has missing rows (missing genes)
  #print(head(exon_length))
  exon_length <- as.data.frame(exon_length[match(row.names(count_df), row.names(exon_length)),])
  row.names(exon_length) <- row.names(count_df)
  colnames(exon_length)[1] <- "exon_length"
  
  grouped_col_mean <- function (counts, metadata, meta_factor) {
    factor_index <-  which(colnames(metadata)==meta_factor)
    group_index <- metadata[, factor_index] 
    #print(group_inddex)
    groups <- unique(group_index)
    #print("groups = ")
    #print(groups)
    #print("group_index =")
    #print(group_index)
    counts_avg <- as.data.frame(lapply(groups, function(x) {rowMeans(counts[,which(group_index == x)])}))
    colnames(counts_avg) <- metadata[match(groups, group_index), factor_index]
    return(counts_avg)
  }

  
count_table2rpm <- function(counts, colData) {  
  # calculate RPM and average 
  count_sums <- colSums(counts/1000000)
  count_RPM <- sweep(counts, MARGIN = 2, count_sums, "/")
  colnames(count_RPM) <- str_c(colnames(count_RPM), "_RPM")
  count_RPM_avg <- grouped_col_mean(count_RPM, colData, "condition")
  colnames(count_RPM_avg) <- str_c(colnames(count_RPM_avg), "_RPM_avg") 
  all_rpm <- merge(x = count_RPM, y = count_RPM_avg, by.x = 0 , by.y = 0, all = F)
  row.names(all_rpm) <- all_rpm$Row.names 
  all_rpm$Row.names <- NULL
  return(all_rpm)
}  
  
  
#################
#Calculate RPKM
#################  
count_table2rpkm <- function(counts, ex_length, colData) {
    
    #library(GenomicFeatures)
    #txdb=makeTranscriptDbFromUCSC(genome='mm9',tablename='ensGene')
    #ex_by_gene=exonsBy(txdb,'gene')
    
    counts <- counts[order(row.names(counts)),]
    counts <- as.matrix(counts)
    geneLengthsInKB <- as.matrix(ex_length[,1]) / 1000
    millionsMapped <- colSums(counts) / 1000000
    millionsMapped <- as.matrix(millionsMapped)
    rpm = sweep(counts, MARGIN = 2, millionsMapped, "/")
    rpkm = sweep(rpm, MARGIN = 1, geneLengthsInKB, "/")
    colnames(rpkm) <- str_c(colnames(rpkm), "_RPKM")
 	avg_rpkm <- grouped_col_mean(rpkm, colData, "condition")
 	colnames(avg_rpkm) <- str_c(colnames(avg_rpkm), "_avg_RPKM")   
    all_rpkm <- merge(x = rpkm, y = avg_rpkm, by.x = 0 , by.y = 0, all = F)
    row.names(all_rpkm) <- all_rpkm$Row.names
    all_rpkm$Row.names <- NULL
    return(all_rpkm)
  }

#################
#Calculate TPM
#################
# TPM = rg X rl x 10^6 / flg x T
counts2tpm <- function(count_table, exon_length, read_length = 101) {
    
  #calculate the total number of transcripts per library 
  Tn <- function(rg, rl, flg) {(rg * rl)/ flg } 
  Tn_sum <- colSums(as.data.frame(lapply(count_table[,1:ncol(count_table)], function(x) {mapply(Tn, x, read_length, exon_length[,"exon_length"])})))
  print(Tn_sum)
  # calculate TPM 
  TPM <- function(rg, rl, flg) {(rg *  rl * 1e6)/(flg)}
  count_TPM = sweep(count_table, MARGIN = 2, Tn_sum, "/")
  #row.names(count_TPM) <- row.names(count_table)
  colnames(count_TPM) <- str_c(colnames(count_TPM), "_TPM")
  count_TPM_avg <- grouped_col_mean(count_TPM, colData, "condition")
  colnames(count_TPM_avg) <- str_c(colnames(count_TPM_avg), "_TPM_avg")  
  all_tpm <- merge(x = count_TPM, y = count_TPM_avg, by.x = 0 , by.y = 0, all = F)
  row.names(all_tpm) <- all_tpm$Row.names
  all_tpm$Row.names <- NULL
  return(all_tpm)
 }

rpm_vals <- count_table2rpm(count_df, colData)
#print(head(rpm_vals))
all_df <- merge(x = count_df, y = rpm_vals, by.x = 0 , by.y = 0, all = F)

if (!is.null(exon_length)) { 
  rpkm_vals <- count_table2rpkm(count_df, exon_length, colData)
  #print(dim(rpkm_vals))
  #print(head(rpkm_vals))
  tpm_vals <- counts2tpm(count_df, exon_length, read_length = 101)
  print(dim(tpm_vals))
  all_df <- merge(x = count_df, y = rpkm_vals, by.x = 0 , by.y = 0, all = F)
  all_df <- merge(x = all_df, y = tpm_vals, by.x = 1 , by.y = 0, all = F)
}

 
# Run DESeq package and create a CountDataSet object
if (!is.null(de_prefix)) {
  #print(colData)
  dds <- DESeqDataSetFromMatrix(countData = count_df, colData = colData, design = ~ condition)
  dds <- estimateSizeFactors(dds)
  idx <- rowSums( counts(dds, normalized=TRUE) >= read_cutoff ) >= 3
  dds <- dds[idx,]
  dds <- DESeq(dds)
  DE_results <- as.data.frame(results(dds))
  DE_results$subset_2F <- "unchanged"
  DE_results[ !is.na(DE_results$padj) & DE_results$padj < 0.05 & DE_results$log2FoldChange < -1, 7] <- str_c(levels(meta)[1], "_high")
  DE_results[ !is.na(DE_results$padj) & DE_results$padj < 0.05 & DE_results$log2FoldChange > 1, 7] <- str_c(levels(meta)[2], "_high")
  DE_results$subset_1.5F <- "unchanged"
  DE_results[ !is.na(DE_results$padj) & DE_results$padj < 0.05 & DE_results$log2FoldChange < -log2(1.5), 8] <- str_c(levels(meta)[1], "_high")
  DE_results[ !is.na(DE_results$padj) & DE_results$padj < 0.05 & DE_results$log2FoldChange > log2(1.5), 8] <- str_c(levels(meta)[2], "_high")
  colnames(DE_results) <- str_c(de_prefix, "_", colnames(DE_results))
  all_df <- merge(x = all_df, y = DE_results, by.x = 1 , by.y = 0, all = F)
}

# create data_frame with all data
row.names(all_df) <- all_df$Row.names
all_df$Row.names <- NULL

return (all_df)
}
