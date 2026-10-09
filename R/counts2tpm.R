#' Convert a count table to TPM
#'
#' Transcripts per million from raw counts and exon lengths, using a fixed read
#' length rather than a per-library mean fragment length.
#'
#' \deqn{T = \sum_g \frac{r_g \cdot r_l}{fl_g}, \quad
#'        TPM_g = \frac{r_g \cdot r_l \cdot 10^6}{fl_g \cdot T}}
#'
#' where \eqn{r_g} is reads on feature \eqn{g}, \eqn{r_l} the read length and
#' \eqn{fl_g} the feature length.
#'
#' @param count_table Features x samples data frame of raw counts, row names
#'   being feature IDs.
#' @param exon_length Data frame or matrix whose first column holds the exonic
#'   length of each feature, in the same row order as `count_table`.
#' @param read_length Read length in bases.
#'
#' @return A features x samples data frame of TPM.
#'
#' @seealso [counts_to_tpm()], which takes a per-library mean fragment length
#'   and excludes features shorter than it; [count_table2rpkm()] for RPKM.
#'
#' @examples
#' \dontrun{
#' tpm <- counts2tpm(counts, exon_lengths, read_length = 101)
#' }
#' @export
# takes count data table and produces a table with TPM values
# need exon length array and what cols of count table to avg
# can only do one average at a time

# input: count_table 
#	row.names = gene names 

# T = total number of transcripts per sequencing library
# T = sum (rg x rl / flg) 
# rg = reads mapped to specific feature 
# flg = feature length 
# rl = read length 

# TPM = rg X rl x 10^6 / flg x T

counts2tpm <- function(count_table, exon_length, read_length = 101) {

# calculate the total number of transcripts per library 
Tn <- function(rg, rl, flg) {(rg * rl)/ flg } 
Tn_sum <- colSums(as.data.frame(lapply(count_table[,1:ncol(count_table)], function(x) {mapply(Tn, x, read_length, exon_length[,1])})))

# calculate TPM 
TPM <- function(rg, rl, flg) {(rg *  rl * 1e6)/(flg)}
count_TPM <- as.data.frame(lapply(count_table[,1:ncol(count_table)], function(x) {mapply(TPM, x, read_length, exon_length[,1])}))
count_TPM = sweep(count_TPM, MARGIN = 2, Tn_sum, "/")
row.names(count_TPM) <- row.names(count_table)



return(count_TPM)

}