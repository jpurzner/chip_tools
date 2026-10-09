#' Read a sample metadata table
#'
#' Reads a tab-delimited metadata file, sorts by the `order` column and drops
#' rows where `order` is zero, which is the convention used to exclude a sample
#' from an analysis without deleting its row.
#'
#' Expected columns: `name` (becomes row names), `order`, `group`, `labels`,
#' then any further factors.
#'
#' @param filename Path to the metadata file.
#'
#' @return A data frame of metadata, row names taken from the first column.
#'
#' @seealso [grouped_col_mean()], [count_table_trim()]
#'
#' @examples
#' \dontrun{
#' meta <- load_metadata("metadata.txt")
#' }
#' @export
# load and process metadata for RNAseq, 
# tab delimited format with the following headings
#	name<\t>order<\t>group<\t>labels<\t>other_factors	

load_metadata <- function (filename) {
	metadata <- utils::read.table(filename, sep="\t", header = T, row.names = 1)
	metadata <- metadata[order(metadata$order),]
	metadata <- metadata[metadata$order != 0,]
	return(metadata)
}