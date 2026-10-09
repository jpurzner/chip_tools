#' Merge several tab-delimited count files on their row names
#'
#' Reads each file with the first column as row names and inner-joins them in
#' sequence, dropping any column name already seen. Intended for the
#' per-sample count files that `htseq-count` and friends leave in a directory.
#'
#' @param files Character vector of file paths, e.g. `dir(pattern = "counts")`.
#' @param verbose If `TRUE` print each file as it is read. The original always
#'   printed, with no way to quiet it.
#'
#' @return A data frame with the intersection of the row names across files.
#'
#' @section Note:
#' Rows are inner-joined (`all = FALSE`), so a feature missing from any one
#' file is dropped from the result entirely.
#'
#' @seealso [load_metadata()]
#'
#' @examples
#' \dontrun{
#' counts <- load_merge(dir(pattern = "[.]counts$"))
#' }
#' @export
load_merge <- function(files, verbose = FALSE) {

# if your files don't have row names then you should use plyr 
# all <- load_merge(dir()) to get all files in a directory

if (verbose) print(files)
#file_id <- file(filename, "r")
#files <- readLines(file_id)

`%ni%` <- Negate(`%in%`)

#print(files[1])
all_merge <- utils::read.table(files[1], header = T, row.names = 1, sep = "\t")
for (f in files[2:length(files)]) {
	if (verbose) print(f)
	# read next df 
	temp <- utils::read.table(f, header = T, row.names = 1, sep = "\t")
	# prune duplicated column names 
	temp <- subset(temp ,select = names(temp) %ni% intersect(colnames(all_merge), colnames(temp)))	
	# add next df to the common df 
	all_merge <- merge(x = all_merge, y = temp, by.x = 0, by.y = 0, all = F)
	row.names(all_merge) <- all_merge$Row.names
	all_merge$Row.names <- NULL		
}	

#close(file_id)	
return (all_merge)	
}