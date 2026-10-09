#' Time at which expression reaches its half-maximum
#'
#' For each gene, finds the interpolated time point at which the trajectory
#' crosses the midpoint between its minimum and maximum, restricted to the
#' monotonic stretch between those two extremes. This is the "t50" used
#' throughout the time-course work as a single-number summary of when a gene
#' turns on.
#'
#' @param expr_align_list Either a genes x timepoints matrix, or a list of such
#'   matrices (one per data set) as returned by [time_expr_interp()].
#'
#' @return For a matrix, a named vector of t50 per gene. For a list, a genes x
#'   data sets matrix.
#'
#' @section Note:
#' Genes whose maximum is below zero return `NA`. The result is capped at the
#' last observed (non-interpolated) time point of the data set.
#'
#' @seealso [expr_half_max_min()] for separate up and down timings,
#'   [time_expr_interp()] to build the input, [t50_to_compressed()] to map the
#'   result back onto sampled time points.
#'
#' @examples
#' \dontrun{
#' interp <- time_expr_interp(expr_list, time_ref)
#' t50 <- expr_half_max(interp)
#' }
#' @export
expr_half_max <- function(expr_align_list) {

# calculates the half max time point between the 
  
  
find_y0 <- function (j) {
  j <- j[!is.na(j)]
  if (max(j) < 0) { 
    return(NA)
  } else { 
    # discard points that are not between the max and min 
    start_j <- min(which.min(j), which.max(j)) 
    end_j <- max(which.min(j), which.max(j)) 
    j_trim <- j[c(start_j:end_j)]
    mid_exp <- (max(j) + min(j))/2
    max_i <- tryCatch(stats::approx(x = j_trim, y = c(1:length(j_trim)), xout = mid_exp)$y, error=function(e) NA)
    max_i <- max_i + start_j
    return(max_i)
  }
}

if (is.list(expr_align_list)) {
  ds_time <- lapply(expr_align_list, function (x) which(!is.na(x[1,])))  
  all_time <- length(expr_align_list[[1]][1,])
  half_max <- sapply(c(1:length(expr_align_list)), 
                     function(i) apply(expr_align_list[[i]], 1, 
                                       function(j) min(find_y0(j), max(ds_time[[i]]), na.rm = TRUE) ))  
} else { 
  ds_time <- which(!is.na(expr_align_list[1,]))
  all_time <- length(expr_align_list[1,])
  half_max <- apply(expr_align_list, 1,  function(j) min(find_y0(j), max(ds_time), na.rm = TRUE) ) 
} # end if 
return(half_max)
}