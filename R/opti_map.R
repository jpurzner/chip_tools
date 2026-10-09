#' Greedily match two cluster labellings
#'
#' Two clusterings of the same objects (say k-means on the expression matrix
#' and k-means on its t-SNE embedding) use arbitrary integer labels, so cluster
#' 3 in one is unrelated to cluster 3 in the other. This greedily pairs them:
#' repeatedly take the largest remaining group in `x`, match it to whichever
#' `y` group it most overlaps, then drop both from consideration.
#'
#' Used by [plot_tsne_kmeans()] to keep the heatmap and the t-SNE panel
#' coloured consistently.
#'
#' @param x,y Two cluster labellings of the same objects, in the same order.
#'
#' @return A two-column data frame mapping each `x` label to a `y` label.
#'
#' @section Note:
#' This is a greedy heuristic, not an optimal assignment. For a guaranteed
#' optimum use the Hungarian algorithm, e.g. `clue::solve_LSAP()`.
#'
#' @examples
#' opti_map(c(1, 1, 2, 2, 3, 3), c(2, 2, 3, 3, 1, 1))
#' @export


opti_map <- function(x,y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  comp <- table(x, y)
  res <- data.frame(x = c(1:length(unique(x))), y = c(1:length(unique(x))))
  #print(comp)
  #print(res)
  for (i in 1:(length(unique(x))-1)) {
    #print(rowSums(comp))
    x_max <- which.max(rowSums(comp))
    #print(x_max)
    x_val <- as.numeric(names(x_max))
    #print(x_val)
    
    x_ind <- which(as.numeric(row.names(comp)) == x_val)
    y_vals <- comp[x_ind,]
    #print(length(y_vals))
    y_ind <- which.max(y_vals)
    y_max <- as.numeric(colnames(comp)[y_ind])
    #print(y_vals)
    #print(y_max)
    
    #print(res)  
    #print(comp)
    
    if ((dim(comp)[1]) > 2) {
      # remove the row and col from table
      res[i,1] <- x_val
      res[i,2] <- y_max
      comp <- comp[-x_max,]
      comp <- comp[, -y_ind]  
      } else {
      res[i,1] <- row.names(comp)[1]
      res[i,2] <- colnames(comp)[1]
      res[i+1,1] <- row.names(comp)[2]
      res[i+1,2] <- colnames(comp)[2]
    }
  }
  return(res)
}