#' Two-way table as a heatmap
#'
#' Works like a two-variable `table()` call but returns a heatmap, with
#' optional row normalisation, a binary presence mask, and filling of
#' combinations that never occur.
#'
#' @param x,y The two variables to cross-tabulate.
#' @param m Minimum count for a cell to be drawn.
#' @param x_title,y_title Axis titles; default to the deparsed `x` and `y`.
#' @param col_choice Fill colour ramp.
#' @param lower_col_lim Count at which the colour ramp starts.
#' @param binary If `TRUE` draw presence/absence rather than counts.
#' @param norm_y Row normalisation: `"none"`, or a mode handled by the internal
#'   `filled_table()`.
#' @param fill_table If `TRUE` include absent combinations as zero cells.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' heatmap_table(calls$H3K4me3_class, calls$H3K27me3_class)
#' }
#' @import ggplot2
#' @export
heatmap_table <- function(x, y, m = 0, x_title = NULL, y_title = NULL,  col_choice = "blue", lower_col_lim = 1, binary = F, norm_y = "none", fill_table = T) {  

  # heatmap_table works like a 2 variable table call but generates a heatmap 
  # binary flag creates a mask of the feature
  
  if (is.null(x_title)) { 
    x_title <- deparse(substitute(x))
    }
  if (is.null(y_title)) {
    y_title <- deparse(substitute(y))
  }
  
  
  filled_table <- function(x, y, m = 0, norm_y) {
    
    # function obtained from http://stackoverflow.com/questions/24519794/r-max-function-ignore-na
    # by coffe
    my.max <- function(x,m) ifelse( !all(is.na(x)), max(x,m, na.rm=T), NA)
    
    # make 2d histogram using table
    mat <- table(x, y)

    # a bunch of code to handle rows or columns with no values arg!
    bmax = max(my.max(x, m),my.max(y, m))
    print(bmax)
    
    missing_col <- c(1:bmax)[!is.element(c(1:bmax), colnames(mat))]
    missing_row <- c(1:bmax)[!is.element(c(1:bmax), row.names(mat))]
    if (length(missing_col) >= 1) {
      col_fill <- matrix(0,dim(mat)[1],length(missing_col))
      colnames(col_fill) <- missing_col
      mat <- cbind(mat, col_fill)
    } # end if 
    if (length(missing_row) >= 1) {
      row_fill <- matrix(0,length(missing_row),dim(mat)[2])
      row.names(row_fill) <- missing_row
      mat <- rbind(mat, row_fill)
    } # end if 
    
    if (norm_y == "sum") {
      y_sums <- table(y)
      mat <- sweep(mat, MARGIN = 2, y_sums, "/")
      mat <- round(mat, digits = 2)
    }
    if (norm_y == "mean") {
      y_mean = colMeans(mat)
      mat <- sweep(mat, MARGIN = 2, y_mean, "/")
      mat <- round(mat, digits = 2)      
      
    }
    
    # reorder the rows and columns
    #mat <- mat[,naturalorder(colnames(mat))]
    #mat <- mat[naturalorder(row.names(mat)),]
    return(mat)
  }  

if (fill_table) { 
  mat <- filled_table(x, y, m, norm_y)  
  
} else {
  mat <- table(x, y)
}
  
mat[mat == 0] <- NA
  
if (binary) {
  mat[!is.na(mat)] <- 1  
  lower_col_lim = 0.0001
}

print(mat)

mat_vec <- c(mat)
mat_vec <- mat_vec[order(mat_vec)]

dd <- melt(mat)
colnames(dd) <- c("Var.1", "Var.2", "value")

matp <- ggplot(dd, aes(as.factor(Var.1), as.factor(Var.2), group=Var.2)) +
  geom_tile(aes(fill = value),  colour = "black") + 
  geom_text(aes(fill = value, label = value), size = 6) +
  scale_fill_gradient(low = "white", high = col_choice, trans = "log", na.value = "white", 
                      # if you want a lengend set guide to T and uncomment the legend settings
                      limits=c(lower_col_lim , max(mat_vec)),  guide = FALSE) +
  scale_y_discrete(name = y_title) +
  scale_x_discrete(name = x_title) +
  theme(axis.text.x = element_text(size = 22, colour = "black"), 
        axis.text.y = element_text(size=22, colour = "black"),  
        axis.title = element_text(size=22, colour = "black"),
        # superseded by the three lines above; axis.text.x was being passed
        # to theme() twice, which is an error, so heatmap_table() could never
        # draw. Left here commented for reference.
        #axis.title.y = element_text(size = 14, angle = 90),
        #axis.text.y = element_text(size = 14),
        #axis.text.x = element_text(size = 14, angle = 0),
        #axis.title.x = element_text(size = 14),
        #legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal", 
        #legend.title = element_text(size=14), legend.text = element_text(size=14), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        # uncomment below to allow a border
        #panel.border = element_rect(colour = "black", fill=NA), 
        panel.background = element_blank()) +
  #plot.margin = unit(c(-0.2, 1, 1, 1), "cm")) +
  coord_fixed()
#end of ggplot2  

return(matp)
}