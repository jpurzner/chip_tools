#' Remove overlaps between cluster outlines
#'
#' Takes the polygons from [outline_tsne()] and iteratively subtracts
#' overlapping regions so no two outlines intersect, re-ordering and
#' re-trimming the boundary after each pass.
#'
#' @param tsne_outl List of per-cluster boundary data frames from
#'   [outline_tsne()].
#' @param point_dist_window Rolling window used to smooth point-to-point
#'   distances before trimming.
#' @param kern Unused; retained for call compatibility.
#' @param strip_path_dist Points whose smoothed neighbour distance exceeds this
#'   are dropped, which removes the long spurs left by polygon subtraction.
#'
#' @return A list of non-overlapping boundary data frames.
#'
#' @section Dependencies:
#' Requires `rgeos`, **archived from CRAN in October 2023**. The function
#' therefore cannot run on a current R installation without installing rgeos
#' from the archive. Porting the three `rgeos` calls
#' (`gIntersects`/`gBuffer`/`gDifference`) to `sf`
#' (`st_intersects`/`st_buffer`/`st_difference`) is the fix; `sf` is listed in
#' `Suggests` in anticipation. Also needs `prevR` for `orderPoints()`, which
#' was never declared.
#'
#' @section Fixes:
#' Point-to-point distance was computed as `sqrt(dy^2 + dy^2)` — the x
#' difference was never used, so the spur filter cut on the y axis alone.
#'
#' @seealso [outline_tsne()], [spline.poly()]
#'
#' @examples
#' \dontrun{
#' clean <- remove_intersect_all(outline_tsne(tsne_res))
#' }
#' @export
remove_intersect_all <- function(tsne_outl, point_dist_window = 2,  kern = FALSE, strip_path_dist = 1) {
  require_pkg("rgeos",
              reason = paste("rgeos was archived from CRAN in October 2023;",
                             "install it from the archive, or port this",
                             "function to sf - see NEWS.md."))
  require_pkg("sp", "raster", "prevR")

  
  
  
  rem_overlap <- function (tsne_outl) {
    poly_list <- lapply(tsne_outl,function(x)   sp::Polygon(as.matrix(x[,c(1,2)])))
    polys_list <-  lapply(c(1:length(poly_list)) ,function(n) sp::Polygons(poly_list[n], n)) 
    polysp_list <- lapply(polys_list, function(x) sp::SpatialPolygons(list(x)))
  
    # test for all relevant intersections 
    comparison <- expand.grid(comp1 = c(1:length(polysp_list)), comp2 = c(1:length(polysp_list))) 
    comparison <- comparison[!(comparison$comp1 == comparison$comp2),] 
    comparison$intersect <- mapply(function(c1, c2) rgeos::gIntersects(polysp_list[[c1]], polysp_list[[c2]]) , comparison$comp1, comparison$comp2)
    comparison <- comparison[comparison$intersect,]
  
    # iterate over the comparisons 
    for (n in c(1:nrow(comparison))) {
      p1 <- polysp_list[[comparison[n,1]]]
      p2 <- polysp_list[[comparison[n,2]]]
      p1 <- rgeos::gBuffer(p1, width=0)
      p2 <- rgeos::gBuffer(p2, width=0)
      polysp_list[[comparison[n,1]]] <- rgeos::gDifference(p1, p2)
    }

    tsne_no_int <- lapply(c(1:length(polysp_list)), function (n) {
      df <- raster::geom(polysp_list[[n]])
      df <- as.data.frame(df)
      df <- df[,c("x", "y")]
      df$k <- n
      return(df)
      })
    
    return(tsne_no_int)
  } # end rem_overlap 
  

  
  strip_long_path <- function(df) {
    # Was sqrt(dy^2 + dy^2): the x difference was never used, so point-to-point
    # distance was |dy| * sqrt(2) and the path-length filter cut on y alone.
    df$p2pd <- c(0, sqrt((diff(df$x)^2) + (diff(df$y)^2))) 
    df$kern_d <- c(rep(0,point_dist_window-1), zoo::rollmean(df$p2pd, point_dist_window, na.pad = FALSE))
    df <- df[df$kern_d  < strip_path_dist ,]
    #plot(df$kern_d)
    return(df)
  }
  
tsne_no_int <- rem_overlap(tsne_outl) 
tsne_no_int <- lapply(tsne_no_int, function(tsne) tsne[prevR::orderPoints(tsne$x, tsne$y,clockwise =  TRUE),])
tsne_no_int <- lapply(tsne_no_int,  strip_long_path)
tsne_no_int <- lapply(tsne_no_int, function(tsne) tsne[prevR::orderPoints(tsne$x, tsne$y,clockwise =  TRUE),])
  for (n in 1:5) { 
    tsne_no_int <- rem_overlap(tsne_no_int) 
    tsne_no_int <- lapply(tsne_no_int, function(tsne) tsne[prevR::orderPoints(tsne$x, tsne$y,clockwise =  TRUE),])
    tsne_no_int <- lapply(tsne_no_int,  strip_long_path)
    tsne_no_int <- lapply(tsne_no_int, function(tsne) tsne[prevR::orderPoints(tsne$x, tsne$y,clockwise =  TRUE),])
  }
  return(tsne_no_int)

}
