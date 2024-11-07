remove_intersect_all <- function(tsne_outl, point_dist_window = 2,  kern = FALSE, strip_path_dist = 1) {
  
  require(sp)
  require(rgeos)
  require(raster)
  require(zoo)
  
  
  rem_overlap <- function (tsne_outl) {
    poly_list <- lapply(tsne_outl,function(x)   Polygon(as.matrix(x[,c(1,2)])))
    polys_list <-  lapply(c(1:length(poly_list)) ,function(n) Polygons(poly_list[n], n)) 
    polysp_list <- lapply(polys_list, function(x) SpatialPolygons(list(x)))
  
    # test for all relevant intersections 
    comparison <- expand.grid(comp1 = c(1:length(polysp_list)), comp2 = c(1:length(polysp_list))) 
    comparison <- comparison[!(comparison$comp1 == comparison$comp2),] 
    comparison$intersect <- mapply(function(c1, c2) gIntersects(polysp_list[[c1]], polysp_list[[c2]]) , comparison$comp1, comparison$comp2)
    comparison <- comparison[comparison$intersect,]
  
    # iterate over the comparisons 
    for (n in c(1:nrow(comparison))) {
      p1 <- polysp_list[[comparison[n,1]]]
      p2 <- polysp_list[[comparison[n,2]]]
      p1 <- gBuffer(p1, width=0)
      p2 <- gBuffer(p2, width=0)
      polysp_list[[comparison[n,1]]] <- gDifference(p1, p2)
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
    df$p2pd <- c(0,sqrt((diff(df$y)^2) + (diff(df$y)^2))) 
    df$kern_d <- c(rep(0,point_dist_window-1), rollmean(df$p2pd, point_dist_window, na.pad = FALSE))
    df <- df[df$kern_d  < strip_path_dist ,]
    #plot(df$kern_d)
    return(df)
  }
  
tsne_no_int <- rem_overlap(tsne_outl) 
tsne_no_int <- lapply(tsne_no_int, function(tsne) tsne[orderPoints(tsne$x, tsne$y,clockwise =  TRUE),])
tsne_no_int <- lapply(tsne_no_int,  strip_long_path)
tsne_no_int <- lapply(tsne_no_int, function(tsne) tsne[orderPoints(tsne$x, tsne$y,clockwise =  TRUE),])
  for (n in 1:5) { 
    tsne_no_int <- rem_overlap(tsne_no_int) 
    tsne_no_int <- lapply(tsne_no_int, function(tsne) tsne[orderPoints(tsne$x, tsne$y,clockwise =  TRUE),])
    tsne_no_int <- lapply(tsne_no_int,  strip_long_path)
    tsne_no_int <- lapply(tsne_no_int, function(tsne) tsne[orderPoints(tsne$x, tsne$y,clockwise =  TRUE),])
  }
  return(tsne_no_int)

}
