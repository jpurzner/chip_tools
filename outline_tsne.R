outline_tsne <- function(df, kern = TRUE, rem_dist = 0,  dec = 0 ,
                         remove_inner_window_size = 5, dbscan = TRUE, 
                         dbscan_eps = NULL, dbscan_minp = 3, area2point_delta_max = 6,
                         dbscan_k = 1, order_k = 1,
                        kernel_width = 5, 
                         debug_outline = "none") {
  # tsne coords listed as x and y columns in df 
  # debug options "outline_rot", "dbscan_parameters", "dbscan_results", "remove_inner"
  
  require(dplyr)
  require(plyr)
  require(tidyr)
  require(contoureR)
  require(KernSmooth)
  require(dbscan)
  require(zoo)
  require(sp)
  require(Biobase)
  
  
  
  rotate <- function(df, degree) {
    
    dfr <- df
    degree <- pi * degree / 180
    l <- sqrt(df$x^2 + df$y^2)
    teta <- atan(df$y / df$x)
    dfr$x <- round(l * cos(teta - degree))
    dfr$y <- round(l * sin(teta - degree))
    return(dfr)
  }
  
  outline_rot <- function (df, degree) {
    df <- df[,c("x", "y")]
    df <- as.data.frame(as.matrix(df)*(10^dec))
    df$x <- as.numeric(df$x)
    df$y <- as.numeric(df$y)
    df <- sweep(df, 2, 100*(10^(dec*2)), `+`)
    #df <- sweep(df, 2, 1000000, `+`)
    dfr <- rotate(df, degree)
    dfry <- dfr %>% 
      mutate( x_r = round(x), y_r = round(y)) %>% 
      group_by(y_r) %>% 
      dplyr::summarize(x_max = max(x_r), 
                       x_min = min(x_r)) %>% 
      tidyr::gather( key, value, -y_r) %>% 
      dplyr::select(x = value, y = y_r) 
    dfry <- as.data.frame(dfry)
    
    dfrx <- dfr %>%
      mutate( x_r = round(x), y_r = round(y)) %>% 
      group_by(x_r) %>% 
      dplyr::summarize(y_max = max(y_r), 
                       y_min = min(y_r)) %>% 
      tidyr::gather( key, value, -x_r) %>% 
      dplyr::select(x = x_r, y =value)
    dfrx <- as.data.frame(dfrx)
    dfr <- rbind(dfry, dfrx)
    
    dfr$x <- as.numeric(dfr$x)
    dfr$y <- as.numeric(dfr$y)
    dfr <- rotate(dfr, -degree)
    #dfr <- sweep(dfr, 2, 100, `-`)
    dfr <- sweep(dfr, 2, 100*(10^(dec*2)), `-`)
    dfr <- as.data.frame(as.matrix(dfr)/(10^dec))
  }
  
  outline_list <- lapply(seq(from = 0, to = 90, by = 10), function(n) outline_rot(df, n))
  
  outline <- bind_rows(outline_list)
  
  center_x  <- mean(outline$x)
  center_y  <- mean(outline$y)
  
  # determine the distance to furthest point, better then center for irregular shapped clusters 
  outline$max_dist <- rowMax(as.matrix(dist(outline[,c("x","y")])))
  outline$dist2center <- sqrt((outline$x - center_x)^2 + (outline$y - center_y)^2)
  
  #hist(outline$dist2center, breaks =100)
  
  outline_nocenter <- subset(outline, dist2center > rem_dist)
  
  outline_nocenter <- outline_nocenter[!duplicated(outline_nocenter),]
  
  #id2 <- boxplot.stats(outline_nocenter$dist2center, coef=2)
  #outline_nocenter <- outline_nocenter %>%
  #  filter(dist2center > id2$stats[1] & dist2center < id2$stats[5])
  
  outline_nocenter <- outline_nocenter[orderPoints(outline_nocenter$x, outline_nocenter$y,clockwise =  TRUE),]
  
  if (debug_outline == "outline_rot") {
    return(outline_nocenter)
  }
  
  
  # from https://chitchatr.wordpress.com/2015/01/23/calculating-the-area-of-a-convex-hull/
  chull_area <- function(df) {
    chull.coords <- chull(df)
    chull.coords <- df[chull.coords, c(1,2)]
    chull.poly <- Polygon(chull.coords, hole=F)
    chull.area <- chull.poly@area
    num.points <- dim(df)[1]
    if ("opt" %in% colnames(df)) {
      res <- c(chull.area, num.points, df$opt[1])
      res <- t(as.data.frame(res))
      res <- as.data.frame(res)
      colnames(res) <- c("hull_area", "num_points", "opts_cl")
    } else {
      res <- c(chull.area, num.points)
      res <- t(as.data.frame(res))
      res <- as.data.frame(res)
      colnames(res) <- c("hull_area", "num_points")
      
    }
    #print(res)
    return(res)
  } # end chull_area 
  
  
  
  if (dbscan) {
    point_cutoff <- 0.5
    
    dbscan_area <- function (df, dbscan_eps, dbscan_minp, k = 1, min_points = 30,  area = FALSE) {
      
      opt <- optics(df[,c(1,2)], eps = dbscan_eps, minPts = dbscan_minp)
      #plot(opt)
      opt <- extractDBSCAN(opt, eps_cl = dbscan_eps)
      df$opt <-  opt$cluster
      
      
      high_clust <- sort(table(df$opt),decreasing=TRUE)[c(1:k)]
      high_clust <- high_clust[high_clust > min_points]
      high_clust <- high_clust[!is.na(high_clust)]
      
      high_clust <- as.numeric(names(high_clust))
      #high_clust <- as.numeric(names(sort(table(df$opt),decreasing=TRUE)[c(1:k)]))
      
      #print(high_clust)
      df <- df %>% 
        dplyr::filter(opt %in% high_clust)
    
      
      df$opt <- as.factor(df$opt) 
      df$opt <- mapvalues(df$opt, from = levels(df$opt), to = c(1:length(levels(df$opt))))
      df <- as.data.frame(df)
      print(head(df))
      
      if (area) {
        area_k <- df %>%
          dplyr::group_by(opt) %>%
          dplyr::do(chull_area(.))
        return(area_k)
      } else { # else (area)
        return(df)
      }  # end (area)
    } # end of dbscan_area function  
    
    
    if (is.null(dbscan_eps) | (debug_outline == "dbscan_parameters")){
      # determine optimal eps setting 
      area_by_eps <- lapply(seq(from = 1.5, to = 20, by = 0.5), 
                            function (eps) dbscan_area(outline_nocenter, eps, dbscan_minp, k = dbscan_k,  area = TRUE))
      #return(area_by_eps)
      all_points_area <-  chull_area(outline_nocenter) 
      
      area_by_eps <- lapply(area_by_eps, as.data.frame)
      # this line merges all clusters
      area_by_eps <- lapply(area_by_eps, function(x) colSums(x[,c("hull_area", "num_points")]))
    
      names(area_by_eps) <- seq(from = 1.5, to = 20, by = 0.5) 
      area_by_eps <- bind_rows(area_by_eps)
      #return(area_by_eps)
      area_by_eps <- as.data.frame(t(area_by_eps))
      colnames(area_by_eps) <- c("hull_area", "num_points")
      area_by_eps$eps <- as.numeric(row.names(area_by_eps))
    
      area_by_eps$num_points <- area_by_eps$num_points / dim(outline_nocenter)[1]
      
      area_by_eps$hull_area <- area_by_eps$hull_area / all_points_area$hull_area
      area_by_eps<- area_by_eps[area_by_eps$num_points > point_cutoff,]
    
      area_by_eps$area_points_d <-  c(0,diff(area_by_eps$hull_area))/c(0,diff(area_by_eps$num_points))
      
      eps_too_high <- min(which(!area_by_eps$area_points_d < area2point_delta_max), dim(area_by_eps)[1])
      
      area_by_eps_trim <- area_by_eps[c(1:max((eps_too_high-1),2)),]
      area_by_eps_trim <- area_by_eps_trim[!is.na(area_by_eps_trim$area_points_d),]
      
      
      optimal_eps <- area_by_eps_trim[dim(area_by_eps_trim)[1],"eps"]
      area_by_eps$optimal_eps <- optimal_eps
      
      outline_nocenter <- dbscan_area(outline_nocenter, optimal_eps, dbscan_minp, k = dbscan_k)
      if (debug_outline == "dbscan_parameters") {
        return(area_by_eps)
      } # end  (debug_dbscan)
    } else { # (is.null(dbscan_eps) | debug_dbscan)
    
    outline_nocenter <- dbscan_area(outline_nocenter, dbscan_eps, dbscan_minp, k = dbscan_k)
    } #end  (is.null(dbscan_eps) | debug_dbscan)
  } # end if(dbscan)
  
  if (debug_outline == "dbscan_results") {
    return(outline_nocenter)
  }
  
  
  #n <- dim(outline_nocenter)[1]
  ##outline_nocenter$pdist <- 0
  #outline_nocenter$pdist[2:n] <- sqrt((outline_nocenter$x[2:n] - outline_nocenter$x[1:n-1]) ^ 2 + (outline_nocenter$y[2:n] - outline_nocenter$y[1:n-1]) ^ 2)
  
  #outline_nocenter$m <- 0
  #outline_nocenter$m[2:n] <- (outline_nocenter$y[2:n] - outline_nocenter$y[1:n-1]) /  (outline_nocenter$x[2:n] - outline_nocenter$x[1:n-1]) 
  #outline_nocenter$dm1 <- 0
  #outline_nocenter$d1m[2:n] <- diff(outline_nocenter$m)
  
  
  #hist(outline_nocenter$pdist, breaks = 40)
  
  
  # slidding window that keeps points furthest from center 
  remove_inner <- function(df,window_size = remove_inner_window_size, filter_knn = TRUE, knn_mean = 2.5, max_dist = FALSE) {
    # remove points that are not crowded 
    if (filter_knn) {
      nn <- kNN(df[,c(1,2)], window_size)
      which(apply(nn$dist, 1, max) < 3)
      keep_far <- which(apply(nn$dist, 1, mean) > knn_mean)
    }
    # must have column dist2center 
    if (max_dist) {
      max_window <- rollapply(df$max_dist, width = window_size, by = 1, align = "left",  function(x) which.max(x))
    } else {
      max_window <- rollapply(df$dist2center, width = window_size, by = 1, align = "left",  function(x) which.max(x))
    }
    ind_list <- lapply(c(1:(dim(df)[1]-window_size)),  function (n) seq(from = n, to= n+(window_size-1), by = 1))
    max_ind <- sapply(c(1:(length(max_window)-1)), function (n) ind_list[[n]][max_window[n]])
    max_ind <- max_ind[!duplicated(max_ind)]
    df <- df[unique(c(max_ind,keep_far)),]
    return(df)
  }
  
  order_clockwise <- function(outline_nocenter, centers = 1) {
    if (centers == 1) {
      outline_nocenter <- outline_nocenter[orderPoints(outline_nocenter$x, outline_nocenter$y,clockwise =  TRUE),]
      return(outline_nocenter)
    } else {
      k_obj <- kmeans(outline_nocenter, centers = centers)
      outline_nocenter$km <- k_obj$cluster 
      outline_nocenter <- split(outline_nocenter, outline_nocenter$km) 
      outline_nocenter <- mapply(function(outline, kc) {
        df <- outline[orderPoints(outline$x, outline$y,clockwise =  TRUE, xm = kc[1], ym = kc[2] ),]
        return(df)
      }, outline_nocenter,split(k_obj$centers, seq(nrow(k_obj$centers))), SIMPLIFY = TRUE)
      
      outline_nocenter <- lapply(c(1, dim(outline_nocenter)[2]), function(n)  as.data.frame(outline_nocenter[,n]))
      
      #poly_list <- lapply(outline_nocenter,function(x)   Polygon(as.matrix(x[,c(1,2)])))
      #polys_list <-  lapply(c(1:length(poly_list)) ,function(n) Polygons(poly_list[n], n)) 
      #polysp_list <- lapply(polys_list, function(x) SpatialPolygons(list(x)))
      
      #poly_merge <- aggregate(rbind(polysp_list[[1]], polysp_list[[2]]))
      
      #poly_merge <- gUnion(polysp_list[[1]], polysp_list[[2], drop_lower_td = TRUE)
      #poly_merge <- gUnaryUnion(poly_merge)
      
      
      #poly_merge_df <- raster::geom(poly_merge)
      #poly_merge_df <- as.data.frame(poly_merge_df)
      #align_check <- expand.grid(c(1:centers),c(1:centexrs))
      #align_check <- align_check[ align_check$Var1 != align_check$Var2, ]
      #dist_align <-lapply(split(align_check, seq(nrow(align_check))),  
      #                          function(n) sqrt(((outline_nocenter[[n[[1]]]][1,"x"] - 
      #                            outline_nocenter[[n[[1]]]][nrow(outline_nocenter[[n[[1]]]]),"x"])^2) + 
      #                      ((outline_nocenter[[n[[1]]]][1,"y"] - 
      #                         outline_nocenter[[n[[1]]]][nrow(outline_nocenter[[n[[1]]]]),"y"])^2)))
      #align_check$dist <- unlist(dist_align)                     
      #align_check <- align_check[order(align_check$dist),c(1,2)]
      #align_order <- unique(unlist(align_check))
      
      outline_nocenter <- bind_rows(outline_nocenter)
      which.max(sqrt(diff(outline_nocenter$y)^2 + diff(outline_nocenter$y)^2))
      #print(head(outline_nocenter))
      return(outline_nocenter)
      
    }
    
  }
  
  if (dbscan_k > 1) {
    outline_nocenter <- split(outline_nocenter, outline_nocenter$opt)
    outline_nocenter <- lapply(outline_nocenter, function(outline) order_clockwise(outline, centers = order_k) )
    outline_nocenter <- lapply(outline_nocenter, function(outline) remove_inner(outline, window_size = remove_inner_window_size))
  } else { 
    
    outline_nocenter <-  order_clockwise(outline_nocenter, centers = order_k)
    outline_nocenter <- remove_inner(outline_nocenter, window_size = remove_inner_window_size)
  
  }
  if (debug_outline == "remove_inner") {
    return(outline_nocenter)
  }
  
  
  if (kern) {
    kern_outline <- function(outline_nocenter) { 
      outline_nocenter <-  order_clockwise(outline_nocenter, centers = order_k)
      fit_y <- locpoly(c(1:length(outline_nocenter$x)), outline_nocenter$y, bandwidth = kernel_width)
      fit_x <- locpoly(c(1:length(outline_nocenter$x)), outline_nocenter$x, bandwidth = kernel_width)
      outline_kern <- data.frame(x = fit_x$y, y =  fit_y$y, opt = unique(outline_nocenter$opt))
      outline_kern <- outline_kern[orderPoints(outline_kern$x, outline_kern$y,clockwise =  TRUE),]
      return(outline_kern)
    }
    if (dbscan_k > 1) {
      outline_kern <- lapply(outline_nocenter, function(outline) kern_outline(outline))
      outline_kern <- bind_rows(outline_kern)
    } else {
      outline_kern <- kern_outline(outline_nocenter)
      
    }
    return(outline_kern)
  } else {
    if (dbscan_k > 1) {
      outline_nocenter <- bind_rows(outline_nocenter)
      return(outline_nocenter)
    } else {
      return(outline_nocenter)  
    }
  }
  
  
  
  
  
  
}