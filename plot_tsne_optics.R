plot_tsne_optics <- function(data_log, tsne, k, meta = NULL) {
  #library(RColorBrewer)
  library(grid)
  library(ggplot2)
  library(pheatmap)
  #library(grDevices)
  #library(gridExtra)
  require(gtable)  
  require(cowplot)
  require(plyr)
  library(RColorBrewer)
  
  tsne_res <- as.data.frame(tsne$Y)
  row.names(tsne_res) <- row.names(data_log)
  colnames(tsne_res) <- c("x", "y")
  
  tsne_km <- kmeans(tsne_res, k, iter.max = 4000, nstart = 50) 
  data_log_k_obj <- kmeans(data_log, k, iter.max = 4000, nstart = 50)
  tsne_res$tsne_k <- as.factor(tsne_km$cluster)
  tsne_res$data_k <- as.factor(data_log_k_obj$cluster)
  
  data_log <- as.data.frame(data_log)
  data_log$tsne_k <- tsne_res$tsne_k
  data_log$data_k <- tsne_res$data_k
  data_log <- data_log[order(data_log$data_k),]
  heat_anno <- as.data.frame(data_log[,c((dim(data_log)[2]-1):(dim(data_log)[2]))])

  
  data2tsne_k_map <- opti_map(tsne_res$tsne_k,tsne_res$data_k)

  # remap the values tsne_k to data_k 
  tsne_res$tsne_k <- mapvalues(tsne_res$tsne_k, from=data2tsne_k_map$x, to=data2tsne_k_map$y)
  data_log$tsne_k <- mapvalues(data_log$tsne_k, from=data2tsne_k_map$x, to=data2tsne_k_map$y)
  
  tsne_res$tsne_k <- factor(tsne_res$tsne_k, levels = c(1:k))
  data_log$tsne_k <- factor(data_log$tsne_k, levels = c(1:k))
  
  if (k > 10) {
    # from http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    k_colors=sample(col_vector, k)
  } else {
    k_colors <- brewer.pal(k,"Paired")
  }
  names(k_colors) <- c(1:k)
  
  tsne_res$tsne_colors <- mapvalues(tsne_res$tsne_k, from = names(k_colors), to = k_colors)
  
  ann_colors = list(
    data_k = k_colors,
    tsne_k = k_colors
  )
  
  breaksList <- seq(0, 10, by = 0.1)
  # set first and last value arbitrarily high and low fill extreme values
  breaksList[1] <- -20
  breaksList[length(breaksList)] <- 50
  
  color_vals <-  colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  
  if (is.null(meta)) {
    h1 <- pheatmap(data_log[c(1:(dim(data_log)[2]-2))], 
                   cluster_rows=FALSE,cluster_cols=FALSE, 
                   show_rownames = FALSE, 
                   color =  color_vals,
                   breaks = breaksList,
                   annotation_row = heat_anno,
                   annotation_colors = ann_colors,
                   fontsize_row = 12,
                   silent = TRUE)  
  } else {
    h1 <- pheatmap(data_log[c(1:(dim(data_log)[2]-2))], 
                   cluster_rows=FALSE,cluster_cols=FALSE, 
                   show_rownames = FALSE,
                   #show_colnames = FALSE,
                   annotation_col = meta,
                   color =  color_vals,
                   breaks = breaksList,
                   labels_col = meta$condition,
                   annotation_row = heat_anno,
                   annotation_colors = ann_colors,
                   fontsize_row = 12,
                   silent = TRUE)      
  }
  
  g1 <- ggplot(tsne_res, aes(x, y, colour = data_k)) + 
    geom_point(size=0.1, alpha = 0.6) +  
    scale_colour_manual(values = k_colors)
  g2 <- ggplot(tsne_res, aes(x, y, colour = tsne_k)) + 
    geom_point(size=0.1, alpha = 0.6) + 
    scale_colour_manual(values = k_colors)

  g1 <- ggplotGrob(g1)
  g2 <- ggplotGrob(g2)
  
  data_log <- data_log[order(data_log$tsne_k),]
  heat_anno <- as.data.frame(data_log[,c((dim(data_log)[2]-1):(dim(data_log)[2]))])
  
  if (is.null(meta)) {
    h2 <- pheatmap(data_log[c(1:(dim(data_log)[2]-2))], 
                 cluster_rows=FALSE,cluster_cols=TRUE, 
                 show_rownames = FALSE,
                 color =  color_vals,
                 breaks = breaksList,
                 annotation_row = heat_anno,
                 annotation_colors = ann_colors,
                 fontsize_row = 12,
                 silent = TRUE)
  } else {
    h2 <- pheatmap(data_log[c(1:(dim(data_log)[2]-2))], 
      cluster_rows=FALSE,cluster_cols=FALSE, 
      show_rownames = FALSE,
      color =  color_vals,
      breaks = breaksList,
      #show_colnames = FALSE,
      annotation_col = meta,
      labels_col = meta$condition,
      annotation_row = heat_anno,
      annotation_colors = ann_colors,
      fontsize_row = 12,
      silent = TRUE)
    
    
    
  }
  h1 <- h1$gtable
  h2 <- h2$gtable

  # this needs to be adjust the following grob indexes 
  h1 <- gtable_remove_grobs(h1, names = c( "annotation_legend", "legend", "col_annotation_names"))
  h2 <- gtable_remove_grobs(h2, names = c( "annotation_legend", "legend", "col_annotation_names"))
  
  #print(h2)
  
  bg <- rectGrob(gp = gpar(alpha = 0))
  
  g_mat <- matrix(list(g1,bg,g2,bg,h1,bg,h2,bg), ncol = 2)
  gt_mat <- gtable_matrix("demo", grobs = g_mat,  
                          widths = unit(c(2,2), "null"),
                          heights = unit(c(1,0.2,1,0.2), "null"))
  grid.newpage()
  grid.draw(gt_mat)
  invisible(tsne_res)
}
