  plot_tsne_kmeans <- function(data_log, tsne, k, meta = NULL, labels=NULL, 
                               tsne_cluster_mode = "kmeans", cluster_order = NULL, sort_rows = FALSE, 
                               clusters = NULL, heatmap_colour = "YlOrRd") {
  
  # labels: the genes that will be txt labeled in TSNE plot
  # meta: data frame the first column being the columns in data_log and the second column 
  #       being the factor they represent
    
  # sort_rows if true will sort within clusters based upon cluster_order
  #           if false will perform hierarchical clustering on the k-means clusters 
  #           will disregard sort_rows if no cluster_order provided
  
  #TO DO handle the case where cluster order has 1,2,3,n entries only works for 2 now   
    
  require(grid)
  require(ggplot2)
  require(pheatmap)
  require(ggrepel)
  #library(grDevices)
  #library(gridExtra)
  require(gtable)  
  require(cowplot)
  require(plyr)
  require(dplyr)
  require(RColorBrewer)
  require(kknn)
  require(reshape2)
  require(mclust)
  source('~/Dropbox/jp_seq/opti_map.R')
  
  set.seed(42)
  

  
  
  # handle output from Rtnse or the python scikit learn version, which is an array
  if (is.data.frame(tsne)) {
    tsne_res <- tsne
  } else {
    tsne_res <- as.data.frame(tsne$Y)
  }
  row.names(tsne_res) <- row.names(data_log)
  colnames(tsne_res) <- c("x", "y")
  
  
  if (tsne_cluster_mode == "kmeans") {
    tsne_km <- kmeans(tsne_res, k, iter.max = 4000, nstart = 50) 
    tsne_res$tsne_k <- as.factor(tsne_km$cluster)
  } else if (tsne_cluster_mode == "spectral") {
    tsne_km <- specClust(tsne_res, centers=k, nn = 50)
    tsne_res$tsne_k <- as.factor(tsne_km$cluster)
  } else if (tsne_cluster_mode == "mclust") {
    tsne_mod <- densityMclust(tsne_res)
    tsne_res$tsne_k <- as.factor(tsne_mod$classification)
    
    
  }
    
  
  # setup the color vals for the heatmap
  breaksList <- quantile(as.matrix(data.matrix(data_log)), seq(4,96)/100)
  # set first and last value arbitrarily high and low fill extreme values
  breaksList[1] <- -20
  breaksList[length(breaksList)] <- 50
  color_vals <-  colorRampPalette(brewer.pal(n = 7, name = heatmap_colour))(length(breaksList))
  
  # k-means clustering on the multidimensional data 
  if (is.null(clusters)){ 
    data_log_k_obj <- kmeans(data_log, k, iter.max = 4000, nstart = 50)
    tsne_res$data_k <- as.factor(data_log_k_obj$cluster)
  } else {
    tsne_res$data_k <- as.factor(clusters)  
  }
  #print(head(tsne_res))  
  data_log <- as.data.frame(data_log)
  data_log$tsne_k <- tsne_res$tsne_k
  data_log$data_k <- tsne_res$data_k
  
  
  # reorder clusters based upon one or label
  # flatten the data df to add the metadata info 
  data_log_m <- data_log 
  data_log_m$gene <- row.names(data_log_m)
  data_log_m <- melt(data_log_m, id.vars = c("gene", "data_k", "tsne_k"))
  data_log_m$condition <- mapvalues(data_log_m$variable, 
                                    from=as.character(row.names(meta)),
                                    to = as.character(meta[,1]))
  data_log_m$condition <- as.factor(data_log_m$condition)
  # take mean of each metric by data_k and condition 
  data_log_m <- plyr::ddply(data_log_m, .(data_k, condition), summarise, mean_val = mean(value))
  #normalize matrix by column
  data_log_m <- dcast(data_log_m, data_k ~ condition, value.var = "mean_val")
  
  # re-order the data_log matrix based upon options 
  # cluster_order 
  # sort_rows 
  if (!is.null(cluster_order)) {
    #order clusters   
    cluster_shuffle <- c(1:k)
    names(cluster_shuffle) <- order(data_log_m[,cluster_order[1]], 1-data_log_m[,cluster_order[2]])
  
    # remap the data clusters to the new order 
    tsne_res$data_k <- mapvalues(tsne_res$data_k, from=names(cluster_shuffle), to=cluster_shuffle)
    data_log$data_k <- mapvalues(data_log$data_k, from=names(cluster_shuffle), to=cluster_shuffle)
    # re level so the order of the plots is correct 
    tsne_res$data_k <- factor(tsne_res$data_k, levels = c(1:k))
    data_log$data_k <- factor(data_log$data_k, levels = c(1:k))  
  
    if (sort_rows) {
      # re-order data_log based on data_k 
      # get the columns that the first clusters are sorted 
      data_log <- data_log[order(data_log$data_k, rowMeans(data_log[,which(meta$condition == cluster_order[1])])),]     
      
    } else { 
      # split by data_k and order by hclust 
      
      hc_order <- function(m) {
        m1 <- m[,c(2:(ncol(m)-2))]
        m1 <- as.matrix(m1)
        d <- dist(m1)
        hc <- hclust(d)
        m <- as.data.frame(m[hc$order,])
        #return(hc$order)
        return(m)
      }
      
      data_log <- data_log[order(data_log$data_k),]
      
      data_log_hc_order <- as.data.frame(data_log) %>% 
        dplyr::add_rownames() %>%
         dplyr::group_by(data_k) %>%
        #dplyr::select(-one_of(c("tsne_k"))) %>%
        #dplyr::do(head(.))  
        dplyr::do(hc_order(.))
      data_log <- as.data.frame(data_log_hc_order[,c(2:ncol(data_log_hc_order))])
      row.names(data_log) <- data_log_hc_order$rowname
      print(head(data_log))
    }
  } else {
    data_log <- data_log[order(data_log$data_k),]
  }
  
  heat_anno <- as.data.frame(data_log[,c((dim(data_log)[2]-1):(dim(data_log)[2]))])
  row.names(heat_anno) <- row.names(data_log)
  
  data2tsne_k_map <- opti_map(tsne_res$tsne_k,tsne_res$data_k)

  # remap the values tsne_k to data_k 
  tsne_res$tsne_k <- mapvalues(tsne_res$tsne_k, from=data2tsne_k_map$x, to=data2tsne_k_map$y)
  data_log$tsne_k <- mapvalues(data_log$tsne_k, from=data2tsne_k_map$x, to=data2tsne_k_map$y)
  
  tsne_res$tsne_k <- factor(tsne_res$tsne_k, levels = c(1:k))
  data_log$tsne_k <- factor(data_log$tsne_k, levels = c(1:k))
  
  if (k > 11) {
    # from http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
    #qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    #col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    #k_colors=sample(col_vector, k)
    k_colors <- colorRampPalette(brewer.pal(11, "Spectral"))(k)
  } else {
    k_colors <- brewer.pal(k,"Spectral")
  }
  names(k_colors) <- c(1:k)
  
  tsne_res$tsne_colors <- mapvalues(tsne_res$tsne_k, from = names(k_colors), to = k_colors)
  tsne_res$data_colors <- mapvalues(tsne_res$data_k, from = names(k_colors), to = k_colors)
  
  ann_colors = list(
    data_k = k_colors,
    tsne_k = k_colors
  )
  

  
  if (is.null(meta)) {
    h1 <- pheatmap(as.matrix(data_log[c(1:(dim(data_log)[2]-2))]), 
                   cluster_rows=FALSE,cluster_cols=FALSE, 
                   show_rownames = FALSE, 
                   color =  color_vals,
                   breaks = breaksList,
                   annotation_row = heat_anno,
                   annotation_colors = ann_colors,
                   fontsize_row = 12,
                   silent = TRUE)  
  } else {
    h1 <- pheatmap(as.matrix(data_log[c(1:(ncol(data_log)-2))]), 
                   cluster_rows=FALSE,cluster_cols=FALSE, 
                   show_rownames = FALSE,
                   #show_colnames = FALSE,
                   #annotation_col = meta,
                   color =  color_vals,
                   breaks = breaksList,
                   labels_col = meta$condition,
                   annotation_row = heat_anno,
                   annotation_colors = ann_colors,
                   fontsize_row = 12,
                   silent = TRUE)      
  }
  
  # generate the TSNE scatter plots
  g1 <- ggplot(tsne_res, aes(x, y, colour = data_k)) + 
    geom_point(size=0.1, alpha = 0.6) +  
    scale_colour_manual(values = k_colors) + 
    guides(colour = guide_legend(override.aes = list(size=4))) +
    theme(legend.position="bottom", 
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.line=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
  
  g2 <- ggplot(tsne_res, aes(x, y, colour = tsne_k)) + 
    geom_point(size=0.1, alpha = 0.6) + 
    scale_colour_manual(values = k_colors) + 
    guides(colour = guide_legend(override.aes = list(size=4))) +
    theme(legend.position="bottom", 
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(), 
          axis.line=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
  
  # add labels to figures if they are provided 
  if (!is.null(labels)) {
    label_only <- tsne_res[row.names(tsne_res) %in% labels,]
    label_only$name <- row.names(label_only)
    g1 <- g1 + geom_text_repel(data = label_only, aes(x=x, y=y, label = name), 
                               point.padding = unit(1, 'lines'), 
                               colour = "black") 
  }
  
  # convert tabel to grob
  g1 <- ggplotGrob(g1)
  g2 <- ggplotGrob(g2)
  
  # move the legend to the empty space between the two 
  legend1 <- gtable_filter(g1,"guide-box")
  #sp <- gtable_filter(g1, "spacer")
  g1 <- gtable_filter(g1,"panel")
  #g1 <- gtable_remove_grobs(g1, names = c("guide-box", "ylab", "xlab"), trim = TRUE)
  
  
  legend2 <- gtable_filter(g2,"guide-box")
  #sp <- gtable_filter(g1, "spacer")
  g2 <- gtable_filter(g2,"panel")
  #return(g1)
  
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
      #annotation_col = meta,
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
  
  g_mat <- matrix(list(g1,legend1,g2,legend2,h1,bg,h2,bg), ncol = 2)
  gt_mat <- gtable_matrix("demo", grobs = g_mat,  
                          widths = unit(c(2,2), "null"),
                          heights = unit(c(1,0.2,1,0.2), "null"))
  grid.newpage()
  grid.draw(gt_mat)
  invisible(tsne_res)
}
