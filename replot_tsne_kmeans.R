  replot_tsne_kmeans <- function(data_log, tsne_res, merge_col, meta = NULL, average_reps = FALSE,   
                                cluster_order = NULL, sort_rows = FALSE, add_anno = NULL, 
                              heatmap_colour = "YlOrRd", plot_what ="heatmap", mask_low = NULL, 
                              label_cols = NULL, label_factors = NULL,  split_horizontal_label= FALSE, label_colours = NULL, gene_labels = NULL, 
                              row_gaps =  NULL, col_gaps = NULL) {
  
  # replots the figures from plot_tsne_kmeans, but lets you get individual pieces for further 
  # modificaiton 
  #
  # pieces are seletec with plot_what 
  #  
  # plot_what: "data_k_tsne", "tsne_k_tsne"  
    
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
  require(gplots)  
  source('~/Dropbox/jp_seq/opti_map.R')
  
  set.seed(42)
  
  
  # setup the color vals for the heatmap
  breaksList <- quantile(as.matrix(data.matrix(data_log)), seq(4,96)/100)
  # set first and last value arbitrarily high and low fill extreme values
  breaksList[1] <- -20
  breaksList[length(breaksList)] <- 50
  color_vals <-  colorRampPalette(brewer.pal(n = 9, name = heatmap_colour))(length(breaksList))
  

  #reord_histone <- match(row.names(data_log), tsne_res[,merge_col])
  reord_histone <- match(tsne_res[,merge_col], row.names(data_log) )
  reord_histone <- reord_histone[!is.na(reord_histone)]
  data_log <- data_log[reord_histone,]
  
  data_log <- as.data.frame(data_log)
  data_log$data_k <- tsne_res$data_k
  
  # to split horizontally the data 
  if (split_horizontal_label) {
    old_k <- tsne_res$data_k
    inital_k <- max(as.numeric(tsne_res$data_k))
    new_k <- inital_k * tsne_res[,label_cols ] %in% label_factors 
    tsne_res$data_k <-tsne_res$data_k +  new_k
  }
  
  
  k <- max(as.numeric(tsne_res$data_k))
  
  # reorder clusters based upon one or label
  # flatten the data df to add the metadata info 
  data_log_m <- data_log 
  data_log_m$gene <- row.names(data_log_m)
  data_log_m <- melt(data_log_m, id.vars = c("gene", "data_k"))
  data_log_m$condition <- mapvalues(data_log_m$variable, 
                                    from=as.character(row.names(meta)),
                                    to = as.character(meta[,1]))
  data_log_m$condition <- as.factor(data_log_m$condition)
  
  if (average_reps) {
    # average the replicates 
    rep_mean <- plyr::ddply(data_log_m, .(gene, condition), summarise, mean_val = mean(value))
    rep_mean <- dcast(rep_mean, gene ~ condition, value.var = "mean_val")
    row.names(rep_mean) <- rep_mean$gene
    rep_mean$gene <- NULL
    data_log <- rep_mean 
    reord_histone <- match(tsne_res[,merge_col], row.names(data_log) )
    reord_histone <- reord_histone[!is.na(reord_histone)]
    data_log <- data_log[reord_histone,]
  
    data_log <- as.data.frame(data_log)
    data_log$data_k <- tsne_res$data_k
    
    # set meta to null otherwise colnames broken
    meta <- NULL
  }
  
  

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
        if (dim(m)[1] > 3) {
          m1 <- m[,c(2:(ncol(m)-1))]
          m1 <- as.matrix(m1)
          d <- dist(m1)
          hc <- hclust(d)
          m <- as.data.frame(m[hc$order,])
        } # end if   
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
  

  #----------------------------------------------------------------------------------------------
  # extract the colours from the tsne_res table
  
  # when spliting the table horizontally by a factor need to reset the previous k values now the array 
  # is ordered 
  if (split_horizontal_label) {
    tsne_res$data_k <- as.factor(old_k)
    data_log$data_k <- plyr::mapvalues(data_log$data_k,
                                       from = c(c(1:initial_k), (c(1:initial_k) + initial_k)), 
                                        to = c(c(1:initial_k), c(1:initial_k)))
    
      
  }  
  
  color_df <- unique(tsne_res[,c("data_k", "data_colors")])
  color_df <- color_df[order(color_df$data_k),]
  k_colors <- factor(color_df$data_colors, levels = color_df$data_colors)
  names(k_colors) <- color_df$data_k
  
  ann_colors = list(
    data_k = k_colors
  )
  
  heat_anno <- data.frame(row.names  = row.names(data_log),  data_k = data_log[,dim(data_log)[2]])
  
  if (!is.null(mask_low )) {
    
    mask_low <- mask_low[match(row.names(data_log) , row.names(mask_low)),]
    data_log[,c(1:(ncol(data_log)-1))] <- as.matrix( mask_low[, colnames(data_log[1:(ncol(data_log)-1)])]) *  as.matrix(data_log[,c(1:(ncol(data_log)-1))])
    color_vals[1] <- col2hex("grey")
  }
  
  
  if (!is.null(label_cols )) {
    label_anno_df <- as.data.frame(sapply(label_factors, function(x) ifelse(tsne_res[,colnames(tsne_res) == label_cols] == x, 1, 0))) 
    row.names(label_anno_df) <-  tsne_res$keep_transcript
    heat_anno <- merge(x = heat_anno,  y = label_anno_df, by.x = 0, by.y = 0, all  = FALSE)
    row.names(heat_anno) <- heat_anno$Row.names
    heat_anno$Row.names <- NULL
    tf_col <- as.factor(col2hex(c("white","black")))
    names(tf_col) <- c(0, 1)
    tf_col_l <- list(tf_col, tf_col)
    names(tf_col_l) <- label_factors
    ann_colors <- c(tf_col_l, ann_colors)
    
  } 
  
  if (!is.null(gene_labels )) { 
    ordered_gene_names <- tsne_res[match(row.names(data_log) , tsne_res$keep_transcript), "mgi_symbol"]
    labels_row = rep("", dim(data_log)[1])
    labels_row[ordered_gene_names %in% gene_labels] = as.character(ordered_gene_names[ordered_gene_names %in% gene_labels])
    labels_row_df <- data.frame(position = c(1:length(labels_row)), gene = labels_row)
    labels_row_df <- labels_row_df[labels_row_df[,2] != "",]  
    
    labels_row_p <- ggplot(labels_row_df, aes(x = 1, y = position, label = gene)) + 
      geom_text_repel(nudge_x = 0.1,
                      nudge_y = 10,
                      direction = "y") + 
    scale_y_reverse( limits=c((dim(data_log)[1]),0), expand = c(0,0)) +
    scale_x_continuous(limits = c(1,1.6), expand = c(0,0))
    
    
  } else {
    labels_row = rep("", dim(data_log)[1])
  }
  
  if (split_horizontal_label) {
    label_i <- min(which(row.names(data_log) %in% row.names(heat_anno)[rowSums(heat_anno[,c(2,3)]) > 0]))
    row_gaps <- label_i
  }  

  if (plot_what == "heatmap") {
    if (is.null(meta)) {
      h1 <- pheatmap(as.matrix(data_log[c(1:(dim(data_log)[2]-1))]), 
                   cluster_rows=FALSE,cluster_cols=FALSE, 
                   #show_rownames = FALSE, 
                   labels_row = labels_row,
                   color =  color_vals,
                   breaks = breaksList,
                   annotation_row = heat_anno,
                   annotation_colors = ann_colors,
                   annotation_names_row = FALSE,
                   fontsize_row = 10,
                   fontsize_col = 12,
                   width =  5, 
                   height = 3,
                   gaps_row = row_gaps,
                   gaps_col = col_gaps,
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
    h2 <- h1$gtable
    # this needs to be adjust the following grob indexes 
    h2 <- gtable_remove_grobs(h2, names = c( "annotation_legend", "legend", "col_annotation_names", "row_names"))
    if (!is.null(label_cols )) {
      labels_row_g <- ggplotGrob(labels_row_p)
      labels_row_g <- gtable_filter(labels_row_g,"panel")
      h2 <- gtable_add_cols(h2, h2$widths[1], pos = -1)
      h2 <- gtable_add_cols(h2, h2$widths[1] *6 , pos = -1)
      h2 <- gtable_add_grob(h2,h2$grobs[[3]] , t = 1, b = 1, l = 3, r = 3, name = "gene_clusters")
      h2 <- gtable_add_grob(h2, labels_row_g$grobs[[1]], t = 1, b = 1, l = 4, r = 4, name = "gene_names")
      h2 <- gtable_remove_grobs(h2, names = c( "row_annotation"))
      h2 <- gtable_add_cols(h2, h2$widths[2], pos = 0)
      bg <- rectGrob(gp = gpar(alpha = 0))
      h2 <- gtable_add_grob(h2, bg, t = 1, b = 1, l = 1, r = 1, name = "fill")
    }  
    
    grid.newpage()
    grid.draw(h2)
    return(h2)
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
  

  
  # add labels to figures if they are provided 
  if (!is.null(gene_labels)) {
    label_only <- tsne_res[row.names(tsne_res) %in% gene_labels,]
    label_only$name <- row.names(label_only)
    g1 <- g1 + geom_text_repel(data = label_only, aes(x=x, y=y, label = name), 
                               point.padding = unit(1, 'lines'), 
                               colour = "black") 
  }
  
  
}
