line_scatter_facet <- function(df, log2_norm_list, meta_list, pseudo_time, 
                               gene_col = 1, facet_var = "cluster_id", 
                               grouping = NULL, all_genes = NULL, FC_size = FALSE, 
                               ncols = 4, ylimits = c(-2,2), grid = FALSE) {
  
  # grouping 2 column table with the linear index of facet var in col1 and the group name in col2 
  
  # to create an overlay of the overall TSNE just pass the histone_cluster output from the initial 
  # TSNE figure to all_genes
  
  source('~/Dropbox/jp_seq/plot_cluster_line_EBseq.R')
  require(grid)
  require(ggplot2)
  require(pheatmap)
  require(ggrepel)
  require(gtable)  
  require(cowplot)
  require(plyr)
  require(RColorBrewer)
  require(Biobase)
  
  df_size <- dim(df)[2]

  if (FC_size) {
  fc <- data.frame(genes = row.names(log2_norm_list[[1]]), 
          ds1 = rowMax(log2_norm_list[[1]]) - rowMin(log2_norm_list[[1]]),
          ds2 = rowMax(log2_norm_list[[2]]) - rowMin(log2_norm_list[[2]]))
  df <- merge(x = df, y =  fc, by.x = gene_col, by.y = 1, all = FALSE)
  print(head(df))
  
  }  
  
  # cross reference meta
  if (is.null(all_genes) & (!FC_size)) {
    sp <- ggplot(df, aes(x= x, y = y, colour = data_k)) +
      geom_point(size = 0.3) + 
      facet_wrap(c(facet_var)) + 
      scale_colour_manual(values = data_colors) + 
      guides(colour = guide_legend(override.aes = list(size=3), ncol = 3)) + 
      theme(legend.position="bottom",
          legend.title=element_blank(),  
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(), 
          axis.line=element_blank())
          #plot.background = element_rect(fill = 'white', colour = 'red'))
          #panel.background=element_blank(),
          #panel.border=element_blank(),
          #panel.grid.major=element_blank(),
          #panel.grid.minor=element_blank())
    
  } else if (!is.null(all_genes) & (!FC_size)) {
    sp <- ggplot(df, aes(x= x, y = y, colour = data_k)) +
      geom_point(size = 0.3) + 
      geom_point(data = all_genes, aes(x = x, y = y, colour = data_k), alpha = 0.05, size = 0.1) +
      facet_wrap(c(facet_var)) + 
      scale_colour_manual(values = data_colors) + 
      guides(colour = guide_legend(override.aes = list(size=3), ncol = 3)) + 
      theme(legend.position="bottom",
            legend.title=element_blank(),  
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(), 
            axis.line=element_blank())
    #plot.background = element_rect(fill = 'white', colour = 'red'))
    #panel.background=element_blank(),
    #panel.border=element_blank(),
    #panel.grid.major=element_blank(),
    #panel.grid.minor=element_blank())
    
  } else if (!is.null(all_genes) & FC_size) {
    sp <- ggplot(df, aes(x= x, y = y, fill = data_k, size = ds1)) +
      geom_point(shape = 21, alpha = 0.8) + 
      scale_size(range = c(0.2,4)) + 
      geom_point(data = all_genes, aes(x = x, y = y, colour = data_k), alpha = 0.05, size = 0.1) +
      facet_wrap(c(facet_var)) + 
      scale_colour_manual(values = data_colors) + 
      scale_fill_manual(values = data_colors) +
      guides(colour = guide_legend(override.aes = list(size=3), ncol = 3)) + 
      theme(legend.position="bottom",
            legend.title=element_blank(),  
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(), 
            axis.line=element_blank())
    #plot.background = element_rect(fill = 'white', colour = 'red'))
    #panel.background=element_blank(),
    #panel.border=element_blank(),
    #panel.grid.major=element_blank(),
    #panel.grid.minor=element_blank())
    
    
  }
  else if (is.null(all_genes) & FC_size) {
    sp <- ggplot(df, aes(x= x, y = y, fill = data_k, size = ds1)) +
      geom_point(shape = 21, alpha = 0.8) + 
      scale_size(range = c(0.2,4)) + 
      facet_wrap(c(facet_var)) + 
      scale_colour_manual(values = data_colors) + 
      scale_fill_manual(values = data_colors) +
      guides(colour = guide_legend(override.aes = list(size=3), ncol = 3)) + 
      theme(legend.position="bottom",
            legend.title=element_blank(),  
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(), 
            axis.line=element_blank())
    #plot.background = element_rect(fill = 'white', colour = 'red'))
    #panel.background=element_blank(),
    #panel.border=element_blank(),
    #panel.grid.major=element_blank(),
    #panel.grid.minor=element_blank())
    
  }
  #sp8 <- sp 
  
  sp <- ggplotGrob(sp)
  
  #all_strips <- names(sp$grobs)[grep("strip", names(sp$grobs))]
  #all_axis <- names(sp$grobs)[grep("axis", names(sp$grobs))]
  #print(c(all_strips, all_axis))
  
  gb <- gtable_filter(sp, "guide-box", trim = TRUE)
  gb <- gb$grobs
  gb <- gtable_remove_grobs(table = gb[[1]], name = "legend.box.background")
  gb1 <- gb[1,1]
  gb1 <- gtable_filter(gb1, "guide", trim = TRUE)
  if (FC_size) {
    gb2 <- gb[1,3]
    gb2 <- gtable_filter(gb2, "guide", trim = TRUE)
  }
  #gb1 <- gb1$grobs
  #gb2 <- gb2$grobs
  
  bg <- rectGrob(gp = gpar(alpha = 0))
  
  sp2 <- gtable_filter(sp, "panel", trim = TRUE)
  
  # need to change to generalize 
 
  #sp2 <- gtable_filter(sp2, "gTree", trim = TRUE)

  sp2 <- sp2$grobs
  # removes the zeroGrob panels 
  sp2 <- sp2[which(lapply(sp2, length) == 5)]
  
  # generate the line graphs
  gene2cluster <- data.frame(gene = df[,gene_col], cluster_id = df[,facet_var])
  print(head(gene2cluster))
  lp <- plot_cluster_line_EBseq(log2_norm_list, meta_list, gene2cluster, p_time = pseudo_time, ylimits = ylimits)
  #lp8 <- lp
  lp <- ggplotGrob(lp)
  lp2 <- gtable_filter(lp, "panel", trim = TRUE)
  lp2 <- lp2$grobs
  # removes the zeroGrob panels 
  lp2 <- lp2[which(lapply(lp2, length) == 5)]
  yaxis <- gtable_filter(lp, "axis-l-1-1", trim = TRUE)
  yaxis <- yaxis$grobs

  # determine dimensions of the plot matrix  
  plot_num <- length(sp2) + 2
  
  empty_needed <- ncols - (plot_num %% ncols)
  nrows <- (plot_num + empty_needed) / ncols 
  
  # create lists of the grobs and merge into matrix 
  if (FC_size) {
    sp2 <- c(sp2, list(gb1), list(gb2), rep(list(bg), empty_needed) )  
  } else {
    sp2 <- c(sp2, list(gb1), rep(list(bg), empty_needed) )  
  }
  sp2_mat <- matrix(sp2, ncol = ncols)
  lp2 <- c(lp2,    rep(list(bg), empty_needed + 2))  
  lp2_mat <- matrix(lp2, ncol = ncols)
  
  plot_mat_list <- list(lp2_mat, sp2_mat)
  all_mat <- do.call(rbind, plot_mat_list)[order(sequence(sapply(plot_mat_list, nrow))), ]
  
  # set up grouping coloring
  if (!is.null(grouping)) {
    group_cols <- brewer.pal(length(unique(grouping[,2])),"Dark2")
    cols_df <- data.frame(groups = unique(grouping[,2]), cols = group_cols)
    grouping <- merge(x = grouping, y = cols_df, by.x = 2, by.y = 1)
    grouping <- grouping[,c(2,1,3)]
    facet2group <- data.frame(index = c(1:length(levels(df[,facet_var]))), facet = levels(df[,facet_var]))
    facet2group <- merge(x = facet2group, y = grouping, by.x = 1, by.y = 1, all.x = TRUE)
    #facet2group[is.na(facet2group$group_name),3] <- "none"   
    # drop the index 
    facet2group <- facet2group[,c(2:4)]
    colnames(facet2group)[3] <- "group_cols"
    df <- merge(x = df, y = facet2group, by.x = which(colnames(df)==facet_var), by.y = 1, all.x = TRUE)
    df <- df[,c(2:df_size, 1,df_size+1, df_size+2)]
  }

  rect <- rectGrob(gp = gpar(col = "black", fill = "white"))
  
  gt_mat <- gtable_matrix("demo", grobs = all_mat,  
                          widths = unit(rep(2,ncols), "null"),
                          heights = unit(rep(c(1,2),nrows), "null"))
  
  # generate a grid around plots 
  if (grid) {
    n = 1 
    for (j in c(1:nrows)) {
      for (i in c(1:ncols)) {
        if ((!is.null(grouping)) & (n %in% grouping[,1] )) {
          curr_col <- grouping[which(grouping[,1] == n),3]
          col_rect <- rectGrob(gp = gpar(col = curr_col, fill = "white", lwd = 5))
          gt_mat <- gtable_add_grob(gt_mat, col_rect, t =(j*2)-1, b = (j*2), l = i, z = -Inf)
        } else {
          gt_mat <- gtable_add_grob(gt_mat, rect, t =(j*2)-1, b = (j*2), l = i, z = -Inf)
        }
        n = n + 1 
      } # end for 
    } # end for 
  } # end if 
  
  # add the y-axis for the expression plots 
  gt_mat <- gtable_add_cols(gt_mat, widths = unit(0.2, "null"), pos = 0 )
  gt_mat <- gtable_add_cols(gt_mat, widths = unit(0.2, "null"), pos = 0 )
  for (i in seq(1,dim(all_mat)[1], by = 2)) {
    gt_mat <- gtable_add_grob(gt_mat, yaxis, t =i, b = i, l = 2, r= 2, z = Inf, clip = "off", name = "ylab-r")
  }
  
  grid.newpage()
  grid.draw(gt_mat)
  invisible(df)
} # end of function 