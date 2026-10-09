#' Mean expression lines per cluster
#'
#' One small-multiple panel per cluster, each showing the mean trajectory of
#' that cluster's genes in every data set.
#'
#' @param expr_list List of genes x timepoints expression matrices, one per
#'   data set.
#' @param ds_names Data set names, in the order of `expr_list`.
#' @param gene2cluster Two-column data frame mapping `gene` to `cluster_id`.
#' @param p_time Pseudotime axis values.
#'
#' @return A `ggplot` object.
#'
#' @seealso [plot_cluster_line_EBseq()] for the version that takes sample
#'   metadata and draws confidence ribbons; [plot_group_timeline()].
#'
#' @examples
#' \dontrun{
#' plot_cluster_line(expr_list, c("scott", "hatten"), gene2cluster, p_time)
#' }
#' @import ggplot2
#' @export
plot_cluster_line <- function (expr_list, ds_names, gene2cluster, p_time) {
  require_pkg("Mfuzz", "ggrepel")

  #require(naturalsort)
  
  # loop through the list of expression objects and extract data, merge with community clusters, 
  #  flatten the data and rbind into a common df for use with ggplot2
  expr <- exprs(expr_list[[1]])
  expr <- merge(x = gene2cluster, y = expr, by.x = 1, by.y = 0)
  expr <- melt(expr, id.vars = c("gene", "cluster_id"))
  colnames(expr)[3] <- "time"
  expr$ds <- ds_names[1]
  print(head(expr))
  expr <- expr[,c(1,2,3,5,4)]
  
  for (i in 2:length(expr_list)) {
    expr_tmp <- exprs(expr_list[[i]])
    expr_tmp <- merge(x = gene2cluster, y = expr_tmp, by.x = 1, by.y = 0)
    expr_tmp <- melt(expr_tmp, id.vars = c("gene","cluster_id"))
    colnames(expr_tmp)[3] <- "time"
    expr_tmp$ds <- ds_names[i]
    expr_tmp <- expr_tmp[,c(1,2,3,5,4)]
    expr <- rbind(expr, expr_tmp)
  }

  
  # to overlap plot the data we need to use a pseudotime 
  # this can be specified using the following df 
  # 
  # ds  time   pseudotime
  #Scott_GNP  E15  0
  #Scott_GNP  P1  0 
  
  # series of re-ordering to make the pseudotime df 
  ### this block of code is defunct unless we need to seperate different ds with same time
  #pseudo_time <- ddply(expr, c("ds", "time"), summarise, N = length(value))
  #pseudo_time$N <- NULL
  #pseudo_time$ds <- factor(pseudo_time$ds, levels = ds_names)
  #pseudo_time <- pseudo_time[order(pseudo_time$ds),]
  #pseudo_time$time <- factor(pseudo_time$time, levels = p_time)
  #pseudo_time <- pseudo_time[order(pseudo_time$time),]
  #pseudo_time$ptime <-  as.numeric(pseudo_time$time)  
  #row.names(pseudo_time) <- seq(length=nrow(pseudo_time)) 
  
  ## TODO automatically try and make the pseudotime table
  # finds the data set with the highest number of time points
  #ds_largest <- table(pseudo_time$ds)
  #ds_largest <- names(ds_largest[order(ds_largest, decreasing  = T)])[1]
  #largest_time <- pseudo_time[pseudo_time$ds == ds_largest, 2]
  
  ##END of defunct code 
  
  # much shorter pseudo time setting 
  ptime_list <- seq(length=length(p_time))
  names(ptime_list) <- p_time
  ptime_list <- factor(ptime_list)
  ptime_list <- as.data.frame(ptime_list)
  print(head(ptime_list))
  print(head(expr))
  expr <- merge(x = expr, y = ptime_list, by.x = 3, by.y = 0, all.x = T)
 
  expr <- expr[,c(2:4,1,6,5)]
  colnames(expr)[5] <- "ptime"
  # add the pseudo to df
  #expr$ptime <- ptime_list[expr$time]
  #expr <- expr[,c(1:5,7,6)]
  #print(head(expr))
  
  # collapse the data frame to get averages
  mean_expr <- ddply(expr, c("cluster_id", "ds", "time", "ptime"), summarise, mean = mean(value), sd = sd(value), sem = sd(value)/sqrt(length(value)))
  #print(head(mean_expr))
  
  p <- ggplot(mean_expr, aes(x=ptime, y=mean, colour = ds, group = ds)) + geom_line() + 
    geom_ribbon(aes(ymax = mean + sd, ymin = mean - sd, fill = ds), alpha = 0.3, colour=NA) + 
    facet_wrap(~ cluster_id) + geom_point() + 
    geom_text_repel(aes(label=time), 
                    nudge_y = ifelse(mean_expr$ds == ds_names[2], -0.5, ifelse(mean_expr$ds == ds_names[1], 0.5, 0.75)), 
                    size = 3, force =1 ) + 
    #scale_y_continuous(limits = c(-1, 1)) + 
    theme(legend.position="bottom", 
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks=element_blank(), 
          strip.background = element_blank(),
          strip.text.x = element_blank())
    
  
  return(p)
  
  # link genes to co_id
  
  
  
  
  
} 