ggroc <- function(roc, showAUC = TRUE, interval = 0.2, breaks = seq(0, 1, interval)){
  require(pROC)
  
  if(class(roc) == "roc") {
    roc_df = data.frame(x = rev(roc$specificities), y = rev(roc$sensitivities), name = "ROC")
  } else if (is.list(roc)) {
    roc_df_list <- lapply(c(1:length(roc)), function (n) { df = data.frame(x = rev(roc[[n]]$specificities), 
                                                            y = rev(roc[[n]]$sensitivities), 
                                                            name = gsub("_", " ", names(roc)[n]))
                                            return(df)
    })
    roc_df = bind_rows(roc_df_list)
  } else {
    simpleError("Please provide roc object from pROC package")
  }
  
  ggroc <- ggplot(roc_df, aes(x = x, y = y, group= name, color = name)) +
    geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0), alpha = 0.5, color = "black") + 
    geom_step() +
    scale_x_reverse(name = "Specificity",limits = c(1,0), breaks = breaks, expand = c(0.001,0.001)) + 
    scale_y_continuous(name = "Sensitivity", limits = c(0,1), breaks = breaks, expand = c(0.001, 0.001)) +
    theme_bw() + 
    theme(axis.ticks = element_line(color = "grey80")) +
    coord_equal() + 
    theme(legend.position = c(0.8, 0.2), 
          legend.title = element_blank())
}