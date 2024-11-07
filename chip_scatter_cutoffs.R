library(dplyr)
library(ggplot2)
library(ggExtra)
library(ggrepel)

chip_scatter_cutoffs <- function(data, x_col, y_col, x_condition, y_condition, annotation) {
  # Extract cutoff values for x and y conditions
  x_cutoffs <- annotation %>%
    filter(replicate_names == x_condition) %>%
    summarise(min = cutoff_50, max = CrossingSmoothedValue) %>%
    collect() %>%
    unlist()
  
  y_cutoffs <- annotation %>%
    filter(replicate_names == y_condition) %>%
    summarise(min = cutoff_50, max = CrossingSmoothedValue) %>%
    collect() %>%
    unlist()
  
  # Prepare the data
  plot_data <- data %>%
    mutate(x_value = Winsorize(get(x_col), maxval = 0.9),
           y_value = Winsorize(get(y_col), maxval = 0.9),
           diff_gene = as.factor(diff_gene)) %>%
    select(x_value, y_value, diff_gene, mgi_symbol)
  
  # Modify labels to replace underscores and remove 'rep' and following
  modified_x_label <- gsub("_", " ", gsub("rep.*", "", x_col))
  modified_y_label <- gsub("_", " ", gsub("rep.*", "", y_col))
  
  
  # Generate the plot
  plot <- ggplot(plot_data, aes(x = x_value, y = y_value, fill = diff_gene, colour = diff_gene)) +
    geom_rect(aes(xmin = x_cutoffs['min'], xmax = x_cutoffs['max'], ymin = y_cutoffs['min'], ymax = y_cutoffs['max']),
              fill = NA, color = "black", inherit.aes = FALSE) +
    geom_point(aes(alpha = diff_gene, size = diff_gene)) +
    scale_colour_manual(values = c("differentiation" = "red", "other" = "grey")) +
    scale_fill_manual(values = c("differentiation" = "red", "other" = "grey")) +
    scale_alpha_manual(values = c("differentiation" = 1, "other" = 0.3)) +
    scale_size_manual(values = c("differentiation" = 1, "other" = 0.5)) +
    geom_text_repel(data = plot_data %>% filter(mgi_symbol %in% label_genes), 
                    aes(label = mgi_symbol),
                    colour = "black", 
                    point.padding = 0.2,
                    box.padding = 0.2,
                    min.segment.length = 0,
                    size = 4, 
                    seed = 42,
                    alpha = 1, 
                    force = 10) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = modified_x_label, y = modified_y_label) +
    
    geom_point(data = subset(data, mgi_symbol %in% label_genes & max_cat == "P14_P56" ), 
               aes(x = Winsorize(get(x_col), maxval = 0.9), 
                   y = Winsorize(get(y_col), maxval = 0.9)), 
               colour = "black", 
               shape = 1,         # Hollow squares
               size = 2, alpha = 1) + 
    
    geom_point(data = subset(data, mgi_symbol %in% c("Lhx9", "Hoxa10")), 
               aes(x = Winsorize(get(x_col), maxval = 0.9), 
                   y = Winsorize(get(y_col), maxval = 0.9)), 
               colour = "black",
               shape = 4,         # Hollow squares
               size = 2, alpha = 1) + 
    
    geom_text_repel(data = subset(data, mgi_symbol %in% c("Lhx9", "Hoxa10")), 
                    aes(x = Winsorize(get(x_col), maxval = 0.9), 
                        y = Winsorize(get(y_col), maxval = 0.9),
                        label = mgi_symbol), 
                    colour = "black",
                    point.padding = 0.2,
                    box.padding = 0.2,
                    min.segment.length = 0,
                    size = 4, 
                    seed = 42,
                    alpha = 1, 
                    force = 10) +   
    theme_classic() + 
    theme(legend.position = "none")
  
  plot_with_marginals <- ggExtra::ggMarginal(plot, type = "density", groupFill = TRUE)
  
  return(plot_with_marginals)
}


