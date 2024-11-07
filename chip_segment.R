chip_segment <- function(data_frame, deriv_cutoff = 3e-4, G_range = 1:5) {
  library(mclust)
  library(dplyr)
  library(ggplot2)
  
 
  process_column <- function(column_name, data_frame) {
    # Rank based on the specified column
    ranked_data <- data_frame %>%
      mutate(Rank = rank(-get(column_name)))
    
    # Smooth the specified column
    smoothed <- loess(get(column_name) ~ Rank, data = ranked_data, span = 0.01, control = loess.control(surface = "direct"))
    ranked_data$Smoothed <- predict(smoothed)
    
    # Calculate first derivative
    ranked_data <- ranked_data %>%
      arrange(desc(Rank)) %>%
      mutate(Deriv1 = c(NA, diff(Smoothed)))
    
    # Smooth the first derivative
    smoothed_deriv1 <- loess(Deriv1 ~ Rank, data = ranked_data, span = 0.1, na.action = na.exclude, control = loess.control(surface = "direct"))
    ranked_data$SmoothedDeriv1 <- predict(smoothed_deriv1, na.action = na.exclude)
    
    return(ranked_data)
  }
  
  
  
  # Main section of code 
  results_list <- list()
  # Assuming the first column is not to be processed (e.g., gene names)
  for (col_name in names(histone_rlog_prune_mean_only)[-1]) {
    results_list[[col_name]] <- process_column(col_name, histone_rlog_prune_mean_only)
  }
  
  
  
  
  # Initialize a list to store crossing points for each column
  crossing_points_list <- list()
  
  # Iterate through each column's results in results_list
  for (col_name in names(results_list)) {
    # Extract the processed data for the column
    processed_data <- results_list[[col_name]]
    
    # Identify the crossing point
    crossing_point <- processed_data %>%
      filter(SmoothedDeriv1 <= deriv_cutoff) %>%
      arrange(Rank) %>%
      slice(1)  # Take the first row if multiple crossings exist
    
    # Extract crossing rank and smoothed value
    crossing_rank <- crossing_point$Rank
    crossing_smoothed_value <- crossing_point$Smoothed
    
    # Store the results in a list
    crossing_points_list[[col_name]] <- data.frame(
      Column = col_name,
      CrossingRank = ifelse(length(crossing_rank) > 0, crossing_rank, NA),
      CrossingSmoothedValue = ifelse(length(crossing_smoothed_value) > 0, crossing_smoothed_value, NA)
    )
  }
  
  # Combine all the results into a single data frame
  crossing_points_df <- do.call(rbind, crossing_points_list)
  
  # Optionally, remove rows with NA (where no crossing point was found)
  crossing_points_df <- na.omit(crossing_points_df)
  
  # Display the combined results
  print(crossing_points_df)
  
  # Convert CrossingRank and CrossingSmoothedValue into numeric if they're not already
  crossing_points_df$CrossingRank <- as.numeric(as.character(crossing_points_df$CrossingRank))
  crossing_points_df$CrossingSmoothedValue <- as.numeric(as.character(crossing_points_df$CrossingSmoothedValue))
  
  # Merge crossing points with combined_data for plotting
  # This assumes that there's a matching mechanism between crossing_points_df$Column and combined_data$Dataset
  # Here, we're creating a mock-up for the illustration
  # In your actual implementation, ensure the data structures are compatible for merging
  combined_data_with_crossings <- merge(combined_data, crossing_points_df, by.x = "Dataset", by.y = "Column", all.x = TRUE)
  
  
  
  # Combine all datasets into one data frame, adding a 'Dataset' column based on the list names
  combined_data <- bind_rows(lapply(names(results_list), function(name) {
    results_list[[name]] %>%
      mutate(Dataset = name)
  }), .id = "DatasetID") %>%
    mutate(Dataset = as.factor(Dataset))  # Convert 'Dataset' to a factor for coloring
  
  # Assuming crossing points have been calculated and added to each dataset in results_list
  # If not, you'll need to perform that step before this
  library(ggplot2)
  library(patchwork)
  
  # Assuming combined_data is correctly defined and available
  
  plot_smoothed <- ggplot(combined_data, aes(x = Rank, y = Smoothed, color = Dataset)) +
    geom_line() +
    geom_hline(data = combined_data %>% distinct(Dataset, .keep_all = TRUE), aes(yintercept = crossing_smoothed_value), linetype = "dashed") +
    geom_point(data = combined_data_with_crossings, aes(x = CrossingRank, y = CrossingSmoothedValue), size = 3, shape = 21, fill = "white") +
    
    geom_vline(data = combined_data %>% distinct(Dataset, .keep_all = TRUE), aes(xintercept = crossing_rank), linetype = "dashed") +
    scale_color_brewer(palette = "Set1", name = "Dataset") +
    scale_x_reverse() +
    theme_minimal() +
    labs(x = "Rank", y = "Smoothed Value") +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  plot_deriv1 <- ggplot(combined_data, aes(x = Rank, y = SmoothedDeriv1, color = Dataset)) +
    geom_line() +
    #geom_vline(aes(xintercept = crossing_rank), linetype = "dashed") +
    geom_hline(yintercept = 3e-4) +
    scale_color_brewer(palette = "Set1", name = "Dataset") +
    scale_x_reverse() +
    theme_minimal() +
    labs(x = "Rank", y = "SmoothedDeriv1") +
    theme(
      plot.title = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"  # Hide legend in the second plot
    )
  
  # Combine the plots with a single common legend
  combined_plot <- plot_smoothed / plot_deriv1 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  # Display the combined plot
  combined_plot
  

  chip_mclust_icl_all <- function(data, G_range = 1:5, tail_cutoff = NULL) {
    results <- list()  # Initialize an empty list to store results
    
    # Iterate through each column in the dataframe
    for (column_name in names(data)) {
      # Apply tail_cutoff if specified
      if (!is.null(tail_cutoff)) {
        data_to_process <- data[[column_name]][data[[column_name]] <= tail_cutoff]
      } else {
        data_to_process <- data[[column_name]]
      }
      
      # Check if data_to_process is not empty
      if(length(data_to_process) > 0) {
        # Ensure data is in the correct format for Mclust
        mclust_data <- data.frame(column_data = data_to_process)
        
        # Perform ICL to determine the best model
        icl_results <- mclustICL(mclust_data, G = G_range)
        
        # Find the best model based on ICL
        best_icl_value <- max(icl_results)
        best_model_info <- which(icl_results == best_icl_value, arr.ind = TRUE)
        best_g <- best_model_info[1]
        best_model <- names(icl_results)[best_model_info[2]]
        
        # Re-run Mclust with the best parameters
        best_mclust_model <- Mclust(mclust_data, G = best_g, model = best_model)
        
        # Store classification results
        results[[column_name]] <- best_mclust_model$classification
      } else {
        results[[column_name]] <- NA  # Assign NA if data was pruned away completely
      }
    }
    
    return(results)
  }
  
  histone_rlog_prune_mean_only$cluster <- NULL
  results_all <- chip_mclust_icl_all(histone_rlog_prune_mean_only)
  return(results_all)
}
