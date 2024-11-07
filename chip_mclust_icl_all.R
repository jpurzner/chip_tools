chip_mclust_icl_all <- function(data, crossing_points_df, G_range = 2:8) {
  results <- list()  # Initialize an empty list to store results for each column
  
  # Iterate through each column in the dataframe
  for (column_name in names(data)) {
    column_data <- data[[column_name]]  # Define column_data here to ensure it's always available
    
    # Determine tail_cutoff for the current column based on crossing_points_df
    tail_cutoff <- crossing_points_df[crossing_points_df$Column == column_name, "CrossingSmoothedValue"]
    
    if (length(tail_cutoff) > 0) {
      tail_cutoff <- as.numeric(tail_cutoff)  # Ensure it's numeric
      above_cutoff_indices <- which(column_data > tail_cutoff)
      
      # Data points for clustering (below cutoff)
      data_to_process <- column_data[-above_cutoff_indices]
    } else {
      data_to_process <- column_data
      above_cutoff_indices <- integer(0)  # No indices above cutoff
    }
    
    if(length(data_to_process) > 0) {
      mclust_data <- data.frame(column_data = data_to_process)
      icl_results <- mclustICL(mclust_data, G = G_range)
      best_icl_value <- max(icl_results)
      best_model_info <- which(icl_results == best_icl_value, arr.ind = TRUE)
      best_g <- best_model_info[1]
      best_model <- names(icl_results)[best_model_info[2]]
      best_mclust_model <- Mclust(mclust_data, G = best_g, model = best_model)
      
      final_classification <- rep(NA, length(column_data))
      final_classification[-above_cutoff_indices] <- best_mclust_model$classification
      
      # Assign a unique group ID for data points above the cutoff
      final_classification[above_cutoff_indices] <- max(final_classification, na.rm = TRUE) + 1
      
      results[[column_name]] <- final_classification
    } else {
      results[[column_name]] <- rep(NA, length(column_data))  # All data was pruned away
    }
  }
  
  # Combine the results into a single dataframe
  results_df <- do.call(cbind, results)
  colnames(results_df) <- paste(colnames(results_df), "class", sep = "_")
  
  # Including original data for reference
  final_results <- cbind(data, results_df)
  
  return(final_results)
}