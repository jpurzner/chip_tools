chip_binarize_tail_removed <- function(data_frame, tail_cutoff_table) {
    results <- list() # Initialize a list to store the results
    
    # Iterate through each column in the data frame
    for (column_name in names(data_frame)) {
      # Retrieve the corresponding cutoff value for the column
      cutoff_value <- tail_cutoff_table[tail_cutoff_table$Column == column_name, "CrossingSmoothedValue"]
      
      if (length(cutoff_value) > 0) {
        # Trim the data based on the cutoff value
        trimmed_values <- data_frame[[column_name]][data_frame[[column_name]] < cutoff_value]
        max_val <- max(trimmed_values)
        # Apply the binarize_counts function
        results[[column_name]] <- binarize_counts(trimmed_values, return_cutoff_mean = TRUE, exclude_extreme = TRUE, label_name = column_name, mu = c(0.5, max_val - 0.2))
      } else {
        # If no cutoff value is found, return NA or an appropriate value
        results[[column_name]] <- NA
      }
    }
    
    # Convert the list  into a dataframe
    results <- do.call(rbind, lapply(results, function(x) {
      df <- as.data.frame(t(x))
      df$Column <- rownames(df)
      rownames(df) <- NULL
      return(df)
    }))
    
    # Ensure 'Column' is a character for both dataframes
    results$Column <- as.character(row.names(results))
    tail_cutoff_table$Column <- as.character(tail_cutoff_table$Column)
    
    # Merge the dataframes on the 'Column' column
    merged_df <- merge(tail_cutoff_table, results, by = "Column")
    
    return(merged_df)
}


  