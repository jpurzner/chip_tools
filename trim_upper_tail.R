require(dplyr)
require(ggplot2)
require(RColorBrewer)
require(patchwork)
require(tidyr)

trim_upper_tail <- function(data_frame, derivative_threshold = 3e-4,  generate_plots = TRUE) {
  
  # TODO add plot only function so users can set their threshold 
  # TODO autoset the derivative_threshold 
  
  
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
  
  results_list <- list()
  
  # Assuming the first column is not to be processed (e.g., gene names)
  #for (col_name in names(data_frame)[-1]) {
  for (col_name in names(data_frame)) {
    results_list[[col_name]] <- process_column(col_name, data_frame)
  }
  
  # Initialize a list to store crossing points for each column
  crossing_points_list <- list()
  
  # Iterate through each column's results in results_list
  for (col_name in names(results_list)) {
    # Extract the processed data for the column
    processed_data <- results_list[[col_name]]
    
    # Identify the crossing point
    crossing_point <- processed_data %>%
      filter(SmoothedDeriv1 <= derivative_threshold) %>%
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
  
  # Convert CrossingRank and CrossingSmoothedValue into numeric if they're not already
  crossing_points_df$CrossingRank <- as.numeric(as.character(crossing_points_df$CrossingRank))
  crossing_points_df$CrossingSmoothedValue <- as.numeric(as.character(crossing_points_df$CrossingSmoothedValue))
  
  # Generate plots if requested
  if (generate_plots) {
    print(create_combined_plots(results_list, crossing_points_df,derivative_threshold))
    print(create_histogram_plot(data_frame, crossing_points_df))
  }
  
  return(crossing_points_df)
}

# generates the rank plots
create_combined_plots <- function(results_list, crossing_points_df, derivative_threshold) {
  combined_data <- bind_rows(lapply(names(results_list), function(name) {
    results_list[[name]] %>%
      mutate(Dataset = name)
  }), .id = "DatasetID") %>%
    mutate(Dataset = as.factor(Dataset))
  
  combined_data_with_crossings <- merge(combined_data, crossing_points_df, by.x = "Dataset", by.y = "Column", all.x = TRUE)
  
  plot_smoothed <- ggplot(combined_data_with_crossings, aes(x = Rank, y = Smoothed, color = Dataset)) +
    geom_line() +
    geom_hline(aes(yintercept = CrossingSmoothedValue), linetype = "dashed") +
    geom_vline(aes(xintercept = CrossingRank), linetype = "dashed") +
    geom_point(aes(x = CrossingRank, y = CrossingSmoothedValue), size = 2, shape = 21, fill = NA) +
    #scale_color_brewer(palette = "Set1") +
    scale_x_reverse() +
    theme_minimal() +
    labs(x = "Rank", y = "Smoothed Value") +
    theme(legend.position = "bottom")
  
  plot_deriv1 <- ggplot(combined_data_with_crossings, aes(x = Rank, y = SmoothedDeriv1, color = Dataset)) +
    geom_line() +
    geom_hline(yintercept = derivative_threshold, linetype = "dashed") +
    #scale_color_brewer(palette = "Set1") +
    scale_x_reverse() +
    theme_minimal() +
    theme(legend.position = "none")
  
  combined_plot <- plot_smoothed / plot_deriv1 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  return(combined_plot)
}

# plotting function for histogram with line at the cutoff 
create_histogram_plot <- function(data_frame, crossing_points_df) {
  
  data_frame$gene <- row.names(data_frame)
  data_frame <- as_tibble(data_frame)
  
  # Convert the data frame to long format using pivot_longer() from tidyr
  data_frame_long <- data_frame %>%
    pivot_longer(
      -gene, 
      names_to = "variable", 
      values_to = "value"
    )
  
  # Merge crossing_points_df with the long format data frame
  data_frame_long <- data_frame_long %>%
    left_join(crossing_points_df, by = c("variable" = "Column"))
  
  # Draw histograms with a vertical line at the CrossingSmoothedValue, all within the pipe
  tail_histogram <- data_frame_long %>%
    ggplot(aes(x = value)) + 
    facet_wrap(~variable, scales = "free_x") + 
    geom_histogram(bins = 100, color = "black", fill = "gray") +
    geom_vline(data = . %>% distinct(variable, CrossingSmoothedValue), aes(xintercept = CrossingSmoothedValue), color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(x = "Value", y = "Frequency") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(tail_histogram)
}
