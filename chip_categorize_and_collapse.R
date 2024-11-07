library(dplyr)
library(tidyr)

# Function to categorize data and collapse categories as per conditions
chip_categorize_and_collapse <- function(data, histone_break_table) {
  # Step 1: Categorize data for each replicate
  categorized_data <- data %>%
    rowwise() %>%
    mutate(across(everything(), ~ {
      replicate_name <- cur_column()
      if(!replicate_name %in% histone_break_table$replicate_names) return(.)
      
      cutoff_info <- histone_break_table %>%
        dplyr::filter(replicate_names == replicate_name) %>%
        dplyr::select(CrossingSmoothedValue, cutoff_50) %>%
        slice(1)
      
      if(is.na(cutoff_info$CrossingSmoothedValue) | is.na(cutoff_info$cutoff_50)) return(NA)
      
      if(. < cutoff_info$cutoff_50) {
        "unbound"
      } else if(. >= cutoff_info$cutoff_50 & . < cutoff_info$CrossingSmoothedValue) {
        "bound"
      } else {
        "tail"
      }
    }, .names = "{.col}_cat")) %>%
    ungroup()
  
  # Step 2: Collapse categories along condition
  condition_mapping <- histone_break_table %>%
    dplyr::select(replicate_names, condition) %>%
    distinct()
  
  collapsed_data <- categorized_data %>%
    pivot_longer(cols = ends_with("_cat"), names_to = "replicate_cat", values_to = "category") %>%
    dplyr::mutate(replicate_name = gsub("_cat$", "", replicate_cat)) %>%
    left_join(condition_mapping, by = c("replicate_name" = "replicate_names")) %>%
    group_by(gene_id, condition) %>%
    dplyr::summarize(category = case_when(
      "tail" %in% category ~ "tail",
      "bound" %in% category ~ "bound",
      TRUE ~ "unbound"
    ), .groups = "drop")
  
  return(collapsed_data)
}

# Example usage:
# Assuming `data` is your dataframe with replicate data and gene_id as rows
# and `histone_break_table` is loaded as shown
# results <- categorize_and_collapse(data, histone_break_table)