chip_gene_group_comp <- function(mark, data, annotation, remove_below_cutoff = FALSE, use_averaged_data = FALSE, t50_cutoff = NULL) {
  # Filter for replicates of the specific mark
  rep_names <- annotation %>%
    dplyr::filter(condition == mark) %>%
    dplyr::select(replicate_names) %>%
    dplyr::pull()
  
  # Extract cutoff_50, mean_upper, and mean_lower values for each replicate
  cutoff_values <- annotation %>%
    dplyr::filter(replicate_names %in% rep_names) %>%
    dplyr::select(replicate_names, cutoff_50, mean_upper, mean_lower)
  
  # Determine column name for averaged data
  averaged_data_column <- mark
  
  # merge the cutoff data 
  data_limited <- data %>%
    #dplyr::filter(diff_gene == "differentiation", H3K27me3_bin == 1) %>%
    dplyr::select(!!sym(mark), all_of(rep_names), diff_gene, t50_cat, avg_t50, mgi_symbol, !!averaged_data_column) %>%
    tidyr::pivot_longer(cols = all_of(rep_names), names_to = "replicate", values_to = "value") %>%
    dplyr::left_join(cutoff_values, by = c("replicate" = "replicate_names")) # Join cutoff values
  
  # Optionally remove data below cutoff_50
  if (remove_below_cutoff) {
    data_limited <- data_limited %>%
      dplyr::filter(value >= cutoff_50)
  }
  

  # Prepare data
  if (use_averaged_data) {
    data_limited <- data_limited %>%
      dplyr::mutate(value = .data[[mark]])  #%>%
    #dplyr::left_join(cutoff_values, by = c("replicate" = "replicate_names")) # Join cutoff values
  } else {

  }
  
  if (!is.null(t50_cutoff)) {
    data_limited <- data_limited %>%
      dplyr::mutate(diff_gene = as.factor(diff_gene),
                    t50_cat = if_else(avg_t50 > t50_cutoff, "late", "early")) 
  } else {
    data_limited <- data_limited %>%
      dplyr::filter(t50_cat %in% c("Late", "Mid"))
  }

  
  if (use_averaged_data) {
    
    # Compute p-values from Wilcoxon rank-sum test for each replicate and t50_cat combination
    if (nrow(data_limited) > 0 && length(unique(data_limited$t50_cat)) > 1) {
      early_values <- data_limited %>% dplyr::filter(t50_cat == "early") %>% dplyr::pull(value)
      late_values <- data_limited %>% dplyr::filter(t50_cat == "late") %>% dplyr::pull(value)
      if (length(early_values) > 1 && length(late_values) > 1) { # Ensure enough data for comparison
        p_values <- wilcox.test(early_values, late_values)$p.value
      } else {
        p_values <- NA
      }
    } else {
      p_values <- NA
    }
    
    
    # Create the average plot
    plot <- data_limited %>%
      ggplot(aes(x = t50_cat, y = value)) + 
      geom_violin() +
      geom_dotplot(stackdir = "centerwhole",stackgroups = TRUE,
                   binpositions = "all", binaxis = "y", binwidth = 0.01, size = 0.01) + 
      
      #geom_hline(aes(yintercept = cutoff_50), linetype = "dashed", color = "red") +
      #geom_hline(aes(yintercept = mean_upper), linetype = "dashed", color = "blue") +
      #geom_hline(aes(yintercept = mean_lower), linetype = "dashed", color = "blue") +
      stat_summary(fun = mean, geom = "point", color = "red", size = 3, aes(group = t50_cat)) +
      geom_text( aes(label = sprintf("p = %e", p_values), x = 1.5, y = Inf), hjust = 1.1, vjust = 2, check_overlap = TRUE) +
      labs(x = "t50 Category", y = "Value")
    
  } else {  
    
    # Compute p-values from Wilcoxon rank-sum test for each replicate and t50_cat combination
    p_values <- list()
    for (rep in unique(data_limited$replicate)) {
      data_subset <- data_limited %>% filter(replicate == rep)
      if (nrow(data_subset) > 0 && length(unique(data_subset$t50_cat)) > 1) {
        early_values <- data_subset %>% filter(t50_cat == "early") %>% pull(value)
        late_values <- data_subset %>% filter(t50_cat == "late") %>% pull(value)
        if (length(early_values) > 1 && length(late_values) > 1) { # Ensure enough data for comparison
          p_values[[rep]] <- wilcox.test(early_values, late_values)$p.value
        } else {
          p_values[[rep]] <- NA
        }
      } else {
        p_values[[rep]] <- NA
      }
    }
    data_pvalues <- data.frame(replicate = names(p_values), p_value = unlist(p_values), stringsAsFactors = FALSE)
    
    
    # Create the plot
    plot <- data_limited %>%
      ggplot(aes(x = t50_cat, y = value)) + 
      
      geom_violin() +
      geom_dotplot(stackdir = "centerwhole",stackgroups = TRUE,
                   binpositions = "all", binaxis = "y", binwidth = 0.01, size = 0.01) + 
      geom_hline(aes(yintercept = cutoff_50), linetype = "dashed", color = "red") +
      geom_hline(aes(yintercept = mean_upper), linetype = "dashed", color = "blue") +
      geom_hline(aes(yintercept = mean_lower), linetype = "dashed", color = "blue") +
      stat_summary(fun = mean, geom = "point", color = "red", size = 3, aes(group = t50_cat)) +
      geom_text(data = data_pvalues, 
                aes(label = sprintf("p = %e", p_value), x = 1.5, y = Inf), hjust = 1.1, vjust = 2, check_overlap = TRUE) +
      facet_wrap(~replicate) +
      
      labs(x = "t50 Category", y = "Value")
    
  }
  
  
  
  return(plot)
}
