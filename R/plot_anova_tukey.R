#' Faceted plots annotated with ANOVA and Tukey HSD results
#'
#' Fits a two-way ANOVA within each combination of two grouping variables,
#' runs Tukey's HSD on `comparison_var`, and returns one annotated plot per
#' level of `group_var1` with significance brackets.
#'
#' @param data A data frame in long form.
#' @param response_var Name of the measured column.
#' @param group_var1,group_var2 Names of the two grouping columns; a separate
#'   model is fitted within each combination.
#' @param comparison_var Name of the column whose levels are compared.
#' @param y_label Y-axis label.
#' @param title_prefix String prepended to each panel title.
#' @param block_var Name of the blocking covariate added to the ANOVA formula.
#'   This was hard-coded to `cell_line`, so the function only worked on data
#'   with a column of that name.
#'
#' @return A list of `ggplot` objects, one per level of `group_var1`.
#'
#' @section Note:
#' Uses `dplyr::do()`, which has been superseded by `group_modify()`/`nest()`.
#' It still works but is a deprecation warning waiting to happen.
#'
#' @examples
#' \dontrun{
#' plots <- plot_anova_tukey(measurements, "intensity", "drug", "timepoint",
#'                           "dose", y_label = "Intensity",
#'                           block_var = "cell_line")
#' }
#' @import ggplot2
#' @export
# Define a function to create plots with ANOVA and Tukey HSD results
plot_anova_tukey <- function(data, response_var, group_var1, group_var2,
                             comparison_var, y_label, title_prefix = "",
                             block_var = NULL) {
  
  # Perform ANOVA for each group
  anova_results <- data %>%
    dplyr::group_by(!!rlang::sym(group_var1), !!rlang::sym(group_var2)) %>%
    dplyr::do(anova = stats::aov(
        stats::as.formula(paste(response_var, "~", comparison_var,
                                if (is.null(block_var)) "" else paste("+", block_var))),
        data = .)) %>%
    dplyr::ungroup()
  
  # Apply Tukey's HSD test to the ANOVA results
  tukey_results <- anova_results %>%
    dplyr::mutate(tukey = lapply(anova, function(x) {
      tukey_df <- as.data.frame(stats::TukeyHSD(x, comparison_var)[[comparison_var]])
      tukey_df$pair <- rownames(tukey_df)  # Add the pair information
      tukey_df
    })) %>%
    dplyr::select(-anova) %>%
    tidyr::unnest(tukey)
  
  
  # Process Tukey results and format for plotting
  tukey_df <- tukey_results %>%
    dplyr::group_by(!!rlang::sym(group_var1), !!rlang::sym(group_var2)) %>%
    dplyr::mutate(
      group1 = as.factor(as.numeric(sub("-.*", "", pair))),
      group2 = as.factor(as.numeric(sub(".*-", "", pair))),
      p.adj = signif(`p adj`, 2),
      p_label = ifelse(p.adj < 0.001, "p < 0.001", paste("p =", p.adj)),
      y.position = max(data[[response_var]]) + (0.05 * seq_along(pair))
    ) %>%
    dplyr::mutate(
      x = (as.numeric(as.character(group1)) + (as.numeric(as.character(group2)) - as.numeric(as.character(group1))) / 2),
      y = y.position
    )
  
  # Get unique levels of the group_var1 (e.g., cell_cycle_drug)
  unique_groups <- unique(data[[group_var1]])
  
  # Create an empty list to store the plots
  plot_list <- list()
  
  # Loop over each group and create the plot dynamically
  for (group in unique_groups) {
    # Filter the data for the current group
    data_filtered <- data %>% dplyr::filter(!!rlang::sym(group_var1) == group)
    tukey_filtered <- tukey_df %>% dplyr::filter(!!rlang::sym(group_var1) == group)
    
    # Create the plot for the current group
    plot <- data_filtered %>%
      ggplot(aes_string(x = comparison_var, y = response_var)) + 
      geom_point(aes(shape = cell_line), position = position_jitter(width = 0.1, height = 0), size = 2, alpha = 0.6) + 
      geom_point(data = data_filtered %>%
                   dplyr::group_by(!!rlang::sym(comparison_var)) %>%
                   dplyr::summarise(mean_response = mean(!!rlang::sym(response_var))), 
                 aes_string(x = comparison_var, y = "mean_response"), 
                 colour = "black", shape = 5, size = 2) +
      geom_line(data = data_filtered %>%
                  dplyr::group_by(!!rlang::sym(comparison_var)) %>%
                  dplyr::summarise(mean_response = mean(!!rlang::sym(response_var))), 
                aes_string(x = comparison_var, y = "mean_response"), 
                size = 1) + 
      # Add Tukey HSD p-values for the current group
      stat_pvalue_manual(tukey_filtered, 
                         label = "p_label",   
                         x = "x",             
                         y.position = "y",    
                         tip.length = 0.02, 
                         bracket.size = 0.5,
                         size = 4, 
                         color = "black") +
      ylim(0, max(data[[response_var]]) * 1.2) +  # Adjust y-limit dynamically
      labs(x = paste("Concentration of", comparison_var), 
           y = y_label, 
           title = paste0(title_prefix, " ", group)) +  # Use the group name dynamically as the title
      theme_classic() +
      theme(text = element_text(size = 18))
    
    # Store the plot in the list
    plot_list[[group]] <- plot
  }
  
  # Combine all the plots using patchwork
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = 2)  # Adjust ncol as needed
  
  # Return the combined plot
  return(combined_plot)
}
