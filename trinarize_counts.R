trinarize_counts <- function (counts, log_data = FALSE, exclude_extreme = FALSE, 
                              return_cutoff = FALSE, return_mean = FALSE, 
                              return_cutoff_mean = FALSE, label_name = NULL) {
  require(mixtools)
  require(ggplot2)
  
  if (log_data) {
    x <- log10(counts + 1)
  } else { 
    x <- counts
  }
  
  if (exclude_extreme) {
    x <- x[x > 0.1 & x < 1]
  }
  
  model <- normalmixEM(x = x, k = 3, arbvar = FALSE, maxit = 10000, maxrestarts = 200, epsilon = 1e-20)
  
  # Assuming ordered by mean for simplicity
  ordered_indices <- order(model$mu)
  cutoff_low <- (model$mu[ordered_indices[1]] + model$mu[ordered_indices[2]]) / 2
  cutoff_high <- (model$mu[ordered_indices[2]] + model$mu[ordered_indices[3]]) / 2
  
  # For plotting
  x <- as.data.frame(x)
  
  norm_funct <- function(x, lambda, mu, sigma) {
    return(lambda * dnorm(x, mean = mu, sd = sigma))
  }
  
  # Adding density functions for each component
  x$nd1 <- norm_funct(x$x, model$lambda[1], model$mu[1], model$sigma[1])
  x$nd2 <- norm_funct(x$x, model$lambda[2], model$mu[2], model$sigma[2])
  x$nd3 <- norm_funct(x$x, model$lambda[3], model$mu[3], model$sigma[3])  # Third component
  

  
  # Adding total mixture density for plotting
  #x$mix_density <- approx(x_range, total_density, x$x)$y
  
  

  
  # Plot
  m <- ggplot(x) + 
    geom_histogram(aes(x = x, y = ..density..), binwidth = 0.01, color = "black", fill = "grey", alpha = 0.1) +
    geom_line(aes(x = x, y = nd1), color = "dodgerblue1", size = 1, alpha = 0.8) +
    geom_line(aes(x = x, y = nd2), color = "seagreen4", size = 1, alpha = 0.8) +
    geom_line(aes(x = x, y = nd3), color = "red", size = 1, alpha = 0.8) +  # Third component line
    geom_vline(xintercept = cutoff_low, linetype = "dotdash", color = "purple", alpha = 0.5) +
    geom_vline(xintercept = cutoff_high, linetype = "dotdash", color = "orange", alpha = 0.5) +  # Second cutoff
    theme(axis.title.y = element_text(size = 16, angle = 90),
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    xlab("log10 + 1 counts") + 
    ylab("density")
  
  if (!is.null(label_name)) {
    m <- m + ggtitle(label_name)
  }
  
  print(m)
  
  return(model)
  # No changes below this line, handle return logic as in your original function...
  
  # This placeholder is where you'd continue with the logic to return cutoffs, binarize the data, etc.,
  # similar to how you've handled it in your original 'binarize_counts' function.
}
