binarize_counts_p <- function (counts, log_data = FALSE, exclude_extreme = FALSE, return_all = FALSE) {
  # take counts and seperate signal from background using mixtools with a mixture of normal distributions
  # bulk of code found at http://stats.stackexchange.com/questions/57993/how-to-explain-how-i-divided-a-bimodal-distribution-based-on-kernel-density-esti 
  # modified to use next gen sequencing count data
  
  
  
  require(mixtools)
  
  # log transform the count data 
  if (log_data) {
    x <- log10(counts + 1)
  } else { 
    x <- counts
  }
  
  if (exclude_extreme) {
    x <- counts[counts > 0 & counts < 1]
  }
  
  model <- gammamixEM(x= scott_max_expr_nona_noz$max_expr_Scott_l , k=2, alpha = c(12, 0.2), beta = c(0.04, 0.05), lambda =c(.7, .3), maxit =10000, maxrestarts = 10, epsilon = 1e-8)
  index.lower <- which.min(c(model$gamma.pars[1,1]*model$gamma.pars[2,1], model$gamma.pars[1,2]*model$gamma.pars[2,2])) # Index of component with lower mean
  
  
  find.cutoff <- function(proba=0.5, i=index.lower) {
    ## Cutoff such that Pr[drawn from bad component] == proba
    f <- function(x) {
      proba - (model$lambda[i]*dgamma(x, model$gamma.pars[1,i], 1/model$gamma.pars[2,i]) /
                 (model$lambda[1]*dgamma(x, model$gamma.pars[1,1], 1/model$gamma.pars[2,1]) + 
                    model$lambda[2]*dgamma(x, model$gamma.pars[1,2], 1/model$gamma.pars[2,2])))
  
    }
    # find minimum of f and set lower limit of uniroot
    low <- x[which.min(f(x))]
    up <- max(x)
    
    cutoff <- try(uniroot(f=f, lower=low, upper=up)$root)
    if (class(cutoff) == "try-error") {
      print("warning: failed to find root")
      cutoff <- (model$gamma.pars[1,1]*model$gamma.pars[2,1] +  model$gamma.pars[1,2]*model$gamma.pars[2,2])/2
    }
    
    return(cutoff)  # Careful with division by zero if changing lower and upper
  } # end of find cutoffs
  
  cutoffs <- c(find.cutoff(proba=0.5, i = index.lower), find.cutoff(proba=0.75, i = index.lower))  # Around c(1.8, 1.5)
  
  # plot using ggplot2 
  x <- as.data.frame(x)
  
  gamma_funct <- function(x, lambda, alpha, beta) {
    return(lambda * dgamma(x, shape = alpha, rate = 1/beta))
  }
  
  x$nd1 <- gamma_funct(x$x, model$lambda[1], model$gamma.pars[1,1], model$gamma.pars[2,1])
  x$nd2 <- gamma_funct(x$x, model$lambda[2], model$gamma.pars[1,2], model$gamma.pars[2,2])
  
  
  m <- ggplot(x) + geom_histogram(aes(x=x, y= ..density..), binwidth = 0.01, color = "black", fill = "grey", alpha = 0.1) 
  m <- m + geom_line(aes(x=x, y=nd1), color = "dodgerblue1", size = 1, alpha = 0.8) 
  m <- m + geom_line(aes(x=x, y=nd2), color = "seagreen4", size = 1, alpha = 0.8) 
  m <- m + geom_vline(xintercept = cutoffs[1], linetype = "dotdash", alpha = 0.5) 
  m <- m + theme(axis.title.y=element_text(size=16, angle = 90), axis.text.y = element_text(size=16), axis.text.x = element_text(size=16), axis.title.x = element_text(size = 16),  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) 
  m <- m + xlab("log10 + 1 counts") + ylab("density")
  print(m)
  
  
  
  
  # determine cutoff for un-transformed counts
  if (log_data) {
    cutoffs <- (10^cutoffs) - 1
  }
  
  
  # binarize the data based upon this cutoff
  results <- as.data.frame(counts)
  results$bin <- 0 
  results[results$counts >= cutoffs[2], 2] <- 1
  
  if (return_all)  {
    return(list(results, model, cutoffs))

  } else {
    return(results)  
  }
  
  
}# end of binarize_counts