enhancer_quant <- function(enhancer_list, promoter_df, dist_cutoff = 1e6, linear_cutoff = 50e3, par = FALSE) {
  
  # TODO sporatic error with parallel still 
  
  require(plyr)
  require(dplyr)
  require(reshape2)
  
  # promoter_df format 
  #    chr   start    stop               name strand
  #    1 chr1 3661578 3661579 ENSMUST00000070533      -
  #    2 chr1 4350472 4350473 ENSMUST00000027032      -
  #    3 chr1 4486493 4486494 ENSMUST00000027035      -
  #    4 chr1 4486443 4486444 ENSMUST00000116652      -
  #    5 chr1 4775790 4775791 ENSMUST00000130201      -
  #    6 chr1 4775819 4775820 ENSMUST00000156816      -
  # 
  #  enhancer_list format 
  #      chr   start    stop        name counts length mid_point      norm
  #  1 chr1 6499300 6501799 epic_peak11    100   2499   6500550 0.8264463
  #  2 chr1 6512100 6518499 epic_peak12    182   6399   6515300 1.5041322
  #  3 chr1 6640800 6643799 epic_peak13     68   2999   6642300 0.5619835
  #  4 chr1 6738800 6747999 epic_peak15    304   9199   6743400 2.5123967
  #  5 chr1 6748700 6753399 epic_peak16    147   4699   6751050 1.2148760
  #  6 chr1 9638100 9641399 epic_peak24    169   3299   9639750 1.3966942
  #
  # par default true, run as parallel, requires package parallel 
  
  
  
  # there are not that many combinations of enhancers to promoters so will test all combinations 
  # not efficient !
  
  # split by chromosome
  
  promoter_df$chr <- as.factor(promoter_df$chr)   
  promoter_df$name <- as.character(promoter_df$name)
  
  
  enh_dist_summary <- function (chr_num, promoter, enhancer, dist_cutoff = 1e6) {
    promoter <- subset(promoter, chr == chr_num)
    enhancer <- subset(enhancer, chr == chr_num)
    # compare each promoter to each enhancer and obtain a distance  
    pe_dist <- outer(promoter$stop, enhancer$mid_point, `-`) 
    row.names(pe_dist) <- promoter$name
    colnames(pe_dist) <- enhancer$name 
    pe_dist <- melt(pe_dist)
    colnames(pe_dist) <- c("promoter", "enhancer", "dist")
    pe_dist$dist <- abs(pe_dist$dist)
    pe_dist <- subset(pe_dist, dist <= dist_cutoff)
    pe_dist$promoter <- as.character(pe_dist$promoter)
    pe_dist$enhancer <- as.character(pe_dist$enhancer)
    return(pe_dist)
  }
  
  scale_quant_bydist <- function(quant, dist, linear_cutoff = 50e3,  dist_cutoff = 1e6) {
    if (dist <= linear_cutoff) {
      return(quant)
    } else if ( dist > linear_cutoff & dist <= dist_cutoff )  {
      #return(quant * ((-dist/(dist_cutoff - linear_cutoff)) + (dist_cutoff/(dist_cutoff-linear_cutoff))))
      return(quant * ((-dist + dist_cutoff)/(dist_cutoff - linear_cutoff)))
    } else {
      return(0)
      
    }
  }
  
  enh_dist_quant <- function(promoter, enhancer, dist_cutoff = 1e6) {  
    pe_sum <- bind_rows(lapply(levels( promoter$chr), function (ch) enh_dist_summary(ch,promoter,enhancer)))
    # map enhancer name back to quantification 
    pe_sum <- as.data.frame(pe_sum)
    enhancer <- as.data.frame(enhancer) 
    pe_sum$quant <- suppressWarnings(plyr::mapvalues(pe_sum$enhancer, enhancer$name, enhancer$norm))
    pe_sum$quant <- as.numeric(pe_sum$quant)
    pe_sum$promoter <- as.factor( pe_sum$promoter)
    # summarize by promoter name
    pe_prom <- pe_sum %>%  
      dplyr::group_by(promoter) %>% 
      dplyr::summarise(enhancer_quant = sum(mapply(function (q,d) 
        scale_quant_bydist(q, d, linear_cutoff = linear_cutoff, dist_cutoff = dist_cutoff), quant,dist)))
    #dplyr::summarise(enhancer_quant = mapply(function (q,d) sum(quant*dist))
    return(pe_prom)
  }
  
  if (par) {
    library(parallel)
    # Calculate the number of cores
    no_cores <- detectCores() - 1
    # Initiate cluster
    cl <- makeCluster(no_cores)
    #clusterExport(cl, c("plyr", "dplyr", "base"))
    clusterEvalQ(cl, {
      library(plyr)
      library(dplyr)
      library(reshape2)
    })
    enh_quant_list <- parLapply(cl,enhancer_list, function(x) enh_dist_quant(promoter_df, x, dist_cutoff))
    stopCluster(cl)
  } else {  
    enh_quant_list <- lapply(enhancer_list, function(x) enh_dist_quant(promoter_df, x, dist_cutoff))
  } 
  
  
  # make all dfs the same length and order 
  
  enh_quant_list_same <- lapply(enh_quant_list, function(d) {
    d$promoter <- as.character(d$promoter)
    df <- left_join( promoter_df[,c(1:4)], d, by = c("name" = "promoter") )
    df <- df[,c("name","enhancer_quant")]
    row.names(df) <- df$name
    df$name <- NULL
    df$enhancer_quant[is.na(df$enhancer_quant)] <- 0  
    return(df)
  })
  
  
  names(enh_quant_list_same) <- names(enhancer_list)
  #return(enh_quant_list_same)
  enh_quant <- dplyr::bind_cols(enh_quant_list_same)
  colnames(enh_quant) <- str_c(names(enhancer_list) , "_enhancer") 
  row.names(enh_quant) <- row.names(enh_quant_list_same[[1]])
  
  return(enh_quant)
}
