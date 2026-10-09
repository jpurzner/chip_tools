#' Summarise a variable by group with SE and confidence interval
#'
#' Counts, mean, standard deviation, standard error and a confidence interval
#' for one measure across any number of grouping variables. Handy for building
#' error bars for ggplot2.
#'
#' @param data A data frame.
#' @param measurevar Name of the column to summarise.
#' @param groupvars Character vector of grouping column names.
#' @param na.rm If `TRUE`, drop `NA` before summarising (and do not count them
#'   in `N`).
#' @param conf.interval Width of the confidence interval.
#' @param .drop Passed to [plyr::ddply()]; drop unused factor combinations.
#'
#' @return A data frame with one row per group and columns `N`, the measure
#'   name (holding the mean), `sd`, `se` and `ci`.
#'
#' @source Winston Chang, *Cookbook for R*,
#'   <http://www.cookbook-r.com/Manipulating_data/Summarizing_data/>.
#'   Included verbatim apart from namespacing the plyr calls.
#'
#' @examples
#' summarySE(ToothGrowth, measurevar = "len", groupvars = c("supp", "dose"))
#' @export


## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = stats::sd(xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- stats::qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}