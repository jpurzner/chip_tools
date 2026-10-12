#' @keywords internal
"_PACKAGE"

# Packages in Imports are imported wholesale so the bare verb names used
# throughout the older functions resolve at call time. Declaring them in
# utils::globalVariables() (see globals.R) silences R CMD check but does NOT
# make them callable -- that distinction cost a release.
#
# Anything in Suggests cannot be imported here and must be called with an
# explicit pkg:: prefix at the call site, behind a require_pkg() guard.
#
# dplyr is imported wholesale but plyr is NOT: they collide on mutate,
# arrange, desc, summarise, summarize, rename and count, and the code here
# wants dplyr's semantics. Only plyr's own two functions are imported.
#
#' @import dplyr
#' @import ggplot2
#' @importFrom plyr ddply mapvalues
#' @importFrom reshape2 melt dcast
#' @importFrom tidyr pivot_longer pivot_wider separate_rows unnest
#' @importFrom grDevices chull colorRampPalette hcl hcl.colors
#' @importFrom graphics hist
#' @importFrom grid gpar grid.draw grid.newpage rectGrob unit
#' @importFrom methods new
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom rlang .data sym :=
#' @importFrom stats anova approx approxfun as.dist complete.cases density
#'   dgamma dist dnorm hclust kmeans loess loess.control median na.exclude
#'   na.omit predict quantile runif sd setNames spline time uniroot var
#'   wilcox.test
#' @importFrom utils head read.table
#' @importFrom zoo na.approx rollapply rollmean
NULL
