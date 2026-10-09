#' ROC curve, or several, as a ggplot
#'
#' Draws one or more pROC `roc` objects as step curves with the chance diagonal.
#' A named list gives one coloured curve per element, with underscores in the
#' names turned into spaces for the legend.
#'
#' @param roc A pROC `roc` object, or a named list of them.
#' @param showAUC Currently unused; retained for call compatibility.
#' @param interval Spacing of axis breaks.
#' @param breaks Explicit axis breaks, overriding `interval`.
#'
#' @return A `ggplot` object.
#'
#' @section Fixes:
#' Tested `class(roc) == "roc"`, which warns (and in older R silently took the
#' first element) on objects with several classes; now uses `inherits()`. The
#' bad-input branch constructed a condition object without signalling it, so
#' invalid input fell through to an error about `roc_df` instead. The plot was
#' also assigned over the function's own name.
#'
#' @examples
#' \dontrun{
#' ggroc(list(H3K4me3 = roc1, H3K27me3 = roc2))
#' }
#' @import ggplot2
#' @export
ggroc <- function(roc, showAUC = TRUE, interval = 0.2, breaks = seq(0, 1, interval)){
  require_pkg("pROC")

  
  if (inherits(roc, "roc")) {
    roc_df = data.frame(x = rev(roc$specificities), y = rev(roc$sensitivities), name = "ROC")
  } else if (is.list(roc)) {
    roc_df_list <- lapply(c(1:length(roc)), function (n) { df = data.frame(x = rev(roc[[n]]$specificities), 
                                                            y = rev(roc[[n]]$sensitivities), 
                                                            name = gsub("_", " ", names(roc)[n]))
                                            return(df)
    })
    roc_df <- dplyr::bind_rows(roc_df_list)
  } else {
    stop("`roc` must be a roc object from pROC, or a named list of them",
         call. = FALSE)
  }
  
  p <- ggplot(roc_df, aes(x = x, y = y, group= name, color = name)) +
    geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0), alpha = 0.5, color = "black") + 
    geom_step() +
    scale_x_reverse(name = "Specificity",limits = c(1,0), breaks = breaks, expand = c(0.001,0.001)) + 
    scale_y_continuous(name = "Sensitivity", limits = c(0,1), breaks = breaks, expand = c(0.001, 0.001)) +
    theme_bw() + 
    theme(axis.ticks = element_line(color = "grey80")) +
    coord_equal() + 
    theme(legend.position = c(0.8, 0.2), 
          legend.title = element_blank())
  p
}