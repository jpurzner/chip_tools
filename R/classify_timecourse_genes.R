#' Classify time-course genes by when they peak and when they turn on
#'
#' The central gene classifier of the GNP/MB time-course analysis. Cuts genes
#' into developmental windows on two axes — the time of their expression
#' maximum (`max_cat`) and their half-maximum time or t50 (`t50_cat`) — then
#' overrides both to `"unchanged"` for genes whose fold change or p-value does
#' not clear the thresholds, and to `"unexpressed"` for genes below the
#' expression floor in both data sets.
#'
#' The override order matters and is deliberate: `"unexpressed"` is applied
#' last, so a gene that is both unexpressed and unchanging is reported as
#' unexpressed.
#'
#' Four copies of this lived in the `transcript_timecourse*.Rmd` notebooks,
#' with the cut points and column names hard-coded. They are arguments here.
#'
#' @param df Data frame with one row per gene and the columns named by the
#'   `*_col` arguments below.
#' @param fc_cutoff Genes with `avg_fc` below this become `"unchanged"`.
#' @param p_val_cutoff Significance threshold applied to both p-value columns.
#' @param expr_cutoffs Named numeric vector of per-data-set expression floors.
#'   A gene below *every* floor becomes `"unexpressed"`. Names must be columns
#'   of `df`.
#' @param pval_both If `TRUE` a gene must clear `p_val_cutoff` in **both**
#'   data sets to be called changing (`|` on the failure condition). If
#'   `FALSE` (default) it need only clear it in one (`&` on the failure
#'   condition). This is the single most consequential argument here and the
#'   notebooks flipped it between runs.
#' @param max_breaks,max_labels Cut points and labels for the time of maximum.
#'   `max_breaks` is passed to [cut()] with `right = FALSE`, so it needs one
#'   more element than `max_labels`.
#' @param t50_breaks,t50_labels Cut points and labels for t50.
#' @param max_col,t50_col Columns holding the time of maximum and the t50.
#' @param fc_col Column holding the average fold change.
#' @param pval_cols Character vector of p-value columns.
#' @param table_only If `TRUE` return the counts per combined category rather
#'   than the classified data frame.
#'
#' @return `df` with factor columns `max_cat`, `t50_cat` and `max_t50_cat`
#'   added — or, with `table_only = TRUE`, a named integer vector of counts
#'   per `max_t50_cat` level including the empty ones.
#'
#' @section Note:
#' Factor levels for `"unchanged"` and `"unexpressed"` are added up front, so
#' `table()` on the result reports zero-count categories rather than dropping
#' them. That is what makes `table_only = TRUE` safe to `rbind` across
#' parameter sweeps, which is how the notebooks chose their cutoffs.
#'
#' @seealso [expr_half_max()] to compute t50, [time_of_change()] for the
#'   underlying rise/fall description, [t50window()] to bin the result
#'
#' @examples
#' df <- data.frame(
#'   avg_max = c(1, 5, 30, 10),
#'   avg_t50 = c(2, 10, 25, 20),
#'   avg_fc = c(3, 4, 5, 0.2),
#'   scott_gnp_min_pval = c(0.001, 0.001, 0.001, 0.5),
#'   hatten_trap_min_pval = c(0.001, 0.4, 0.001, 0.5),
#'   max_expr_Scott = c(100, 100, 100, 100),
#'   max_expr_Hatten = c(100, 100, 100, 100)
#' )
#' classify_timecourse_genes(
#'   df,
#'   fc_cutoff = 2, p_val_cutoff = 0.05,
#'   expr_cutoffs = c(max_expr_Scott = 5, max_expr_Hatten = 5)
#' )[, c("max_cat", "t50_cat", "max_t50_cat")]
#' @export
classify_timecourse_genes <- function(df,
                                      fc_cutoff,
                                      p_val_cutoff,
                                      expr_cutoffs,
                                      pval_both = FALSE,
                                      max_breaks = c(0, 2, 20, 64),
                                      max_labels = c("E15", "P0_7", "P14_P56"),
                                      t50_breaks = c(0, 6, 18, 64),
                                      t50_labels = c("Early", "Mid", "Late"),
                                      max_col = "avg_max",
                                      t50_col = "avg_t50",
                                      fc_col = "avg_fc",
                                      pval_cols = c("hatten_trap_min_pval",
                                                    "scott_gnp_min_pval"),
                                      table_only = FALSE) {

  needed <- c(max_col, t50_col, fc_col, pval_cols, names(expr_cutoffs))
  absent <- setdiff(needed, names(df))
  if (length(absent)) {
    stop("df is missing column(s): ", paste(absent, collapse = ", "),
         call. = FALSE)
  }
  if (is.null(names(expr_cutoffs)) || any(names(expr_cutoffs) == "")) {
    stop("expr_cutoffs must be a named vector, names being columns of df",
         call. = FALSE)
  }
  stopifnot(length(max_breaks) == length(max_labels) + 1L,
            length(t50_breaks) == length(t50_labels) + 1L)

  new_df <- df
  extra <- c("unchanged", "unexpressed")

  new_df$max_cat <- cut(df[[max_col]], max_breaks, right = FALSE,
                        labels = max_labels)
  new_df$t50_cat <- cut(df[[t50_col]], t50_breaks, right = FALSE,
                        labels = t50_labels)
  new_df$max_t50_cat <- paste0("Max_", new_df$max_cat, " t50_", new_df$t50_cat)

  # Declare the override levels before assigning them, so table() on the
  # result reports zero-count categories instead of silently dropping them.
  new_df$max_cat <- factor(new_df$max_cat, levels = c(max_labels, extra))
  new_df$t50_cat <- factor(new_df$t50_cat, levels = c(t50_labels, extra))
  all_combos <- paste0("Max_", rep(max_labels, each = length(t50_labels)),
                       " t50_", rep(t50_labels, times = length(max_labels)))
  new_df$max_t50_cat <- factor(new_df$max_t50_cat,
                               levels = c(all_combos, extra))

  # A gene fails the significance screen if its fold change is too small, or
  # if its p-values miss the threshold. `pval_both` picks which.
  pvals <- as.matrix(df[, pval_cols, drop = FALSE])
  pval_fail <- if (pval_both) {
    apply(pvals > p_val_cutoff, 1, any)
  } else {
    apply(pvals > p_val_cutoff, 1, all)
  }
  unchanged <- df[[fc_col]] < fc_cutoff | pval_fail

  # Unexpressed only if below the floor in EVERY data set.
  below <- vapply(names(expr_cutoffs),
                  function(col) df[[col]] < expr_cutoffs[[col]],
                  logical(nrow(df)))
  if (nrow(df) == 1L) below <- matrix(below, nrow = 1L)
  unexpressed <- apply(below, 1, all)

  for (col in c("max_cat", "t50_cat", "max_t50_cat")) {
    new_df[[col]][unchanged]   <- "unchanged"
    new_df[[col]][unexpressed] <- "unexpressed"
  }

  if (table_only) {
    counts <- table(new_df$max_t50_cat)
    return(stats::setNames(as.vector(counts), names(counts)))
  }
  new_df
}


#' Compress part of a continuous axis
#'
#' A scales transformation that divides the span between `from` and `to` by
#' `factor`, leaving everything outside it linear. Use it when most of the
#' data sits at one end of the range but a few points far out still have to be
#' shown — the alternative, a log axis, distorts the dense region too.
#'
#' @param from,to Bounds of the region to compress.
#' @param factor How much to compress by.
#'
#' @return A transformation object for the `trans` argument of a ggplot2
#'   continuous scale.
#'
#' @source The widely circulated `squish_trans` recipe; used in
#'   `3_gene_classification.Rmd`.
#'
#' @examples
#' \dontrun{
#' ggplot(df, aes(x, y)) +
#'   geom_point() +
#'   scale_y_continuous(trans = squish_trans(10, 100, 20))
#' }
#' @export
squish_trans <- function(from, to, factor) {
  require_pkg("scales")

  trans <- function(x) {
    if (any(is.na(x))) return(x)
    isq <- x > from & x < to
    ito <- x >= to
    x[isq] <- from + (x[isq] - from) / factor
    x[ito] <- from + (to - from) / factor + (x[ito] - to)
    x
  }

  inv <- function(x) {
    if (any(is.na(x))) return(x)
    isq <- x > from & x < from + (to - from) / factor
    ito <- x >= from + (to - from) / factor
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from) / factor))
    x
  }

  scales::trans_new("squished", trans, inv)
}


#' Map mouse gene symbols to human orthologs
#'
#' Queries Ensembl BioMart for the human orthologs of a set of MGI symbols.
#' Used in the scRNA-seq and human-MB notebooks to carry a mouse-derived gene
#' set onto human data.
#'
#' @param x Character vector of MGI symbols.
#' @param unique_rows Passed to the BioMart query as `uniqueRows`.
#'
#' @return A data frame of `mgi_symbol` and `hgnc_symbol`.
#'
#' @section Note:
#' This needs network access to Ensembl and is slow; cache the result rather
#' than calling it in a loop. `biomaRt::getLDS()` has been removed from recent
#' biomaRt, so on a current install this errors — the replacement is a
#' `getBM()` call against the `hsapiens_homolog_associated_gene_name`
#' attribute. Kept as the notebooks had it, with the failure made explicit
#' rather than silent.
#'
#' The notebook version printed `head()` of the result and returned the full
#' two-column frame despite computing a unique human-symbol vector first.
#'
#' @examples
#' \dontrun{
#' orthologs <- convert_mouse_gene_list(c("Gli1", "Ccnd1", "Mycn"))
#' }
#' @export
convert_mouse_gene_list <- function(x, unique_rows = FALSE) {
  require_pkg("biomaRt")

  if (!"getLDS" %in% getNamespaceExports("biomaRt")) {
    stop("biomaRt::getLDS() is not available in biomaRt ",
         as.character(utils::packageVersion("biomaRt")),
         ". Use biomaRt::getBM() against the ",
         "'hsapiens_homolog_associated_gene_name' attribute instead.",
         call. = FALSE)
  }

  human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  getLDS <- getExportedValue("biomaRt", "getLDS")
  getLDS(attributes = c("mgi_symbol"),
         filters = "mgi_symbol",
         values = x,
         mart = mouse,
         attributesL = c("hgnc_symbol"),
         martL = human,
         uniqueRows = unique_rows)
}
