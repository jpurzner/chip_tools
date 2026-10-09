#' Average columns within levels of a factor
#'
#' Collapses a features x samples matrix to features x groups by taking the
#' row mean within each level of `meta`.
#'
#' This is the lighter-weight sibling of [grouped_col_mean()]: it takes the
#' grouping factor directly rather than a metadata data frame plus a column
#' name. Copy-pasted into four of the `1_transcript_timecourse*.Rmd`
#' notebooks before being lifted here.
#'
#' @param df Features x samples data frame.
#' @param meta Factor of length `ncol(df)` giving each column's group. Columns
#'   come back in `levels(meta)` order, so `meta` is also how you set the
#'   output column order — relevel it rather than reordering afterwards.
#' @param na.rm Passed to [rowMeans()].
#'
#' @return A data frame with one column per level of `meta`, carrying the row
#'   names of `df`.
#'
#' @seealso [grouped_col_mean()] when the grouping lives in a metadata frame.
#'
#' @examples
#' df <- data.frame(a1 = c(1, 2), a2 = c(3, 4), b1 = c(10, 20))
#' average_df(df, factor(c("a", "a", "b")))
#' @export
average_df <- function(df, meta, na.rm = TRUE) {
  meta <- as.factor(meta)
  if (length(meta) != ncol(df)) {
    stop("meta has length ", length(meta), " but df has ", ncol(df),
         " columns; they must correspond one to one", call. = FALSE)
  }

  # The notebook version used rowMeans() with no na.rm and no empty-level
  # guard, so a dropped factor level produced a zero-column subset and
  # rowMeans() errored.
  out <- lapply(levels(meta), function(lv) {
    cols <- which(meta == lv)
    if (length(cols) == 0L) {
      return(rep(NA_real_, nrow(df)))
    }
    rowMeans(as.matrix(df[, cols, drop = FALSE]), na.rm = na.rm)
  })

  out <- as.data.frame(out, stringsAsFactors = FALSE)
  colnames(out) <- levels(meta)
  rownames(out) <- rownames(df)
  out
}


#' The default ggplot2 discrete colour palette
#'
#' Reproduces the hues ggplot2 assigns to `n` discrete levels, for when a
#' base-graphics panel or a `pheatmap` annotation has to sit beside a ggplot
#' and match its colours.
#'
#' @param n Number of colours.
#'
#' @return Character vector of `n` hex colours.
#'
#' @source The widely circulated `gg_color_hue` recipe; appears in both
#'   `5_MB_histones*.Rmd` notebooks.
#'
#' @examples
#' gg_color_hue(3)
#' @export
gg_color_hue <- function(n) {
  stopifnot(length(n) == 1L, n >= 1)
  hues <- seq(15, 375, length.out = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}


#' Binarize the value column of a long-format data frame
#'
#' Convenience wrapper that runs [binarize_counts()] on `df$value` and writes
#' the 0/1 call back as `df$group`. Written for use inside a
#' `split()`/`lapply()` over a long-format signal table, which is how the
#' `4_promoter_histones*.Rmd` notebooks call it — six copies of this three-line
#' function across them.
#'
#' @param df Data frame with a numeric `value` column.
#' @param ... Passed to [binarize_counts()], e.g. `log_data`, `mu`, `plot`.
#'
#' @return `df` with an integer `group` column added.
#'
#' @seealso [binarize_counts()]
#'
#' @examples
#' \dontrun{
#' long_signal |>
#'   split(~replicate) |>
#'   lapply(call_mix_bin, plot = FALSE) |>
#'   dplyr::bind_rows()
#' }
#' @export
call_mix_bin <- function(df, ...) {
  if (!"value" %in% names(df)) {
    stop("df needs a 'value' column", call. = FALSE)
  }
  bin_res <- binarize_counts(df$value, ...)
  df$group <- bin_res$bin
  df
}


#' Percentile rank of each integer value in a vector
#'
#' For every integer from `from` to `to`, the percentage of `vector` at or
#' below it. Used in `6_Human_MB.Rmd` to place a mouse-derived score on the
#' human expression distribution.
#'
#' @param vector Numeric vector.
#' @param from,to Integer range to evaluate over; default to the range of
#'   `vector`.
#' @param round Decimal places for the returned percentage.
#'
#' @return A data frame of `value` and `perc.rank`.
#'
#' @examples
#' percentile_ranks(c(1, 2, 2, 3, 8))
#' @export
percentile_ranks <- function(vector, from = NA, to = NA, round = 1) {
  vector <- stats::na.omit(vector)
  if (length(vector) == 0L) {
    stop("vector has no non-missing values", call. = FALSE)
  }
  if (is.na(from)) from <- floor(min(vector))
  if (is.na(to))   to   <- ceiling(max(vector))

  vals <- seq(from, to)
  data.frame(
    value = vals,
    perc.rank = round(100 * vapply(vals,
                                   function(x) mean(vector <= x),
                                   numeric(1)), round)
  )
}


#' A z-scored metagene score across samples
#'
#' Z-scores each gene across samples, then averages over a gene set to give
#' one score per sample. Scoring per gene first stops highly expressed genes
#' dominating the set. From `survival_analysis.Rmd`, where the metagene is the
#' covariate in a Cox model.
#'
#' @param gene_ids Gene identifiers to include.
#' @param expr_df Data frame with an identifier column and one column per
#'   sample.
#' @param sample_ids Column names of `expr_df` holding expression.
#' @param id_col Name of the identifier column in `expr_df`.
#'
#' @return Named numeric vector, one metagene score per entry of `sample_ids`.
#'
#' @section Note:
#' The notebook version hard-coded the identifier column as
#' `EnsemblGeneGeneID_from_ensemblv77`; it is the `id_col` argument here.
#' Genes with zero variance across samples z-score to `NaN` and are dropped
#' from the mean by `na.rm`.
#'
#' @examples
#' expr <- data.frame(gene = c("g1", "g2", "g3"),
#'                    s1 = c(1, 5, 9), s2 = c(2, 4, 8), s3 = c(9, 1, 2))
#' compute_metagene(c("g1", "g2"), expr, c("s1", "s2", "s3"), id_col = "gene")
#' @export
compute_metagene <- function(gene_ids, expr_df, sample_ids, id_col) {
  if (!id_col %in% names(expr_df)) {
    stop(sQuote(id_col), " is not a column of expr_df", call. = FALSE)
  }
  missing_cols <- setdiff(sample_ids, names(expr_df))
  if (length(missing_cols)) {
    stop("expr_df is missing sample column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  keep <- expr_df[[id_col]] %in% gene_ids
  if (!any(keep)) {
    stop("none of gene_ids are present in expr_df[[id_col]]", call. = FALSE)
  }

  mat <- as.matrix(expr_df[keep, sample_ids, drop = FALSE])
  mat_z <- t(scale(t(mat)))
  colMeans(mat_z, na.rm = TRUE)
}


#' Flag GO rows containing any gene from a set
#'
#' Adds a logical column marking GO terms whose comma-separated `genes` string
#' includes at least one gene of interest — the quick way to ask which
#' enriched terms are driven by your differentiation set. From
#' `7_GNP_diff_go_terms*.Rmd`.
#'
#' @param go_data GO result data frame with a `genes` column of
#'   comma-separated symbols.
#' @param gene_set Character vector of genes to look for.
#' @param sep Separator splitting the `genes` column.
#'
#' @return `go_data` with a logical `is_diff_gene` column, plus an integer
#'   `n_diff_gene` column counting the matches.
#'
#' @seealso [go_to_gene_table()], [topGO_wrap()]
#'
#' @examples
#' go <- data.frame(GO.ID = c("GO:1", "GO:2"),
#'                  genes = c("Gli1, Ccnd1", "Actb, Gapdh"))
#' check_diff_genes(go, c("Gli1", "Mycn"))
#' @export
check_diff_genes <- function(go_data, gene_set, sep = ",\\s*") {
  if (!"genes" %in% names(go_data)) {
    stop("go_data needs a 'genes' column", call. = FALSE)
  }

  # The notebook version used rowwise() + any(), which is slow on a long GO
  # table and returned only the logical. vapply over the split list is both
  # faster and lets us report the count as well.
  split_genes <- strsplit(as.character(go_data$genes), sep)
  hits <- vapply(split_genes,
                 function(g) sum(trimws(g) %in% gene_set),
                 integer(1))

  go_data$is_diff_gene <- hits > 0L
  go_data$n_diff_gene <- hits
  go_data
}


#' Look up one cell of a per-condition count table
#'
#' Returns `table[table$condition == condition, category]`, or `NA` if that
#' condition is absent or duplicated. Appears in three of the
#' `4_promoter_histones*.Rmd` notebooks as a closure over a global
#' `category_counts`.
#'
#' @param condition Condition to look up.
#' @param category Column name to return.
#' @param counts_table Data frame with a `condition` column.
#'
#' @return A single value, or `NA` if the lookup is not unique.
#'
#' @section Note:
#' The notebook version read `category_counts` from the global environment,
#' so it only worked where that object happened to be bound. It is the
#' `counts_table` argument here.
#'
#' @examples
#' tab <- data.frame(condition = c("P7", "P56"), bound = c(10, 20))
#' lookup_count("P7", "bound", tab)
#' lookup_count("P14", "bound", tab)
#' @export
lookup_count <- function(condition, category, counts_table) {
  if (!"condition" %in% names(counts_table)) {
    stop("counts_table needs a 'condition' column", call. = FALSE)
  }
  if (!category %in% names(counts_table)) {
    stop(sQuote(category), " is not a column of counts_table", call. = FALSE)
  }
  row <- counts_table[counts_table$condition == condition, , drop = FALSE]
  if (nrow(row) != 1L) {
    return(NA)
  }
  row[[category]]
}
