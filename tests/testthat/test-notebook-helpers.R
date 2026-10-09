# Helpers lifted out of the Ezh2_2022 notebooks. These tests pin the
# behaviour the notebooks relied on, plus the guards added on extraction.

test_that("average_df averages within factor levels and respects level order", {
  df <- data.frame(a1 = c(1, 2), a2 = c(3, 4), b1 = c(10, 20))

  out <- average_df(df, factor(c("a", "a", "b")))
  expect_equal(colnames(out), c("a", "b"))
  expect_equal(out$a, c(2, 3))
  expect_equal(out$b, c(10, 20))

  # column order follows levels(meta), so releveling reorders the output
  reordered <- average_df(df, factor(c("a", "a", "b"), levels = c("b", "a")))
  expect_equal(colnames(reordered), c("b", "a"))
})

test_that("average_df rejects a mismatched factor and survives an empty level", {
  df <- data.frame(a1 = 1:2, a2 = 3:4)
  expect_error(average_df(df, factor("a")), "correspond one to one")

  # a dropped level made the notebook version error inside rowMeans()
  meta <- factor(c("a", "a"), levels = c("a", "b"))
  out <- average_df(df, meta)
  expect_equal(colnames(out), c("a", "b"))
  expect_true(all(is.na(out$b)))
})

test_that("average_df and grouped_col_mean agree", {
  df <- data.frame(a1 = c(1, 2), a2 = c(3, 4), b1 = c(10, 20))
  meta <- c("a", "a", "b")
  expect_equal(
    as.matrix(average_df(df, factor(meta))),
    as.matrix(grouped_col_mean(df, data.frame(condition = meta), "condition"))
  )
})

test_that("gg_color_hue returns n distinct hex colours", {
  cols <- gg_color_hue(3)
  expect_length(cols, 3)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", cols)))
  expect_length(unique(cols), 3)
  expect_length(gg_color_hue(1), 1)
})

test_that("call_mix_bin writes a group column and demands a value column", {
  skip_if_not_installed("mixtools")
  set.seed(1)
  df <- data.frame(value = c(stats::rnorm(300, 0.3, 0.08),
                             stats::rnorm(200, 0.7, 0.08)))

  out <- call_mix_bin(df, plot = FALSE)
  expect_true("group" %in% names(out))
  expect_equal(nrow(out), nrow(df))
  expect_true(all(out$group %in% c(0L, 1L)))

  expect_error(call_mix_bin(data.frame(x = 1)), "needs a 'value' column")
})

test_that("percentile_ranks reports the cumulative percentage", {
  out <- percentile_ranks(c(1, 2, 2, 3, 8))
  expect_equal(out$value, 1:8)
  expect_equal(out$perc.rank[out$value == 2], 60)   # 3 of 5 values <= 2
  expect_equal(out$perc.rank[out$value == 8], 100)
  expect_error(percentile_ranks(c(NA, NA)), "no non-missing values")
})

test_that("compute_metagene z-scores per gene before averaging", {
  expr <- data.frame(gene = c("g1", "g2", "g3"),
                     s1 = c(1, 10, 9), s2 = c(2, 20, 8), s3 = c(3, 30, 2))

  out <- compute_metagene(c("g1", "g2"), expr, c("s1", "s2", "s3"),
                          id_col = "gene")
  expect_named(out, c("s1", "s2", "s3"))
  # g1 and g2 both rise monotonically, so the metagene must rise too
  expect_true(all(diff(out) > 0))
  # per-gene z-scoring means g2's 10x larger scale does not dominate
  expect_equal(mean(out), 0, tolerance = 1e-8)
})

test_that("compute_metagene errors informatively on bad input", {
  expr <- data.frame(gene = "g1", s1 = 1)
  expect_error(compute_metagene("g1", expr, "s1", id_col = "nope"),
               "not a column of expr_df")
  expect_error(compute_metagene("g1", expr, "s9", id_col = "gene"),
               "missing sample column")
  expect_error(compute_metagene("zz", expr, "s1", id_col = "gene"),
               "none of gene_ids are present")
})

test_that("check_diff_genes flags and counts hits", {
  go <- data.frame(GO.ID = c("GO:1", "GO:2", "GO:3"),
                   genes = c("Gli1, Ccnd1", "Actb, Gapdh", "Mycn"))

  out <- check_diff_genes(go, c("Gli1", "Mycn"))
  expect_equal(out$is_diff_gene, c(TRUE, FALSE, TRUE))
  expect_equal(out$n_diff_gene, c(1L, 0L, 1L))

  expect_error(check_diff_genes(data.frame(x = 1), "Gli1"),
               "needs a 'genes' column")
})

test_that("lookup_count returns NA rather than erroring on a bad key", {
  tab <- data.frame(condition = c("P7", "P56"), bound = c(10, 20))
  expect_equal(lookup_count("P7", "bound", tab), 10)
  expect_true(is.na(lookup_count("P14", "bound", tab)))

  # duplicated condition is also not a unique lookup
  dup <- data.frame(condition = c("P7", "P7"), bound = c(1, 2))
  expect_true(is.na(lookup_count("P7", "bound", dup)))

  expect_error(lookup_count("P7", "nope", tab), "not a column of counts_table")
})
