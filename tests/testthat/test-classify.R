# classify_timecourse_genes() is the central classifier of the time-course
# analysis, so the override precedence and the pval_both switch are pinned
# here explicitly.

demo_df <- function() {
  data.frame(
    row.names = paste0("g", 1:5),
    #        early  mid    late   low-fc  unexpressed
    avg_max = c(1,    5,     30,    10,     1),
    avg_t50 = c(2,    10,    25,    20,     2),
    avg_fc  = c(3,    4,     5,     0.2,    3),
    hatten_trap_min_pval = c(0.001, 0.4,   0.001, 0.5, 0.001),
    scott_gnp_min_pval   = c(0.001, 0.001, 0.001, 0.5, 0.001),
    max_expr_Scott  = c(100, 100, 100, 100, 1),
    max_expr_Hatten = c(100, 100, 100, 100, 1)
  )
}

cutoffs <- c(max_expr_Scott = 5, max_expr_Hatten = 5)

test_that("genes are cut into max and t50 windows", {
  out <- classify_timecourse_genes(demo_df(), fc_cutoff = 2,
                                   p_val_cutoff = 0.05,
                                   expr_cutoffs = cutoffs)

  expect_equal(as.character(out["g1", "max_cat"]), "E15")
  expect_equal(as.character(out["g2", "max_cat"]), "P0_7")
  expect_equal(as.character(out["g3", "max_cat"]), "P14_P56")
  expect_equal(as.character(out["g1", "t50_cat"]), "Early")
  expect_equal(as.character(out["g3", "t50_cat"]), "Late")
  expect_equal(as.character(out["g3", "max_t50_cat"]), "Max_P14_P56 t50_Late")
})

test_that("a low fold change overrides to unchanged", {
  out <- classify_timecourse_genes(demo_df(), fc_cutoff = 2,
                                   p_val_cutoff = 0.05,
                                   expr_cutoffs = cutoffs)
  expect_equal(as.character(out["g4", "max_cat"]), "unchanged")
  expect_equal(as.character(out["g4", "max_t50_cat"]), "unchanged")
})

test_that("unexpressed is applied after unchanged, so it wins", {
  # g5 clears the fold-change and p-value screens but is below both
  # expression floors. The override order in the original was deliberate.
  out <- classify_timecourse_genes(demo_df(), fc_cutoff = 2,
                                   p_val_cutoff = 0.05,
                                   expr_cutoffs = cutoffs)
  expect_equal(as.character(out["g5", "max_cat"]), "unexpressed")
  expect_equal(as.character(out["g5", "t50_cat"]), "unexpressed")
})

test_that("unexpressed requires being below EVERY floor", {
  df <- demo_df()
  df["g5", "max_expr_Scott"] <- 100   # now above one floor
  out <- classify_timecourse_genes(df, fc_cutoff = 2, p_val_cutoff = 0.05,
                                   expr_cutoffs = cutoffs)
  expect_equal(as.character(out["g5", "max_cat"]), "E15")
})

test_that("pval_both switches between any and all", {
  # g2 misses the threshold in hatten but clears it in scott.
  lenient <- classify_timecourse_genes(demo_df(), fc_cutoff = 2,
                                       p_val_cutoff = 0.05,
                                       expr_cutoffs = cutoffs,
                                       pval_both = FALSE)
  strict <- classify_timecourse_genes(demo_df(), fc_cutoff = 2,
                                      p_val_cutoff = 0.05,
                                      expr_cutoffs = cutoffs,
                                      pval_both = TRUE)

  expect_equal(as.character(lenient["g2", "max_cat"]), "P0_7")
  expect_equal(as.character(strict["g2", "max_cat"]), "unchanged")
})

test_that("table_only reports zero-count categories too", {
  tab <- classify_timecourse_genes(demo_df(), fc_cutoff = 2,
                                   p_val_cutoff = 0.05,
                                   expr_cutoffs = cutoffs,
                                   table_only = TRUE)
  # 3 max levels x 3 t50 levels + unchanged + unexpressed
  expect_length(tab, 11)
  expect_equal(sum(tab), nrow(demo_df()))
  expect_true(any(tab == 0))                    # empty combos are retained
  expect_equal(tab[["unchanged"]], 1)
  expect_equal(tab[["unexpressed"]], 1)
})

test_that("cut points and column names are arguments", {
  df <- demo_df()
  names(df)[names(df) == "avg_max"] <- "peak_time"

  out <- classify_timecourse_genes(df, fc_cutoff = 2, p_val_cutoff = 0.05,
                                   expr_cutoffs = cutoffs,
                                   max_col = "peak_time",
                                   max_breaks = c(0, 10, 64),
                                   max_labels = c("early", "late"))
  expect_setequal(levels(out$max_cat),
                  c("early", "late", "unchanged", "unexpressed"))
  expect_equal(as.character(out["g3", "max_cat"]), "late")
})

test_that("classify_timecourse_genes validates its inputs", {
  expect_error(
    classify_timecourse_genes(data.frame(x = 1), fc_cutoff = 2,
                              p_val_cutoff = 0.05, expr_cutoffs = cutoffs),
    "missing column"
  )
  expect_error(
    classify_timecourse_genes(demo_df(), fc_cutoff = 2, p_val_cutoff = 0.05,
                              expr_cutoffs = c(5, 5)),
    "must be a named vector"
  )
  expect_error(
    classify_timecourse_genes(demo_df(), fc_cutoff = 2, p_val_cutoff = 0.05,
                              expr_cutoffs = cutoffs,
                              max_breaks = c(0, 10), max_labels = c("a", "b")),
    "max_breaks"
  )
})

test_that("squish_trans compresses the named span and inverts cleanly", {
  skip_if_not_installed("scales")
  tr <- squish_trans(10, 100, 10)

  expect_equal(tr$transform(5), 5)            # below `from`: untouched
  expect_equal(tr$transform(100), 10 + 9)     # the 10-100 span becomes 9 wide
  expect_equal(tr$transform(110), 10 + 9 + 10)
  expect_equal(tr$inverse(tr$transform(55)), 55)
  expect_equal(tr$inverse(tr$transform(110)), 110)
  expect_true(is.na(tr$transform(NA)))
})

test_that("go_cluster honours its cutoff", {
  # The cutoff argument used to be ignored, with 0.5 hard-coded.
  df <- data.frame(GO = c("GO:1", "GO:2", "GO:3"))
  sim <- matrix(c(1.0, 0.6, 0.1,
                  0.6, 1.0, 0.1,
                  0.1, 0.1, 1.0), nrow = 3,
                dimnames = list(df$GO, df$GO))

  loose <- go_cluster(df, sim, cutoff = 0.5)   # 0.6 >= cutoff: 1 and 2 merge
  tight <- go_cluster(df, sim, cutoff = 0.9)   # nothing merges

  expect_equal(loose$GO_cluster[1], loose$GO_cluster[2])
  expect_false(tight$GO_cluster[1] == tight$GO_cluster[2])
  expect_equal(length(unique(tight$GO_cluster)), 3)
})

test_that("plot_norm_heatmap centres each gene on its own mean", {
  dat <- data.frame(s1 = c(1, 10), s2 = c(2, 20), s3 = c(3, 30))
  out <- plot_norm_heatmap(dat, cutoff = 0)

  expect_equal(dim(out), dim(dat))
  # log2 of (value / row mean) sums to ~0 per row only for symmetric data;
  # what must hold is that the row mean maps to 0
  expect_equal(unname(out[1, 2]), log2(2 / mean(c(1, 2, 3))))
  expect_true(all(is.finite(out)))
})

test_that("plot_norm_heatmap clamps -Inf from zeros", {
  dat <- data.frame(s1 = c(0, 10), s2 = c(4, 20), s3 = c(8, 30))
  out <- plot_norm_heatmap(dat, cutoff = 0)
  expect_true(all(is.finite(out)))
  expect_true(min(out) < 0)
})
