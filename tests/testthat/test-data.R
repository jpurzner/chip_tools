# The example datasets are a contract: the docs, examples and vignettes all
# assert things about their shape and scale. These tests pin that contract so
# a regenerated dataset cannot quietly break them.

test_that("gnp_histone_promoter has the documented shape and scale", {
  expect_s3_class(gnp_histone_promoter, "data.frame")
  expect_equal(dim(gnp_histone_promoter), c(3000L, 21L))
  expect_true(all(vapply(gnp_histone_promoter, is.numeric, logical(1))))

  # already 0-1 normalised, i.e. wins_norm_histones() output
  vals <- unlist(gnp_histone_promoter)
  expect_gte(min(vals, na.rm = TRUE), 0)
  expect_lte(max(vals, na.rm = TRUE), 1)

  # row names are transcript IDs, mixed Ensembl and UCSC as in the source
  rn <- rownames(gnp_histone_promoter)
  expect_true(all(grepl("^(ENSMUST|uc)", rn)))
  expect_gt(sum(grepl("^ENSMUST", rn)), sum(grepl("^uc", rn)))
  expect_false(anyDuplicated(rn) > 0)

  expect_true(all(c("GNP_H3K4me3_rep1", "GNP_H3K27me3_rep1", "GNP_ATACseq")
                  %in% names(gnp_histone_promoter)))
})

test_that("gnp_histone_promoter is bimodal enough to threshold", {
  # the whole segmentation vignette rests on this
  skip_if_not_installed("mixtools")
  fit <- binarize_counts(gnp_histone_promoter$GNP_H3K4me3_rep1,
                         plot = FALSE, return_what = "summary")
  expect_gt(fit[["mean_upper"]], fit[["mean_lower"]])
  expect_gt(fit[["cutoff_50"]], fit[["mean_lower"]])
  expect_lt(fit[["cutoff_50"]], fit[["mean_upper"]])
})

test_that("gnp_timecourse has the columns classify_timecourse_genes needs", {
  expect_s3_class(gnp_timecourse, "data.frame")
  expect_equal(nrow(gnp_timecourse), 3000L)

  needed <- c("avg_t50", "avg_max", "avg_fc",
              "hatten_trap_min_pval", "scott_gnp_min_pval",
              "max_expr_Scott", "max_expr_Hatten",
              "ensembl_transcript_id", "mgi_symbol")
  expect_true(all(needed %in% names(gnp_timecourse)))
})

test_that("the expression floors are on the TPM-like scale, not the 0-1 one", {
  # Using the 0-1.8 scaled max.val.* columns with a cutoff of 100 classifies
  # every transcript "unexpressed". This asserts we shipped the right family.
  expect_gt(stats::median(gnp_timecourse$max_expr_Scott, na.rm = TRUE), 10)
  expect_gt(stats::median(gnp_timecourse$max_expr_Hatten, na.rm = TRUE), 10)
  expect_gt(max(gnp_timecourse$max_expr_Scott, na.rm = TRUE), 1000)
})

test_that("classify_timecourse_genes gives a usable spread on the real sample", {
  cl <- classify_timecourse_genes(
    gnp_timecourse, fc_cutoff = 2, p_val_cutoff = 0.001,
    expr_cutoffs = c(max_expr_Scott = 100, max_expr_Hatten = 100)
  )
  counts <- table(cl$max_cat)

  # every developmental category should be populated -- if one is empty the
  # cutoffs or the columns are wrong
  for (lv in c("E15", "P0_7", "P14_P56")) {
    expect_gt(counts[[lv]], 0)
  }
  expect_gt(counts[["unchanged"]], 0)

  # 427 rows have NA time-course summaries; the classifier propagates the NA
  # rather than inventing a category, so table() covers fewer rows than df.
  expect_equal(sum(counts) + sum(is.na(cl$max_cat)), nrow(gnp_timecourse))
  expect_equal(sum(is.na(cl$max_cat)), sum(is.na(gnp_timecourse$avg_max)))
})

test_that("gnp_timecourse and gnp_histone_promoter can be joined", {
  overlap <- intersect(gnp_timecourse$ensembl_transcript_id,
                       rownames(gnp_histone_promoter))
  expect_gt(length(overlap), 2000)   # documented as 2,946
})

test_that("gnp_enhancers is formatted as enhancer_quant expects", {
  expect_type(gnp_enhancers, "list")
  expect_length(gnp_enhancers, 5L)
  expect_named(gnp_enhancers)

  for (nm in names(gnp_enhancers)) {
    d <- gnp_enhancers[[nm]]
    expect_true(all(c("chr", "start", "stop", "name", "counts",
                      "length", "mid_point", "norm") %in% names(d)),
                info = nm)
    expect_s3_class(d$chr, "factor")
    expect_equal(levels(d$chr), "chr19", info = nm)
    expect_gt(nrow(d), 100, label = nm)
    # mid_point must sit inside the peak or distances are wrong
    expect_true(all(d$mid_point >= d$start & d$mid_point <= d$stop), info = nm)
    expect_equal(d$length, d$stop - d$start, info = nm)
    # median-normalised, so the median of norm is ~1
    expect_gt(stats::median(d$norm), 0)
  }
})

test_that("gnp_promoters is formatted as enhancer_quant expects", {
  expect_s3_class(gnp_promoters, "data.frame")
  expect_equal(nrow(gnp_promoters), 870L)
  expect_true(all(c("chr", "start", "stop", "name", "strand")
                  %in% names(gnp_promoters)))
  expect_equal(levels(gnp_promoters$chr), "chr19")
  expect_type(gnp_promoters$name, "character")
  expect_setequal(unique(gnp_promoters$strand), c("+", "-"))
})

test_that("enhancer_quant runs on the shipped data and honours dist_cutoff", {
  wide <- enhancer_quant(gnp_enhancers, gnp_promoters,
                         dist_cutoff = 1e6, linear_cutoff = 20e3)
  expect_s3_class(wide, "data.frame")
  expect_equal(nrow(wide), nrow(gnp_promoters))
  expect_equal(ncol(wide), length(gnp_enhancers))
  expect_true(all(grepl("_enhancer$", colnames(wide))))
  expect_setequal(rownames(wide), gnp_promoters$name)
  expect_true(all(wide >= 0))

  # dist_cutoff was silently ignored before 0.1.0; a tighter window must give
  # strictly less total signal
  tight <- enhancer_quant(gnp_enhancers, gnp_promoters,
                          dist_cutoff = 5e4, linear_cutoff = 5e4)
  expect_lt(sum(tight), sum(wide))
})

test_that("enhancer_quant emits no message about duplicate names", {
  # the old bind_cols() construction printed a "New names:" block per sample
  expect_silent(
    enhancer_quant(gnp_enhancers[1:2], gnp_promoters,
                   dist_cutoff = 1e5, linear_cutoff = 10e3)
  )
})

test_that("t50window drops unused factor levels and undersized groups", {
  # An unused level previously produced a zero-row group that failed with
  # "length of 'dimnames' [2] not equal to array extent".
  d <- data.frame(gene = paste0("g", 1:40),
                  t50 = seq(1, 80, length.out = 40),
                  group = factor(rep("keep", 40),
                                 levels = c("keep", "never_used")))
  w <- t50window(d, window_number = 5)
  expect_named(w, "keep")
  expect_equal(dim(w[[1]]), c(8L, 5L))
  # column names are each window's mean t50, increasing
  expect_true(all(diff(as.numeric(colnames(w[[1]]))) > 0))

  d2 <- rbind(d, data.frame(gene = "x", t50 = 5, group = "tiny"))
  d2$group <- factor(d2$group)
  expect_warning(w2 <- t50window(d2, window_number = 5), "fewer than 5 rows")
  expect_named(w2, "keep")

  # the undersized-group warning fires on the way to the error
  expect_error(suppressWarnings(t50window(d[1:3, ], window_number = 5)),
               "no group has at least")
})
