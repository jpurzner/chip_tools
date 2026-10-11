# Build the package's example datasets from the source analysis files.
#
# Not run at install or check time -- data-raw/ is in .Rbuildignore. Re-run by
# hand when the source tables change:
#
#   Rscript data-raw/make-datasets.R
#
# The sources are the real analysis inputs under ~/Dropbox/GNP_MB_integrate.
# They total ~60 MB and do not belong in a package, so each dataset here is a
# deliberate subset, chosen to keep the structure the functions depend on:
#
#   * the promoter matrix is a RANDOM row sample -- the mixture and rank-curve
#     machinery only cares about the shape of each column's distribution, so a
#     random sample of transcripts preserves what matters
#   * the enhancer and promoter coordinate tables are restricted to ONE
#     CHROMOSOME, not randomly sampled -- enhancer_quant() works per chromosome
#     and weights by promoter-to-enhancer distance, so random sampling would
#     destroy the spatial structure and make the result meaningless
#
# Everything is saved with compress = "xz", which on these tables beats the
# gzip default by roughly 40%.

SRC <- "~/Dropbox/GNP_MB_integrate"
src <- function(...) path.expand(file.path(SRC, ...))

set.seed(42)          # row sampling must be reproducible
N_PROMOTER <- 3000    # transcripts kept in the promoter matrix
CHR        <- "chr19" # chromosome kept in the coordinate tables

stopifnot(dir.exists(path.expand(SRC)))


## gnp_histone_promoter --------------------------------------------------
## 0-1 normalised ChIP/ATAC/RNA signal over transcript promoters.
## Source: GNP_combo_quant_2019_02_28/gnp_histone_normalized_2019_02_28.txt
##         24,209 transcripts x 21 columns, 9.1 MB

message("reading promoter signal matrix ...")
full <- read.table(
  src("GNP_combo_quant_2019_02_28", "gnp_histone_normalized_2019_02_28.txt"),
  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
)

keep <- sort(sample(nrow(full), N_PROMOTER))
gnp_histone_promoter <- full[keep, , drop = FALSE]

message(sprintf("  %d x %d  (from %d x %d)",
                nrow(gnp_histone_promoter), ncol(gnp_histone_promoter),
                nrow(full), ncol(full)))


## gnp_timecourse -------------------------------------------------------
## Time-course summary per transcript: t50, time of maximum, fold change,
## per-study p-values and expression maxima.
## Source: Ezh2_2022/toc_t50_2023_05_17.txt
##         24,631 transcripts x 115 columns, 20.6 MB
##
## Only the columns classify_timecourse_genes() and the t50 functions need
## are kept; the other ~95 are per-study intermediates.

message("reading time-course table ...")
toc_full <- read.delim(src("Ezh2_2022", "toc_t50_2023_05_17.txt"),
                       check.names = FALSE)

toc_cols <- c(
  "ensembl_transcript_id", "ensembl_gene_id", "mgi_symbol",
  "ensembl_chr", "ensembl_chromStart", "ensembl_chromEnd", "ensembl_strand",
  "avg_t50", "avg_max", "avg_fc",
  "hatten_trap_min_pval", "scott_gnp_min_pval",
  # the expression floors classify_timecourse_genes() screens on. These are
  # on a TPM-like scale (median ~170); do NOT confuse them with
  # max.val.Scott / max.val.Hatten, which are time_of_change() outputs on a
  # 0-1.8 normalised scale and will classify everything "unexpressed" if
  # used with a cutoff of 100.
  "max_expr_Scott", "max_expr_Hatten",
  "ensembl_exon_length"
)
missing_cols <- setdiff(toc_cols, names(toc_full))
if (length(missing_cols)) {
  stop("time-course source is missing: ", paste(missing_cols, collapse = ", "))
}

# Match the promoter matrix where possible, so the two datasets can be joined
# in an example. Top up with a random sample if the overlap is short.
in_both <- intersect(rownames(gnp_histone_promoter),
                     toc_full$ensembl_transcript_id)
message(sprintf("  %d transcripts in both tables", length(in_both)))

rows <- which(toc_full$ensembl_transcript_id %in% in_both)
if (length(rows) < N_PROMOTER) {
  top_up <- sample(setdiff(seq_len(nrow(toc_full)), rows),
                   N_PROMOTER - length(rows))
  rows <- sort(c(rows, top_up))
}

gnp_timecourse <- toc_full[rows, toc_cols, drop = FALSE]
rownames(gnp_timecourse) <- NULL

message(sprintf("  %d x %d  (from %d x %d)",
                nrow(gnp_timecourse), ncol(gnp_timecourse),
                nrow(toc_full), ncol(toc_full)))


## gnp_enhancers + gnp_promoters ----------------------------------------
## Enhancer peak quantification and promoter TSS coordinates, both on one
## chromosome so the distance weighting in enhancer_quant() is meaningful.
## Sources: 2016_06_18_integration/enhancers/*_no1KTSS_bl_noK4me3_quant.bed
##          2016_06_18_integration/enhancers/mm_ensembl67_pc_gnp_v1_sort_TSS_2.bed

enh_dir <- src("2016_06_18_integration", "enhancers")

enh_files <- list.files(enh_dir, pattern = "no1KTSS_bl_noK4me3_quant[.]bed$",
                        full.names = TRUE)
stopifnot(length(enh_files) > 0)
enh_names <- gsub("_no1KTSS_bl_noK4me3_quant[.]bed$", "", basename(enh_files))

# This is format_enh() from enhancers/enhancer_metri.Rmd, plus the median
# normalisation applied in the next chunk of that notebook.
format_enh <- function(df) {
  df <- df[, c(1:4, 7)]
  colnames(df) <- c("chr", "start", "stop", "name", "counts")
  df$length    <- df$stop - df$start
  df$mid_point <- round(df$start + df$length / 2)
  df$norm      <- df$counts / stats::median(df$counts)
  df$chr       <- as.factor(df$chr)
  df
}

message("reading enhancer peaks ...")
gnp_enhancers <- lapply(enh_files, function(f) {
  df <- read.table(f, sep = "\t", header = FALSE)
  # Normalise on the FULL set, then subset -- doing it the other way round
  # would compute the median from one chromosome and shift every value.
  df <- format_enh(df)
  df <- df[df$chr == CHR, , drop = FALSE]
  df$chr <- factor(df$chr, levels = CHR)
  rownames(df) <- NULL
  df
})
names(gnp_enhancers) <- enh_names

for (nm in names(gnp_enhancers)) {
  message(sprintf("  %-32s %5d peaks on %s", nm, nrow(gnp_enhancers[[nm]]), CHR))
}

message("reading promoter TSS ...")
prom <- read.table(file.path(enh_dir, "mm_ensembl67_pc_gnp_v1_sort_TSS_2.bed"),
                   sep = "\t", header = FALSE)
colnames(prom) <- c("chr", "start", "stop", "name", "score", "strand")
gnp_promoters <- prom[prom$chr == CHR, c("chr", "start", "stop", "name", "strand")]
gnp_promoters$chr  <- factor(gnp_promoters$chr, levels = CHR)
gnp_promoters$name <- as.character(gnp_promoters$name)
rownames(gnp_promoters) <- NULL

message(sprintf("  %d promoters on %s", nrow(gnp_promoters), CHR))


## save ----------------------------------------------------------------

usethis_save <- function(obj, name) {
  f <- file.path("data", paste0(name, ".rda"))
  e <- new.env(); assign(name, obj, envir = e)
  save(list = name, file = f, envir = e, compress = "xz", version = 2)
  message(sprintf("  %-24s %6.1f KB", basename(f), file.size(f) / 1024))
}

message("saving ...")
usethis_save(gnp_histone_promoter, "gnp_histone_promoter")
usethis_save(gnp_timecourse,       "gnp_timecourse")
usethis_save(gnp_enhancers,        "gnp_enhancers")
usethis_save(gnp_promoters,        "gnp_promoters")

message(sprintf("total data/ = %.1f KB",
                sum(file.size(list.files("data", full.names = TRUE))) / 1024))
