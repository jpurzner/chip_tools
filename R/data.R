#' Example datasets
#'
#' Four subsets of the real GNP/MB analysis inputs, included so the examples
#' and vignettes run on actual data rather than simulated frames. Together
#' they are about 640 KB, down from roughly 60 MB of source tables.
#'
#' They are subsets, not the analysis. Anything computed from them is an
#' illustration of the call sequence, not a published result — the figures in
#' the paper were made from the full tables. See `data-raw/make-datasets.R`
#' in the repository for exactly how each was derived.
#'
#' @section How they were subset:
#' The sampling scheme differs by dataset, deliberately:
#'
#' * [gnp_histone_promoter] and [gnp_timecourse] are **random row samples**.
#'   The mixture-model and rank-curve machinery only depends on the shape of
#'   each column's distribution, which a random sample of transcripts
#'   preserves.
#' * [gnp_enhancers] and [gnp_promoters] are restricted to **one chromosome**
#'   (chr19), not randomly sampled. [enhancer_quant()] works per chromosome
#'   and weights by promoter-to-enhancer distance, so a random sample would
#'   destroy the spatial structure and make the output meaningless.
#'
#' @name chiptools-data
#' @keywords datasets
NULL


#' Histone-mark signal over GNP promoters
#'
#' Normalised ChIP-seq, ATAC-seq and RNA-seq signal at protein-coding
#' transcript promoters in P7 cerebellar granule neuron precursors. Values are
#' already scaled to 0–1 per column, i.e. they are [wins_norm_histones()]
#' output, so they are ready for [trim_upper_tail()] and [binarize_counts()].
#'
#' The three `MIRA_low` columns are DNA methylation; `ATACseq` is chromatin
#' accessibility; `RNAseq_rep1..3` are expression. The rest are histone marks
#' and chromatin proteins, several with two replicates.
#'
#' @format A data frame of 3,000 rows and 21 numeric columns, one per assay.
#'   Row names are transcript identifiers, mixed: roughly seven in ten are
#'   Ensembl (`ENSMUST...`) and the rest UCSC (`uc007aia.2`), as in the
#'   source matrix.
#'   Join on them with that in mind.
#'   \describe{
#'     \item{GNP_MIRA_low_rep1}{DNA methylation (MIRA), replicate 1}
#'     \item{GNP_MIRA_low_rep2}{DNA methylation (MIRA), replicate 2}
#'     \item{GNP_MIRA_low_rep3}{DNA methylation (MIRA), replicate 3}
#'     \item{GNP_H3K4me3_rep1}{active promoter mark}
#'     \item{GNP_H3K4me1_rep1}{enhancer mark, replicate 1}
#'     \item{GNP_H3K4me1_rep2}{enhancer mark, replicate 2}
#'     \item{GNP_H3K27me3_rep1}{repressive PRC2 mark, replicate 1}
#'     \item{GNP_H3K27me3_rep2}{repressive PRC2 mark, replicate 2}
#'     \item{GNP_H3K27ac_rep1}{active enhancer/promoter mark, replicate 1}
#'     \item{GNP_H3K27ac_rep2}{active enhancer/promoter mark, replicate 2}
#'     \item{GNP_Ring1_rep1}{PRC1 component RING1}
#'     \item{GNP_Ring1b_rep2}{PRC1 component RING1B}
#'     \item{GNP_H2Aubi119_rep1}{PRC1 catalytic product, replicate 1}
#'     \item{GNP_H2Aubi119_rep2}{PRC1 catalytic product, replicate 2}
#'     \item{GNP_Pol2_ser5_abcam}{initiating RNA Pol II (Abcam antibody)}
#'     \item{GNP_Pol2_ser5_am}{initiating RNA Pol II (Active Motif antibody)}
#'     \item{GNP_ATACseq}{chromatin accessibility}
#'     \item{GNP_H3K36me3}{transcribed-body mark}
#'     \item{RNAseq_rep1}{expression, replicate 1}
#'     \item{RNAseq_rep2}{expression, replicate 2}
#'     \item{RNAseq_rep3}{expression, replicate 3}
#'   }
#'
#' @source A 3,000-transcript random sample (seed 42) of
#'   `GNP_combo_quant_2019_02_28/gnp_histone_normalized_2019_02_28.txt`
#'   (24,209 x 21).
#'
#' @seealso [trim_upper_tail()], [chip_binarize_tail_removed()],
#'   [chip_segment()], [plot_tsne_kmeans()],
#'   `vignette("histone-segmentation", package = "chiptools")`
#'
#' @examples
#' dim(gnp_histone_promoter)
#' summary(gnp_histone_promoter[, c("GNP_H3K4me3_rep1", "GNP_H3K27me3_rep1")])
#'
#' # H3K4me3 is bimodal -- background and real signal
#' hist(gnp_histone_promoter$GNP_H3K4me3_rep1, breaks = 60,
#'      main = "H3K4me3 at promoters", xlab = "normalised signal")
"gnp_histone_promoter"


#' Transcript time-course summary
#'
#' Per-transcript summary of cerebellar development, pooled across the two
#' time-course studies: when each transcript peaks (`avg_max`), when it
#' reaches half its maximum (`avg_t50`), its fold change, and the per-study
#' significance and expression maxima that [classify_timecourse_genes()] uses
#' to screen out genes that do not change or are not expressed.
#'
#' Of the 3,000 rows, 2,946 are also in [gnp_histone_promoter], so the two can
#' be joined on transcript ID.
#'
#' 427 rows have `NA` for the time-course summaries, which is what the source
#' table carries for transcripts that could not be fitted.
#' [classify_timecourse_genes()] propagates that `NA` rather than inventing a
#' category, so `table()` on its output covers 2,573 rows, not 3,000.
#'
#' @format A data frame of 3,000 rows and 15 columns:
#'   \describe{
#'     \item{ensembl_transcript_id, ensembl_gene_id, mgi_symbol}{identifiers}
#'     \item{ensembl_chr, ensembl_chromStart, ensembl_chromEnd, ensembl_strand}{coordinates (mm9)}
#'     \item{avg_t50}{half-maximum time on the interpolated axis — the "when does it turn on" number}
#'     \item{avg_max}{index of the expression maximum on the interpolated axis}
#'     \item{avg_fc}{average fold change across the time course}
#'     \item{hatten_trap_min_pval, scott_gnp_min_pval}{smallest per-study p-value across time points}
#'     \item{max_expr_Scott, max_expr_Hatten}{per-study expression maximum on a TPM-like scale (median ~170), for the expression floor}
#'     \item{ensembl_exon_length}{total exonic length, for TPM/RPKM}
#'   }
#'
#' @source A 3,000-transcript sample of `Ezh2_2022/toc_t50_2023_05_17.txt`
#'   (24,631 x 115), keeping only the 15 columns the package functions use.
#'   Rows were chosen to overlap [gnp_histone_promoter] where possible.
#'
#' @section Careful with the expression columns:
#' The source table carries two similarly named families and they are on
#' different scales. `max_expr_Scott` / `max_expr_Hatten`, kept here, are
#' expression maxima on a TPM-like scale with a median around 170 — these are
#' the ones the `expr_cutoffs = 100` screen is for. The source also has
#' `max.val.Scott` / `max.val.Hatten`, which are [time_of_change()] outputs on
#' a 0–1.8 normalised scale; using those with a cutoff of 100 classifies every
#' transcript `"unexpressed"`. This is why [classify_timecourse_genes()] takes
#' `expr_cutoffs` as a *named* vector rather than two positional arguments.
#'
#' @seealso [classify_timecourse_genes()], [expr_half_max()],
#'   [t50window()], `vignette("timecourse", package = "chiptools")`
#'
#' @examples
#' str(gnp_timecourse, max.level = 1)
#'
#' classified <- classify_timecourse_genes(
#'   gnp_timecourse,
#'   fc_cutoff    = 2,
#'   p_val_cutoff = 0.001,
#'   expr_cutoffs = c(max_expr_Scott = 100, max_expr_Hatten = 100)
#' )
#' table(classified$max_cat)
"gnp_timecourse"


#' Enhancer peak quantification on chromosome 19
#'
#' ChIP-seq and ATAC-seq read counts over enhancer peaks called with EPIC,
#' excluding peaks within 1 kb of a TSS, ENCODE blacklist regions and
#' H3K4me3-marked regions — so what is left is distal regulatory sequence.
#' Counts are median-normalised per sample.
#'
#' Formatted exactly as [enhancer_quant()] expects, and restricted to chr19
#' so that distance weighting against [gnp_promoters] is meaningful.
#'
#' @format A named list of five data frames, one per sample, each with:
#'   \describe{
#'     \item{chr}{factor, `"chr19"` only}
#'     \item{start, stop}{peak bounds (mm9)}
#'     \item{name}{EPIC peak identifier}
#'     \item{counts}{raw reads in the peak}
#'     \item{length}{`stop - start`}
#'     \item{mid_point}{peak centre, which is what distances are measured from}
#'     \item{norm}{`counts` divided by the sample's median count}
#'   }
#'   Peaks per sample: `GNP_H3K27ac_rep1` 252, `GNP_H3K27ac_rep2` 172,
#'   `GNP_H3K4me1_rep1` 697, `GNP_H3K4me1_rep2` 1150,
#'   `GNP_P7_ATACseq` 765.
#'
#' @source chr19 subset of
#'   `2016_06_18_integration/enhancers/*_no1KTSS_bl_noK4me3_quant.bed`
#'   (9,158–48,960 peaks genome-wide). Median normalisation was computed on
#'   the full genome-wide set *before* subsetting, so the values match what
#'   the whole-genome analysis would give.
#'
#' @seealso [enhancer_quant()], [gnp_promoters]
#'
#' @examples
#' sapply(gnp_enhancers, nrow)
#' head(gnp_enhancers$GNP_H3K27ac_rep1)
#'
#' \dontrun{
#' # distance-weighted enhancer signal per promoter
#' eq <- enhancer_quant(gnp_enhancers, gnp_promoters,
#'                      dist_cutoff = 1e6, linear_cutoff = 20e3)
#' head(eq)
#' }
"gnp_enhancers"


#' Promoter transcription start sites on chromosome 19
#'
#' One-base TSS coordinates for protein-coding transcripts expressed in GNPs,
#' formatted as [enhancer_quant()]'s `promoter_df` argument and restricted to
#' chr19 to match [gnp_enhancers].
#'
#' @format A data frame of 870 rows and 5 columns:
#'   \describe{
#'     \item{chr}{factor, `"chr19"` only}
#'     \item{start, stop}{the one-base TSS position (mm9)}
#'     \item{name}{Ensembl transcript ID}
#'     \item{strand}{`"+"` or `"-"`}
#'   }
#'
#' @source chr19 subset of
#'   `2016_06_18_integration/enhancers/mm_ensembl67_pc_gnp_v1_sort_TSS_2.bed`
#'   (Ensembl 67, mm9).
#'
#' @seealso [enhancer_quant()], [gnp_enhancers]
#'
#' @examples
#' head(gnp_promoters)
#' nrow(gnp_promoters)
"gnp_promoters"
