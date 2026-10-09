# chiptools

Analysis helpers accumulated across the GNP/MB cerebellar ChIP-seq and RNA-seq
work. Previously a flat directory of 38 loose `.R` files meant to be `source()`d
one at a time; now an installable R package of 79 documented functions.

Start with `vignette("overview", package = "chiptools")` — it maps the seven
stages of the GNP/MB analysis onto the functions that implement them, and
shows how to migrate notebook code off the old `source()` calls.

```r
# install.packages("remotes")
remotes::install_github("jpurzner/chip_tools", build_vignettes = TRUE)
library(chiptools)
?binarize_counts
vignette("overview", package = "chiptools")
```

Most functions take a **features × samples** data frame — genes, promoters or
peaks down the rows, one column per sample — and are meant to be composed at
the console rather than run as a fixed pipeline.

## What's here

### Thresholding binding signal

Separate bound from unbound by fitting a mixture to the signal distribution,
rather than picking a cutoff by eye.

| Function | Use |
| --- | --- |
| `binarize_counts()` | 2-component normal mixture → 0/1 call per feature |
| `binarize_counts_p()` | the same with a gamma mixture, for un-logged right-skewed counts |
| `trinarize_counts()` | 3-component mixture → low/mid/high call |
| `split_counts()` | fits k = 2…n and keeps the lowest-BIC model |
| `soften_zero()` | thins a dominant zero bin so it stops swamping the fit |

### Segmenting histone-mark intensity

The sequence used throughout the GNP/MB work: find where the long right tail
starts, hold it out, cluster what's left.

| Function | Use |
| --- | --- |
| `trim_upper_tail()` | per-column upper-tail cutoff from the rank/value knee |
| `chip_mclust_icl_all()` | per-column mclust clustering, k chosen by ICL, tail as its own class |
| `chip_segment()` | both of the above in one call |
| `chip_binarize_tail_removed()` | per-column `binarize_counts()` with the tail excluded |
| `chip_categorize_and_collapse()` | unbound/bound/tail per replicate, collapsed per condition |
| `wins_norm_histones()` | winsorise by point count, then scale 0–1 |
| `vsd_norm()` | DESeq2 VST, averaged by condition |

### Differential binding and expression

| Function | Use |
| --- | --- |
| `chip_de_rpm()` | DESeq2 on ChIP counts + RPM columns + up/down calls |
| `RNA_de_tpm_3()` | DESeq2 on RNA counts + RPM/RPKM/TPM columns (current version) |
| `RNA_de_tpm_2()` | the previous version, kept because figures were made with it |
| `counts_to_tpm()` | TPM with per-library mean fragment length |
| `counts2tpm()` | TPM with a fixed read length |
| `count_table2rpkm()` | RPKM averaged over chosen columns |
| `count_table_trim()`, `load_metadata()`, `load_merge()`, `grouped_col_mean()` | input wrangling |

### Time courses

| Function | Use |
| --- | --- |
| `time_expr_interp()` | interpolate studies sampled at different ages onto one axis |
| `time_expr_align()` | remove the between-study offset |
| `time_expr_avg()` | pool aligned time courses |
| `expr_half_max()` | t50, the time a gene reaches half its maximum |
| `expr_half_max_min()` | separate t50 for the largest rise and the largest fall |
| `time_of_change()` | where each trajectory rises and falls, and for how long |
| `t50_to_compressed()` | map t50 back onto sampled time points |
| `t50window()` | bin genes into equal-sized windows along t50 |

### GO enrichment

| Function | Use |
| --- | --- |
| `topGO_wrap()` | topGO with semantic ordering and p-value filtering (current) |
| `topGO_split()` | one enrichment per level of a grouping column |
| `topGO_timeseries()` | one enrichment per t50 window |
| `go_to_gene_table()`, `collapse_GO_columns()`, `organize_genes_and_terms()` | reshape results for printing |

### Figures

| Function | Use |
| --- | --- |
| `plot_tsne_kmeans()` | t-SNE panel + expression heatmap, cluster colours matched via `opti_map()` |
| `replot_tsne_kmeans()` | redraw one panel without re-clustering |
| `plot_tsne_optics()` | the same with OPTICS density clustering |
| `outline_tsne()` | trace polygons around clusters for shading |
| `remove_intersect_all()` | resolve overlaps between those polygons — **see Known limitations** |
| `spline.poly()` | smooth the polygons |
| `plot_cluster_line()`, `plot_cluster_line_EBseq()`, `plot_cluster_line_interp()` | mean trajectory per cluster |
| `plot_group_timeline()` | cluster trajectories ordered along a shared timeline |
| `line_scatter_facet()` | cluster trajectory beside a per-gene scatter |
| `chip_scatter_cutoffs()` | two samples scattered with their binding cutoffs boxed |
| `chip_gene_group_comp()` | early vs late gene groups, with a Wilcoxon p-value |
| `heatmap_table()` | two-way table as a heatmap |
| `plot_anova_tukey()` | faceted plots with ANOVA/Tukey brackets |
| `ggroc()` | ROC curves from pROC objects |

### Clustering and membership graphs

Fuzzy c-means membership across several data sets, turned into a graph that
can be collapsed to a chosen number of clusters.

| Function | Use |
| --- | --- |
| `multi_memb_prod()` | membership product across every cross-data-set cluster combination |
| `memb2graph_edgeeat()` | build the directed co-cluster graph from those products |
| `edgeeat()` | delete a vertex, transferring its gene membership to its neighbours |
| `recursive_edgeeater()` | collapse the graph to a target vertex count |
| `get_co_clusters()` | the hard-assignment equivalent |
| `graph_cluster_plot()`, `graph_cluster_plot_EBseq()` | plot cluster trajectories from a graph |
| `opti_map()` | greedily match two cluster labellings so colours agree |

### Helpers lifted from the notebooks

Repeated inline definitions from the `Ezh2_2022` notebooks, now with
arguments in place of the globals they used to read.

| Function | Use |
| --- | --- |
| `classify_timecourse_genes()` | cut genes into max/t50 windows, override unchanged and unexpressed |
| `average_df()` | average columns within levels of a factor |
| `call_mix_bin()` | `binarize_counts()` over a long-format `value` column |
| `check_diff_genes()` | flag GO rows containing a gene of interest |
| `compute_metagene()` | z-scored metagene score per sample |
| `convert_mouse_gene_list()` | mouse to human orthologs via BioMart |
| `percentile_ranks()` | percentile rank of each integer value |
| `lookup_count()` | one cell of a per-condition count table |
| `gg_color_hue()` | the default ggplot2 discrete palette |
| `squish_trans()` | compress part of a continuous axis |

### Other utilities

`prune_close()`, `enhancer_quant()`, `df_rowfilt()`, `df_rowmatch()`,
`min_finite()`, `summarySE()`, `plot_norm_heatmap()`, `go2sym()`,
`go_cluster()`.

## Vignettes

| Vignette | Covers |
| --- | --- |
| `overview` | the seven pipeline stages, the parameter values the analysis settled on, and how to migrate notebook code |
| `histone-segmentation` | the three-threshold segmentation, why the tail is excluded before fitting, and the clustering alternative |
| `timecourse` | the time reference, interpolate/align/pool, t50, and gene classification |

## Typical session

```r
library(chiptools)

# 1. where does the upper tail of each sample start?
cutoffs <- trim_upper_tail(histone_rlog)

# 2. mixture threshold per sample, tail excluded from the fit
breaks <- chip_binarize_tail_removed(histone_rlog, cutoffs)

# 3. call every gene unbound / bound / tail, collapse replicates
calls <- chip_categorize_and_collapse(signal, breaks)

# 4. when does each gene turn on?
interp <- time_expr_interp(expr_list, time_ref)
t50 <- expr_half_max(time_expr_avg(interp))

# 5. compare binding in early vs late genes
chip_gene_group_comp("H3K27me3", signal, breaks, t50_cutoff = 20)
```

## Dependencies

Only `dplyr`, `ggplot2`, `tidyr`, `reshape2`, `plyr`, `RColorBrewer`, `rlang`
and `zoo` are required. Everything heavier is in `Suggests` and checked at call
time, so you get one actionable message naming the package rather than a
`could not find function` error several frames deep:

```r
#> Error: chip_mclust_icl_all needs 'mclust', which is not installed.
#>   install.packages(c("mclust"))
```

Bioconductor packages (`DESeq2`, `topGO`, `GO.db`, `org.Mm.eg.db`, `GOSemSim`,
`Biobase`, `Mfuzz`) install with `BiocManager::install()`. The GO functions are
hard-wired to mouse annotation (`org.Mm.eg.db`).

## Known limitations

- **`remove_intersect_all()` cannot run on current R.** It needs `rgeos`, which
  was archived from CRAN in October 2023. The three calls
  (`gIntersects`/`gBuffer`/`gDifference`) map onto
  `sf::st_intersects`/`st_buffer`/`st_difference`; `sf` is already in
  `Suggests` in anticipation. Until then, install `rgeos` from the CRAN archive.
- **`expr_half_max_min()` only works on a matrix**, not a list. The list branch
  collapses its six-element result to a single number. Use `lapply()` over your
  data sets instead.
- **`chip_gene_group_comp(t50_cutoff = NULL)`** filters `t50_cat` to
  `"Late"`/`"Mid"` and then tests for `"early"`/`"late"`, so the p-value comes
  back `NA`. Pass a `t50_cutoff`.
- **`wins_norm_histones(median_subtract = TRUE)`** has no effect on the default
  output: min-max rescaling cancels a constant shift. It is observable only
  with `normalize = FALSE`.
- **`enhancer_quant(par = TRUE)`** was flagged as intermittently failing in the
  original source. The serial path is the reliable one.
- The GO and mouse-annotation functions are single-species.

## Development

```r
devtools::load_all()
devtools::test()      # 139 assertions
devtools::document()
devtools::build_vignettes()
devtools::check()
```

[`REVIEW.md`](REVIEW.md) catalogues the defects found when this was
restructured from loose scripts, with file and line references into the
pre-restructure history. Worth reading before changing any of the
mixture-model or segmentation code.

## License

MIT. `counts_to_tpm()`, `summarySE()` and `spline.poly()` are third-party; see
their `@source` tags.
