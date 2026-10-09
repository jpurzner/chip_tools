# chiptools 0.2.0

Brings in the code the `Ezh2_2022` analysis notebooks were carrying inline,
and documents the pipeline they implement.

## Vignettes

Three, all built from code that runs:

- **`overview`** — the seven pipeline stages, the parameter values the
  analysis actually settled on, and how to migrate notebook code off the old
  `source()` calls and the old mixture return flags.
- **`histone-segmentation`** — the three-threshold segmentation, with a worked
  demonstration of why the upper tail has to be excluded before fitting.
- **`timecourse`** — the time reference, interpolate/align/pool, t50, and
  gene classification, including what `pval_both` changes.

## The notebooks now run on the package alone

All 28 files the notebooks `source()`d from `~/Dropbox/jp_seq` are present, so

```r
source("~/Dropbox/jp_seq/binarize_counts.R")   # x 28
```

becomes `library(chiptools)`. The seven that 0.1.0 left out are now in:
`edgeeat()`, `memb2graph_edgeeat()`, `multi_memb_prod()`,
`recursive_edgeeater()` (the fuzzy-membership graph lineage), plus
`go_cluster()`, `go2sym()` and `plot_norm_heatmap()`.

## New functions from the notebooks

Repeated inline definitions, lifted out with their globals turned into
arguments. Renamed to match the rest of the package where needed.

| New | Was | Copies in the notebooks |
| --- | --- | --- |
| `classify_timecourse_genes()` | `group_genes_txchange()` | 4 |
| `call_mix_bin()` | same | 6 |
| `average_df()` | same | 4 |
| `lookup_count()` | same | 3 |
| `gg_color_hue()` | same | 2 |
| `convert_mouse_gene_list()` | `convertMouseGeneList()` | 5 |
| `check_diff_genes()` | same | 1 |
| `compute_metagene()` | same | 1 |
| `percentile_ranks()` | `percentile.ranks()` | 1 |
| `squish_trans()` | same | 1 |

`classify_timecourse_genes()` is the one worth reading the docs for: it is
the central gene classifier, and its cut points, column names and the
`pval_both` switch were all hard-coded in the notebook copies.

79 exports, 139 testthat assertions.

## Bug fixes

- **`topGO_split()` worked only by coincidence.** It split on
  `gnp_multi_split[, split_cat]` — a global — rather than its own `t50_df`
  argument. The notebooks happened to call it as
  `topGO_split(gnp_multi_split, ...)`, so the global and the argument were
  the same object. Any other input was silently split on the wrong data.
- **`edgeeat(memb_cutoff > 0)` transferred membership it should not have.**
  The cutoff mask *replaced* the next-rank mask instead of being combined
  with it, so genes that ranked higher on a surviving vertex had their
  membership moved anyway. Now the intersection of both conditions.
- **`recursive_edgeeater()` returned `metrif_df`**, so `$metric_df` was
  `NULL`. It also accepted `memb_cutoff` and then hard-coded `0` on the
  `edgeeat()` call.
- **`go_cluster()` ignored `cutoff`** — `0.5` was hard-coded in the
  comparison.
- **`plot_norm_heatmap()` could not screen out `-Inf`.** It stayed on a data
  frame, where `max()` works but `is.finite()` errors with "default method
  not implemented for type 'list'". Now works on a matrix throughout, and
  takes the maximum over finite values only.
- **`go2sym()` indexed its result column by position** (`, 5`), so a schema
  change in `org.Mm.eg.db` would have silently returned the wrong column.
  Now named.
- **`average_df()`** (as `average_df` in the notebooks) errored on a dropped
  factor level, because the zero-column subset broke `rowMeans()`.
- **`lookup_count()`, `compute_metagene()`** read globals or hard-coded an
  identifier column name; both are arguments.
- igraph calls moved off the deprecated `get.edge.attribute()` /
  `set.*.attribute()` / `graph.data.frame()` spellings.

## Other changes

- `memb2graph_edgeeat()` and `plot_norm_heatmap()` no longer plot or print
  unconditionally; both gained arguments for it.
- `plot_norm_heatmap()` errors rather than returning an empty matrix when no
  gene clears `cutoff`.
- New `Suggests`: `AnnotationDbi`, `biomaRt`, `knitr`, `pbapply`, `rmarkdown`,
  `scales`.

# chiptools 0.1.0

First packaged release. Previously 38 loose `.R` files at the repository root,
meant to be `source()`d individually.

[`REVIEW.md`](REVIEW.md) is the full defect catalogue with line references;
this is the summary.

## Structure

- Installable R package: `DESCRIPTION`, `NAMESPACE`, `R/`, roxygen2
  documentation for all 62 exported functions, `man/` pages, 69 `testthat`
  assertions, `.gitignore`, MIT `LICENSE`.
- `remotes::install_github("jpurzner/chip_tools")` then `library(chiptools)`.
  No more `source()`.
- Dependencies declared. Only eight packages are required; everything heavier
  is in `Suggests` and checked at call time with one actionable message naming
  the missing package.
- 25 ChIP/RNA-seq sibling functions added from `~/Dropbox/jp_seq`, including
  `opti_map()`, which three files already depended on but which was not in the
  repository at all.

## Breaking changes

- **`binarize_counts()`, `split_counts()`, `trinarize_counts()`** replace the
  mutually exclusive `return_cutoff` / `return_mean` / `return_cutoff_mean`
  flags with one `return_what` argument. The old flag combinations were
  resolved by a cascading `if`/`else if` chain, so setting two did something
  arbitrary. `return_cutoff_mean = TRUE` is now `return_what = "summary"`.
- **`trinarize_counts()`** gains `proba` and `plot` and drops the unreachable
  flags; it now returns calls rather than the raw `mixEM` object.
- **`split_counts()`** replaces the hard-coded `k = 2:4` loop with `min_k` /
  `max_k`, and returns the selected model, the per-`k` comparison table, or
  both.
- **`chip_segment()`** now uses its own arguments and returns a list of
  `cutoffs` and `classes`. It previously read a global and could not run.
- **`chip_scatter_cutoffs()`** takes `label_genes`, `mark_genes`,
  `mark_category` as arguments; `label_genes` was read from the global
  environment.
- **`plot_anova_tukey()`** takes `block_var`; the blocking covariate was
  hard-coded as `cell_line`.
- **`plot_tsne_kmeans()`** gains `palette_n` (default 9). The deleted
  `replot_plot_tsne_kmeans.R` was the same function with 7.
- **`wins_norm_histones()`** gains `normalize`.
- **`soften_zero()`** no longer plots by default (`plot = FALSE`).
- **`load_merge()`, `chip_de_rpm()`** print only with `verbose = TRUE`.
- **`chip_mclust_icl_all()`** has `crossing_points_df = NULL` as a default.
- **`ggroc()`** now errors on invalid input instead of falling through.
- **`chip_de_tpm()` is gone** — it was byte-identical to `chip_de_rpm()` apart
  from the name, and computed RPM. Use `chip_de_rpm()`.
- **`prune_gene.R` is gone.** It defined a second, older copy of
  `prune_close()`; there was never a `prune_gene()` function.

## Bug fixes

Wrong results:

- **BIC was always equal to AIC.** `n <- length(data)` took the length of the
  base R `data` *function*, i.e. 1, so `log(n)` was 0. The parameter count was
  also wrong (double-counted the variance), and `split_counts()` penalised BIC
  by the component count while penalising AIC by the parameter count.
- **`binarize_counts_p()` fitted a hard-coded global**
  (`scott_max_expr_nona_noz$max_expr_Scott_l`) instead of its own input, so
  every caller got one particular analysis's answer.
- **`grouped_col_mean()` misaligned columns** whenever any condition had a
  single replicate, and dropped those conditions from the output. This fed
  `chip_de_rpm()` and both `RNA_de_tpm_*()`.
- **`exclude_extreme` discarded the log transform** in all four mixture
  functions — it subset the raw counts and overwrote the transformed values.
- **`chip_mclust_icl_all()` returned all-`NA`** for any column with nothing
  above the tail cutoff, because `x[-integer(0)]` is empty rather than `x`. It
  also passed `model =`, which is not an `Mclust()` argument, so the model
  family was re-selected by BIC instead of fixed to the ICL choice.
- **`enhancer_quant()` ignored `dist_cutoff`** — it was not forwarded to the
  inner distance call.
- **`remove_intersect_all()` measured distance on the y axis only** —
  `sqrt(dy^2 + dy^2)`; `diff(x)` never appeared.
- **`df_rowfilt()` returned its input unchanged** — built `df_new`, returned
  `df_filt`.
- **The topGO parallel workers loaded `GOstats`** while the worker function
  calls topGO, so every worker failed. They now load `topGO` and
  `org.Mm.eg.db`.
- **`wins_norm_histones(median_subtract =)` was applied to the wrong matrix.**
  Now wired correctly, though min-max rescaling cancels a constant shift, so
  it is only observable with `normalize = FALSE` — see
  `?wins_norm_histones`.
- **`topGO_split()` split by an undefined global** (`gnp_multi_split`) rather
  than its own argument.
- **`get_co_clusters()` read an undefined global** (`fuzz_list`) in six places.
- **`chip_segment()` could not run** — ignored its arguments, read a global,
  and used `combined_data` fifty lines before creating it.
- **`soften_zero()` errored** whenever the zero bin was already below target,
  because `sample(x, n)` rejects a negative `n`.
- **Root bracketing.** `find.cutoff()` took its lower bound from
  `x[which.min(f(x))]`, which need not bracket a sign change — the usual cause
  of the `"failed to find root"` fallback. Now bracketed between the component
  means.

Dead code removed:

- `split_counts()`: ~100 unreachable lines after an early `return()`.
- `trinarize_counts()`: the whole return-option body, unreachable after
  `return(model)`.
- `rna_de_tpm2.R`: a TPM computation immediately overwritten by a `sweep()` of
  the raw counts.
- `time_of_change.R`, `wins_norm_histones.R`, `grouped_col_mean.R`: variables
  computed and never read.

Name collisions resolved — four function names were each defined in two files,
and since every file was `source()`d, whichever loaded last silently won.

Nine `source()` calls to hard-coded `~/Dropbox/jp_seq` paths removed; one file
sourced itself.

API modernisation: `inherits()` for `class(x) ==`, `after_stat(density)` for
`..density..`, `linewidth` for `size` on lines, `prevR`/`DescTools`/`stringr`
calls namespaced or replaced with base equivalents.

Attribution added for the three third-party functions (`counts_to_tpm()`,
`summarySE()`, `spline.poly()`).

## Known limitations

- **`remove_intersect_all()` cannot run on current R.** `rgeos` was archived
  from CRAN in October 2023. The three calls map onto
  `sf::st_intersects`/`st_buffer`/`st_difference`; `sf` is in `Suggests` in
  anticipation. This is the top item in REVIEW.md's "worth doing next".
- `expr_half_max_min()` works on a matrix, not a list.
- `chip_gene_group_comp(t50_cutoff = NULL)` yields an `NA` p-value by
  construction.
- The GO functions are hard-wired to mouse annotation.
