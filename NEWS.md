# chiptools 0.4.0

Clears the last blocking limitation: every function in the package now runs on
a current R installation.

## remove_intersect_all() ported from rgeos to sf

`rgeos` was archived from CRAN in October 2023, so this function could not run
at all. The calls map across directly:

| was | now |
| --- | --- |
| `rgeos::gIntersects()` | `sf::st_intersects()` |
| `rgeos::gBuffer(width = 0)` | `sf::st_make_valid()` |
| `rgeos::gDifference()` | `sf::st_difference()` |
| `sp::Polygon()` etc. | `sf::st_polygon()` |
| `raster::geom()` | `sf::st_coordinates()` |

The zero-width-buffer trick for repairing self-intersecting rings becomes
`st_make_valid()`, which is what it was always emulating.

Writing the tests established what the function actually guarantees, which was
never documented: the returned outlines have **zero pairwise overlapping
area**, and their total area equals the area of the union of the inputs.
Nothing is lost and nothing double-counted. Subtraction is *sequential*, not
symmetric — contested area goes to the higher-indexed cluster — which the
documentation now says instead of the reverse.

Four bugs surfaced in the course of the port, none of which could ever have
been hit while the function was unrunnable:

- With no intersecting pairs, `for (n in 1:nrow(comparison))` ran once with
  `n = 1` against a zero-row frame and errored. Non-overlapping input now
  passes straight through.
- `gDifference()` returns `NULL` when one polygon is wholly inside another,
  and that `NULL` went unguarded into `raster::geom()`. Fully consumed
  clusters now return zero rows.
- Intersections were tested once against the *original* polygons and then
  subtracted in a loop that mutated them, so every test after the first was
  stale. Each pair is now re-tested against current geometry.
- Rings were not closed before polygon construction. `sp` tolerated it; `sf`
  does not.

`passes` is now an argument rather than a hard-coded 5.

## orderPoints() is from contoureR, not prevR

A note added in 0.1.0 attributed `orderPoints()` to `prevR`. That was wrong —
`prevR` has never exported a function of that name. It is `contoureR`, which
`outline_tsne()` already depended on for its contour step.

The consequence was not cosmetic. The call sat inside a `try()`, so with
`prevR` listed it failed silently and **the point reordering did nothing at
all**. Making it work exposed the next problem: ordering and spur-trimming
points *across* rings mixes a hole's boundary into its container's and
produces a self-intersecting path. Both steps are now applied per ring, and
the returned frames carry a `ring` column so callers can keep them apart.

## Suggests pruned

`rgeos`, `prevR`, `sp` and `raster` are gone, as above. An audit also found
`cluster`, `naturalsort` and `stringr` declared but never called —
`naturalsort` only ever appeared in commented-out `require()` lines, and the
last `str_c()` went when 0.3.0 replaced it with `paste0()`. Every remaining
`Suggests` entry is now actually reachable from the code, which a check in
the build asserts.

## sp and raster dropped

`outline_tsne()` used `sp::Polygon(...)@area` for one number: the area of a
convex hull. That is the shoelace formula, which is exact for a simple polygon
and needs no dependency — verified identical to `sp` to ten decimal places.

With that and the `rgeos` port, `sp`, `raster` and `rgeos` are gone from
`Suggests` entirely. `sf` is the only spatial dependency left.

## License detection

GitHub reported the license as `NOASSERTION`, because its detector reads
`LICENSE` while R's `MIT + file LICENSE` convention wants that file to be
valid DCF — just `YEAR` and `COPYRIGHT HOLDER`.

These turn out to be mutually exclusive. Pointing `DESCRIPTION` at a
differently named file is rejected (`Non-standard license specification`; the
filename is not negotiable), and GitHub ignores `LICENSE.md` whenever
`LICENSE` exists. So `LICENSE` now carries the full MIT text, GitHub detects
it, and `R CMD check` reports `License stub is invalid DCF` — a NOTE, not a
warning, and a deliberate trade on a package that is not going to CRAN. It is
documented in the README so it is not mistaken for an oversight.

252 testthat assertions, up from 222. `R CMD check` reports two NOTEs, both
benign and both documented in the README: `Mfuzz` not installed on the check
machine, and the deliberate license stub above. The dependency NOTE is down
from four packages to one, since `rgeos`, `prevR`, `sp` and `raster` are all
gone.

# chiptools 0.3.0

Adds four example datasets, and fixes a dependency-resolution regression
introduced in 0.1.0 that left 23 functions unable to run.

## The regression (read this if you used 0.1.0 or 0.2.0)

0.1.0 removed the in-function `require()` calls as part of making this a
proper package. Those calls were what made bare verb names like `left_join()`
and `melt()` resolve at run time. They were replaced with a
`utils::globalVariables()` declaration, which silences `R CMD check` but
**does not make a name callable**. The result: 23 of the 79 exported functions
failed with `could not find function` on their first real call.

It went undetected because the test suite exercised the mixture and utility
functions, not these, and the vignette chunks that touch them were
`eval = FALSE`. It surfaced the moment `enhancer_quant()` was run on real data.

Fixed by importing the `Imports` packages wholesale in the package
documentation, and giving every `Suggests` call an explicit `pkg::` prefix
behind a `require_pkg()` guard. All 79 functions now have fully resolvable
calls, asserted by a `codetools` audit over every export.

Two things the fix turned up:

- **`gtable_remove_grobs()` is in cowplot, not gtable.** Four functions were
  calling it, and `gtable::gtable_remove_grobs` does not exist.
- **`stat_pvalue_manual()` is from ggpubr**, which was not in `DESCRIPTION`
  at all. `plot_anova_tukey()` could never have run.

`str_c()` calls were replaced with base `paste0()`, dropping stringr.

## Example datasets

Four subsets of the real analysis inputs, ~690 KB in total, down from roughly
60 MB of source tables. `data-raw/make-datasets.R` builds them and is not
shipped.

| Dataset | Contents |
| --- | --- |
| `gnp_histone_promoter` | 3,000 x 21 normalised ChIP/ATAC/RNA signal at GNP promoters |
| `gnp_timecourse` | 3,000 x 15 per-transcript t50, time of maximum, fold change, p-values, expression maxima |
| `gnp_enhancers` | five enhancer peak tables (172–1,150 peaks), chr19 |
| `gnp_promoters` | 870 promoter TSS coordinates, chr19 |

The sampling differs by dataset on purpose. The two matrices are random row
samples, which preserves each column's distribution — all the mixture and
rank-curve machinery needs. The two coordinate tables are restricted to one
chromosome instead, because `enhancer_quant()` works per chromosome and
weights by distance, so a random sample would destroy the spatial structure.

The vignettes now run on these rather than on simulated frames. The
segmentation vignette shows the real H3K4me3 distribution and measures how
much the upper tail moves the threshold.

## Bug fixes

- **`t50window()` failed on an unused factor level.** A filtered grouping
  column that still carried its original levels produced a zero-row group and
  then died several steps later with `length of 'dimnames' [2] not equal to
  array extent`. It now drops unused levels, warns about and skips groups
  smaller than `window_number`, and errors clearly when nothing is left.
- **`enhancer_quant()` printed a block of noise per sample.** It `bind_cols()`ed
  single-column frames that were all named `enhancer_quant`, so dplyr renamed
  them and said so, and it used `plyr::mapvalues()` where the lookup table was
  deliberately larger than the keys. Both replaced with direct construction.

## Documentation

- `classify_timecourse_genes()` now documents that it propagates `NA` rather
  than assigning a category, so `table()` on its output can cover fewer rows
  than the input.
- `gnp_timecourse` carries a warning about the two similarly named expression
  column families in the source table: `max_expr_*` are TPM-like with a median
  around 170, `max.val.*` are `time_of_change()` outputs on a 0–1.8 scale.
  Pass the latter with a cutoff of 100 and every transcript comes back
  `"unexpressed"`.

222 testthat assertions, up from 139. `R CMD check` clean with one NOTE
(four `Suggests` not installed on the check machine).

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
