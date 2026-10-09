# Code review: chip_tools before restructuring

Findings from reading all 38 files of the original flat script collection,
plus the ~25 sibling files pulled in from `~/Dropbox/jp_seq`. Line references
are to the pre-restructure versions — `git log --follow` on any file in `R/`
reaches them.

Three automated audits backed the reading:

```r
# 1. every file parses
for (f in list.files("R", "[.]R$", full.names = TRUE)) parse(f)

# 2. undefined variables, per function
codetools::findGlobals(fn, merge = FALSE)$variables

# 3. unresolvable function calls, per function
codetools::findGlobals(fn, merge = FALSE)$functions
```

Audit 2 is what caught the hard-coded globals; audit 3 caught the undeclared
packages. Both are worth re-running after any edit.

---

## 1. Name collisions — silent, and the worst of it

Four function names were each defined in two files. Because the workflow was
`source()` every file in the directory, **whichever file loaded last won**, with
no warning. The outcome depended on filesystem order.

| Name | Defined in | Resolution |
| --- | --- | --- |
| `plot_tsne_kmeans` | `plot_tsne_kmeans.R`, `replot_plot_tsne_kmeans.R` | duplicate deleted |
| `prune_close` | `prune_close.R`, `prune_gene.R` | duplicate deleted |
| `time_expr_avg` | `time_expr_avg.R`, `expr_list_avg.R` | broken copy deleted |
| `trinarize_counts` | `trinarize_counts.R`, `trinarize.R` | merged |

Details:

- **`replot_plot_tsne_kmeans.R`** was a 307-line copy of `plot_tsne_kmeans.R`
  differing in **one character** — `brewer.pal(n = 7)` instead of `n = 9` — and
  defining the same function name. Both behaviours are preserved as the
  `palette_n` argument.
- **`prune_gene.R`** defined `prune_close()`, not `prune_gene()`. There was
  never a `prune_gene()` function. It was the older of the two copies, missing
  the `complete.cases()` guard the newer one has.
- **`expr_list_avg.R`** defined `time_expr_avg()` and read an undefined global
  `gene_num`; `time_expr_avg.R` is the same function with `gene_num` computed
  from the input. The broken copy loaded first alphabetically, so it happened
  to be overwritten by the working one. Luck, not design.
- **`chip_mclust_icl_all`** also existed twice: once as a top-level function and
  once as a stale nested copy inside `chip_segment()`.

## 2. Functions that could not run at all

### `chip_segment.R` — three independent faults

```r
for (col_name in names(histone_rlog_prune_mean_only)[-1])   # line 33
```

- Ignored both of its own arguments (`data_frame`, `G_range`) and read an
  undefined global `histone_rlog_prune_mean_only`.
- Used `combined_data` at line 83, fifty lines **before** creating it at line 88.
- Carried the stale nested copy of `chip_mclust_icl_all()` noted above.

Rewritten as the two-step pipeline it was reaching for:
`trim_upper_tail()` → `chip_mclust_icl_all()`.

### `binarize_counts_p.R` — fitted the wrong data

```r
model <- gammamixEM(x = scott_max_expr_nona_noz$max_expr_Scott_l, ...)   # line 21
```

The function's own `x` was computed and then discarded; the model was fitted to
a hard-coded global from one particular analysis. Every caller got that
analysis's answer regardless of what they passed in. (The newer copy in
`jp_seq` had already fixed this; that version is what the package carries.)

### `trinarize.R` — truncated, does not parse

The file ends mid-expression at line 64:

```r
axis.title.x = element_text(size = 16
```

`parse()` fails with `unexpected end of input`. It was never loadable.

### `chip_scatter_cutoffs.R` — undefined global, undeclared package

- `label_genes` (lines 41, 55) was read from the global environment; it is now
  a parameter.
- Calls `DescTools::Winsorize()` but only loaded `dplyr`, `ggplot2`, `ggExtra`
  and `ggrepel`.

### `get_co_clusters.R` — undefined global, four times

Read `fuzz_list` for the data-set names, the loop length and the row names
(lines 9, 11, 12, 15, 17, 19) while taking `memb_list` as its argument.
Everything now derives from `memb_list`.

### `topGO_split.R` — split by an undefined global

```r
gene_list <- split(t50_df, f = gnp_multi_split[, split_cat])   # line 23
```

Grouping came from `gnp_multi_split`, not the function's own `t50_df`.

### `outline_tsne.R`, `remove_intersect_all.R` — `orderPoints()` never declared

Both call `prevR::orderPoints()` without requiring `prevR`. This is the usual
reason the t-SNE outlining failed on a machine other than the one it was
written on.

## 3. Dead code after an early `return()`

| File | Dead region | What was lost |
| --- | --- | --- |
| `split_counts.R` | lines 52–151, after `return(best_model)` at line 50 | the entire cutoff/binarize/return-option body — ~100 lines, a verbatim copy of `binarize_counts.R` |
| `trinarize_counts.R` | everything after `return(model)` at line 69 | all four documented return options; only the raw `mixEM` object was reachable |

Both carried comments describing behaviour that could never execute
(`"handle return logic as in your original function"`).

Smaller instances:

- `rna_de_tpm2.R:97` computed `count_TPM` with `mapply()` over the TPM formula,
  then **line 101 overwrote it** with a `sweep()` of the raw counts. The first
  computation was dead; the returned TPM came from the sweep.
- `time_of_change.R:7` set `d1sum <- rowSums(diff1)`, overwritten at line 26 by
  `rowSums(df)`.
- `wins_norm_histones.R:58` computed `histone_rlog_medsubt`, never used it
  (see §5).
- `grouped_col_mean.R:15` computed `unique_counts`, never used it (see §4).

## 4. Wrong results

These ran to completion and returned plausible-looking wrong answers — the
most dangerous category.

### BIC was always wrong

`binarize_counts.R:108`, `split_counts.R:137`:

```r
n <- length(data)
```

`data` is the base R function. `length(data)` is **1**, so `log(n)` is 0 and
`bic == aic` for every model, always. Any model selection that compared them
was comparing a number to itself.

The parameter count was also wrong: `length(lambda) + 2 * length(mu)` double
counts the variance and ignores that mixing proportions have one fewer free
parameter than components. A 2-component mixture with `arbvar = FALSE` has
2 means + 1 shared variance + 1 free proportion = 4, not 6.

And in `split_counts.R:32`, BIC was penalised by `k` (the component count)
while AIC on the line above used `kv` (the parameter count) — so the two
criteria were on different scales and selection was biased toward larger `k`.

Fixed once in `mixture_fit_stats()`; pinned by a test.

### `grouped_col_mean()` misaligned columns

```r
grouped_counts <- counts[, which(!(group_index %in% unique_groups))]
...
rowMeans(grouped_counts[, which(group_index == x)])
```

`grouped_counts` is a **column subset**, but it is then indexed with positions
taken from the **unsubset** `group_index`. Whenever any condition had exactly
one replicate, the means were computed over the wrong columns. Single-replicate
conditions were separately computed into the unused `unique_counts` and
silently dropped from the output, and the output column names were built from
a vector of a different length than the number of output columns.

This function feeds `chip_de_rpm()` and both `RNA_de_tpm_*()`, so the error
propagated into the per-condition mean RPM columns of every differential
analysis run on a design with an unreplicated condition.

### `exclude_extreme` discarded the log transform

`binarize_counts.R:18`, `binarize_counts_p.R:18`, `split_counts.R:20`,
`trinarize.R:15`:

```r
x <- log10(counts + 1)            # transform
if (exclude_extreme) {
  x <- counts[counts > 0.1 & counts < 1]   # subset the RAW counts, overwriting x
}
```

With both flags set the mixture was fitted to raw counts on a `0.1 < x < 1`
window chosen for the log scale. On real count data that window is usually
empty or near-empty.

### `chip_mclust_icl_all()` returned all-NA columns

```r
above_cutoff_indices <- which(column_data > tail_cutoff)   # integer(0) if none
data_to_process <- column_data[-above_cutoff_indices]      # EMPTY, not all of x
```

`x[-integer(0)]` returns an empty vector, not `x`. Any column whose values all
fell below the tail cutoff was clustered on nothing, failed the
`length(...) > 0` check, and came back entirely `NA`. The same pattern
recurred on the assignment two lines later.

Also: `Mclust(mclust_data, G = best_g, model = best_model)` — `model` is not an
`Mclust()` argument. It was swallowed by `...` and the model family was
re-selected by BIC rather than fixed to the one ICL had chosen. Correct
argument is `modelNames`.

### `enhancer_quant()` ignored its distance cutoff

```r
enh_dist_summary(ch, promoter, enhancer)   # dist_cutoff not forwarded
```

The inner function fell back to its own `1e6` default, so the caller's
`dist_cutoff` had no effect on which pairs were considered.

### `remove_intersect_all()` measured distance on one axis

```r
df$p2pd <- c(0, sqrt((diff(df$y)^2) + (diff(df$y)^2)))
```

`diff(df$x)` never appears. Point-to-point distance was `|dy| * sqrt(2)`, so
the spur-trimming filter cut on the y axis alone.

### `df_rowfilt()` was a no-op

Built the aligned frame as `df_new`, then `return(df_filt)` — returned its own
input unchanged, having printed three dimensions on the way.

### topGO workers loaded the wrong package

All three topGO functions did:

```r
clusterEvalQ(cl, library(GOstats))
```

while the worker function `gostarter()` calls topGO. Every worker failed with a
missing-function error. Workers now load `topGO` and `org.Mm.eg.db`.

## 5. Arguments that did nothing

- **`wins_norm_histones(median_subtract =)`** — the median-subtracted matrix
  was computed into `histone_rlog_medsubt` and then the **original**
  `histone_rlog` was winsorised. Now wired correctly, but note it still cannot
  change the default output: min-max rescaling is shift-invariant, so a
  constant shift cancels exactly. A `normalize = FALSE` argument was added so
  the effect is reachable. Documented in `?wins_norm_histones`.
- **`chip_gene_group_comp(t50_cutoff = NULL)`** — filters `t50_cat` to
  `"Late"`/`"Mid"` and then tests for `"early"`/`"late"`, so the Wilcoxon test
  has one group and the p-value is `NA`. Left as-is, documented.
- **`time_expr_interp(meta =)`** — accepted and never read.
- **`remove_intersect_all(kern =)`**, **`ggroc(showAUC =)`** — same.
- **`expr_half_max_min()`** list branch — wraps the six-element result in
  `min(find_y0(j), ...)`, collapsing it to one number. Documented as
  matrix-only.

## 6. Machine-specific paths

Nine `source()` calls hard-coded `~/Dropbox/jp_seq`, so the repository was not
usable by anyone else — including the author on another machine:

```
line_scatter_facet.R:11   plot_cluster_line_EBseq.R
plot_tsne_kmeans.R:29     opti_map.R          # not even in the repo
replot_tsne_kmeans.R:41   opti_map.R
plot_group_timeline.R:9   plot_group_timeline.R   # sourced itself
plot_group_timeline.R:10  t50_to_compressed.R
rna_de_tpm3.R:17          counts_to_tpm.R
vsd_norm.R:3              grouped_col_mean.R
topGO_wrap.R:64           collapse_GO_columns.R
topGO_wrap.R:65           organize_genes_and_terms.R
```

`opti_map.R` was required by three files and was **not in the repository at
all**, so a fresh clone could not run the t-SNE figures. `plot_group_timeline.R`
sourced itself.

All nine targets are now package functions and all nine `source()` calls are
gone.

## 7. Undeclared packages

Audit 3 found calls with no matching `require()`. Beyond `orderPoints` and
`Winsorize` already noted:

- `str_c()` (stringr) in `enhancer_quant.R`, `chip_de_rpm.R` — replaced with
  base `paste0()`.
- `bind_rows()` (dplyr) in `ggroc.R`.
- `opti_map()` in `plot_tsne_optics.R` — called but, unlike
  `plot_tsne_kmeans.R`, never sourced, so it only worked if the other file had
  been loaded first.
- `rowMax`/`rowMin` (Biobase) in `line_scatter_facet.R`.

## 8. API and style

- **`class(x) == "try-error"`** in `binarize_counts.R:38`,
  `binarize_counts_p.R:38`, `split_counts.R:67`; **`class(roc) == "roc"`** in
  `ggroc.R:4`. `class()` can return a vector, making the comparison warn in
  R ≥ 4.2 (and silently use only the first element before that). Now
  `inherits()` / `tryCatch()`.
- **`ggroc.R:15`** constructed `simpleError(...)` without signalling it, so bad
  input fell through to an error about `roc_df` instead. Now `stop()`.
- **`ggroc.R:18`** assigned the plot over the function's own name.
- **`..density..`** (deprecated in ggplot2 3.4) → `after_stat(density)`;
  **`size =`** on lines → `linewidth =`.
- **Root bracketing** in `find.cutoff()`: `low <- x[which.min(f(x))]` is not
  guaranteed to bracket a sign change, which is why the
  `"warning: failed to find root"` fallback fired so often. Now bracketed
  between the two component means, where the crossing must lie.
- **Unconditional `print()`** in `grouped_col_mean()`, `min_finite()`,
  `load_merge()`, `chip_de_rpm()`, `soften_zero()` (which always drew a
  histogram). Now behind `verbose`/`plot` arguments.
- **`soften_zero()`** called `sample(x, n)` with a negative `n` whenever the
  zero bin was already below target, which errors. Now returns the input
  unchanged.
- **`plot_anova_tukey()`** hard-coded `+ cell_line` into the ANOVA formula, so
  it only worked on data with a column of that name. Now the `block_var`
  argument.
- **`plot_tsne_kmeans.R`, `replot_tsne_kmeans.R`** had their entire bodies
  indented two spaces at top level.

## 9. Attribution

Three functions are third-party and were unattributed or under-attributed.
`@source` tags added:

| Function | Origin |
| --- | --- |
| `counts_to_tpm()` | Kamil Slowikowski, gist `c6ab0348747f86e2748b` |
| `summarySE()` | Winston Chang, *Cookbook for R* |
| `spline.poly()` | whuber, GIS StackExchange 24827 |

`binarize_counts()` already credited its StackExchange origin in a comment;
that is now in the roxygen block.

## 10. Not carried over

| File | Why |
| --- | --- |
| `replot_plot_tsne_kmeans.R` | 1-char duplicate of `plot_tsne_kmeans.R`, same function name → `palette_n` |
| `prune_gene.R` | older duplicate of `prune_close()` |
| `expr_list_avg.R` | broken duplicate of `time_expr_avg()` |
| `trinarize.R` | truncated, does not parse; merged into `trinarize_counts()` |
| `chip_de_tpm.R` | byte-identical to `chip_de_rpm.R` but for the function name — and it computes RPM, so the name was wrong |
| `gostarter.R` | older GOstats version; all three topGO functions define their own `gostarter()` internally |
| `flatten_width.R` | 3-line stub with a body of `col` |
| `log2_norm_df.R` | a script, not a function |
| `egl_remove_high.R`, `elg_igl_normalize.R` | cerebellar imaging layer quantification, not ChIP/RNA |
| `memb2graph.R`, `multi_overlap*.R`, `prune_group_clusters.R` | fuzzy c-means membership lineage; belongs with the separate `fuzzycMerge` package |

One regression was found going the other way: the newer `rna_de_tpm3.R` in
`jp_seq` had `as.data.frdame()` for `as.data.frame()`, which errors on the
first call. The repository version was correct and is what the package carries.

## Worth doing next

Not done here, since this pass was restructure-and-fix rather than rewrite:

1. **Port `remove_intersect_all()` off `rgeos`** to `sf`. It is the only
   function that cannot run on current R.
2. **`chip_categorize_and_collapse()`** uses `rowwise()` with
   `across(everything())` — the obvious bottleneck on a large gene table.
3. **`enhancer_quant()`** compares every promoter against every enhancer per
   chromosome (its own comment says `not efficient !`). An interval join
   (`data.table::foverlaps()`, `GenomicRanges::findOverlaps()`) is the fix.
4. **`plot_cluster_line_EBseq()` / `_interp()` / `plot_group_timeline()`** are
   three ~260-line files whose pairwise diffs are under 90 lines. They want to
   be one function with arguments.
5. **`plot_anova_tukey()`** uses `dplyr::do()`, superseded by
   `group_modify()`/`nest()`.
6. **Single-species GO.** `org.Mm.eg.db` is hard-wired in all three topGO
   functions.
7. **Tidy evaluation.** The plotting functions use `get(column_name)` and bare
   column names throughout; `.data[[ ]]` and `{{ }}` would remove the remaining
   `codetools` noise and the `R CMD check` NOTEs about undefined globals.
