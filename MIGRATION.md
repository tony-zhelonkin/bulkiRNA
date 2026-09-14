# Migration: RNAseq-toolkit scripts → bulkiRNA package

For projects that `source()`d the toolkit from a submodule.

## API stability contract

`bulkirna_api()` returns the machine-readable public surface. Every export has a `lifecycle`:

- `stable` follows semantic versioning; incompatible changes require a major release after
  any applicable deprecation cycle.
- `experimental` may change or disappear without a major version bump.
- `deprecated` was a `0.x` tier. It is empty as of **1.0.0**: the 21 names were removed from
  the public API and their implementations kept as non-exported golden-baseline fixtures.

The separate `frozen` column recorded the 24 signatures inherited from the script library.
Twenty-one of those left the public API in 1.0.0, so it now lists the remaining 3. `superseded_by` may describe a sequence or technique,
and says plainly when part of an old behaviour has no replacement.

## 1. Replace the sourcing

```r
# before
rtk <- "01_modules/RNAseq-toolkit/scripts"
source(file.path(rtk, "GSEA/GSEA_processing/run_gsea.R"))
source(file.path(rtk, "GSEA/GSEA_plotting/gsea_dotplot.R"))
library(clusterProfiler)

# after
library(bulkiRNA)
```

GSEA now runs through `fgsea` directly, so `clusterProfiler` leaves the dependency list.

## 2. Every old name worked through 0.x, and none works in 1.0.0

Through the `0.x` series each old name was exported, warned once through `.Deprecated()`, and
named its migration path, so a project ran after step 1 alone. **In 1.0.0 they are gone from
the public API**, so step 1 is no longer sufficient and the table below is mandatory. A migration path may use several functions, and two deliberately note behaviour
with no successor. Step 4 covers the two changes that reach further.

| Old | New |
|---|---|
| `download_gatom_references(dest_dir = )` | `gatom_download_refs(dir = )` |
| `run_gsea()` | `gs_ranks()` + `gs_test()` |
| `run_gsea_analysis()` | `gs_ranks()` + `gsdb_msigdb()` + `gs_test()` + `gs_plot_*()` |
| `normalize_gsea_results()` | drop it — `gs_test()` already returns a tibble |
| `empty_gsea_tibble()` | `gs_test()`; filter its result to zero rows for an empty `gs_result`; there is no exported constructor, by design |
| `load_reference_db()` | `gsdb_load()` |
| `list_reference_dbs()` | `gsdb_list()` |
| `parse_gmx()` | `gsdb_from_file()` |
| `parse_mitoxplorer()` | `gsdb_load("mitoxplorer")`, or `gsdb_from_file()` for an arbitrary file |
| `list_to_term2gene()` | `gsdb_register()` |
| `filter_by_size()` | `min_size=` / `max_size=` on `gsdb_msigdb()`, `gsdb_load()` or `gsdb_from_file()` |
| `convert_human_to_mouse()` | `gsdb_msigdb(species = "Mus musculus", db_species = "HS")` |
| `gsea_dotplot()`, `gsea_dotplot_facet()` | `gs_plot_dot()` |
| `gsea_barplot()` | `gs_plot_bar()` |
| `gsea_running_sum_plot()` | `gs_plot_running()` |
| `plot_all_gsea_results()` | `gs_plot_dot()`, `gs_plot_bar()`, `gs_plot_running()` and `gs_save()` |
| `create_standard_volcano()` | `de_volcano()` |
| `create_MD_plot()` | `de_md_plot()` |
| `custom_minimal_theme_with_grid()` | `theme_bulki()` |
| `save_gsea_log()` | `gs_save()` for the plot/table artifact; the free-text log has no replacement |

## 3. The typical rewrite

```r
# before
res <- run_gsea(de_table, rank_metric = "t", species = "Mus musculus",
                collection = "C2", subcollection = "CP:REACTOME")
keep <- res@result$setSize >= 15 & res@result$setSize <= 500
res@result <- res@result[keep, ]
tbl <- normalize_gsea_results(res, database = "Reactome", contrast = co)

# after
ranks <- gs_ranks(de_table, metric = "t")
db    <- gsdb_msigdb("Mus musculus", collection = "C2", subcollection = "CP:REACTOME")
tbl   <- gs_test(ranks, db, min_size = 15, max_size = 500)
```

Three lines collapse into one because the size bound became an argument and the result
arrives already normalized.

## 4. What changes in the output

**Results are tibbles.** `gs_test()` and the `run_gsea()` shim both return a `gs_result`
tibble, so `@result` raises an error, and so do `@geneSets` and `@geneList`. Read the
columns directly.

| old | new |
|---|---|
| `NES` | `stat` |
| `pvalue` | `p_value` |
| `setSize` | `n_genes_tested` |
| `core_enrichment` (slash-joined string) | `leading_edge` (list column) |

**`direction` reads `"up"` and `"down"`.** A filter written as `direction == "Up"` matches
zero rows and reports nothing. Search for that comparison before trusting a re-run — it is
the one change that stays quiet.

**Size filtering happens before the test, which moves `padj`.** The old path adjusted across
the whole collection and filtered by `setSize` afterwards. `gs_test(min_size=, max_size=)`
filters first, so Benjamini-Hochberg runs over the family it actually tested. Expect `padj`
to shift wherever the size bound bites.

**The leading edge is fgsea's.** In every case measured so far it is a subset of the old
clusterProfiler one, which lifts `gene_ratio`.

**`pathway_name` reads better.** `format_pathway_name()` is a real function with a large
acronym dictionary. Names improve; ids stay identical, so joins keep working.

**`genes_full_set` moves to the caller.** It used to come from `gseaResult@geneSets`.
Recompute it from the `gs_db` you tested against:

```r
paste(intersect(db[[pathway_id]], universe), collapse = "/")
```

**Give migrated scripts a new cache filename.** A staleness check that compares names only
(`load_or_compute(expected_keys=)`) accepts an old cache full of `gseaResult`s and loads it
straight into new code. A fresh filename keeps the two eras apart.

## 5. Retiring the submodule

Once a project runs on the package and `bulkiRNA` is installed in the image:

```bash
git submodule deinit -f 01_modules/RNAseq-toolkit
git rm -f 01_modules/RNAseq-toolkit
```

Keep the submodule until then. Its pinned commit is how you read the old implementations.
