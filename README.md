# bulkiRNA

Bulk RNA-seq differential expression and gene-set analysis, as an installed R package.

**Version:** 0.5.0 · **License:** MIT · **Author:** Anton Zhelonkin

This package succeeds `RNAseq-toolkit`, a folder of scripts each project `source()`d from a
submodule. There is now one installed copy and one line to reach it: `library(bulkiRNA)`.
The old function names still work — see [Legacy names](#legacy-names).

---

## Install

```r
remotes::install_github("tony-zhelonkin/bulkiRNA@v0.5.0")
```

Pin a ref. The package lives on tagged releases and on the `feat/bulkirna-package` branch,
so `@v0.5.0` or `@feat/bulkirna-package` both resolve.

The hard dependencies are small: `msigdbr`, `fgsea`, `ggplot2`, `dplyr` and a handful of
base-adjacent packages. `limma`, `edgeR`, `GSVA`, `gatom`, `org.*.eg.db` and `plotly` sit
in `Suggests`, so install each one when you reach the feature that uses it.

To see the whole optional set at once — what is present, which version, and the exact command
for anything missing:

```r
bulkirna_check_deps("all")                 # or "de", "annotation", "scoring", "network", ...
bulkirna_check_deps("de", error = TRUE)    # stops when something is missing; for CI
```

---

## Quick start

```r
library(bulkiRNA)

# 1. rank genes from a limma topTable (rownames = symbols)
ranks <- gs_ranks(de_table, metric = "t")

# 2. pick a gene-set database
db <- gsdb_msigdb("Mus musculus", collection = "H")

# 3. test
res <- gs_test(ranks, db, min_size = 15, max_size = 500)

# 4. plot
gs_plot_dot(res, top = 20)
```

`res` is a tibble with class `gs_result`. Read its columns directly, pipe it through `dplyr`
verbs, and combine results with `rbind()`.

```r
res |> gs_filter(padj = 0.05) |> gs_top(10, by = "stat")
gs_write(res, "results/gsea", name = "hallmark")   # tables + provenance
gs_read("results/gsea", name = "hallmark")
```

To hand results to a downstream master table, `gs_to_master()` serializes them to a versioned
schema and `gs_validate_master()` checks one that already exists. Both work on tables, not
files, so where the table lives stays your decision:

```r
master <- gs_to_master(res, db = db, universe = names(ranks))
gs_validate_master(master)     # every problem at once, or invisible(df)
```

The validator earns its keep on the failure that is otherwise invisible: a derived column
NA-filled beside a finite `padj`, which is what a column allowlist plus `rbind()` produces.

---

## How it fits together

Four layers, in dependency order. Each layer has one job: providers supply gene sets,
compute functions return data, and renderers turn that data into plots.

| Layer | Functions | Produces |
|---|---|---|
| **Providers** | `gsdb_msigdb()` `gsdb_load()` `gsdb_from_file()` `gsdb_register()` `gsdb_list()` `gsdb_info()` | `gs_db` |
| **Compute** | `gs_test()` `gs_score()` `gs_coregulation()` | `gs_result`, `gs_matrix` |
| **Result ops** | `gs_filter()` `gs_top()` `gs_split()` `gs_leading_edge()` `gs_read()` `gs_write()` `gs_to_master()` `gs_validate_master()` | `gs_result`, `tibble` |
| **Renderers** | `gs_plot_dot()` `gs_plot_bar()` `gs_plot_heatmap()` `gs_plot_running()` | `ggplot` |

`gs_test()` reads its input and picks the matching method: a named numeric vector of ranks
runs preranked GSEA through fgsea, a character vector of genes runs over-representation
through `fora()`, and a `gs_matrix` runs the matrix path. `gs_stat_types()` names the
statistic each method returns.

`gs_coregulation()` runs GESECA on a genes x samples matrix. It returns an
unsigned `pct_var` result, with `direction = NA` rather than labelling every
positive percentage as up.

`gs_score()` is the per-sample counterpart (GSVA, ssGSEA). It returns a `gs_matrix`, which
`gs_plot_heatmap()` consumes.

### Gene-set sources

```r
gsdb_msigdb("Mus musculus", collection = "C2", subcollection = "CP:REACTOME")
gsdb_from_file("custom.gmt", database = "MyDB")   # .gmt / .gmx
gsdb_load("mitoxplorer")                          # registered reference DBs
```

Every provider returns a `gs_db`, so the compute and plotting layers treat all sources
alike.

### Differential expression and QC

| Function | Purpose |
|---|---|
| `de_volcano()` `de_volcano_grid()` | volcano; FDR decides significance, raw p sets the y-axis |
| `de_pca()` `de_pca_3d()` | PCA from a `DGEList` |
| `de_md_plot()` `de_bfc_plot()` | mean-difference from a fit; logFC vs B-statistic |
| `build_dge()` `read_counts_matrix()` `read_metadata()` `annotate_genes()` | inputs |

### Utilities

`theme_bulki()` (publication ggplot2 theme) · `gs_save()` (a plot with the table behind it)
· `format_pathway_name()` (biological capitalisation from a ~400-term dictionary) ·
`write_session_provenance()` · `ensure_dir()` · `bulkirna_check_deps()` ·
`bulkirna_api()` (machine-readable lifecycle and signature-freeze registry) ·
`bulkirna_stochastic()` (which functions consume randomness, and each one's seed)

`write_session_provenance()` records the `bulkiRNA` version, every hard-dependency version,
the bundled registry version, `RNGkind()` and the stochastic seed defaults, and any shared
reference-data snapshot resolved this session —
the package version is the unit of reproducibility, so it is stated outright rather than
inferred from `sessionInfo()`.

### Network modules

`gatom_de()` `gatom_download_refs()` `gatom_genes()` `gatom_module()`
`gatom_refs()` `gatom_save_html()` — these need `gatom` from Bioconductor and `mwcsr`
from CRAN.

---

## Legacy names

**Removed from the public API in 1.0.0.** These 21 names were shims through the `0.x`
series. They are no longer exported, so `library(bulkiRNA)` will not find them:

`run_gsea()` `run_gsea_analysis()` `normalize_gsea_results()` `gsea_dotplot()`
`gsea_dotplot_facet()` `gsea_barplot()` `gsea_running_sum_plot()` `plot_all_gsea_results()`
`create_standard_volcano()` `create_MD_plot()` `custom_minimal_theme_with_grid()`
`load_reference_db()` `list_reference_dbs()` `filter_by_size()` `parse_gmx()`
`parse_mitoxplorer()` `list_to_term2gene()` `convert_human_to_mouse()`
`empty_gsea_tibble()` `save_gsea_log()` `download_gatom_references()`

Their implementations survive as non-exported fixtures. The golden harness captured its
baseline at `752481f`, before any restructuring, and 17 of its 20 cases run through these
functions — so they are what proves today's numbers still match the original script
library. [MIGRATION.md](MIGRATION.md) names the replacement for each.

Two changes reach past the shims, because the return type itself changed:

- **Results are tibbles.** `run_gsea()` returns a `gs_result`, so `@result` raises an error.
  Pass size bounds to `gs_test(min_size=, max_size=)`, and read `stat`, `p_value` and
  `n_genes_tested` in place of `NES`, `pvalue` and `setSize`.
- **`direction` reads `"up"` and `"down"`.** A filter written as `direction == "Up"`
  matches zero rows and says nothing about it. Search for that comparison when migrating.

[MIGRATION.md](MIGRATION.md) carries the full mapping.

---

## A caveat worth knowing

`msigdbr` (as of 26.1.0) keys its ortholog cache on the species alone. When it maps human
sets to another species, the first collection of a session fixes the gene space, and every
later collection is trimmed to fit it. Reactome queried after Hallmark keeps 3,688 of its
10,762 mapped genes; GO:BP keeps 4,313 of 15,988. The run completes normally and reports
under-tested sets with an inflated `padj`.

`gsdb_msigdb()` handles this in two steps: it clears the stale cache before each query, then
measures the result's gene coverage and raises an error when the numbers look trimmed. The
second step stands on its own, so the guard holds even after upstream changes its internals.
Reach for `gsdb_msigdb()` to get that protection.

---

## Development

R runs in a container.

```bash
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/cache \
  -v /path/to/.msigdb-cache:/cache -v "$PWD":/pkg -w /pkg scdock-r-dev:v0.5.13 \
  Rscript -e 'devtools::test(".")'
```

`--user` grants write permissions and `HOME` gives msigdbr its runtime cache. Add
`--network host` to fill a cold MSigDB cache.

Two gates guard every commit:

```r
devtools::test(".")
```
```bash
Rscript tests/golden/verify_golden.R    # 20/20, exit 0
```

`tests/golden/` holds rendered-output baselines that catch silent visual regressions.
Refresh a single baseline with `capture_golden.R --cases=<name>`; a bare run rewrites all
20.

`NAMESPACE` is generated from roxygen comments — refresh it with `devtools::document()`.

[CONVENTIONS.md](CONVENTIONS.md) holds the house style. [docs/](docs/) holds the extended
documentation.

---

## License

MIT — Copyright (c) 2025 Anton Zhelonkin. See [LICENSE.md](LICENSE.md).
