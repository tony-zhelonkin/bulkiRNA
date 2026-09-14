# API reference

The installed help pages carry the per-function detail, generated from the source and always
current:

```r
help(package = "bulkiRNA")
?gs_test
```

This file is the index: what exists, and which layer it belongs to.

## Gene-set providers → `gs_db`

| Function | Purpose |
|---|---|
| `gsdb_msigdb()` | MSigDB collection, with ortholog mapping on request |
| `gsdb_from_file()` | `.gmt` / `.gmx` file |
| `gsdb_load()` | bundled reference database by name |
| `gsdb_list()` | the registered databases as a table |
| `gsdb_info(name)` | metadata for one registered database, by name |
| `gsdb_register()` | add a database to the registry |

```r
gsdb_list()                    # mitopathways, mitoxplorer, mito_unified, ...
gsdb_info("mitoxplorer")       # name, description, source_url, species
gsdb_load("mitoxplorer")       # the gs_db itself
```

## Compute → `gs_result` / `gs_matrix`

| Function | Purpose |
|---|---|
| `gs_ranks()` | DE table → named ranking vector |
| `gs_test()` | preranked GSEA (fgsea), over-representation (`fora`), or matrix input |
| `gs_score()` | per-sample scores (GSVA, ssGSEA) → `gs_matrix` |
| `gs_stat_types()` | the statistic each method returns |

`gs_test()` reads the class of its input and picks the method: a named numeric vector runs
preranked GSEA, a character vector runs over-representation, and a `gs_matrix` runs the
matrix path.

## Result operations

| Function | Purpose |
|---|---|
| `gs_filter()` | by `padj`, `p_value`, `stat`, `direction`, `database`, `contrast`, size, or name pattern |
| `gs_top()` | top n, per group or per direction on request |
| `gs_split()` | split by a column |
| `gs_leading_edge()` | leading-edge genes |
| `gs_write()` `gs_read()` | round-trip tables together with provenance |
| `gs_to_master()` | serialize to the versioned master-table schema; returns a tibble, writes nothing |
| `gs_validate_master()` | check a master table against the schema, reporting every problem at once |

A `gs_result` is a tibble. `dplyr` verbs, `rbind()`, `summary()` and `as_tibble()` all work
on it, and `[` carries its attributes through.

## Renderers → `ggplot`

| Function | Purpose |
|---|---|
| `gs_plot_dot()` | dotplot; faceting, up/down split, comparison, highlighting |
| `gs_plot_bar()` | statistic barplot |
| `gs_plot_running()` | running enrichment score |
| `gs_plot_heatmap()` | pathway × sample or pathway × contrast heatmap |

Each renderer takes a `gs_result` or `gs_matrix` and returns a `ggplot`.

## Differential expression

| Function | Purpose |
|---|---|
| `de_volcano()` | volcano; FDR decides significance, raw p sets the y-axis |
| `de_volcano_grid()` | several volcanoes in a row, sharing one legend |
| `de_pca()` `de_pca_3d()` | PCA from a `DGEList` |
| `de_md_plot()` | mean-difference (MA), from a limma fit and a `coef` |
| `de_bfc_plot()` | logFC vs B-statistic, from a DE table |

## Inputs and annotation

`read_counts_matrix()` · `read_metadata()` · `build_dge()` · `annotate_genes()`

## Utilities

`theme_bulki()` · `gs_save()` (a plot with the table behind it) · `format_pathway_name()` ·
`write_session_provenance()` · `ensure_dir()` · `bulkirna_check_deps()`

`bulkirna_check_deps()` reports the optional-dependency closure — presence, version, and the
correct install command — grouped by feature area (`de`, `annotation`, `scoring`, `network`,
`plots`, `io`). `error = TRUE` makes it a preflight suitable for a CI or image smoke test.

## Network modules (need `gatom` and `mwcsr`)

`gatom_refs()` · `gatom_de()` · `gatom_download_refs()` · `gatom_genes()` ·
`gatom_module()` · `gatom_save_html()`

## Deprecated

21 shims for superseded API names. They work and they warn.
[../MIGRATION.md](../MIGRATION.md) maps each one to its replacement.
