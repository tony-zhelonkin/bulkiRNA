# Workflows

Recipes. `?fn` carries the argument detail, and
[API_REFERENCE.md](API_REFERENCE.md) carries the function index.

## GSEA over several collections and contrasts

The pattern behind most pipelines: rank once per contrast, then loop the collections.

```r
library(bulkiRNA)

DBS <- list(
  Hallmark = list(collection = "H",  subcollection = NULL),
  Reactome = list(collection = "C2", subcollection = "CP:REACTOME"),
  GO_BP    = list(collection = "C5", subcollection = "GO:BP")
)

for (co in names(de_results)) {
  ranks <- gs_ranks(de_results[[co]], metric = "t")

  per_db <- lapply(names(DBS), function(nm) {
    db <- gsdb_msigdb("Mus musculus", collection = DBS[[nm]]$collection,
                      subcollection = DBS[[nm]]$subcollection)
    gs_test(ranks, db, contrast = co, min_size = 15, max_size = 500)
  })

  res <- do.call(rbind, per_db)      # rbind() works on gs_result
  gs_write(res, "results/gsea", name = co)
}
```

`contrast =` labels the rows, so results from many contrasts stack with `rbind()` and then
facet or compare inside one plot.

## Custom gene sets

```r
db  <- gsdb_from_file("pathways.gmt", database = "MyDB", species = "Mus musculus")
res <- gs_test(gs_ranks(de_table, metric = "t"), db, min_size = 10)
```

Every provider returns a `gs_db`, so custom sets travel the same route as MSigDB ones.

## Over-representation instead of preranked

Pass genes in place of ranks:

```r
sig <- rownames(de_table)[de_table$adj.P.Val < 0.05]
ora <- gs_test(sig, db)             # dispatches to fgsea::fora()
gs_plot_dot(ora, top_n = 20)
```

## Per-sample scores and a heatmap

```r
mat <- gs_score(expr, db, method = "gsva")   # gs_matrix; needs GSVA
gs_plot_heatmap(mat)
```

## Plotting

```r
res |> gs_filter(padj = 0.05) |> gs_top(20, by = "stat") |> gs_plot_dot()

gs_plot_dot(res, top_n = 10, facet = "database")       # one panel per database
gs_plot_dot(res, top_n = 10, direction = "up")         # up-regulated only
gs_plot_dot(ab, top_n = 10, compare = "contrast")      # contrasts side by side
gs_plot_running(res, ranks = ranks, db = db, top_n = 3) # takes ranks and sets
```

Save a plot together with the table behind it, so any figure traces back to its numbers:

```r
gs_save(p, "figures/hallmark_dot.pdf", width = 10, height = 8, data = res)
```

## Differential expression

```r
de_volcano(de_table, decision_by = "fdr", p_cutoff = 0.05, fc_cutoff = 1,
           highlight_gene = c("Il6", "Tnf"), title = "Treatment vs Control")

de_pca(dge, colour_by = "condition")
de_bfc_plot(de_table)                              # takes the DE table
de_md_plot(fit, coef = "Treatment_vs_Control")     # takes the limma fit
```

`de_volcano_grid()` combines several volcanoes into a row under one shared legend:

```r
de_volcano_grid(list(p1, p2, p3), labels = c("A", "B", "C"))
```

## Reproducibility

```r
write_session_provenance("results/provenance.txt")
```

`gs_write()` records provenance alongside its tables already.

## One gotcha

Reach for `gsdb_msigdb()` when mapping human MSigDB sets onto another species. The provider
carries a guard against an upstream caching bug that otherwise trims every collection after
the first in a session; the [README](../README.md) describes it. A direct
`msigdbr::msigdbr()` call runs without that guard.
