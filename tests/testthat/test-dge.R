test_that("build_dge returns a normalized DGEList", {
  skip_if_not_installed("edgeR")
  dge <- fake_dge()
  expect_true(methods::is(dge, "DGEList"))
  expect_true(all(dge$samples$norm.factors > 0))
  expect_equal(nrow(dge), 40L)
})

test_that("non-integer counts are rounded with a message, or refused", {
  skip_if_not_installed("edgeR")
  counts <- matrix(c(1.4, 2.6, 3.1, 4.9), nrow = 2,
                   dimnames = list(c("g1", "g2"), c("S1", "S2")))
  samples <- data.frame(group = c("a", "b"), row.names = colnames(counts))
  genes <- data.frame(gene = rownames(counts))

  expect_message(dge <- build_dge(counts, samples, genes), "Rounding")
  expect_equal(unname(dge$counts[, 1]), c(1L, 3L))
  expect_error(build_dge(counts, samples, genes, round_nonint = FALSE),
               "must be integers")
})

test_that("mismatched sample or gene tables are an error", {
  skip_if_not_installed("edgeR")
  counts <- matrix(1:4, nrow = 2,
                   dimnames = list(c("g1", "g2"), c("S1", "S2")))
  bad_samples <- data.frame(group = c("a", "b"), row.names = c("S2", "S1"))
  genes <- data.frame(gene = rownames(counts))

  expect_error(build_dge(counts, bad_samples, genes), "same order")
  expect_error(
    build_dge(counts, data.frame(group = c("a", "b"),
                                 row.names = c("S1", "S2")),
              genes[1, , drop = FALSE]),
    "`genes_df` has 1 row")
})

test_that("annotate_genes strips versions and never returns an NA Symbol", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("AnnotationDbi")
  ids <- c("ENSMUSG00000000001.5", "ENSMUSG00000000003", "ENSMUSG99999999999")
  ann <- annotate_genes(ids, use_biomart = FALSE)

  expect_equal(nrow(ann), 3L)
  expect_equal(ann$Ensembl, sub("\\..*$", "", ids))
  expect_false(anyNA(ann$Symbol))
  # An unmappable ID falls back to its own stable ID rather than NA.
  expect_equal(ann$Symbol[3], "ENSMUSG99999999999")
  expect_equal(names(ann), c("Symbol", "Ensembl", "ENTREZID", "gene_biotype",
                             "input_gene_name"))
})

test_that("annotate_genes threads a named input_gene_name through", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("AnnotationDbi")
  ids <- c("ENSMUSG00000000001.5", "ENSMUSG00000000003")
  ign <- c(ENSMUSG00000000001 = "Quant1", ENSMUSG00000000003 = "Quant2")
  ann <- annotate_genes(ids, use_biomart = FALSE, input_gene_name = ign)
  expect_equal(ann$input_gene_name, c("Quant1", "Quant2"))

  ann2 <- annotate_genes(ids, use_biomart = FALSE)
  expect_true(all(is.na(ann2$input_gene_name)))
})

test_that("annotate_genes validates species", {
  expect_error(
    annotate_genes("ENSMUSG1", species = "rat"),
    "`species` must be one of"
  )
})

test_that("annotate_genes uses human annotation for aliases and partial matches", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("AnnotationDbi")

  common <- annotate_genes("ENSG00000141510", species = "human",
                           use_biomart = FALSE)
  partial <- annotate_genes("ENSG00000141510", species = "Homo",
                            use_biomart = FALSE)

  expect_identical(common$Symbol, "TP53")
  expect_identical(common$ENTREZID, "7157")
  expect_identical(partial$Symbol, "TP53")
  expect_identical(partial$ENTREZID, "7157")
  expect_identical(common, partial)
})

test_that("annotate_genes rejects an absent species", {
  expect_error(
    annotate_genes("ENSG00000141510", species = NULL, use_biomart = FALSE),
    "`species` must be one of"
  )
})
