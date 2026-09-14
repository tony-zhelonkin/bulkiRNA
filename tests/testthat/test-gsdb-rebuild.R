test_that("human-to-mouse conversion preserves provenance for retained sets", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  skip_if_not_installed("homologene")
  testthat::local_mocked_bindings(
    homologene = function(genes, inTax, outTax) {
      data.frame(human = "ACTB", mouse = "Actb")
    },
    .package = "homologene"
  )
  provenance <- list(source = "test", snapshot = "snapshot-1")
  set_provenance <- tibble::tibble(
    set_name = c("A", "B"),
    source_row = c(1L, 2L)
  )
  db <- bulkiRNA:::gs_db(
    list(A = "ACTB", B = "UNMAPPED"),
    database = "human_sets",
    species = "Homo sapiens",
    pathway_names = c(A = "Set A", B = "Set B"),
    pathway_descriptions = c(A = "Description A", B = "Description B"),
    set_provenance = set_provenance,
    provenance = provenance
  )

  mouse <- bulkiRNA:::.gsdb_human_to_mouse(db)

  expect_identical(names(mouse), "A")
  expect_identical(attr(mouse, "provenance"), provenance)
  expect_identical(
    attr(mouse, "set_provenance"),
    set_provenance[1L, , drop = FALSE]
  )
  expect_identical(attr(mouse, "pathway_descriptions"),
                   c(A = "Description A"))
})
