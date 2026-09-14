write_tsv <- function(df, ext = ".tsv", sep = "\t") {
  f <- tempfile(fileext = ext)
  utils::write.table(df, f, sep = sep, row.names = FALSE, quote = FALSE)
  f
}

test_that("featureCounts-processed shape loads", {
  f <- write_tsv(data.frame(Geneid = c("g1", "g2"), S1 = c(1, 2), S2 = c(3, 4)))
  m <- read_counts_matrix(f)
  expect_equal(dim(m), c(2L, 2L))
  expect_equal(rownames(m), c("g1", "g2"))
  expect_true(is.na(attr(m, "input_gene_name")))
})

test_that("raw featureCounts metadata columns are dropped", {
  f <- write_tsv(data.frame(Geneid = c("g1", "g2"), Chr = "1", Start = 1,
                            End = 2, Strand = "+", Length = 2,
                            S1 = c(1, 2), S2 = c(3, 4)))
  m <- read_counts_matrix(f)
  expect_equal(colnames(m), c("S1", "S2"))
})

test_that("raw featureCounts output with its real leading `#` comment line loads", {
  f <- tempfile(fileext = ".txt")
  writeLines(c(
    '# Program:featureCounts v2.0.1; Command:"featureCounts" "-a" "g.gtf" "-o" "counts.txt" "s1.bam" "s2.bam"',
    paste(c("Geneid", "Chr", "Start", "End", "Strand", "Length",
            "s1.bam", "s2.bam"), collapse = "\t"),
    paste(c("g1", "1", "1", "100", "+", "100", "3", "5"), collapse = "\t"),
    paste(c("g2", "1", "200", "300", "-", "100", "4", "6"), collapse = "\t")
  ), f)
  m <- read_counts_matrix(f)
  expect_equal(dim(m), c(2L, 2L))
  expect_equal(rownames(m), c("g1", "g2"))
  expect_equal(colnames(m), c("s1", "s2"))
  expect_equal(unname(m["g1", ]), c(3, 5))
})

test_that("Salmon gene-level input keeps gene_name as an attribute", {
  f <- write_tsv(data.frame(gene_id = c("ENSG1.3", "ENSG2.1"),
                            gene_name = c("A", "B"),
                            S1 = c(1, 2), S2 = c(3, 4)))
  m <- read_counts_matrix(f)
  ign <- attr(m, "input_gene_name")
  expect_equal(ign[["ENSG1"]], "A")
  expect_equal(colnames(m), c("S1", "S2"))
})

test_that("sample columns lose paths and aligner suffixes", {
  f <- write_tsv(data.frame(Geneid = "g1", check.names = FALSE))
  df <- data.frame(Geneid = "g1", a = 1, b = 2)
  names(df) <- c("Geneid", "/tmp/aln/S1.bam", "S2.sorted")
  f <- write_tsv(df)
  expect_equal(colnames(read_counts_matrix(f)), c("S1", "S2"))
})

test_that("an unsupported shape errors", {
  f <- write_tsv(data.frame(foo = 1, bar = 2))
  expect_error(read_counts_matrix(f), "Unsupported counts file format")
  expect_error(read_counts_matrix(tempfile()), "does not exist")
})

test_that("the separator is sniffed, not guessed from the extension", {
  # comma-separated data in a .tsv-named file
  f <- write_tsv(data.frame(Geneid = c("g1", "g2"), S1 = 1:2, S2 = 3:4),
                 ext = ".tsv", sep = ",")
  expect_equal(dim(read_counts_matrix(f)), c(2L, 2L))
})

test_that("read_metadata canonicalises the sample column", {
  df <- data.frame(a = c("S1", "S2"), group = c("WT", "KO"))
  names(df)[1] <- "Sample ID"
  f <- write_tsv(df, ext = ".csv", sep = ",")
  md <- read_metadata(f)
  expect_true("Sample_ID" %in% names(md))
  expect_equal(md$Sample_ID, c("S1", "S2"))
})

test_that("read_metadata reports what it needed and what it found", {
  f <- write_tsv(data.frame(x = 1, y = 2), ext = ".csv", sep = ",")
  expect_error(read_metadata(f), "sample-ID column")

  f2 <- write_tsv(data.frame(Sample_ID = "S1", group = "WT"),
                  ext = ".csv", sep = ",")
  expect_error(read_metadata(f2, required_cols = c("batch", "sex")), "batch")
})

test_that(".aggregate_duplicate_ids sums duplicate rows and keeps gene names", {
  m <- matrix(1:8, nrow = 4, dimnames = list(c("b", "a", "b", "c"), NULL))
  attr(m, "input_gene_name") <- c(b = "B1", a = "A", b2 = "B2", c = NA)
  out <- .aggregate_duplicate_ids(m)

  expect_equal(rownames(out), c("a", "b", "c"))   # split() order, documented
  expect_equal(unname(out["b", ]), c(1L + 3L, 5L + 7L))
  expect_equal(unname(attr(out, "input_gene_name")[["b"]]), "B1")
})

test_that(".aggregate_duplicate_ids is a no-op without duplicates", {
  m <- matrix(1:4, nrow = 2, dimnames = list(c("a", "b"), NULL))
  expect_identical(.aggregate_duplicate_ids(m), m)
})

test_that("write_session_provenance writes sessionInfo and creates parents", {
  .clear_ref_resolutions()
  f <- file.path(tempdir(), "prov-subdir", "provenance.txt")
  unlink(dirname(f), recursive = TRUE)
  expect_equal(write_session_provenance(f, genome_build = "mm10"), f)

  lines <- readLines(f)
  expect_true(any(grepl("Genome build: mm10", lines)))
  expect_true(any(grepl("^bulkiRNA version: ", lines)))
  expect_true(any(grepl("^  dplyr: ", lines)))
  expect_true(any(grepl(
    "^Bundled reference registry version: 1\\.0\\.0$", lines
  )))
  expect_true(any(grepl(
    "^Reference data: none resolved this session$", lines
  )))
  expect_true(any(grepl("--- sessionInfo ---", lines)))
})

test_that("write_session_provenance records resolved snapshot tags", {
  .clear_ref_resolutions()
  on.exit(.clear_ref_resolutions(), add = TRUE)
  root <- withr::local_tempdir()
  source_dir <- file.path(root, "coresh")
  snapshot <- file.path(source_dir, "coresh_snapshot_test")
  dir.create(snapshot, recursive = TRUE)
  linked <- suppressWarnings(
    file.symlink(snapshot, file.path(source_dir, "current"))
  )
  if (!isTRUE(linked)) {
    skip("This filesystem does not support symbolic links")
  }
  withr::local_envvar(c(REFCACHE_ROOT = root))

  .ref_path("coresh")
  f <- tempfile(fileext = ".txt")
  write_session_provenance(f)

  lines <- readLines(f)
  expect_true(any(grepl("coresh: coresh_snapshot_test", lines, fixed = TRUE)))
})
