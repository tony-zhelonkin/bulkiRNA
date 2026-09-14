#' Build a normalized DGEList
#'
#' Wraps `edgeR::DGEList()` + `edgeR::calcNormFactors()` with the two checks
#' that catch the mistakes that actually happen: non-integer counts (edgeR
#' silently misbehaves on them, so they are rounded with a message or refused)
#' and sample/gene tables that are out of order with the count matrix.
#'
#' @param count_mat Numeric matrix or data frame of counts, genes in rows.
#' @param samples_df Sample metadata; `rownames()` must equal
#'   `colnames(count_mat)` in the same order.
#' @param genes_df Gene annotation with one row per row of `count_mat`, in the
#'   same order.
#' @param round_nonint Logical. Round non-integer counts (with a message) or
#'   stop.
#' @param norm_method Normalisation method passed to
#'   `edgeR::calcNormFactors()`.
#' @return A `DGEList` with normalisation factors computed.
#' @export
#' @examples
#' if (requireNamespace("edgeR", quietly = TRUE)) {
#'   counts <- matrix(1:12, nrow = 3,
#'     dimnames = list(paste0("G", 1:3), paste0("S", 1:4)))
#'   samples <- data.frame(group = c("a", "a", "b", "b"),
#'                         row.names = colnames(counts))
#'   build_dge(counts, samples, data.frame(gene = rownames(counts)))
#' }
build_dge <- function(count_mat, samples_df, genes_df,
                      round_nonint = TRUE, norm_method = "TMM") {
  .require_pkg("edgeR", "`build_dge()`", 'BiocManager::install("edgeR")')

  cm <- as.matrix(count_mat)
  storage.mode(cm) <- "numeric"

  if (any(abs(cm - round(cm)) > .Machine$double.eps^0.5)) {
    if (round_nonint) {
      message("Rounding non-integer counts for edgeR.")
      cm <- round(cm)
    } else {
      stop("Counts must be integers for edgeR; set `round_nonint = TRUE` to ",
           "round them.", call. = FALSE)
    }
  }
  storage.mode(cm) <- "integer"

  if (!identical(colnames(cm), rownames(samples_df))) {
    stop("`colnames(count_mat)` and `rownames(samples_df)` must match, in the ",
         "same order.", call. = FALSE)
  }
  if (nrow(cm) != nrow(genes_df)) {
    stop(sprintf("`genes_df` has %d row(s) but `count_mat` has %d.",
                 nrow(genes_df), nrow(cm)), call. = FALSE)
  }

  dge <- edgeR::DGEList(counts = cm, genes = genes_df, samples = samples_df)
  edgeR::calcNormFactors(dge, method = norm_method)
}

#' Annotate Ensembl gene IDs with symbols, Entrez IDs and biotype
#'
#' Maps version-stripped Ensembl gene IDs through the species `org.*.eg.db`
#' package, then -- when `use_biomart = TRUE` and the network is reachable --
#' refines the symbol and adds `gene_biotype` from Ensembl via biomaRt. The
#' biomaRt step is best effort: a failure leaves the org.db result intact
#' rather than aborting the annotation.
#'
#' Genes with no symbol keep their stable Ensembl ID as `Symbol`, so the column
#' is never `NA` and downstream joins and plots do not silently lose rows.
#'
#' @param ens_ids Character vector of Ensembl gene IDs, with or without
#'   version suffixes.
#' @param species A human or mouse alias accepted by [.species()]. Partial
#'   scientific names accepted by the historical `match.arg()` call remain
#'   supported.
#' @param use_biomart Logical. Attempt the biomaRt refinement.
#' @param input_gene_name Optional character vector of symbols supplied by the
#'   quantifier, either named by version-stripped ID (as produced by
#'   [read_counts_matrix()]) or parallel to `ens_ids`.
#' @param biomart_version Optional Ensembl release for `useEnsembl()`.
#' @param biomart_host Optional biomaRt host.
#' @return A tibble with one row per element of `ens_ids` and columns `Symbol`,
#'   `Ensembl`, `ENTREZID`, `gene_biotype`, `input_gene_name`.
#' @export
#' @examples
#' if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
#'   annotate_genes(c("ENSMUSG00000000001.5", "ENSMUSG00000000003"),
#'                  use_biomart = FALSE)
#' }
annotate_genes <- function(ens_ids,
                           species = c("Mus musculus", "Homo sapiens"),
                           use_biomart = TRUE,
                           input_gene_name = NULL,
                           biomart_version = NULL,
                           biomart_host = NULL) {
  if (identical(species, c("Mus musculus", "Homo sapiens"))) {
    species <- species[[1L]]
  }
  sp <- .species(species)
  species <- sp$scientific
  orgdb <- sp$orgdb
  dataset <- sp$biomart_dataset

  .require_pkg("AnnotationDbi", "`annotate_genes()`",
              'BiocManager::install("AnnotationDbi")')
  .require_pkg(orgdb, sprintf("`annotate_genes(species = \"%s\")`", species),
              sprintf('BiocManager::install("%s")', orgdb))

  stable <- sub("\\..*$", "", ens_ids)
  db <- getExportedValue(orgdb, orgdb)

  pull <- function(column) {
    unname(suppressMessages(AnnotationDbi::mapIds(
      db, keys = stable, column = column, keytype = "ENSEMBL",
      multiVals = "first")))
  }
  symbol <- pull("SYMBOL")
  entrez <- pull("ENTREZID")

  ann <- tibble::tibble(
    Symbol       = ifelse(!is.na(symbol) & nzchar(symbol), symbol, stable),
    Ensembl      = stable,
    ENTREZID     = entrez,
    gene_biotype = NA_character_
  )

  if (use_biomart) {
    .require_pkg("biomaRt", "`annotate_genes(use_biomart = TRUE)`",
                 'BiocManager::install("biomaRt")')
    bt <- try({
      margs <- list("genes", dataset = dataset)
      if (!is.null(biomart_version)) margs[["version"]] <- biomart_version
      if (!is.null(biomart_host))    margs[["host"]]    <- biomart_host
      mart <- do.call(biomaRt::useEnsembl, margs)
      biomaRt::getBM(
        attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name"),
        filters = "ensembl_gene_id", values = unique(stable), mart = mart)
    }, silent = TRUE)

    if (!inherits(bt, "try-error") && nrow(bt)) {
      bt <- bt[!duplicated(bt$ensembl_gene_id), , drop = FALSE]
      idx <- match(ann$Ensembl, bt$ensembl_gene_id)
      bm_symbol  <- bt$external_gene_name[idx]
      bm_biotype <- bt$gene_biotype[idx]
      ann$Symbol <- ifelse(!is.na(bm_symbol) & nzchar(bm_symbol),
                           bm_symbol, ann$Symbol)
      ann$gene_biotype <- bm_biotype
    }
  }

  ann$input_gene_name <- if (is.null(input_gene_name)) {
    NA_character_
  } else if (is.null(names(input_gene_name))) {
    as.character(input_gene_name)
  } else {
    unname(input_gene_name[ann$Ensembl])
  }

  ann[, c("Symbol", "Ensembl", "ENTREZID", "gene_biotype", "input_gene_name")]
}
