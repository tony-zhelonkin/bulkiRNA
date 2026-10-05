#' GATOM metabolic-network modules
#'
#' @description
#' GATOM is three calls: `gatom::makeMetabolicGraph()` builds an
#' atom-resolved metabolic graph, `gatom::scoreGraph()` scores it, and
#' `mwcsr::solve_mwcsp()` finds the maximum-weight connected subgraph. The
#' wrapper is the same three layers, each a pure function that checks its
#' inputs and records what it did as attributes on its result:
#'
#' ```
#' gatom_refs() + gatom_de() -> gatom_graph() -> gatom_score() -> gatom_solve()
#'                                                                gatom_solver()
#' ```
#'
#' [gatom_module()] runs the three in order. [gatom_genes()] and
#' [gatom_pathways()] read a module; [gatom_save_html()] writes one.
#'
#' The traps the raw API only documents are enforced here:
#'
#' 1. `pval` must be **raw** p-values -- BUM scoring degrades silently on FDR
#'    values. [gatom_de()] refuses an adjusted-p column.
#' 2. gatom keeps the 12,000 genes ranked highest by `baseMean` and drops the
#'    rest before mapping. [gatom_graph()] reports how many genes that
#'    removed.
#' 3. The genes live on the graph **edges**, not the vertices.
#'    [gatom_genes()] reads the edges.
#' 4. `met.db` is required even when `met.de` is `NULL`. [gatom_refs()] carries
#'    it, so [gatom_graph()] cannot omit it.
#' 5. A p-value distribution that does not fit a beta-uniform mixture makes
#'    `scoreGraph()` set every score to 0 with only a warning. [gatom_score()]
#'    stops instead.
#'
#' @name gatom
#' @keywords internal
NULL

#' Load the GATOM reference files
#'
#' GATOM needs a network, its metabolite database, a species annotation and,
#' for the Rhea, combined and lipid networks, a supplementary gene-to-reaction
#' table. They are read from one directory, `dir`, under the file names
#' [gatom_download_refs()] writes. When `dir` holds a `SHA256SUMS`, each file
#' is verified against it.
#'
#' The returned object carries `met.db` alongside the network, because
#' `gatom::makeMetabolicGraph()` needs `met.db` even when `met.de = NULL`, and
#' the `gene2reaction` table, because the vignette passes it with every
#' network that has one.
#'
#' @param species Character(1) human or mouse species alias.
#' @param dir Character(1) directory holding the reference files.
#' @param download Deprecated. Call [gatom_download_refs()] first instead.
#' @param network Character(1) network: `"kegg"`, `"rhea"`, `"combined"`
#'   (KEGG + Rhea + BiGG transport) or `"lipids"` (the Rhea lipid
#'   subnetwork).
#' @return An object of class `gatom_refs`: a list with `network`, `met_db`,
#'   `org_anno`, `gene2reaction_extra` (`NULL` for KEGG), `species`,
#'   `network_name`, `files` (paths), `sha256` and `checksums` (`"verified"`
#'   or `"no SHA256SUMS"`).
#' @examples
#' \dontrun{
#' gatom_download_refs(dir = "00_data/references/gatom",
#'                     species = "Mus musculus", networks = "combined")
#' refs <- gatom_refs("Mus musculus", network = "combined")
#' }
#' @export
gatom_refs <- function(species = "Homo sapiens",
                       dir = "00_data/references/gatom",
                       download = FALSE, network = "kegg") {
  sp <- .species(species)
  .gatom_check_network(network)
  if (!is.character(dir) || length(dir) != 1L || is.na(dir) || !nzchar(dir)) {
    stop("`dir` must be a single directory path.", call. = FALSE)
  }
  if (!is.logical(download) || length(download) != 1L || is.na(download)) {
    stop("`download` must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(download)) {
    warning("`download` is deprecated and will be removed in bulkiRNA 2.0.0. ",
            "Call gatom_download_refs(dir = \"", dir, "\") once, then ",
            "gatom_refs(dir = \"", dir, "\").", call. = FALSE)
    gatom_download_refs(dir = dir, species = sp$gatom_download,
                        networks = network)
  }

  wanted <- .gatom_ref_files(network, sp)
  found <- file.path(dir, wanted)
  names(found) <- names(wanted)
  if (!all(file.exists(found))) {
    stop("GATOM reference file(s) not found in `dir` (", dir, "): ",
         paste(sprintf("`%s`", wanted[!file.exists(found)]), collapse = ", "),
         ".\nFetch them with gatom_download_refs(dir = \"", dir,
         "\", species = \"", sp$gatom_download, "\", networks = \"", network,
         "\").", call. = FALSE)
  }

  sha <- .gatom_sha256(found)
  names(sha) <- names(found)
  checksums <- .gatom_verify_sums(dir, wanted, sha)

  extra <- NULL
  if ("gene2reaction_extra" %in% names(found)) {
    extra <- utils::read.delim(found[["gene2reaction_extra"]],
                               colClasses = "character")
    if (!all(c("gene", "reaction") %in% names(extra))) {
      stop("`dir` holds a malformed `", wanted[["gene2reaction_extra"]],
           "`: it needs `gene` and `reaction` columns. Refetch it with ",
           "gatom_download_refs(overwrite = TRUE).", call. = FALSE)
    }
  }

  structure(
    list(
      network             = readRDS(found[["network"]]),
      met_db              = readRDS(found[["met_db"]]),
      org_anno            = readRDS(found[["org_anno"]]),
      gene2reaction_extra = extra,
      species             = sp$scientific,
      network_name        = network,
      files               = found,
      sha256              = sha,
      checksums           = checksums
    ),
    class = "gatom_refs"
  )
}

#' Verify loaded GATOM files against the directory's `SHA256SUMS`
#'
#' @param dir Character(1) reference directory.
#' @param wanted Named character vector of file names.
#' @param sha Named character vector of their hashes, in `wanted` order.
#' @return `"verified"` when `SHA256SUMS` exists and lists every file with
#'   the same hash; `"no SHA256SUMS"` when it does not exist. Stops on a
#'   mismatch or an unlisted file.
#' @keywords internal
.gatom_verify_sums <- function(dir, wanted, sha) {
  sums_file <- file.path(dir, "SHA256SUMS")
  if (!file.exists(sums_file)) return("no SHA256SUMS")
  if (anyNA(sha)) {
    stop("`dir` has a SHA256SUMS that this R cannot check: ",
         "tools::sha256sum() needs R >= 4.5.", call. = FALSE)
  }
  lines <- readLines(sums_file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  listed <- stats::setNames(sub("^([0-9a-fA-F]+)\\s+\\*?.*$", "\\1", lines),
                            sub("^[0-9a-fA-F]+\\s+\\*?", "", lines))
  unlisted <- wanted[!wanted %in% names(listed)]
  if (length(unlisted)) {
    stop("`dir`'s SHA256SUMS does not list ",
         paste(sprintf("`%s`", unlisted), collapse = ", "),
         ". Rewrite it with gatom_download_refs(dir = \"", dir, "\").",
         call. = FALSE)
  }
  bad <- wanted[tolower(listed[wanted]) != tolower(sha)]
  if (length(bad)) {
    stop("`dir` holds file(s) that do not match its SHA256SUMS: ",
         paste(sprintf("`%s`", bad), collapse = ", "),
         ". Refetch them with gatom_download_refs(overwrite = TRUE).",
         call. = FALSE)
  }
  "verified"
}

#' Print a `gatom_refs` object
#'
#' @param x A `gatom_refs` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @examples
#' \dontrun{
#' print(gatom_refs("Homo sapiens"))
#' }
#' @export
print.gatom_refs <- function(x, ...) {
  cat("<gatom_refs>", x$species, x$network_name %||% "",
      paste0("(", x$checksums %||% "unchecked", ")"), "\n")
  for (nm in names(x$files)) cat(" ", nm, ":", x$files[[nm]], "\n")
  invisible(x)
}

.gatom_col <- function(expr, data, env, arg) {
  val <- tryCatch(eval(expr, data, env), error = function(e) {
    stop("`", arg, "` could not be evaluated in `x`: ", conditionMessage(e),
         call. = FALSE)
  })
  if (length(val) == 1L) val <- rep(val, nrow(data))
  if (length(val) != nrow(data)) {
    stop("`", arg, "` must be length 1 or nrow(x) (", nrow(data), "); got ",
         length(val), ".", call. = FALSE)
  }
  val
}

#' Build a validated GATOM `gene.de` table
#'
#' Assembles the four columns `gatom::makeMetabolicGraph()` needs -- `ID`,
#' `pval`, `log2FC`, `baseMean` -- sorted by p-value and deduplicated on `ID`
#' (lowest p-value wins).
#'
#' **`pval` must be raw.** A column whose *name* looks adjusted (`adj`, `fdr`,
#' `padj`, `qval`/`q.val`) is an error, as is any value outside `[0, 1]`, as
#' is an all-`NA` column. GATOM's BUM model is fitted to a raw p-value
#' distribution; FDR values fit it to nonsense without complaining.
#'
#' **`baseMean` is used only as a rank.** gatom keeps the 12,000 genes with
#' the highest `baseMean` and drops the rest, so any monotone expression
#' measure gives the same graph: `AveExpr`, `2^AveExpr` and DESeq2's
#' `baseMean` are interchangeable. [gatom_graph()] reports how many genes the
#' cut-off removed.
#'
#' Arguments are evaluated inside `x`, so bare column names, expressions
#' (`2^AveExpr`) and constants (`1`) all work.
#'
#' @param x A data frame of differential-expression results.
#' @param id Gene identifier column -- symbols, Entrez, Ensembl or RefSeq,
#'   matching one of the annotation's `mapFrom` types.
#' @param pval **Raw** p-value column.
#' @param log2FC Effect-size column.
#' @param baseMean Mean-expression column on any monotone scale, or a
#'   constant placeholder.
#' @return A data frame of class `gatom_de` with columns `ID`, `pval`,
#'   `log2FC`, `baseMean`.
#' @examples
#' tt <- data.frame(symbol = c("IDO1", "KMO", "KYNU"),
#'                  P.Value = c(1e-6, 2e-4, 0.03),
#'                  logFC = c(2.1, -1.4, 0.8),
#'                  AveExpr = c(8, 7, 6))
#' gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC,
#'          baseMean = AveExpr)
#' @export
gatom_de <- function(x, id, pval, log2FC, baseMean) {
  if (!is.data.frame(x)) {
    stop("`x` must be a data frame of DE results.", call. = FALSE)
  }
  if (!nrow(x)) stop("`x` has zero rows.", call. = FALSE)
  env <- parent.frame()

  pval_expr <- substitute(pval)
  pval_label <- paste(deparse(pval_expr), collapse = "")
  if (grepl("adj|fdr|padj|q\\.?val", pval_label, ignore.case = TRUE)) {
    stop("`pval` looks like an adjusted p-value (`", pval_label, "`). GATOM ",
         "scores with a BUM model fitted to RAW p-values; FDR values fit it ",
         "silently and wrongly. Pass the raw p-value column.", call. = FALSE)
  }

  ID <- as.character(.gatom_col(substitute(id), x, env, "id"))
  p  <- .gatom_col(pval_expr, x, env, "pval")
  fc <- .gatom_col(substitute(log2FC), x, env, "log2FC")
  bm <- .gatom_col(substitute(baseMean), x, env, "baseMean")

  if (!is.numeric(p)) {
    stop("`pval` must be numeric; got ", class(p)[[1L]], ".", call. = FALSE)
  }
  if (all(is.na(p))) {
    stop("`pval` is entirely NA; GATOM cannot score an empty p-value column.",
         call. = FALSE)
  }
  bad <- !is.na(p) & (p < 0 | p > 1)
  if (any(bad)) {
    stop("`pval` must lie in [0, 1]; ", sum(bad), " value(s) do not (range ",
         signif(min(p, na.rm = TRUE), 3), " to ",
         signif(max(p, na.rm = TRUE), 3),
         "). Raw p-values, not scores or -log10 p.", call. = FALSE)
  }
  if (!is.numeric(fc)) {
    stop("`log2FC` must be numeric; got ", class(fc)[[1L]], ".", call. = FALSE)
  }
  if (!is.numeric(bm)) {
    stop("`baseMean` must be numeric; got ", class(bm)[[1L]], ".",
         call. = FALSE)
  }
  if (all(is.na(bm))) {
    stop("`baseMean` is entirely NA. Pass a mean-expression column such as ",
         "`AveExpr` or DESeq2's `baseMean`; a constant such as 1 is an ",
         "acceptable placeholder.", call. = FALSE)
  }

  out <- data.frame(ID = ID, pval = p, log2FC = fc, baseMean = bm,
                    stringsAsFactors = FALSE)
  keep <- !is.na(out$ID) & nzchar(out$ID) & !is.na(out$pval)
  if (!any(keep)) {
    stop("No rows left after dropping missing `id`/`pval`.", call. = FALSE)
  }
  if (!all(keep)) {
    message(sum(!keep), " row(s) dropped: missing `id` or `pval`.")
  }
  out <- out[keep, , drop = FALSE]
  out <- out[order(out$pval), , drop = FALSE]
  out <- out[!duplicated(out$ID), , drop = FALSE]
  rownames(out) <- NULL
  structure(out, class = c("gatom_de", "data.frame"))
}

# ---- layer 1: the graph ------------------------------------------------------

# Attributes each layer adds, carried forward so a module alone describes
# the run that produced it.
.gatom_graph_attrs <- c(
  "species", "network", "topology", "keep_reactions_without_enzymes",
  "id_type", "genes_in_de", "genes_kept", "gene_keep_top", "graph_genes",
  "graph_nodes", "graph_edges", "refs_files"
)
.gatom_score_attrs <- c(
  "k_gene", "k_met", "score_seed", "gene_threshold", "gene_bum_alpha",
  "gene_fdr", "met_threshold", "met_bum_alpha", "met_fdr"
)

.gatom_carry <- function(to, from, names) {
  for (nm in names) attr(to, nm) <- attr(from, nm)
  to
}

.gatom_check_refs <- function(refs) {
  if (!inherits(refs, "gatom_refs")) {
    stop("`refs` must be a `gatom_refs` object from gatom_refs().",
         call. = FALSE)
  }
  invisible(refs)
}

.gatom_check_de <- function(de) {
  if (is.null(de)) return(invisible(NULL))
  if (!is.data.frame(de)) {
    stop("`de` must be a data frame; build it with gatom_de().",
         call. = FALSE)
  }
  needed <- c("ID", "pval", "log2FC", "baseMean")
  if (!all(needed %in% names(de))) {
    stop("`de` is missing column(s): ",
         paste(setdiff(needed, names(de)), collapse = ", "),
         ". Build it with gatom_de().", call. = FALSE)
  }
  if (!nrow(de)) stop("`de` has zero rows.", call. = FALSE)
  invisible(de)
}

.gatom_check_met_de <- function(met_de) {
  if (!is.null(met_de) && (!is.data.frame(met_de) || !nrow(met_de))) {
    stop("`met_de` must be NULL or a non-empty data frame of metabolite ",
         "DE results (an ID column, `pval` and `log2FC`).", call. = FALSE)
  }
  invisible(met_de)
}

.gatom_check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(x)
}

#' Build a GATOM metabolic graph
#'
#' Calls `gatom::makeMetabolicGraph()` on a DE table and a reference set, and
#' reports what happened to the genes on the way: how many the table held,
#' how many survived gatom's `gene.keep.top` cut-off (the 12,000 highest by
#' `baseMean`), and how many reached the graph. The graph's genes are the
#' universe for [gatom_pathways()].
#'
#' The vignette uses `topology = "atoms"` for the KEGG, Rhea and combined
#' networks and recommends `"metabolites"` for the lipid network.
#'
#' gatom detects the id type of a table longer than 1,000 rows from a random
#' sample of 1,000 rows. That sample is drawn under a fixed seed, so the graph
#' depends only on the inputs and the caller's random stream is untouched.
#'
#' @param de A `gatom_de` table from [gatom_de()], or `NULL` when only
#'   metabolite data is given.
#' @param refs A `gatom_refs` object from [gatom_refs()].
#' @param topology `"atoms"` (default) or `"metabolites"`: what a vertex is.
#' @param met_de Optional metabolite DE table, or `NULL`.
#' @param keep_reactions_without_enzymes Logical(1); keep reactions no gene
#'   maps to. gatom's default, `FALSE`.
#' @return The graph as an `igraph`, with attributes `species`, `network`,
#'   `topology`, `keep_reactions_without_enzymes`, `id_type` (as gatom
#'   detected it), `genes_in_de`, `genes_kept`, `gene_keep_top`,
#'   `graph_genes` (Entrez ids on the edges), `graph_nodes`, `graph_edges`
#'   and `refs_files`.
#' @examples
#' \dontrun{
#' g <- gatom_graph(de, refs)
#' attr(g, "genes_kept")
#' }
#' @export
gatom_graph <- function(de, refs, topology = c("atoms", "metabolites"),
                        met_de = NULL, keep_reactions_without_enzymes = FALSE) {
  .require_pkg("gatom", "gatom_graph()", 'BiocManager::install("gatom")')
  .require_pkg("igraph", "gatom_graph()")
  .gatom_check_refs(refs)
  .gatom_check_de(de)
  .gatom_check_met_de(met_de)
  topology <- match.arg(topology)
  .gatom_check_flag(keep_reactions_without_enzymes,
                    "keep_reactions_without_enzymes")
  if (is.null(de) && is.null(met_de)) {
    stop("`de` and `met_de` are both NULL; GATOM needs gene data, ",
         "metabolite data, or both.", call. = FALSE)
  }

  keep_top <- eval(formals(gatom::makeMetabolicGraph)$gene.keep.top)
  id_type <- NA_character_
  genes_in_de <- 0L
  genes_kept <- 0L
  gene_meta <- NULL
  # gatom detects id types from sample(): pinned, so the graph is a pure
  # function of its inputs and the caller's stream is left as it was.
  detect_seed <- 42L
  if (!is.null(de)) {
    de <- as.data.frame(de)
    gene_meta <- .with_pinned_seed(
      detect_seed, gatom::getGeneDEMeta(de, refs$org_anno)
    )
    id_type <- gene_meta$idType
    prepared <- gatom::prepareDE(de, gene_meta)
    genes_in_de <- nrow(de)
    genes_kept <- sum(prepared$signalRank <= keep_top)
  }

  g <- withCallingHandlers(
    .with_pinned_seed(detect_seed, gatom::makeMetabolicGraph(
      network = refs$network,
      topology = topology,
      org.gatom.anno = refs$org_anno,
      gene.de = de,
      gene.de.meta = gene_meta,
      met.db = refs$met_db,          # required even when met.de is NULL
      met.de = met_de,
      gene2reaction.extra = refs$gene2reaction_extra,
      keepReactionsWithoutEnzymes = keep_reactions_without_enzymes
    )),
    # Any other error propagates unchanged.
    error = function(e) {
      if (!grepl("No edges", conditionMessage(e), fixed = TRUE)) return()
      stop("The metabolic graph has no edges: gatom read `de$ID` as ",
           id_type, " ids and none mapped to an enzyme of the ",
           refs$network_name %||% "chosen", " network. Check that `id` in ",
           "gatom_de() holds one of the annotation's id types (",
           paste(names(refs$org_anno$mapFrom), collapse = ", "),
           ") or Entrez ids.", call. = FALSE)
    }
  )

  graph_genes <- unique(as.character(igraph::E(g)$gene))
  graph_genes <- graph_genes[!is.na(graph_genes) & nzchar(graph_genes)]
  if (!is.null(de)) {
    message(sprintf(paste0(
      "gatom_graph(): %d genes in `de`, %d kept by gatom's gene.keep.top ",
      "(%d), %d in the graph"
    ), genes_in_de, genes_kept, keep_top, length(graph_genes)))
  }

  attr(g, "species") <- refs$species
  attr(g, "network") <- refs$network_name
  attr(g, "topology") <- topology
  attr(g, "keep_reactions_without_enzymes") <- keep_reactions_without_enzymes
  attr(g, "id_type") <- id_type
  attr(g, "genes_in_de") <- genes_in_de
  attr(g, "genes_kept") <- genes_kept
  attr(g, "gene_keep_top") <- keep_top
  attr(g, "graph_genes") <- graph_genes
  attr(g, "graph_nodes") <- igraph::vcount(g)
  attr(g, "graph_edges") <- igraph::ecount(g)
  attr(g, "refs_files") <- refs$files
  attr(g, "gatom_layer") <- "graph"
  g
}

# ---- layer 2: the score ------------------------------------------------------

.gatom_check_k <- function(k, arg) {
  if (!is.null(k) &&
        (!is.numeric(k) || length(k) != 1L || is.na(k) || k <= 0)) {
    stop("`", arg, "` must be a single positive number or NULL (50 is the ",
         "GATOM default; larger k gives a larger module).", call. = FALSE)
  }
  invisible(k)
}

.gatom_check_seed <- function(seed) {
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    stop("`seed` must be a single number; the BUM fit and the MWCS solver ",
         "are stochastic and unseeded runs are not reproducible.",
         call. = FALSE)
  }
  invisible(seed)
}

.gatom_check_layer <- function(x, layer, arg, from) {
  if (!inherits(x, "igraph") || !identical(attr(x, "gatom_layer"), layer)) {
    stop("`", arg, "` must be the result of ", from, ".", call. = FALSE)
  }
  invisible(x)
}

# The BUM fit is reported only as console messages; these are their prefixes.
.gatom_score_messages <- c(
  gene_threshold = "Gene p-value threshold: ",
  gene_bum_alpha = "Gene BU alpha: ",
  gene_fdr       = "FDR for genes: ",
  met_threshold  = "Metabolite p-value threshold: ",
  met_bum_alpha  = "Metabolite BU alpha: ",
  met_fdr        = "FDR for metabolites: "
)

#' Score a GATOM graph
#'
#' Calls `gatom::scoreGraph()`, which fits a beta-uniform mixture (BUM) to the
#' gene p-values (and to the metabolite p-values when `k_met` is set) and
#' turns each into a signed weight around a threshold. The threshold is the
#' `k`-th smallest p-value, capped at an FDR of 0.1, so **larger `k` gives a
#' larger module**. 50 is the GATOM default; 25 and 75 are the usual
#' sensitivity branches.
#'
#' The BUM fit uses random starts, so it runs under `seed`. gatom only
#' prints the fit; it is recorded here as attributes. The threshold and BUM
#' alpha are exact, recovered from the scores; the FDR is as gatom prints
#' it, to six decimals, so `0` means below `5e-7`.
#'
#' When the p-values do not fit a BUM (shape `a > 0.5`), gatom sets every
#' score to 0 and warns; a module solved from that graph is noise. This
#' function stops instead.
#'
#' @param g A graph from [gatom_graph()].
#' @param k_gene Numeric(1) gene scoring parameter, or `NULL` to leave genes
#'   unscored (metabolite-only analysis).
#' @param k_met Numeric(1) metabolite scoring parameter, or `NULL`. Needs
#'   `met_de` in [gatom_graph()].
#' @param seed Integer(1) seed for the BUM fit.
#' @return The scored SGMWCS instance as an `igraph`, carrying the graph's
#'   attributes plus `k_gene`, `k_met`, `score_seed`, `gene_threshold`,
#'   `gene_bum_alpha`, `gene_fdr`, `met_threshold`, `met_bum_alpha` and
#'   `met_fdr` (`NA` for an unscored side).
#' @examples
#' \dontrun{
#' gs <- gatom_score(gatom_graph(de, refs), k_gene = 50)
#' attr(gs, "gene_fdr")
#' }
#' @export
gatom_score <- function(g, k_gene = 50, k_met = NULL, seed = 42) {
  .require_pkg("gatom", "gatom_score()", 'BiocManager::install("gatom")')
  .require_pkg("igraph", "gatom_score()")
  .gatom_check_layer(g, "graph", "g", "gatom_graph()")
  .gatom_check_k(k_gene, "k_gene")
  .gatom_check_k(k_met, "k_met")
  .gatom_check_seed(seed)
  if (is.null(k_gene) && is.null(k_met)) {
    stop("`k_gene` and `k_met` are both NULL; nothing would be scored.",
         call. = FALSE)
  }
  if (!is.null(k_gene) && !any(!is.na(igraph::E(g)$pval))) {
    stop("`k_gene` was supplied but the graph carries no gene p-values; ",
         "build it with `de` or leave `k_gene = NULL`.", call. = FALSE)
  }
  if (!is.null(k_met) && !any(!is.na(igraph::V(g)$pval))) {
    stop("`k_met` was supplied but the graph carries no metabolite ",
         "p-values; build it with `met_de` or leave `k_met = NULL`.",
         call. = FALSE)
  }

  fit <- stats::setNames(rep(NA_real_, length(.gatom_score_messages)),
                         names(.gatom_score_messages))
  gs <- .with_pinned_seed(seed, withCallingHandlers(
    gatom::scoreGraph(g, k.gene = k_gene, k.met = k_met),
    message = function(m) {
      txt <- trimws(conditionMessage(m))
      for (nm in names(.gatom_score_messages)) {
        prefix <- .gatom_score_messages[[nm]]
        if (startsWith(txt, prefix)) {
          fit[[nm]] <<- as.numeric(substring(txt, nchar(prefix) + 1L))
        }
      }
    },
    warning = function(w) {
      txt <- conditionMessage(w)
      if (grepl("scores have been assigned to 0", txt, fixed = TRUE)) {
        side <- if (startsWith(txt, "Edge")) "gene" else "metabolite"
        stop("The ", side, " p-values do not fit a beta-uniform mixture ",
             "(BUM shape a > 0.5), so gatom set every ", side, " score to 0 ",
             "and a module would be noise. Check that `",
             if (side == "gene") "de" else "met_de",
             "` holds raw p-values for the whole tested universe, not a ",
             "pre-filtered list.", call. = FALSE)
      }
    }
  ))

  # Guard: the fit is read from gatom's messages, so a gatom release that
  # rewords them would leave these NA without this check.
  scored <- c(if (!is.null(k_gene)) c("gene_threshold", "gene_bum_alpha",
                                      "gene_fdr"),
              if (!is.null(k_met)) c("met_threshold", "met_bum_alpha",
                                     "met_fdr"))
  if (anyNA(fit[scored])) {
    stop("gatom::scoreGraph() did not report ",
         paste(sprintf("`%s`", scored[is.na(fit[scored])]), collapse = ", "),
         "; its messages no longer match `.gatom_score_messages`.",
         call. = FALSE)
  }
  # gatom prints six decimals, so a threshold of 6e-7 prints as 0. The exact
  # threshold and alpha are recovered from the scores themselves and must
  # agree with the printed values.
  if (!is.null(k_gene)) {
    fit[c("gene_bum_alpha", "gene_threshold")] <- .gatom_exact_fit(
      igraph::E(gs)$pval, igraph::E(gs)$score,
      fit[c("gene_bum_alpha", "gene_threshold")], "gene"
    )
  }
  if (!is.null(k_met)) {
    fit[c("met_bum_alpha", "met_threshold")] <- .gatom_exact_fit(
      igraph::V(gs)$pval, igraph::V(gs)$score,
      fit[c("met_bum_alpha", "met_threshold")], "metabolite"
    )
  }

  gs <- .gatom_carry(gs, g, .gatom_graph_attrs)
  attr(gs, "k_gene") <- k_gene
  attr(gs, "k_met") <- k_met
  attr(gs, "score_seed") <- seed
  for (nm in names(fit)) attr(gs, nm) <- fit[[nm]]
  attr(gs, "gatom_layer") <- "scored"
  gs
}

#' Recover the exact BUM alpha and threshold from gatom's scores
#'
#' `scoreGraph()` scores each p-value as `(a - 1) * (log(p) - log(t))`, so
#' the scores are linear in `log(p)` with slope `a - 1`. Two distinct
#' p-values determine `a` and `t` exactly; the rest must lie on the line.
#'
#' @param pval,score Numeric vectors from the scored graph's edges or
#'   vertices.
#' @param printed Numeric(2): the alpha and threshold gatom printed.
#' @param side `"gene"` or `"metabolite"`, for the error message.
#' @return Numeric(2): the exact alpha and threshold, or `printed` when
#'   fewer than two distinct p-values were scored.
#' @keywords internal
.gatom_exact_fit <- function(pval, score, printed, side) {
  ok <- !is.na(pval) & pval > 0 & !is.na(score)
  pts <- unique(data.frame(lp = log(pval[ok]), s = score[ok]))
  if (length(unique(pts$lp)) < 2L) return(printed)
  ends <- pts[c(which.min(pts$lp), which.max(pts$lp)), ]
  slope <- diff(ends$s) / diff(ends$lp)
  log_t <- ends$lp[[1L]] - ends$s[[1L]] / slope
  off_line <- max(abs(pts$s - slope * (pts$lp - log_t)))
  exact <- c(slope + 1, exp(log_t))
  # Printed with %f: agreement to the sixth decimal, plus rounding.
  if (off_line > 1e-6 * max(1, abs(pts$s)) ||
        any(abs(exact - printed) > 1e-6)) {
    stop("gatom::scoreGraph()'s ", side, " scores no longer follow ",
         "(a - 1) * (log(p) - log(threshold)); `.gatom_exact_fit()` needs ",
         "updating for this gatom version.", call. = FALSE)
  }
  exact
}

# ---- layer 3: the solve ------------------------------------------------------

#' Build a GATOM MWCS solver
#'
#' The two solvers the GATOM vignette uses, with the vignette's parameters:
#'
#' - `"rnc"`: `mwcsr::rnc_solver()`, the heuristic relax-and-cut solver the
#'   vignette uses "for simplicity". It needs nothing beyond mwcsr.
#' - `"virgo"`: `mwcsr::virgo_solver(cplex_dir, threads = 4,
#'   penalty = 0.001, log = 1)`, the exact solver the vignette recommends
#'   "for proper analysis quality". It needs Java and IBM CPLEX (12.7 or
#'   later). The edge penalty removes redundant genes from the module.
#'
#' CPLEX is licensed by IBM and is not redistributable, so it is never part
#' of an image: bind an installation read-only into the container and point
#' `CPLEX_HOME` at it (by convention `/opt/cplex`). Without `cplex_dir`,
#' `virgo_solver()` falls back to an approximate mode; this function stops
#' instead.
#'
#' @param type `"rnc"` or `"virgo"`.
#' @param cplex_dir Character(1) CPLEX installation directory, holding
#'   `cplex.jar` and the CPLEX libraries. Virgo only.
#' @param threads Integer(1) CPLEX threads. Virgo only.
#' @param penalty Numeric(1) per-edge penalty. Virgo only.
#' @param timelimit Numeric(1) seconds, or `NULL` for none. Virgo only.
#' @return An `mwcs_solver` object for [gatom_solve()], with attribute
#'   `gatom_solver`: the type and its parameters.
#' @examples
#' if (requireNamespace("mwcsr", quietly = TRUE)) {
#'   solver <- gatom_solver("rnc")
#' }
#' \dontrun{
#' solver <- gatom_solver("virgo", cplex_dir = "/opt/cplex")
#' }
#' @export
gatom_solver <- function(type = c("rnc", "virgo"),
                         cplex_dir = Sys.getenv("CPLEX_HOME"),
                         threads = 4, penalty = 0.001, timelimit = NULL) {
  .require_pkg("mwcsr", "gatom_solver()")
  type <- match.arg(type)

  if (type == "rnc") {
    given <- c(cplex_dir = !missing(cplex_dir), threads = !missing(threads),
               penalty = !missing(penalty), timelimit = !missing(timelimit))
    if (any(given)) {
      stop(paste(sprintf("`%s`", names(given)[given]), collapse = ", "),
           if (sum(given) == 1L) " is a virgo setting" else
             " are virgo settings",
           "; drop it for `type = \"rnc\"`.", call. = FALSE)
    }
    solver <- mwcsr::rnc_solver()
    attr(solver, "gatom_solver") <- list(type = "rnc")
    return(solver)
  }

  if (!is.character(cplex_dir) || length(cplex_dir) != 1L ||
        is.na(cplex_dir) || !nzchar(cplex_dir)) {
    stop("`cplex_dir` is empty: set `CPLEX_HOME` to a CPLEX installation ",
         "(in the container, bind it read-only at /opt/cplex) or pass ",
         "`cplex_dir`. The virgo solver needs CPLEX; use type = \"rnc\" ",
         "without it.", call. = FALSE)
  }
  if (!dir.exists(cplex_dir)) {
    stop("`cplex_dir` does not exist: ", cplex_dir, ". Check the CPLEX ",
         "bind mount and `CPLEX_HOME`.", call. = FALSE)
  }
  if (!is.numeric(threads) || length(threads) != 1L || is.na(threads) ||
        threads < 1) {
    stop("`threads` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(penalty) || length(penalty) != 1L || is.na(penalty) ||
        penalty < 0) {
    stop("`penalty` must be a single non-negative number.", call. = FALSE)
  }
  if (!is.null(timelimit) &&
        (!is.numeric(timelimit) || length(timelimit) != 1L ||
         is.na(timelimit) || timelimit <= 0)) {
    stop("`timelimit` must be NULL or a single positive number of seconds.",
         call. = FALSE)
  }
  solver <- tryCatch(
    mwcsr::virgo_solver(cplex_dir = cplex_dir, threads = threads,
                        penalty = penalty, log = 1, timelimit = timelimit),
    error = function(e) {
      stop("`cplex_dir` (", cplex_dir, ") is not a usable CPLEX ",
           "installation: ", conditionMessage(e), call. = FALSE)
    }
  )
  attr(solver, "gatom_solver") <- list(
    type = "virgo", cplex_dir = cplex_dir, threads = threads,
    penalty = penalty, timelimit = timelimit
  )
  solver
}

.gatom_resolve_solver <- function(solver) {
  if (inherits(solver, "mwcs_solver")) return(solver)
  if (!is.character(solver) || length(solver) != 1L || is.na(solver)) {
    stop("`solver` must be \"rnc\", \"virgo\" or an `mwcs_solver` from ",
         "gatom_solver().", call. = FALSE)
  }
  if (solver %in% c("rnc", "virgo")) return(gatom_solver(solver))
  if (solver %in% c("rmwcs", "annealing")) {
    warning("`solver = \"", solver, "\"` is deprecated and will be removed ",
            "in bulkiRNA 2.0.0: the GATOM authors use \"rnc\" and recommend ",
            "\"virgo\".", call. = FALSE)
    out <- switch(solver,
                  rmwcs = mwcsr::rmwcs_solver(),
                  annealing = mwcsr::annealing_solver())
    attr(out, "gatom_solver") <- list(type = solver)
    return(out)
  }
  stop("`solver` must be \"rnc\", \"virgo\" or an `mwcs_solver` from ",
       "gatom_solver(); got \"", solver, "\".", call. = FALSE)
}

.gatom_solution_weight <- function(solution) {
  weight <- solution$weight
  if (is.null(weight)) NA_real_ else as.numeric(weight)
}

#' Solve a scored GATOM graph for its module
#'
#' Calls `mwcsr::solve_mwcsp()` under `seed`, as the vignette does with
#' `set.seed(42)`, and returns the module with everything that produced it
#' recorded on it.
#'
#' `seed` fixes R-side randomness, which is all `"rnc"` uses. Virgo runs in
#' Java and CPLEX, outside R's generator: it always returns an optimum, but
#' where several subgraphs tie it can return a different one on each call.
#' On a KathleenM module, five solves gave the same weight, the same 50 genes
#' and the same 56 reactions, with 50 or 51 metabolites as the route through
#' unscored metabolites changed. Compare virgo modules by
#' `solution_weight` and genes, not by their exact node and edge sets.
#'
#' @param gs A scored graph from [gatom_score()].
#' @param solver `"rnc"`, `"virgo"`, or an `mwcs_solver` object -- from
#'   [gatom_solver()] or built directly with mwcsr.
#' @param seed Integer(1) seed for the solver.
#' @return The module as an `igraph`, carrying the graph and score attributes
#'   plus `solver` (its type), `solver_params`, `seed`, `solution_weight`,
#'   `solved_to_optimality` (`TRUE` only when the solver proved the optimum),
#'   `n_nodes` and `n_edges`.
#' @examples
#' \dontrun{
#' m <- gatom_solve(gatom_score(gatom_graph(de, refs)), solver = "rnc")
#' attr(m, "solved_to_optimality")
#' }
#' @export
gatom_solve <- function(gs, solver = "rnc", seed = 42) {
  .require_pkg("mwcsr", "gatom_solve()")
  .require_pkg("igraph", "gatom_solve()")
  .gatom_check_layer(gs, "scored", "gs", "gatom_score()")
  .gatom_check_seed(seed)
  solver_obj <- .gatom_resolve_solver(solver)
  params <- attr(solver_obj, "gatom_solver") %||%
    list(type = class(solver_obj)[[1L]])

  solution <- .with_pinned_seed(seed, mwcsr::solve_mwcsp(solver_obj, gs))
  m <- solution$graph

  m <- .gatom_carry(m, gs, c(.gatom_graph_attrs, .gatom_score_attrs))
  attr(m, "solver") <- params$type
  attr(m, "solver_params") <- params
  attr(m, "seed") <- seed
  # This name matches the downstream results-table column.
  attr(m, "solution_weight") <- .gatom_solution_weight(solution)
  attr(m, "solved_to_optimality") <- isTRUE(solution$solved_to_optimality)
  attr(m, "n_nodes") <- igraph::vcount(m)
  attr(m, "n_edges") <- igraph::ecount(m)
  attr(m, "gatom_layer") <- "module"
  m
}

# ---- the one-call recipe -----------------------------------------------------

#' Build, score and solve a GATOM module
#'
#' [gatom_graph()], [gatom_score()] and [gatom_solve()] in order, with one
#' seed for the BUM fit and the solver. Use the layers directly to reuse a
#' graph across `k_gene` values or to compare solvers on one scored graph.
#'
#' `k_gene` is the module-size dial: **larger `k` gives a larger module.**
#' 50 is the GATOM default; 25 and 75 are the standard sensitivity branches.
#'
#' @param de A `gatom_de` table from [gatom_de()], or `NULL` with `met_de`.
#' @param refs A `gatom_refs` object from [gatom_refs()].
#' @param k_gene Numeric(1) gene scoring parameter, or `NULL`; see
#'   [gatom_score()].
#' @param k_met Numeric(1) metabolite scoring parameter, or `NULL` when
#'   `met_de` is `NULL`.
#' @param met_de Optional metabolite DE table.
#' @param seed Integer(1) seed, applied separately to the BUM fit and the
#'   solver.
#' @param solver `"rnc"` (default), `"virgo"`, or an `mwcs_solver` from
#'   [gatom_solver()]. `"rmwcs"` and `"annealing"` are deprecated.
#' @param verbose Logical(1); also report the graph and module sizes.
#' @param gene2reaction_extra Deprecated: [gatom_refs()] now loads the
#'   network's `gene2reaction` table. A data frame given here replaces it.
#' @param topology `"atoms"` (default) or `"metabolites"`; see
#'   [gatom_graph()].
#' @param keep_reactions_without_enzymes Logical(1); see [gatom_graph()].
#' @return The module from [gatom_solve()].
#' @examples
#' \dontrun{
#' refs <- gatom_refs("Mus musculus", network = "combined")
#' de <- gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC,
#'                baseMean = AveExpr)
#' m <- gatom_module(de, refs, k_gene = 50, seed = 42)
#' }
#' @export
gatom_module <- function(de, refs, k_gene = 50, k_met = NULL, met_de = NULL,
                         seed = 42, solver = "rnc", verbose = FALSE,
                         gene2reaction_extra = NULL,
                         topology = c("atoms", "metabolites"),
                         keep_reactions_without_enzymes = FALSE) {
  .require_pkg("gatom", "gatom_module()", 'BiocManager::install("gatom")')
  .require_pkg("mwcsr", "gatom_module()")
  .require_pkg("igraph", "gatom_module()")

  # Everything is checked before the graph is built, which takes seconds.
  .gatom_check_refs(refs)
  if (is.null(de) && is.null(met_de)) {
    stop("`de` and `met_de` are both NULL; GATOM needs gene data, ",
         "metabolite data, or both.", call. = FALSE)
  }
  .gatom_check_de(de)
  .gatom_check_met_de(met_de)
  .gatom_check_k(k_gene, "k_gene")
  .gatom_check_k(k_met, "k_met")
  if (is.null(met_de) && !is.null(k_met)) {
    stop("`k_met` was supplied but `met_de` is NULL; metabolite scoring ",
         "needs metabolite data. Pass `met_de` or leave `k_met = NULL`.",
         call. = FALSE)
  }
  .gatom_check_seed(seed)
  .gatom_check_flag(verbose, "verbose")
  topology <- match.arg(topology)
  .gatom_check_flag(keep_reactions_without_enzymes,
                    "keep_reactions_without_enzymes")
  if (!is.null(gene2reaction_extra)) {
    if (!is.data.frame(gene2reaction_extra) ||
          !all(c("gene", "reaction") %in% names(gene2reaction_extra))) {
      stop("`gene2reaction_extra` must be NULL or a data frame with `gene` ",
           "and `reaction` columns.", call. = FALSE)
    }
    warning("`gene2reaction_extra` is deprecated and will be removed in ",
            "bulkiRNA 2.0.0: gatom_refs() now loads the network's ",
            "gene2reaction table. The table given here replaces it.",
            call. = FALSE)
    refs$gene2reaction_extra <- gene2reaction_extra
  }
  solver_obj <- .gatom_resolve_solver(solver)

  g <- gatom_graph(de, refs, topology = topology, met_de = met_de,
                   keep_reactions_without_enzymes =
                     keep_reactions_without_enzymes)
  if (verbose) {
    message("Graph: ", igraph::vcount(g), " nodes, ", igraph::ecount(g),
            " edges")
  }
  gs <- gatom_score(g, k_gene = k_gene, k_met = k_met, seed = seed)
  m <- gatom_solve(gs, solver = solver_obj, seed = seed)
  if (verbose) {
    message("Module: ", igraph::vcount(m), " nodes, ", igraph::ecount(m),
            " edges")
  }
  m
}

# ---- readers -----------------------------------------------------------------

#' Genes in a GATOM module
#'
#' The genes live on the graph's **edges** -- enzyme-catalysed reactions --
#' in both topologies, so `igraph::V(m)$Symbol` returns `NULL` and looks like
#' an empty result. This reads the edge `label`, which gatom fills with the
#' annotation's gene symbol whatever id type the DE table used.
#'
#' @param m A module from [gatom_module()] or [gatom_solve()].
#' @return A character vector of unique gene symbols, in edge order.
#' @examples
#' \dontrun{
#' gatom_genes(m)
#' }
#' @export
gatom_genes <- function(m) {
  .require_pkg("igraph", "gatom_genes()")
  if (!inherits(m, "igraph")) {
    stop("`m` must be an igraph module from gatom_module(); got ",
         class(m)[[1L]], ".", call. = FALSE)
  }
  edges <- igraph::as_data_frame(m, "edges")
  if (!"label" %in% names(edges)) {
    stop("The edges of `m` carry no `label` column. GATOM modules label ",
         "genes on edges; this graph was either not built by gatom or lost ",
         "its edge attributes.", call. = FALSE)
  }
  sym <- as.character(edges$label)
  unique(sym[!is.na(sym) & nzchar(sym)])
}

#' Pathway over-representation in a GATOM module
#'
#' The vignette's annotation step: `fgsea::fora()` tests the module's genes
#' against the KEGG and Reactome metabolic pathways carried in the species
#' annotation, with the graph's genes as the universe. With `collapse = TRUE`,
#' `fgsea::collapsePathwaysORA()` then marks the non-redundant pathways among
#' those with `padj < 0.05`, the vignette's threshold. Every tested pathway
#' is returned.
#'
#' @param m A module from [gatom_module()] or [gatom_solve()].
#' @param refs The `gatom_refs` object the module was built from.
#' @param universe Character vector of Entrez ids: the genes that could have
#'   entered the module. Defaults to the graph's genes, recorded on `m`.
#' @param min_size Integer(1) smallest pathway tested, after restricting to
#'   `universe`. The vignette's value, 5.
#' @param collapse Logical(1); mark the main pathways among the significant.
#' @return A tibble with one row per tested pathway: `fgsea::fora()`'s
#'   columns (`pathway`, `pval`, `padj`, `foldEnrichment`, `overlap`, `size`,
#'   `overlapGenes`, a list column of Entrez ids)
#'   and `main` (`TRUE` for a main pathway; `NA` throughout when
#'   `collapse = FALSE`). Attributes `n_module_genes`, `n_universe` and
#'   `collapse_padj` record the test.
#' @examples
#' \dontrun{
#' gatom_pathways(m, refs)
#' }
#' @export
gatom_pathways <- function(m, refs, universe = attr(m, "graph_genes"),
                           min_size = 5, collapse = TRUE) {
  .require_pkg("igraph", "gatom_pathways()")
  if (!inherits(m, "igraph")) {
    stop("`m` must be an igraph module from gatom_module(); got ",
         class(m)[[1L]], ".", call. = FALSE)
  }
  .gatom_check_refs(refs)
  pathways <- refs$org_anno$pathways
  if (!is.list(pathways) || !length(pathways)) {
    stop("`refs` carries no pathways in its annotation ",
         "(`refs$org_anno$pathways`).", call. = FALSE)
  }
  if (!is.character(universe) || !length(universe)) {
    stop("`universe` must be a non-empty character vector of Entrez ids; ",
         "a module from gatom_module() records it as attr(m, ",
         "\"graph_genes\").", call. = FALSE)
  }
  if (!is.numeric(min_size) || length(min_size) != 1L || is.na(min_size) ||
        min_size < 1) {
    stop("`min_size` must be a single positive number.", call. = FALSE)
  }
  .gatom_check_flag(collapse, "collapse")

  genes <- unique(as.character(igraph::E(m)$gene))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) {
    stop("`m` has no genes on its edges; there is nothing to annotate.",
         call. = FALSE)
  }

  collapse_padj <- 0.05
  res <- fgsea::fora(pathways = pathways, genes = genes, universe = universe,
                     minSize = min_size)
  main <- rep(NA, nrow(res))
  if (collapse) {
    main <- rep(FALSE, nrow(res))
    sig <- res[res$padj < collapse_padj, ]
    if (nrow(sig)) {
      collapsed <- fgsea::collapsePathwaysORA(sig, pathways = pathways,
                                              genes = genes,
                                              universe = universe)
      main <- res$pathway %in% collapsed$mainPathways
    }
  }
  out <- tibble::as_tibble(as.data.frame(res))
  out$main <- main
  attr(out, "n_module_genes") <- length(genes)
  attr(out, "n_universe") <- length(unique(universe))
  attr(out, "collapse_padj") <- collapse_padj
  out
}

#' Save a GATOM module as a self-contained HTML view
#'
#' Wraps `gatom::saveModuleToHtml()`. `htmlwidgets` needs pandoc to build a
#' self-contained file, and Quarto's bundled pandoc is not exported to child
#' R sessions -- so when pandoc is not on the PATH and `RSTUDIO_PANDOC` is
#' unset, this points `RSTUDIO_PANDOC` at a known location for the duration
#' of the call, and fails with an actionable message when none is found. The
#' parent directory is created if missing.
#'
#' @param m A module `igraph` from [gatom_module()].
#' @param path Character(1) output `.html` path.
#' @param name Character(1) title shown in the view.
#' @return `path`, invisibly.
#' @examples
#' \dontrun{
#' gatom_save_html(m, "03_results/kyn_module.html", name = "Kynurenine")
#' }
#' @export
gatom_save_html <- function(m, path, name = "") {
  .require_pkg("gatom", "gatom_save_html()", 'BiocManager::install("gatom")')
  .require_pkg("igraph", "gatom_save_html()")
  .gatom_check_save(m, path, name)

  if (!nzchar(Sys.getenv("RSTUDIO_PANDOC")) &&
        !nzchar(Sys.which("pandoc"))) {
    candidates <- c("/opt/quarto/bin/tools/x86_64",
                    "/opt/quarto/bin/tools",
                    "/usr/lib/rstudio-server/bin/quarto/bin/tools")
    hit <- candidates[dir.exists(candidates) &
                        file.exists(file.path(candidates, "pandoc"))]
    if (!length(hit)) {
      stop("saveModuleToHtml() needs pandoc, which is not on the PATH and ",
           "`RSTUDIO_PANDOC` is unset. Install pandoc, or point at Quarto's ",
           "copy: Sys.setenv(RSTUDIO_PANDOC = \"/opt/quarto/bin/tools/",
           "x86_64\").", call. = FALSE)
    }
    # Scoped to this call: the session's environment is left as it was.
    Sys.setenv(RSTUDIO_PANDOC = hit[[1L]])
    on.exit(Sys.unsetenv("RSTUDIO_PANDOC"), add = TRUE)
  }

  gatom::saveModuleToHtml(m, path, name = name)
  invisible(path)
}

.gatom_check_save <- function(m, path, name) {
  if (!inherits(m, "igraph")) {
    stop("`m` must be an igraph module from gatom_module(); got ",
         class(m)[[1L]], ".", call. = FALSE)
  }
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be a single non-empty file path.", call. = FALSE)
  }
  if (!is.character(name) || length(name) != 1L) {
    stop("`name` must be a single string.", call. = FALSE)
  }
  parent <- dirname(path)
  if (nzchar(parent) && !identical(parent, ".")) ensure_dir(parent)
  invisible(NULL)
}

#' Save a GATOM module as a PDF with a repelled layout
#'
#' Wraps `gatom::saveModuleToPdf()` with the vignette's call: `n_iter = 100`,
#' `force = 1e-5`, after `set.seed(42)`. The label layout is a stochastic
#' repel, so `seed` fixes it. The vignette's advice: the larger the module,
#' the softer the `force`. The parent directory is created if missing.
#'
#' gatom's layout can place every node on one line, and the drawing then
#' fails. This is common in modules of a few nodes and depends on the seed.
#' The function then stops, closes the devices gatom left open and removes
#' the partial file, so a failed call leaves no state behind.
#'
#' @param m A module `igraph` from [gatom_module()].
#' @param path Character(1) output `.pdf` path.
#' @param name Character(1) title printed on the page.
#' @param n_iter Integer(1) iterations of the label-repel layout.
#' @param force Numeric(1) repel force.
#' @param seed Integer(1) seed for the layout.
#' @return `path`, invisibly.
#' @examples
#' \dontrun{
#' gatom_save_pdf(m, "03_results/module.pdf", name = "M0.vs.M1")
#' }
#' @export
gatom_save_pdf <- function(m, path, name = "", n_iter = 100, force = 1e-5,
                           seed = 42) {
  .require_pkg("gatom", "gatom_save_pdf()", 'BiocManager::install("gatom")')
  .require_pkg("igraph", "gatom_save_pdf()")
  if (!is.numeric(n_iter) || length(n_iter) != 1L || is.na(n_iter) ||
        n_iter < 1) {
    stop("`n_iter` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(force) || length(force) != 1L || is.na(force) ||
        force <= 0) {
    stop("`force` must be a single positive number.", call. = FALSE)
  }
  .gatom_check_seed(seed)
  .gatom_check_save(m, path, name)
  # saveModuleToPdf() opens its devices before it can fail, and leaves them
  # open when it does. Close whatever it opened and drop the partial file.
  devices <- grDevices::dev.list()
  on.exit({
    for (d in setdiff(grDevices::dev.list(), devices)) grDevices::dev.off(d)
  }, add = TRUE)
  tryCatch(
    .with_pinned_seed(seed, gatom::saveModuleToPdf(
      m, file = path, name = name, n_iter = n_iter, force = force
    )),
    error = function(e) {
      for (d in setdiff(grDevices::dev.list(), devices)) grDevices::dev.off(d)
      unlink(path)
      stop("gatom::saveModuleToPdf() could not draw `m`: ",
           conditionMessage(e), ". Its label layout can put every node on ",
           "one line, most often in small modules; another `seed`, `force` ",
           "or `n_iter` may avoid it, or use gatom_save_html().",
           call. = FALSE)
    }
  )
  invisible(path)
}
