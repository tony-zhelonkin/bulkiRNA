#' GATOM metabolic-network modules
#'
#' @description
#' A five-function wrapper around the GATOM workflow (build an atom-resolved
#' metabolic graph, score it, solve for the maximum-weight connected subgraph).
#' The wrapper exists to *enforce* the four traps that the raw API only
#' documents:
#'
#' 1. `pval` must be **raw** p-values -- BUM scoring degrades silently on FDR
#'    values. [gatom_de()] refuses an adjusted-p column.
#' 2. `baseMean` must be on the **linear** scale (`2^AveExpr` from limma).
#'    [gatom_de()] warns on the log-scale smell.
#' 3. With `topology = "atoms"` the genes live on the graph **edges**, not the
#'    vertices. [gatom_genes()] is the only correct extractor.
#' 4. `met.db` is required even when `met.de` is `NULL`. [gatom_refs()] carries
#'    it, so [gatom_module()] cannot omit it.
#'
#' @name gatom
#' @keywords internal
NULL

#' Locate the GATOM reference files
#'
#' GATOM needs three files: an atom-transition network, a metabolite database
#' and a species-specific organism annotation. They are too large to bundle,
#' so they are resolved from disk. Search order is an explicit `dir` first,
#' then the staged container path `/opt/gatom-refs`, then
#' [gatom_download_refs()]'s default `dir`
#' (`"00_data/references/gatom"`) -- note that the staged copy therefore wins
#' over that default. With `download = TRUE` the download destination is
#' searched first, so a fresh fetch is never shadowed by the staged copy.
#'
#' The returned object carries `met.db` alongside the network. That is
#' deliberate: `gatom::makeMetabolicGraph()` needs `met.db` even when
#' `met.de = NULL`, and omitting it breaks topology loading silently.
#'
#' @param species Character(1) human or mouse species alias.
#' @param dir Character(1) directory to search first, or `NULL`.
#' @param download Logical(1); if `TRUE`, missing files are fetched with
#'   [gatom_download_refs()] into `dir` (default
#'   `"00_data/references/gatom"`) before resolving.
#' @param network Character(1) network flavour, either `"kegg"` or
#'   `"combined"`.
#' @return An object of class `gatom_refs`: a list with `network`, `met_db`,
#'   `org_anno`, `species` and `files`.
#' @examples
#' \dontrun{
#' refs <- gatom_refs("Homo sapiens")
#' refs <- gatom_refs("Mus musculus", dir = "00_data/references/gatom",
#'                    download = TRUE)
#' }
#' @export
gatom_refs <- function(species = "Homo sapiens", dir = NULL,
                       download = FALSE, network = "kegg") {
  sp <- .species(species)
  if (!is.character(network) || length(network) != 1L || is.na(network) ||
        !network %in% c("kegg", "combined")) {
    stop("`network` must be one of \"kegg\", \"combined\"; got ",
         paste0("\"", paste(network, collapse = "\", \""), "\"."),
         call. = FALSE)
  }
  if (!is.null(dir) && (!is.character(dir) || length(dir) != 1L)) {
    stop("`dir` must be a single directory path or NULL.", call. = FALSE)
  }
  if (!is.logical(download) || length(download) != 1L || is.na(download)) {
    stop("`download` must be TRUE or FALSE.", call. = FALSE)
  }

  default_dir <- "00_data/references/gatom"
  wanted <- c(
    network  = sprintf("network.%s.rds", network),
    met_db   = sprintf("met.%s.db.rds", network),
    org_anno = sprintf("org.%s.eg.gatom.anno.rds", sp$gatom_short)
  )

  dl_dir <- dir %||% default_dir
  if (isTRUE(download)) {
    gatom_download_refs(
      dir = dl_dir,
      species = sp$gatom_download,
      networks = network
    )
  }

  search_dirs <- .gatom_search_dirs(dir, download, default_dir)
  found <- vapply(wanted, function(fname) {
    hit <- file.path(search_dirs, fname)
    hit <- hit[file.exists(hit)]
    if (length(hit)) hit[[1L]] else NA_character_
  }, character(1))

  if (anyNA(found)) {
    missing <- wanted[is.na(found)]
    stop(
      "GATOM reference file(s) not found: ",
      paste(sprintf("`%s`", missing), collapse = ", "), ".\n",
      "Searched: ", paste(search_dirs, collapse = ", "), ".\n",
      "Fetch them with gatom_download_refs(species = \"",
      sp$gatom_download, "\", networks = \"", network, "\", dir = \"",
      dir %||% default_dir,
      "\"), then pass that directory as `dir`.",
      call. = FALSE
    )
  }

  structure(
    list(
      network  = readRDS(found[["network"]]),
      met_db   = readRDS(found[["met_db"]]),
      org_anno = readRDS(found[["org_anno"]]),
      species  = sp$scientific,
      files    = found
    ),
    class = "gatom_refs"
  )
}

#' Directories to resolve GATOM reference files from, in order
#'
#' The staged container copy is searched before `default_dir`, so with
#' `dir = NULL` a `download = TRUE` run fetched fresh references into
#' `00_data/references/gatom` and then loaded the older `/opt/gatom-refs` ones:
#' the download had no visible effect, and `print(refs)` was the only clue.
#' What was just downloaded therefore goes first.
#'
#' @param dir An explicit directory, or `NULL`.
#' @param download Logical; whether this call fetched the references.
#' @param default_dir The downloader's default destination.
#' @return Character vector of directories, most-preferred first.
#' @keywords internal
.gatom_search_dirs <- function(dir, download, default_dir) {
  unique(c(
    if (isTRUE(download)) dir %||% default_dir,
    dir, "/opt/gatom-refs", default_dir
  ))
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
  cat("<gatom_refs>", x$species, "\n")
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
#' (lowest p-value wins), and validates the two inputs that fail silently:
#'
#' - **`pval` must be raw.** A column whose *name* looks adjusted
#'   (`adj`, `fdr`, `padj`, `qval`/`q.val`) is an error, as is any value
#'   outside `[0, 1]`, as is an all-`NA` column. GATOM's BUM model is fitted
#'   to a raw p-value distribution; FDR values fit it to nonsense without
#'   complaining.
#' - **`baseMean` must be linear.** Negative values, or a maximum below 30,
#'   are the signature of a log-scale column (use `2^AveExpr` from limma) and
#'   raise a warning.
#'
#' Arguments are evaluated inside `x`, so bare column names, expressions
#' (`2^AveExpr`) and constants (`1`) all work.
#'
#' @param x A data frame of differential-expression results.
#' @param id Gene identifier column -- symbols, Entrez or RefSeq, matching the
#'   annotation's `mapFrom`.
#' @param pval **Raw** p-value column.
#' @param log2FC Effect-size column.
#' @param baseMean Linear-scale mean expression (e.g. `2^AveExpr`), or a
#'   constant placeholder.
#' @return A data frame of class `gatom_de` with columns `ID`, `pval`,
#'   `log2FC`, `baseMean`.
#' @examples
#' tt <- data.frame(symbol = c("IDO1", "KMO", "KYNU"),
#'                  P.Value = c(1e-6, 2e-4, 0.03),
#'                  logFC = c(2.1, -1.4, 0.8),
#'                  AveExpr = c(8, 7, 6))
#' gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC,
#'          baseMean = 2^AveExpr)
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
    stop("`baseMean` is entirely NA. Use `2^AveExpr` (limma) or DESeq2's ",
         "`baseMean`; a constant such as 1 is an acceptable placeholder.",
         call. = FALSE)
  }
  if (any(!is.na(bm) & bm < 0)) {
    warning("`baseMean` has negative values -- it looks log-scaled. GATOM ",
            "expects LINEAR expression; use `2^AveExpr`.", call. = FALSE)
  } else if (max(bm, na.rm = TRUE) < 30) {
    warning("`baseMean` maxes out at ", signif(max(bm, na.rm = TRUE), 3),
            " -- that is the log-scale smell. GATOM expects LINEAR ",
            "expression; use `2^AveExpr`.", call. = FALSE)
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

.gatom_solver <- function(solver) {
  if (!is.character(solver) || length(solver) != 1L) {
    stop("`solver` must be one of \"rnc\", \"rmwcs\", \"annealing\".",
         call. = FALSE)
  }
  switch(
    solver,
    rnc       = mwcsr::rnc_solver(),
    rmwcs     = mwcsr::rmwcs_solver(),
    annealing = mwcsr::annealing_solver(),
    stop("`solver` must be one of \"rnc\", \"rmwcs\", \"annealing\"; got \"",
         solver, "\".", call. = FALSE)
  )
}

.gatom_solution_weight <- function(solution) {
  weight <- solution$weight
  if (is.null(weight)) NA_real_ else as.numeric(weight)
}

#' Build, score and solve a GATOM module
#'
#' Runs the three-call GATOM pipeline -- `makeMetabolicGraph()` (with
#' `topology = "atoms"` and `keepReactionsWithoutEnzymes = FALSE`),
#' `scoreGraph()`, then `mwcsr::solve_mwcsp()` -- and returns the module
#' subgraph.
#'
#' `k_gene` is the module-size dial: **smaller `k` gives a larger module.**
#' 50 is the GATOM default; 25 and 75 are the standard sensitivity branches.
#' Scoring and solving are stochastic, so `seed` is reset independently before
#' each call and recorded on the result.
#'
#' @param de A `gatom_de` table from [gatom_de()] (or a data frame with the
#'   same four columns).
#' @param refs A `gatom_refs` object from [gatom_refs()].
#' @param k_gene Numeric(1) gene-score parameter; smaller means larger module.
#' @param k_met Numeric(1) metabolite-score parameter, or `NULL` when
#'   `met_de` is `NULL`.
#' @param met_de Optional metabolite DE table.
#' @param seed Integer(1) RNG seed reset separately before scoring and solving.
#'   For callers using R's default generator, each stochastic call starts from
#'   the same state as before this helper was introduced, so results do not
#'   move. Results do move to the reproducible value when callers had changed
#'   `RNGkind()`.
#' @param solver Character(1): `"rnc"` (default), `"rmwcs"` or `"annealing"`.
#' @param verbose Logical(1); report graph and module sizes.
#' @param gene2reaction_extra Optional data frame of additional gene-to-reaction
#'   mappings with columns `gene` and `reaction`. Passed to
#'   `gatom::makeMetabolicGraph()` as `gene2reaction.extra`.
#' @return The module as an `igraph`, with attributes `solution_weight`,
#'   `k_gene`, `k_met`, `seed`, `solver`, `species`, `graph_nodes`,
#'   `graph_edges`, `n_nodes` and `n_edges`.
#' @examples
#' \dontrun{
#' refs <- gatom_refs("Homo sapiens")
#' de <- gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC,
#'                baseMean = 2^AveExpr)
#' m <- gatom_module(de, refs, k_gene = 50, seed = 42)
#' }
#' @export
gatom_module <- function(de, refs, k_gene = 50, k_met = NULL, met_de = NULL,
                         seed = 42, solver = "rnc", verbose = FALSE,
                         gene2reaction_extra = NULL) {
  .require_pkg("gatom", "gatom_module()", 'BiocManager::install("gatom")')
  .require_pkg("mwcsr", "gatom_module()")
  .require_pkg("igraph", "gatom_module()")

  if (!inherits(refs, "gatom_refs")) {
    stop("`refs` must be a `gatom_refs` object from gatom_refs().",
         call. = FALSE)
  }
  if (!is.data.frame(de)) {
    stop("`de` must be a data frame; build it with gatom_de().", call. = FALSE)
  }
  needed <- c("ID", "pval", "log2FC", "baseMean")
  if (!all(needed %in% names(de))) {
    stop("`de` is missing column(s): ",
         paste(setdiff(needed, names(de)), collapse = ", "),
         ". Build it with gatom_de().", call. = FALSE)
  }
  if (!is.numeric(k_gene) || length(k_gene) != 1L || is.na(k_gene) ||
        k_gene <= 0) {
    stop("`k_gene` must be a single positive number (50 is the GATOM ",
         "default; smaller k gives a larger module).", call. = FALSE)
  }
  if (!is.null(k_met) &&
        (!is.numeric(k_met) || length(k_met) != 1L || is.na(k_met))) {
    stop("`k_met` must be a single number or NULL.", call. = FALSE)
  }
  if (is.null(met_de) && !is.null(k_met)) {
    stop("`k_met` was supplied but `met_de` is NULL; metabolite scoring ",
         "needs metabolite data. Pass `met_de` or leave `k_met = NULL`.",
         call. = FALSE)
  }
  if (!is.null(gene2reaction_extra) &&
        (!is.data.frame(gene2reaction_extra) ||
         !all(c("gene", "reaction") %in% names(gene2reaction_extra)))) {
    stop("`gene2reaction_extra` must be NULL or a data frame with `gene` and ",
         "`reaction` columns. Read the downloaded gene2reaction TSV with ",
         "data.table::fread() and pass the resulting table.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    stop("`seed` must be a single number; the MWCS solver is a heuristic ",
         "and unseeded runs are not reproducible.", call. = FALSE)
  }
  solver_obj <- .gatom_solver(solver)

  g <- gatom::makeMetabolicGraph(
    network = refs$network,
    topology = "atoms",
    org.gatom.anno = refs$org_anno,
    gene.de = as.data.frame(de),
    met.db = refs$met_db,          # required even when met.de is NULL
    met.de = met_de,
    gene2reaction.extra = gene2reaction_extra,
    keepReactionsWithoutEnzymes = FALSE
  )
  n_g <- igraph::vcount(g)
  e_g <- igraph::ecount(g)
  if (isTRUE(verbose)) {
    message("Graph: ", n_g, " nodes, ", e_g, " edges")
  }
  if (n_g == 0L) {
    stop("The metabolic graph is empty -- no `de$ID` matched the annotation. ",
         "Check that `id` uses the identifier type the annotation maps from ",
         "(symbols for org.*.eg.gatom.anno).", call. = FALSE)
  }

  # Seed twice: BioNet's BUM fit inside scoreGraph() uses random starts, so
  # seeding only the solver leaves the *scores* irreproducible (observed:
  # 40- vs 53-node modules from one seed at realistic scale).
  gs <- .with_pinned_seed(
    seed,
    gatom::scoreGraph(g, k.gene = k_gene, k.met = k_met)
  )
  solution <- .with_pinned_seed(
    seed,
    mwcsr::solve_mwcsp(solver_obj, gs)
  )
  m <- solution$graph

  # This name matches the downstream results-table column and avoids ambiguity.
  attr(m, "solution_weight") <- .gatom_solution_weight(solution)
  attr(m, "k_gene") <- k_gene
  attr(m, "k_met") <- k_met
  attr(m, "seed") <- seed
  attr(m, "solver") <- solver
  attr(m, "species") <- refs$species
  attr(m, "graph_nodes") <- n_g
  attr(m, "graph_edges") <- e_g
  attr(m, "n_nodes") <- igraph::vcount(m)
  attr(m, "n_edges") <- igraph::ecount(m)
  if (isTRUE(verbose)) {
    message("Module: ", igraph::vcount(m), " nodes, ", igraph::ecount(m),
            " edges")
  }
  m
}

#' Genes in a GATOM module
#'
#' In an atom-topology graph the vertices are metabolite atoms and the
#' **edges** are enzyme-catalysed reactions, so the gene labels live on the
#' edges. Reading `igraph::V(m)$Symbol` returns `NULL` and looks like an empty
#' result; this function reads the edges.
#'
#' @param m A module `igraph` from [gatom_module()].
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
  if (!"Symbol" %in% names(edges)) {
    stop("The module's edges carry no `Symbol` column. GATOM modules built ",
         "with topology = \"atoms\" label genes on edges; this graph was ",
         "either not built by gatom_module() or lost its edge attributes.",
         call. = FALSE)
  }
  sym <- as.character(edges$Symbol)
  unique(sym[!is.na(sym) & nzchar(sym)])
}

#' Save a GATOM module as a self-contained HTML view
#'
#' Wraps `gatom::saveModuleToHtml()`. `htmlwidgets` needs pandoc on the PATH
#' to build a self-contained file, and Quarto's bundled pandoc is not exported
#' to child R sessions -- so this sets `RSTUDIO_PANDOC` from a known location
#' when it is unset, and fails with an actionable message when none is found.
#' The parent directory is created if missing. This is the only `gatom_*`
#' function that writes to disk.
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
    Sys.setenv(RSTUDIO_PANDOC = hit[[1L]])
  }

  parent <- dirname(path)
  if (nzchar(parent) && !identical(parent, ".")) ensure_dir(parent)
  gatom::saveModuleToHtml(m, path, name = name)
  invisible(path)
}
