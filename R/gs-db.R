#' The `gs_db` gene-set database object
#'
#' Every `gsdb_*` provider returns a `gs_db`: a named list of character
#' vectors mapping set id to gene symbols, which is exactly the shape
#' [fgsea::fgseaMultilevel()], [fgsea::fora()] and `GSVA::gsvaParam()` take.
#'
#' Attributes carried alongside the list:
#'
#' - `pathway_names` -- named character, id -> human-readable label, with
#'   `names()` identical to `names(sets)`.
#' - `database` -- the **machine-typeable registry key** (`"mitopathways"`,
#'   `"msigdb_H"`); lands in `gs_result$database`, where it is a join and
#'   filter key, so it must be stable rather than pretty.
#' - `database_label` -- the display string (`"MitoPathways 3.0"`). Renderers
#'   prettify from here; the data layer never depends on it.
#' - `species` -- `"Mus musculus"` or `"Homo sapiens"`.
#' - `gene_id_type` -- `"symbol"` (the only value today).
#' - `pathway_descriptions` -- optional named character, id -> full-sentence
#'   description (e.g. MSigDB's `gs_description`), or `NULL` when a provider
#'   has none. Unlike `pathway_names` this is not guaranteed to cover every
#'   id. Carried through `[.gs_db` like every other attribute.
#' - `set_provenance` -- optional data frame with exactly one row per set and a
#'   unique `set_name` key. It is subset with the database.
#' - `provenance` -- optional one-row data frame or named list of scalar values
#'   describing the database as a whole. It is carried through subsetting.
#'
#' Empty sets are dropped at construction, genes are de-duplicated within a
#' set, and set ids must be unique.
#'
#' @param sets Named list of character vectors: set id -> gene symbols.
#' @param database Character(1) registry key, e.g. `"mito_unified"`.
#' @param species Character(1) species, e.g. `"Mus musculus"`.
#' @param pathway_names Named character of human-readable labels, or `NULL`
#'   to derive them from the ids.
#' @param database_label Character(1) display name, or `NULL` to reuse
#'   `database`.
#' @param gene_id_type Character(1); only `"symbol"` is supported.
#' @param pathway_descriptions Named character of longer descriptions, or
#'   `NULL` when the provider has none. Not required to cover every id.
#' @param set_provenance A data frame with one row per set and a unique
#'   `set_name` column, or `NULL`.
#' @param provenance A one-row data frame or named list of scalar values, or
#'   `NULL`.
#' @return A `gs_db` object.
#' @keywords internal
gs_db <- function(sets,
                  database,
                  species,
                  pathway_names = NULL,
                  database_label = NULL,
                  gene_id_type = "symbol",
                  pathway_descriptions = NULL,
                  set_provenance = NULL,
                  provenance = NULL) {
  if (!is.list(sets) || (length(sets) > 0L && is.null(names(sets)))) {
    stop("`sets` must be a named list of character vectors.", call. = FALSE)
  }
  if (!is.character(database) || length(database) != 1L || is.na(database) ||
        !nzchar(database)) {
    stop("`database` must be a single non-empty string.", call. = FALSE)
  }
  if (!is.null(database_label) &&
        (!is.character(database_label) || length(database_label) != 1L ||
           is.na(database_label) || !nzchar(database_label))) {
    stop("`database_label` must be a single non-empty string or NULL.",
         call. = FALSE)
  }
  species <- .gsdb_species_label(species)
  if (!identical(gene_id_type, "symbol")) {
    stop("`gene_id_type` must be \"symbol\"; got ",
         paste0("\"", gene_id_type, "\""), ".", call. = FALSE)
  }

  ids <- names(sets)
  if (length(ids) && (anyNA(ids) || any(!nzchar(ids)))) {
    stop("`sets` must have non-missing, non-empty names.", call. = FALSE)
  }
  if (anyDuplicated(ids)) {
    dup <- unique(ids[duplicated(ids)])
    stop("`sets` names must be unique; duplicated: ",
         paste0("`", utils::head(dup, 5L), "`", collapse = ", "), ".",
         call. = FALSE)
  }

  sets <- lapply(sets, function(g) {
    g <- as.character(g)
    unique(g[!is.na(g) & nzchar(g)])
  })

  keep <- vapply(sets, length, integer(1L)) > 0L
  sets <- sets[keep]

  pathway_names <- .gsdb_resolve_names(pathway_names, names(sets))

  if (!is.null(pathway_descriptions)) {
    if (length(pathway_descriptions) && is.null(names(pathway_descriptions))) {
      stop("`pathway_descriptions` must be a *named* character vector ",
           "(id -> description).", call. = FALSE)
    }
    hit <- intersect(names(sets), names(pathway_descriptions))
    pathway_descriptions <- if (length(hit)) {
      stats::setNames(as.character(pathway_descriptions[hit]), hit)
    } else {
      NULL
    }
  }

  provenance <- .gsdb_validate_provenance(provenance)
  set_provenance <- .gsdb_validate_set_provenance(
    set_provenance,
    names(sets)
  )

  structure(
    sets,
    pathway_names  = pathway_names,
    database       = database,
    database_label = database_label %||% database,
    species        = species,
    gene_id_type   = gene_id_type,
    pathway_descriptions = pathway_descriptions,
    set_provenance = set_provenance,
    provenance     = provenance,
    class          = "gs_db"
  )
}

#' Validate database-level provenance
#'
#' @param provenance A one-row data frame, named list, or `NULL`.
#' @return `provenance`, unchanged.
#' @keywords internal
.gsdb_validate_provenance <- function(provenance) {
  if (is.null(provenance)) return(NULL)

  if (is.data.frame(provenance)) {
    if (nrow(provenance) != 1L || !ncol(provenance) ||
        anyNA(names(provenance)) || any(!nzchar(names(provenance))) ||
        anyDuplicated(names(provenance))) {
      stop("`provenance` must be a one-row data frame with unique, ",
           "non-empty column names.", call. = FALSE)
    }
    if (any(!vapply(provenance, function(x) {
      is.atomic(x) && length(x) == 1L
    }, logical(1L)))) {
      stop("`provenance` data-frame columns must each contain one atomic ",
           "scalar value.", call. = FALSE)
    }
    return(provenance)
  }

  if (!is.list(provenance) || !length(provenance) ||
      is.null(names(provenance)) || anyNA(names(provenance)) ||
      any(!nzchar(names(provenance))) || anyDuplicated(names(provenance)) ||
      any(!vapply(provenance, function(x) is.atomic(x) && length(x) == 1L,
                  logical(1L)))) {
    stop("`provenance` must be a one-row data frame or a named list of ",
         "scalar values.", call. = FALSE)
  }
  provenance
}

#' Validate set-level provenance
#'
#' @param set_provenance A set-keyed data frame or `NULL`.
#' @param set_names Character vector of retained set names.
#' @return `set_provenance`, ordered to match `set_names`.
#' @keywords internal
.gsdb_validate_set_provenance <- function(set_provenance, set_names) {
  if (is.null(set_provenance)) return(NULL)
  if (!is.data.frame(set_provenance) ||
      !"set_name" %in% names(set_provenance)) {
    stop("`set_provenance` must be a data frame with a `set_name` column.",
         call. = FALSE)
  }
  keys <- set_provenance$set_name
  if (!is.character(keys) || anyNA(keys) || any(!nzchar(keys)) ||
      anyDuplicated(keys)) {
    stop("`set_provenance$set_name` must contain unique, non-missing, ",
         "non-empty strings.", call. = FALSE)
  }
  missing <- setdiff(set_names, keys)
  extra <- setdiff(keys, set_names)
  if (length(missing) || length(extra)) {
    stop("`set_provenance$set_name` must match the retained set names ",
         "exactly; missing: ",
         if (length(missing)) paste(missing, collapse = ", ") else "none",
         "; extra: ",
         if (length(extra)) paste(extra, collapse = ", ") else "none", ".",
         call. = FALSE)
  }
  set_provenance[match(set_names, keys), , drop = FALSE]
}

#' Is this a `gs_db`?
#'
#' @param x Any object.
#' @return `TRUE` or `FALSE`.
#' @keywords internal
is_gs_db <- function(x) inherits(x, "gs_db")

#' Normalise a species label
#'
#' Known aliases are normalized through [.species()]. Unknown non-empty labels
#' remain supported for user-supplied databases and have underscores replaced
#' by spaces.
#'
#' @param species Character(1).
#' @return Character(1) with underscores replaced by spaces.
#' @keywords internal
.gsdb_species_label <- function(species) {
  .species(species, allow_custom = TRUE)$scientific
}

#' Species label in directory form
#'
#' @param species Character(1).
#' @return Character(1) with spaces replaced by underscores.
#' @keywords internal
.gsdb_species_dir <- function(species) {
  gsub(" ", "_", .gsdb_species_label(species), fixed = TRUE)
}

#' Resolve pathway labels against a set of ids
#'
#' @param pathway_names Named character or `NULL`.
#' @param ids Character vector of set ids.
#' @return Named character with `names()` equal to `ids`.
#' @keywords internal
.gsdb_resolve_names <- function(pathway_names, ids) {
  out <- stats::setNames(ids, ids)
  if (!is.null(pathway_names) && length(pathway_names)) {
    if (is.null(names(pathway_names))) {
      stop("`pathway_names` must be a *named* character vector ",
           "(id -> label).", call. = FALSE)
    }
    hit <- intersect(ids, names(pathway_names))
    out[hit] <- as.character(pathway_names[hit])
  }
  out[is.na(out) | !nzchar(out)] <- ids[is.na(out) | !nzchar(out)]
  out
}

#' Print a `gs_db`
#'
#' Database provenance wraps to the current console width, with a 40-character
#' minimum, and continuation lines are indented under the label.
#'
#' @param x A `gs_db`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.gs_db <- function(x, ...) {
  sizes <- vapply(x, length, integer(1L))
  label <- attr(x, "database_label")
  key <- attr(x, "database")
  cat(sprintf("<gs_db> %s  (%s, %s)\n",
              if (identical(label, key)) key else paste0(label, " [", key, "]"),
              attr(x, "species"), attr(x, "gene_id_type")))
  cat(sprintf("%d sets, %d unique genes",
              length(x), length(unique(unlist(x, use.names = FALSE)))))
  if (length(sizes)) {
    cat(sprintf(", size %d-%d (median %g)",
                min(sizes), max(sizes), stats::median(sizes)))
  }
  cat("\n")
  provenance <- attr(x, "provenance")
  if (!is.null(provenance)) {
    values <- if (is.data.frame(provenance)) {
      as.list(provenance[1L, , drop = FALSE])
    } else {
      provenance
    }
    rendered <- vapply(values, function(value) {
      if (is.na(value)) "NA" else as.character(value)
    }, character(1L))
    width <- getOption("width", 80L)
    if (!is.numeric(width) || length(width) != 1L || is.na(width) ||
        !is.finite(width)) {
      width <- 80L
    }
    provenance_lines <- strwrap(
      paste0(names(rendered), "=", rendered, collapse = "; "),
      width = max(40L, as.integer(width)),
      prefix = "            ",
      initial = "Provenance: "
    )
    cat(paste(provenance_lines, collapse = "\n"), "\n", sep = "")
  }
  if (length(x)) {
    show <- utils::head(names(x), 3L)
    for (nm in show) {
      cat(sprintf("  %s (%d)\n", nm, length(x[[nm]])))
    }
    if (length(x) > 3L) cat(sprintf("  ... %d more\n", length(x) - 3L))
  }
  invisible(x)
}

#' Summarise a `gs_db` as one row per set
#'
#' @param object A `gs_db`.
#' @param ... Ignored.
#' @return A data frame with one row per set: `pathway_id`, `pathway_name`,
#'   `n_genes`.
#' @export
summary.gs_db <- function(object, ...) {
  data.frame(
    pathway_id   = names(object),
    pathway_name = unname(attr(object, "pathway_names")[names(object)]),
    n_genes      = unname(vapply(object, length, integer(1L))),
    row.names    = NULL,
    stringsAsFactors = FALSE
  )
}

#' Subset a `gs_db`, keeping its attributes
#'
#' @param x A `gs_db`.
#' @param i Index (name, position, or logical) of sets to keep.
#' @return A `gs_db` with the selected sets.
#' @export
`[.gs_db` <- function(x, i) {
  if (is.character(i)) {
    bad <- setdiff(i, names(x))
    if (length(bad)) {
      stop("`i` selects sets not in this database: ",
           paste0("\"", utils::head(bad, 5L), "\"", collapse = ", "), ".",
           call. = FALSE)
    }
  }
  sets <- unclass(x)[i]
  set_provenance <- attr(x, "set_provenance")
  if (!is.null(set_provenance)) {
    set_provenance <- set_provenance[
      match(names(sets), set_provenance$set_name),
      ,
      drop = FALSE
    ]
  }
  gs_db(
    sets,
    database       = attr(x, "database"),
    species        = attr(x, "species"),
    pathway_names  = attr(x, "pathway_names"),
    database_label = attr(x, "database_label"),
    gene_id_type   = attr(x, "gene_id_type"),
    pathway_descriptions = attr(x, "pathway_descriptions"),
    set_provenance = set_provenance,
    provenance = attr(x, "provenance")
  )
}

#' Filter a `gs_db` by set size
#'
#' The size filter that used to live in `build_reference_databases.R`. `NULL`
#' bounds mean "no bound".
#'
#' Named with a leading dot so it cannot collide with the exported, deprecated
#' `filter_by_size(result, min_size = 5, max_size = 500)` shim. Two top-level
#' bindings of one name silently resolve by collation order, which is how the
#' shim was being shadowed before this rename.
#'
#' @param db A `gs_db`.
#' @param min_size Integer(1) or `NULL`; smallest set size to keep.
#' @param max_size Integer(1) or `NULL`; largest set size to keep.
#' @param verbose Logical(1); message how many sets were dropped.
#' @return A `gs_db` with out-of-range sets removed.
#' @keywords internal
#' @name dot-gs_filter_size
.gs_filter_size <- function(db, min_size = NULL, max_size = NULL,
                            verbose = FALSE) {
  if (!is_gs_db(db)) {
    stop("`db` must be a `gs_db` object; see `gsdb_load()`.", call. = FALSE)
  }
  if (is.null(min_size) && is.null(max_size)) {
    return(db)
  }
  lo <- if (is.null(min_size)) 1L else as.integer(min_size)
  hi <- if (is.null(max_size)) Inf else as.numeric(max_size)
  if (!is.na(lo) && !is.na(hi) && lo > hi) {
    stop("`min_size` (", lo, ") must not exceed `max_size` (", hi, ").",
         call. = FALSE)
  }
  sizes <- vapply(db, length, integer(1L))
  keep <- sizes >= lo & sizes <= hi
  if (verbose && any(!keep)) {
    message(sprintf("Size filter [%s, %s]: dropped %d of %d sets.",
                    lo, hi, sum(!keep), length(keep)))
  }
  db[keep]
}

#' Convert a `gs_db` to the legacy T2G/T2N pair
#'
#' Exists **solely** so the deprecation shim in `R/deprecated.R` can keep
#' `load_reference_db()` returning its old `clusterProfiler` shape. Nothing
#' new calls it.
#'
#' @param db A `gs_db`.
#' @return `list(T2G = data.frame(gs_name, gene_symbol), T2N =
#'   data.frame(gs_name, description))`.
#' @keywords internal
.gsdb_as_t2g <- function(db) {
  if (!is_gs_db(db)) {
    stop("`db` must be a `gs_db` object.", call. = FALSE)
  }
  ids <- names(db)
  n <- vapply(db, length, integer(1L))
  list(
    T2G = data.frame(
      gs_name     = rep(ids, times = n),
      gene_symbol = unlist(db, use.names = FALSE),
      row.names   = NULL,
      stringsAsFactors = FALSE
    ),
    T2N = data.frame(
      gs_name     = ids,
      description = unname(attr(db, "pathway_names")[ids]),
      row.names   = NULL,
      stringsAsFactors = FALSE
    )
  )
}

#' Build a `gs_db` from a legacy T2G/T2N pair
#'
#' Reads the shipped processed RDS shape (and the parser output of the
#' internal rebuild driver) into the frozen `gs_db` contract.
#'
#' @param x List with `T2G` (`gs_name`, `gene_symbol`) and optional `T2N`
#'   (`gs_name`, `description`).
#' @param database Character(1) registry key.
#' @param species Character(1) species.
#' @param database_label Character(1) display name, or `NULL`.
#' @return A `gs_db`.
#' @keywords internal
.gsdb_from_t2g <- function(x, database, species, database_label = NULL) {
  if (!is.list(x) || is.null(x$T2G)) {
    stop("`x` must be a list with a `T2G` data frame.", call. = FALSE)
  }
  t2g <- x$T2G
  if (!all(c("gs_name", "gene_symbol") %in% names(t2g))) {
    stop("`T2G` must have columns `gs_name` and `gene_symbol`; got ",
         paste0("`", names(t2g), "`", collapse = ", "), ".", call. = FALSE)
  }
  sets <- split(as.character(t2g$gene_symbol), as.character(t2g$gs_name))

  labels <- NULL
  if (!is.null(x$T2N) && all(c("gs_name", "description") %in% names(x$T2N))) {
    t2n <- x$T2N[!duplicated(x$T2N$gs_name), , drop = FALSE]
    labels <- stats::setNames(as.character(t2n$description),
                              as.character(t2n$gs_name))
  }

  gs_db(sets, database = database, species = species,
        pathway_names = labels, database_label = database_label)
}
