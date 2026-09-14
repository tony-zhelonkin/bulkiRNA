#' Read a master-table schema
#'
#' Uses base CSV support so loading the core schema never requires the optional
#' `yaml` package.
#'
#' @param schema_version Character(1) schema version.
#' @return A data frame describing the requested schema.
#' @keywords internal
.gs_master_schema <- function(schema_version) {
  if (!is.character(schema_version) || length(schema_version) != 1L ||
      is.na(schema_version) || !nzchar(schema_version)) {
    shown <- paste0(schema_version, collapse = ", ")
    stop("`schema_version` must be a non-empty character scalar; got ",
         encodeString(shown, quote = "\""), ".", call. = FALSE)
  }

  file <- paste0("master-schema-v", schema_version, ".csv")
  path <- system.file("extdata", file, package = "bulkiRNA")
  if (!nzchar(path)) {
    stop("Unknown `schema_version` ",
         encodeString(schema_version, quote = "\""),
         "; expected a packaged file named ",
         encodeString(file, quote = "\""), ".",
         call. = FALSE)
  }

  schema <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  expected_fields <- c("column", "type", "required", "description")
  if (!identical(names(schema), expected_fields)) {
    stop("Master schema ", encodeString(file, quote = "\""),
         " is malformed: expected fields ",
         paste(expected_fields, collapse = ", "), ".", call. = FALSE)
  }
  if (anyDuplicated(schema$column) ||
      any(!schema$type %in% c("character", "numeric", "integer")) ||
      any(!schema$required %in% c("yes", "no"))) {
    stop("Master schema ", encodeString(file, quote = "\""),
         " has invalid or duplicated declarations.", call. = FALSE)
  }
  schema$required <- schema$required == "yes"
  schema
}

#' List master-table columns in emission order
#'
#' Returns the columns declared by the versioned master-table schema in the
#' order [gs_to_master()] emits them and [gs_validate_master()] expects.
#' Optional columns come first, followed by the required columns.
#'
#' The shipped schema CSV lists `entity_type` last while `gs_to_master()` emits
#' it first, so take the order from here rather than from the file.
#'
#' @param schema_version Character(1) master-table schema version.
#' @param optional Whether to include the optional columns. `gs_to_master()`
#'   emits them only when given an `entity_type`, so the default describes a
#'   table built without one.
#' @return A character vector of master-table column names in emission order.
#' @examples
#' gs_master_columns()
#' gs_master_columns(optional = TRUE)
#' @export
gs_master_columns <- function(schema_version = "1", optional = FALSE) {
  if (!is.logical(optional) || length(optional) != 1L || is.na(optional)) {
    stop("`optional` must be TRUE or FALSE.", call. = FALSE)
  }
  schema <- .gs_master_schema(schema_version)
  required <- schema$column[schema$required]
  if (!optional) {
    return(required)
  }
  c(schema$column[!schema$required], required)
}

#' Test whether a master-table column is coercible to its declared type
#'
#' Coercion is considered safe when it does not turn a non-missing input into a
#' missing value. Integer coercion additionally rejects fractional values.
#'
#' @param x A data-frame column.
#' @param type One of `"character"`, `"numeric"`, or `"integer"`.
#' @return A list with logical `ok` and the coerced `value`.
#' @keywords internal
.gs_master_coerce <- function(x, type) {
  original_missing <- is.na(x)
  if (is.matrix(original_missing) || is.data.frame(original_missing)) {
    return(list(ok = FALSE, value = NULL))
  }

  converted <- tryCatch(
    suppressWarnings(
      switch(
        type,
        character = as.character(x),
        numeric = as.numeric(if (is.factor(x)) as.character(x) else x),
        integer = {
          numeric_value <- as.numeric(
            if (is.factor(x)) as.character(x) else x
          )
          integer_value <- as.integer(numeric_value)
          non_missing <- !is.na(numeric_value)
          if (any(non_missing & numeric_value != integer_value)) {
            stop("fractional value", call. = FALSE)
          }
          integer_value
        }
      )
    ),
    error = function(e) NULL
  )

  if (is.null(converted) || length(converted) != length(x)) {
    return(list(ok = FALSE, value = NULL))
  }
  newly_missing <- !original_missing & is.na(converted)
  list(ok = !any(newly_missing), value = converted)
}

#' Convert a gene-set result to the master-table contract
#'
#' Serializes the normalized columns of a [gs_result-class] into versioned
#' master-table columns. This function only computes and returns data; it never
#' writes a file. Pathway names pass through unchanged.
#'
#' @param res A `gs_result`.
#' @param db A named list (including a `gs_db`) mapping pathway IDs to gene
#'   vectors, or `NULL`. When `NULL`, `genes_full_set` is `NA_character_`.
#' @param universe Character vector of genes tested. Required when `db` is not
#'   `NULL` and used to derive `genes_full_set`.
#' @param entity_type Optional character scalar or one value per result row.
#'   When supplied, it is emitted as the first column.
#' @param schema_version Character(1) master-table schema version.
#' @param stat_as_nes Logical(1). Deliberate override allowing a non-NES
#'   `res$stat` to be written to `nes`.
#' @return A tibble in master-table column order, carrying a `schema_version`
#'   attribute.
#' @examples
#' ranks <- c(A = 2.5, B = 1.8, C = 0.9, D = -0.4, E = -1.2, F = -2.1)
#' db <- gsdb_register(list(SET_A = c("A", "B", "C")),
#'                     database = "demo", species = "Homo sapiens")
#' res <- gs_test(ranks, db, min_size = 2)
#'
#' # `db` and `universe` are what make `genes_full_set` computable (MADR-008).
#' gs_to_master(res, db = db, universe = names(ranks))
#' @export
gs_to_master <- function(res, db = NULL, universe = NULL, entity_type = NULL,
                         schema_version = "1", stat_as_nes = FALSE) {
  if (!inherits(res, "gs_result")) {
    stop("`res` must be a gs_result; got class ",
         encodeString(paste(class(res), collapse = "/"), quote = "\""),
         ".", call. = FALSE)
  }
  validate_gs_result(res)
  if (!is.logical(stat_as_nes) || length(stat_as_nes) != 1L ||
      is.na(stat_as_nes)) {
    stop("`stat_as_nes` must be TRUE or FALSE; got ",
         encodeString(paste0(stat_as_nes, collapse = ", "), quote = "\""),
         ".", call. = FALSE)
  }

  schema <- .gs_master_schema(schema_version)
  stat_types <- unique(res$stat_type)
  bad_stat_types <- stat_types[is.na(stat_types) | stat_types != "NES"]
  if (length(bad_stat_types) && !stat_as_nes) {
    shown <- ifelse(is.na(bad_stat_types), "<NA>", bad_stat_types)
    stop("`res$stat_type` must be \"NES\" before `stat` can populate `nes`; ",
         "got ", paste(encodeString(shown, quote = "\""), collapse = ", "),
         ". Set `stat_as_nes = TRUE` for a deliberate override.",
         call. = FALSE)
  }

  n <- nrow(res)
  if (!is.null(entity_type)) {
    if (!is.character(entity_type) ||
        !(length(entity_type) %in% c(1L, n))) {
      stop("`entity_type` must be a character scalar or have one value per ",
           "result row; got length ", length(entity_type), ".", call. = FALSE)
    }
    entity_type <- rep(entity_type, length.out = n)
  }

  if (!is.null(db)) {
    if (!is.list(db) || is.null(names(db)) || anyNA(names(db)) ||
        any(!nzchar(names(db))) || anyDuplicated(names(db)) ||
        !all(vapply(db, is.character, logical(1L)))) {
      stop("`db` must be a named list of gene vectors or `NULL`; got class ",
           encodeString(paste(class(db), collapse = "/"), quote = "\""),
           ".", call. = FALSE)
    }
    if (!is.character(universe)) {
      stop("`universe` must be a character vector when `db` is supplied; got ",
           "class ",
           encodeString(paste(class(universe), collapse = "/"), quote = "\""),
           ".",
           call. = FALSE)
    }
    missing_sets <- setdiff(unique(res$pathway_id), names(db))
    if (length(missing_sets)) {
      stop("`db` is missing `res$pathway_id` value(s): ",
           paste(encodeString(missing_sets, quote = "\""), collapse = ", "),
           ".", call. = FALSE)
    }
  }

  has_leading_edge <- "leading_edge" %in% names(res)
  if (has_leading_edge) {
    leading_edge_size <- as.integer(lengths(res$leading_edge))
    core_enrichment <- vapply(
      res$leading_edge,
      paste,
      collapse = "/",
      FUN.VALUE = character(1L),
      USE.NAMES = FALSE
    )
  } else {
    leading_edge_size <- rep(NA_integer_, n)
    core_enrichment <- rep(NA_character_, n)
  }

  if (is.null(db)) {
    genes_full_set <- rep(NA_character_, n)
  } else {
    genes_full_set <- vapply(
      res$pathway_id,
      function(pathway_id) {
        paste(intersect(db[[pathway_id]], universe), collapse = "/")
      },
      FUN.VALUE = character(1L),
      USE.NAMES = FALSE
    )
  }

  set_size <- as.integer(res$n_genes_tested)
  out <- tibble::tibble(
    pathway_id = as.character(res$pathway_id),
    pathway_name = as.character(res$pathway_name),
    database = as.character(res$database),
    contrast = as.character(res$contrast),
    nes = as.numeric(res$stat),
    pvalue = as.numeric(res$p_value),
    padj = as.numeric(res$padj),
    set_size = set_size,
    leading_edge_size = leading_edge_size,
    gene_ratio = as.numeric(leading_edge_size / set_size),
    core_enrichment = core_enrichment,
    genes_full_set = genes_full_set,
    direction = unname(
      c(up = "Up", down = "Down", ns = "NS")[res$direction]
    ),
    neg_log_padj = -log10(pmax(as.numeric(res$padj),
                               .Machine$double.xmin))
  )

  required <- schema$column[schema$required]
  if (!is.null(entity_type)) {
    out <- tibble::add_column(out, entity_type = entity_type, .before = 1L)
  }
  expected <- c(if (!is.null(entity_type)) "entity_type", required)
  out <- out[, expected, drop = FALSE]
  attr(out, "schema_version") <- schema_version

  gs_validate_master(out, schema_version = schema_version)
  out
}

#' Validate a master GSEA table
#'
#' Checks the versioned schema, exact column order and declared types, then
#' verifies master-table semantic invariants. With `error = TRUE`, every
#' detected problem is included in one error.
#'
#' @param df A data frame intended to satisfy the master-table schema.
#' @param schema_version Character(1) master-table schema version.
#' @param error Logical(1). Stop with all problems when `TRUE`; return the
#'   problem tibble when `FALSE`.
#' @return With `error = FALSE`, a tibble with `check`, `column`, and `message`
#'   columns (zero rows when valid). With `error = TRUE`, `df` invisibly when
#'   valid, otherwise an error.
#' @examples
#' ranks <- c(A = 2.5, B = 1.8, C = 0.9, D = -0.4, E = -1.2, F = -2.1)
#' db <- gsdb_register(list(SET_A = c("A", "B", "C")),
#'                     database = "demo", species = "Homo sapiens")
#' master <- gs_to_master(gs_test(ranks, db, min_size = 2))
#'
#' gs_validate_master(master, error = FALSE)   # zero rows means clean
#'
#' # The defect this exists to catch: a derived column NA-filled beside a
#' # finite `padj`, which is what a column allowlist plus rbind() produces.
#' broken <- master
#' broken$neg_log_padj <- NA_real_
#' gs_validate_master(broken, error = FALSE)
#' @export
gs_validate_master <- function(df, schema_version = "1", error = TRUE) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame; got class ",
         encodeString(paste(class(df), collapse = "/"), quote = "\""),
         ".", call. = FALSE)
  }
  if (!is.logical(error) || length(error) != 1L || is.na(error)) {
    stop("`error` must be TRUE or FALSE; got ",
         encodeString(paste0(error, collapse = ", "), quote = "\""),
         ".", call. = FALSE)
  }

  schema <- .gs_master_schema(schema_version)
  required <- schema$column[schema$required]
  optional <- schema$column[!schema$required]
  present_optional <- optional[optional %in% names(df)]
  expected_order <- c(present_optional, required)

  checks <- character()
  columns <- character()
  messages <- character()
  add_problem <- function(check, column, message) {
    checks <<- c(checks, check)
    columns <<- c(columns, column)
    messages <<- c(messages, message)
  }

  add_missing_problem <- function(check, column, rows, condition) {
    row_indices <- which(rows)
    if (!length(row_indices)) {
      return(invisible(NULL))
    }
    example_indices <- utils::head(row_indices, 3L)
    if ("pathway_id" %in% names(coerced)) {
      example_ids <- coerced$pathway_id[example_indices]
      example_ids[is.na(example_ids)] <- "<NA>"
    } else {
      example_ids <- rep("<unavailable>", length(example_indices))
    }
    row_label <- if (length(row_indices) == 1L) "row" else "rows"
    add_problem(
      check,
      column,
      paste0(
        "Column `", column, "` is NA on ", length(row_indices),
        " ", row_label, " ", condition, "; example pathway_id(s): ",
        paste(encodeString(example_ids, quote = "\""), collapse = ", "), "."
      )
    )
  }

  missing <- setdiff(required, names(df))
  if (length(missing)) {
    add_problem(
      "required_columns",
      paste(missing, collapse = ", "),
      paste0("Missing required column(s): ", paste(missing, collapse = ", "),
             ".")
    )
  }

  unexpected <- setdiff(names(df), schema$column)
  if (length(unexpected)) {
    add_problem(
      "declared_columns",
      paste(unexpected, collapse = ", "),
      paste0("Unexpected column(s) not declared by schema v", schema_version,
             ": ", paste(unexpected, collapse = ", "), ".")
    )
  }

  if (!identical(names(df), expected_order)) {
    add_problem(
      "column_order",
      NA_character_,
      paste0(
        "Column order does not match schema v", schema_version, "; expected ",
        paste(expected_order, collapse = ", "), "; got ",
        paste(names(df), collapse = ", "), "."
      )
    )
  }

  coerced <- list()
  declared_present <- intersect(schema$column, names(df))
  for (column in declared_present) {
    declared_type <- schema$type[match(column, schema$column)]
    result <- .gs_master_coerce(df[[column]], declared_type)
    if (!result$ok) {
      shown <- tryCatch(
        paste(utils::head(df[[column]], 3L), collapse = ", "),
        error = function(e) "<unprintable>"
      )
      add_problem(
        "column_type",
        column,
        paste0("Column `", column, "` is not coercible to ",
               declared_type, "; got ",
               encodeString(shown, quote = "\""), ".")
      )
    } else {
      coerced[[column]] <- result$value
    }
  }

  if (all(c("padj", "neg_log_padj") %in% names(coerced))) {
    rows <- !is.na(coerced$padj) & is.na(coerced$neg_log_padj)
    add_missing_problem(
      "neg_log_padj_missing",
      "neg_log_padj",
      rows,
      "where `padj` is non-NA"
    )
  }

  if (all(c("core_enrichment", "leading_edge_size") %in% names(coerced))) {
    rows <- !is.na(coerced$core_enrichment) &
      nzchar(coerced$core_enrichment) &
      is.na(coerced$leading_edge_size)
    add_missing_problem(
      "leading_edge_size_missing",
      "leading_edge_size",
      rows,
      "where `core_enrichment` is non-NA and non-empty"
    )
  }

  if (all(c("leading_edge_size", "set_size", "gene_ratio") %in%
          names(coerced))) {
    rows <- !is.na(coerced$leading_edge_size) &
      !is.na(coerced$set_size) & coerced$set_size > 0 &
      is.na(coerced$gene_ratio)
    add_missing_problem(
      "gene_ratio_missing",
      "gene_ratio",
      rows,
      paste0("where `leading_edge_size` and `set_size` are non-NA and ",
             "`set_size` is greater than zero")
    )
  }

  for (column in intersect(c("pathway_id", "database", "contrast"),
                           names(coerced))) {
    rows <- is.na(coerced[[column]])
    add_missing_problem(
      "identity_missing",
      column,
      rows,
      "although identity columns must be non-NA on every row"
    )
  }

  if ("direction" %in% names(coerced)) {
    direction <- coerced$direction
    bad_direction <- unique(direction[
      !is.na(direction) & !direction %in% c("Up", "Down", "NS")
    ])
    if (length(bad_direction)) {
      add_problem(
        "direction_values",
        "direction",
        paste0("Column `direction` must contain only \"Up\", \"Down\", or ",
               "\"NS\"; got ",
               paste(encodeString(bad_direction, quote = "\""), collapse = ", "),
               ".")
      )
    }
  }

  if (all(c("padj", "neg_log_padj") %in% names(coerced))) {
    padj <- coerced$padj
    actual <- coerced$neg_log_padj
    compare <- !is.na(padj) & !is.na(actual)
    if (any(compare)) {
      expected <- -log10(pmax(padj[compare], .Machine$double.xmin))
      observed <- actual[compare]
      tolerance <- sqrt(.Machine$double.eps)
      agrees <- observed == expected
      finite <- is.finite(observed) & is.finite(expected)
      agrees[finite] <- abs(observed[finite] - expected[finite]) <=
        tolerance * pmax(1, abs(expected[finite]))
      if (any(!agrees)) {
        add_problem(
          "neg_log_padj_values",
          "neg_log_padj",
          paste0(
            "Column `neg_log_padj` does not equal ",
            "-log10(pmax(padj, .Machine$double.xmin)); the retired ",
            "cap-at-16 convention is not allowed."
          )
        )
      }
    }
  }

  problems <- tibble::tibble(
    check = checks,
    column = columns,
    message = messages
  )
  if (!error) {
    return(problems)
  }
  if (nrow(problems)) {
    stop(
      "`df` failed master schema v", schema_version, " validation:\n- ",
      paste(problems$message, collapse = "\n- "),
      call. = FALSE
    )
  }
  invisible(df)
}
