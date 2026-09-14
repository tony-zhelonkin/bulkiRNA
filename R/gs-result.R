#' The `gs_result` contract
#'
#' `gs_result` is the single normalized shape every gene-set test returns,
#' whatever the method. It is a [tibble][tibble::tibble] subclass with one row
#' per pathway per contrast.
#'
#' # Core columns
#'
#' Always present, always these types:
#'
#' \describe{
#'   \item{`pathway_id`}{character. Stable identifier, e.g. `"HALLMARK_HYPOXIA"`.}
#'   \item{`pathway_name`}{character. Human-readable label for axes.}
#'   \item{`database`}{character. Which `gsdb_*` provider the set came from.}
#'   \item{`contrast`}{character. Comparison label, e.g. `"KO-WT"`. Pooling
#'     across contrasts is `rbind()`, not a separate code path.}
#'   \item{`method`}{character. `"fgsea"`, `"ora"`, `"geseca"`, or -- for a
#'     [gs_matrix] tested with limma -- the scoring method that produced the
#'     matrix (`"gsva"`, `"ssgsea"`, `"zscore"`, `"plage"`).}
#'   \item{`n_genes`}{integer. Set size in the database.}
#'   \item{`n_genes_tested`}{integer. Set size after intersecting with the data.}
#'   \item{`stat`}{numeric. The effect size. What it *is* is named by `stat_type`.}
#'   \item{`stat_type`}{character. One of [gs_stat_types()].}
#'   \item{`direction`}{character. `"up"`, `"down"`, or `"ns"`; `NA` when a
#'     method's statistic is unsigned.}
#'   \item{`p_value`}{numeric. Raw p-value.}
#'   \item{`padj`}{numeric. Multiplicity-adjusted p-value.}
#' }
#'
#' # Optional columns
#'
#' Method-specific, dropped when not applicable: `es` (raw enrichment score),
#' `log2err` (the estimator's own bound on the log2 error of the p-value),
#' `leading_edge` (list of character vectors), `fold_enrichment`, `overlap`.
#' `log2err = Inf` means the estimate is past reliable resolution, not that the
#' computation failed.
#'
#' Renderers label their statistic axis from `stat_type`, so a GSVA
#' t-statistic is never mislabelled "NES". Never write a `stat` column without
#' the matching `stat_type`.
#'
#' @name gs_result-class
NULL

#' Recognised `stat_type` values
#'
#' The controlled vocabulary for the `stat_type` column of a [gs_result].
#' Names are the values; elements are the axis labels renderers use.
#'
#' @return A named character vector of axis labels.
#' @examples
#' gs_stat_types()
#' @export
gs_stat_types <- function() {
  c(
    NES                  = "NES",
    t                    = "t statistic",
    log2_fold_enrichment = "log2 fold enrichment",
    signed_log10p        = "signed -log10 p",
    pct_var              = "% variance explained"
  )
}

#' Strip a tibble subclass down to a plain tibble
#'
#' @param x A data frame.
#' @return A plain tibble with the same columns.
#' @keywords internal
.as_plain_tibble <- function(x) {
  class(x) <- c("tbl_df", "tbl", "data.frame")
  x
}

# Core column names and their coercion functions, in canonical order.
.gs_core_cols <- c(
  "pathway_id", "pathway_name", "database", "contrast", "method",
  "n_genes", "n_genes_tested", "stat", "stat_type", "direction",
  "p_value", "padj"
)

.gs_optional_cols <- c(
  "es", "log2err", "leading_edge", "fold_enrichment", "overlap"
)

.gs_col_coerce <- list(
  pathway_id     = as.character,
  pathway_name   = as.character,
  database       = as.character,
  contrast       = as.character,
  method         = as.character,
  n_genes        = as.integer,
  n_genes_tested = as.integer,
  stat           = as.numeric,
  stat_type      = as.character,
  direction      = as.character,
  p_value        = as.numeric,
  padj           = as.numeric,
  es             = as.numeric,
  log2err        = as.numeric,
  fold_enrichment = as.numeric,
  overlap        = as.integer
)

#' Construct a `gs_result`
#'
#' The constructor every compute adapter funnels through. Coerces core columns
#' to their contracted types, derives `direction` from the sign of `stat` when
#' absent, orders columns canonically, and validates.
#'
#' @param x A data frame with at least the core columns of [gs_result-class],
#'   minus those supplied via `...` defaults below.
#' @param database,contrast,method,stat_type Scalar defaults, recycled to fill
#'   the corresponding column when `x` does not already carry it. Supplying a
#'   scalar here is the common case: one call to one adapter produces one
#'   database/contrast/method/stat_type.
#' @return A `gs_result`.
#' @keywords internal
gs_result <- function(x,
                      database = NULL,
                      contrast = NULL,
                      method = NULL,
                      stat_type = NULL) {
  if (!is.data.frame(x)) {
    stop("`x` must be a data frame.", call. = FALSE)
  }
  x <- tibble::as_tibble(x)

  defaults <- list(
    database = database, contrast = contrast,
    method = method, stat_type = stat_type
  )
  for (nm in names(defaults)) {
    val <- defaults[[nm]]
    if (is.null(val)) next
    if (length(val) != 1L) {
      stop("`", nm, "` must be a length-1 value.", call. = FALSE)
    }
    if (is.null(x[[nm]])) x[[nm]] <- val
  }

  missing <- setdiff(setdiff(.gs_core_cols, "direction"), names(x))
  if (length(missing)) {
    stop(
      "`x` is missing required gs_result column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(x[["direction"]])) {
    x[["direction"]] <- gs_direction(x[["stat"]])
  }

  for (nm in intersect(names(.gs_col_coerce), names(x))) {
    x[[nm]] <- .gs_col_coerce[[nm]](x[[nm]])
  }

  keep <- c(
    .gs_core_cols,
    intersect(.gs_optional_cols, names(x)),
    setdiff(names(x), c(.gs_core_cols, .gs_optional_cols))
  )
  x <- x[, keep, drop = FALSE]

  validate_gs_result(new_gs_result(x))
}

#' Low-level `gs_result` constructor
#'
#' Attaches the class without checking anything. Use [gs_result()] unless you
#' have just validated `x` yourself.
#'
#' @param x A tibble.
#' @return A `gs_result`.
#' @keywords internal
new_gs_result <- function(x) {
  stopifnot(is.data.frame(x))
  tibble::new_tibble(x, class = "gs_result", nrow = nrow(x))
}

#' Validate a `gs_result`
#'
#' @param x Object to check.
#' @return `x`, invisibly on success; errors otherwise.
#' @keywords internal
validate_gs_result <- function(x) {
  if (!inherits(x, "gs_result")) {
    stop("`x` must be a gs_result; construct it with `gs_test()` or ",
         "`gs_coregulation()`.", call. = FALSE)
  }
  missing <- setdiff(.gs_core_cols, names(x))
  if (length(missing)) {
    stop(
      "`x` is missing core gs_result column(s): ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }

  bad_type <- character()
  for (nm in .gs_core_cols) {
    coerce <- .gs_col_coerce[[nm]]
    expected <- if (identical(coerce, as.character)) {
      is.character(x[[nm]])
    } else if (identical(coerce, as.integer)) {
      is.integer(x[[nm]])
    } else {
      is.numeric(x[[nm]])
    }
    if (!expected) bad_type <- c(bad_type, nm)
  }
  if (length(bad_type)) {
    stop(
      "`x` has gs_result column(s) of the wrong type: ",
      paste(bad_type, collapse = ", "), call. = FALSE
    )
  }

  st <- unique(stats::na.omit(x[["stat_type"]]))
  unknown <- setdiff(st, names(gs_stat_types()))
  if (length(unknown)) {
    stop(
      "Unknown stat_type in `x$stat_type`: ",
      paste(encodeString(unknown, quote = "\""), collapse = ", "),
      ". Allowed: ", paste(names(gs_stat_types()), collapse = ", "),
      call. = FALSE
    )
  }

  dir <- unique(stats::na.omit(x[["direction"]]))
  bad_dir <- setdiff(dir, c("up", "down", "ns"))
  if (length(bad_dir)) {
    stop(
      "`direction` must be \"up\", \"down\" or \"ns\"; got: ",
      paste(encodeString(bad_dir, quote = "\""), collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(x[["leading_edge"]]) && !is.list(x[["leading_edge"]])) {
    stop("`leading_edge` must be a list column.", call. = FALSE)
  }

  if (!is.null(x[["log2err"]])) {
    if (!is.numeric(x[["log2err"]])) {
      stop("`log2err` must be numeric.", call. = FALSE)
    }
    negative <- !is.na(x[["log2err"]]) & x[["log2err"]] < 0
    if (any(negative)) {
      stop("`log2err` must be non-negative when present.", call. = FALSE)
    }
  }

  invisible(x)
}

#' Classify a statistic's sign
#'
#' @param stat Numeric vector.
#' @return Character vector of `"up"`, `"down"`, `"ns"` (for zero or `NA`).
#' @keywords internal
gs_direction <- function(stat) {
  out <- rep("ns", length(stat))
  out[!is.na(stat) & stat > 0] <- "up"
  out[!is.na(stat) & stat < 0] <- "down"
  out
}

#' Axis label for a `gs_result`'s statistic
#'
#' Renderers call this instead of hard-coding `"NES"`.
#'
#' @param x A `gs_result`, or a character vector of `stat_type` values.
#' @return A length-1 character label.
#' @keywords internal
gs_stat_label <- function(x) {
  st <- if (inherits(x, "gs_result")) unique(stats::na.omit(x[["stat_type"]])) else unique(x)
  labs <- gs_stat_types()
  if (length(st) != 1L) {
    return("statistic")
  }
  unname(labs[[st]])
}

# ---- S3 methods -------------------------------------------------------------

#' @importFrom tibble tbl_sum
#' @exportS3Method tibble::tbl_sum
tbl_sum.gs_result <- function(x, ...) {
  dbs <- unique(x[["database"]])
  cons <- unique(x[["contrast"]])
  meth <- unique(x[["method"]])
  c(
    "gs_result" = paste0(
      nrow(x), " x ", ncol(x),
      " [", paste(meth, collapse = "/"), "]"
    ),
    "Databases" = paste(dbs, collapse = ", "),
    "Contrasts" = paste(cons, collapse = ", ")
  )
}

#' Summarise a `gs_result`
#'
#' Counts pathways per database and contrast, and how many pass an FDR cutoff.
#'
#' @param object A `gs_result`.
#' @param padj_cutoff FDR threshold for the significance counts.
#' @param ... Ignored.
#' @return A tibble, one row per database x contrast.
#' @export
summary.gs_result <- function(object, padj_cutoff = 0.05, ...) {
  sig <- !is.na(object[["padj"]]) & object[["padj"]] < padj_cutoff
  key <- paste(object[["database"]], object[["contrast"]], sep = "\r")
  split_key <- unique(key)
  tibble::tibble(
    database = vapply(strsplit(split_key, "\r", fixed = TRUE), `[`, "", 1L),
    contrast = vapply(strsplit(split_key, "\r", fixed = TRUE), `[`, "", 2L),
    n_pathways = vapply(split_key, function(k) sum(key == k), integer(1L),
                        USE.NAMES = FALSE),
    n_sig = vapply(split_key, function(k) sum(key == k & sig), integer(1L),
                   USE.NAMES = FALSE),
    n_up = vapply(split_key, function(k) {
      sum(key == k & sig & object[["direction"]] == "up")
    }, integer(1L), USE.NAMES = FALSE),
    n_down = vapply(split_key, function(k) {
      sum(key == k & sig & object[["direction"]] == "down")
    }, integer(1L), USE.NAMES = FALSE)
  )
}

#' Drop the `gs_result` class
#'
#' @param x A `gs_result`.
#' @param ... Passed to [tibble::as_tibble()].
#' @return A plain tibble.
#' @importFrom tibble as_tibble
#' @exportS3Method tibble::as_tibble
as_tibble.gs_result <- function(x, ...) {
  .as_plain_tibble(x)
}

#' Combine `gs_result` objects
#'
#' Pooling across contrasts or databases is row-binding, because `contrast` and
#' `database` are columns. Optional columns present in only some inputs are
#' filled with `NA`.
#'
#' @param ... `gs_result` objects (or data frames with the core columns).
#' @param deparse.level Ignored; present for `rbind()` compatibility.
#' @return A `gs_result`.
#' @export
rbind.gs_result <- function(..., deparse.level = 1) {
  parts <- list(...)
  parts <- parts[!vapply(parts, is.null, logical(1L))]
  if (!length(parts)) {
    return(NULL)
  }
  flat <- lapply(parts, function(p) .as_plain_tibble(p))
  gs_result(dplyr::bind_rows(flat))
}

#' Subset a `gs_result`
#'
#' Downgrades to a plain tibble when the subset no longer satisfies the
#' contract -- a dropped core column, a wrong type, an out-of-vocabulary
#' `direction` -- so a broken object never survives a `[` or a `select()`.
#'
#' @param x A `gs_result`.
#' @param ... Row/column subscripts, passed to the tibble method.
#' @return A `gs_result` if the core columns survived, else a tibble (or a
#'   vector, if the tibble method returned one).
#' @export
`[.gs_result` <- function(x, ...) {
  out <- NextMethod()
  if (!is.data.frame(out)) {
    return(out)
  }
  .gs_reclass_result(out)
}

#' Re-attach the `gs_result` class only if the contract still holds
#'
#' `[` and every dplyr verb used to re-attach the class after checking column
#' *presence* alone, so `validate_gs_result()` was only ever reached from
#' `gs_result()` itself. `mutate(res, padj = "oops")` or a restored `"Up"`/
#' `"Down"` capitalisation therefore produced an object that still claimed to be
#' a `gs_result`, and `gs_filter(direction = "up")` then returned zero rows
#' silently. Failing the contract downgrades to a plain tibble -- with a
#' warning, because a silent downgrade is the same class of surprise.
#'
#' @param x A data frame produced by a subset or a dplyr verb.
#' @return A `gs_result` if valid, else a plain tibble.
#' @keywords internal
.gs_reclass_result <- function(x) {
  if (!all(.gs_core_cols %in% names(x))) {
    return(.as_plain_tibble(x))
  }
  out <- new_gs_result(x)
  ok <- tryCatch(
    {
      validate_gs_result(out)
      TRUE
    },
    error = function(e) {
      warning("Result no longer satisfies the gs_result contract (",
              conditionMessage(e), "); returning a plain tibble.",
              call. = FALSE)
      FALSE
    }
  )
  if (ok) out else .as_plain_tibble(x)
}

#' Keep `gs_result` through dplyr verbs
#'
#' Downgrades to a plain tibble (with a warning) when a verb leaves the object
#' outside the contract, so a broken object never survives a `select()` or a
#' `mutate()`.
#'
#' @param data The reconstructed data frame.
#' @param template The original `gs_result`.
#' @return A `gs_result` if the core columns survived, else a tibble.
#' @importFrom dplyr dplyr_reconstruct
#' @exportS3Method dplyr::dplyr_reconstruct
dplyr_reconstruct.gs_result <- function(data, template) {
  .gs_reclass_result(data)
}
