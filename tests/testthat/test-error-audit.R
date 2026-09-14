error_audit_r_dir <- function() {
  for (path in c("../../R", "R", testthat::test_path("..", "..", "R"))) {
    if (dir.exists(path) && length(list.files(path, "[.]R$"))) return(path)
  }
  NULL
}

error_audit_stop_calls <- function(r_dir) {
  out <- list()
  add_expression <- function(file, owner, expression) {
    ordinal <- 0L
    walk <- function(node) {
      if (is.call(node) && identical(node[[1L]], as.name("stop"))) {
        ordinal <<- ordinal + 1L
        key <- paste0(basename(file), ":", owner, "#", ordinal)
        out[[key]] <<- node
      }
      if (!is.call(node) && !is.pairlist(node) && !is.expression(node) &&
          !is.list(node)) return(invisible(NULL))
      children <- as.list(node)
      if (is.call(node)) children <- children[-1L]
      # A formal declared without a default, and a skipped index such as the
      # first slot of `x[, 1]`, are both the empty symbol. Binding one to a
      # variable gives that variable missing-argument semantics, so it must be
      # tested by index rather than after assignment. Nearly every function
      # here has such a formal, so the walker aborted before asserting
      # anything: this audit reported an error, never a result.
      for (i in seq_along(children)) {
        if (identical(children[[i]], quote(expr = ))) next
        walk(children[[i]])
      }
      invisible(NULL)
    }
    walk(expression)
  }

  for (file in sort(list.files(r_dir, "[.]R$", full.names = TRUE))) {
    parsed <- parse(file, keep.source = TRUE)
    for (expr in parsed) {
      is_function_assignment <- is.call(expr) &&
        identical(expr[[1L]], as.name("<-")) && is.call(expr[[3L]]) &&
        identical(expr[[3L]][[1L]], as.name("function"))
      if (is_function_assignment) {
        add_expression(file, as.character(expr[[2L]]), expr[[3L]])
      } else {
        add_expression(file, "<top-level>", expr)
      }
    }
  }
  out
}

error_audit_character_literals <- function(x) {
  if (is.character(x)) return(x)
  if (!is.call(x) && !is.pairlist(x) && !is.expression(x)) {
    return(character(0L))
  }
  unlist(lapply(as.list(x), error_audit_character_literals), use.names = FALSE)
}

error_audit_expect_allowlist <- function(observed, allowlist, info = NULL) {
  expected <- names(allowlist)
  if (is.null(expected)) expected <- character(0L)
  expect_identical(sort(as.character(observed)), sort(expected), info = info)
}

test_that("every stop call suppresses calls and names an argument or reason", {
  r_dir <- error_audit_r_dir()
  if (is.null(r_dir)) {
    skip(paste(
      "Source-tree-only error audit: R/ is unavailable when tests run",
      "against the installed package."
    ))
  }
  stops <- error_audit_stop_calls(r_dir)
  # expect_gt() takes no `info`, so the reason rides on expect_true().
  expect_true(length(stops) > 0L,
              info = "The source parser must find stop() calls.")

  no_argument_exceptions <- c(
    "coresh.R:.coresh_validate_object#3" = paste(
      "Validates a deserialized dataset payload; the dynamic context names",
      "the corrupt E1024 field, not a caller argument."
    ),
    "coresh.R:.coresh_validate_object#4" = paste(
      "Validates a deserialized dataset payload; the dynamic context names",
      "the corrupt rownames field, not a caller argument."
    ),
    "coresh.R:.coresh_validate_object#5" = paste(
      "Validates a deserialized dataset payload; the dynamic context names",
      "the corrupt totalVar field, not a caller argument."
    ),
    "coresh.R:.coresh_validate_object#6" = paste(
      "Validates scalar fields in a deserialized dataset payload and names",
      "the offending field dynamically."
    ),
    "coresh.R:coresh_match#4" = paste(
      "Rejects a malformed third-party GESECA result and reports the dataset",
      "identity and observed result shape."
    ),
    "gatom-download.R:gatom_download_refs#2" = paste(
      "This caught transfer failure becomes a warning naming the downloaded",
      "file; no caller argument is invalid."
    ),
    "gatom-download.R:gatom_download_refs#3" = paste(
      "This caught filesystem failure becomes a warning naming the target",
      "file; no caller argument is invalid."
    ),
    "gs-coregulation.R:.gs_geseca#1" = paste(
      "Rejects a malformed third-party GESECA return object, not a caller",
      "argument."
    ),
    "gs-coregulation.R:.gs_geseca#2" = paste(
      "Rejects missing fields in a third-party GESECA return object, not a",
      "caller argument."
    ),
    "gs-master.R:.gs_master_schema#3" = paste(
      "Detects a malformed schema file shipped inside bulkiRNA; reinstalling",
      "the package is the remedy, and no caller argument caused it."
    ),
    "gs-master.R:.gs_master_schema#4" = paste(
      "Detects invalid declarations in a schema shipped inside bulkiRNA; no",
      "caller argument caused the package-data corruption."
    ),
    "gs-master.R:.gs_master_coerce#1" = paste(
      "Private sentinel caught immediately by the same function to report an",
      "integer-coercion failure; it never reaches a caller."
    ),
    "gsdb-file.R:.gsdb_parse_gmt#1" = paste(
      "The private parser receives already-read lines; it reports the malformed",
      "GMT record shape because no public argument name exists at this layer."
    ),
    "gsdb-file.R:.gsdb_parse_gmx#1" = paste(
      "The private parser receives already-read lines; it reports the malformed",
      "GMX row count because no public argument name exists at this layer."
    )
  )
  expect_true(
    all(nzchar(no_argument_exceptions)) &&
      !any(grepl("[\r\n]", no_argument_exceptions)),
    info = "Every no-argument stop exception needs its own one-line reason."
  )

  missing_call_false <- names(stops)[!vapply(stops, function(call) {
    args <- as.list(call)[-1L]
    named <- names(args)
    index <- which(!is.na(named) & named == "call.")
    length(index) == 1L && identical(args[[index]], FALSE)
  }, logical(1L))]
  missing_argument <- names(stops)[!vapply(stops, function(call) {
    any(grepl("`", error_audit_character_literals(call), fixed = TRUE))
  }, logical(1L))]
  option_dependent_quotes <- names(stops)[vapply(stops, function(call) {
    "sQuote" %in% all.names(call, functions = TRUE)
  }, logical(1L))]

  error_audit_expect_allowlist(
    missing_call_false,
    NULL,
    info = "Every stop() call must set call. = FALSE."
  )
  error_audit_expect_allowlist(
    missing_argument,
    no_argument_exceptions,
    info = paste(
      "Every stop() must name an argument in backticks or have an exact,",
      "per-call reason explaining why no caller argument exists."
    )
  )
  error_audit_expect_allowlist(
    option_dependent_quotes,
    NULL,
    info = paste(
      "Error values must use a stable quoted representation; sQuote()",
      "changes with options(useFancyQuotes)."
    )
  )
  error_audit_expect_allowlist(
    character(0L),
    NULL,
    info = "The error exception comparison must stay shape-stable when empty."
  )
})
