#' Ensure a directory exists
#'
#' Creates each directory in `path` if it is missing, including parents.
#' Idempotent. `path` is a **directory** path, never a file path -- if you have a
#' file path, pass `dirname(path)` (or use the internal [ensure_parent_dir()]).
#'
#' @param path Character vector of directory paths.
#' @return `path`, invisibly, so the call can be inlined. Errors if a
#'   directory could not be created (e.g. a read-only parent).
#' @examples
#' d <- file.path(tempdir(), "figures")
#' ensure_dir(d)
#' dir.exists(d)
#' @export
ensure_dir <- function(path) {
  if (is.null(path) || length(path) == 0L) {
    stop("`path` must be a non-empty character vector of directory paths.", call. = FALSE)
  }
  if (!is.character(path) || anyNA(path) || any(!nzchar(path))) {
    stop("`path` must contain non-missing, non-empty directory paths.", call. = FALSE)
  }
  for (p in path) {
    if (!dir.exists(p)) {
      dir.create(p, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(p)) {
        stop("Failed to create directory `", p, "`. ",
             "Check that the parent directory is writable (file.access(): ",
             file.access(dirname(p), mode = 2), ").", call. = FALSE)
      }
    }
  }
  invisible(path)
}

#' Ensure the parent directory of a file path exists
#'
#' File-path counterpart to [ensure_dir()]. Internal; kept separate so the two
#' semantics never live under one name again.
#'
#' @param path Character vector of **file** paths.
#' @return `path`, invisibly.
#' @keywords internal
ensure_parent_dir <- function(path) {
  ensure_dir(unique(dirname(path)))
  invisible(path)
}

#' Require a Suggests package, with a caller-supplied install hint
#'
#' Shared body for the "this needs an optional package" check. Two other
#' near-duplicates of this existed (`.de_require()` in `R/de-utils.R` and
#' `.gatom_require()` in `R/gatom.R`, the latter hardcoding a
#' `BiocManager::install("gatom")` special case) plus several bare
#' `requireNamespace()` calls scattered around `R/`. This is the union of
#' their behaviour: `install` is a caller-supplied hint string (as in
#' `.de_require()`), so call sites needing a `BiocManager::install(...)`
#' message (as `.gatom_require()` special-cased for `pkg == "gatom"`) just
#' pass that string instead of relying on a hardcoded branch.
#'
#' @param pkg Character(1) package name.
#' @param what Character(1) description of what needed it.
#' @param install Character(1) install-command hint shown in the error.
#'   Defaults to `install.packages("<pkg>")`.
#' @return `TRUE`, invisibly; errors if `pkg` is not installed.
#' @keywords internal
.require_pkg <- function(pkg, what,
                         install = sprintf('install.packages("%s")', pkg)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("%s requires the `%s` package. Install it with %s.",
                 what, pkg, install), call. = FALSE)
  }
  invisible(TRUE)
}

#' Default value for NULL
#'
#' @param x,y `x` unless it is `NULL`, in which case `y`.
#' @return `x` or `y`.
#' @name null-default
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x
