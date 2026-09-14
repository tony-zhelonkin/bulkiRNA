# Package-local record of reference data resolved during this session and
# whether each source has emitted its invalid-layout provenance warning.
.ref_resolution_state <- new.env(parent = emptyenv())

#' Resolve a path in the shared reference cache
#'
#' An explicit `path` takes precedence over the shared cache. Otherwise the
#' path is resolved below `$REFCACHE_ROOT/<source>/current`, and the snapshot
#' selected by `current` is recorded for session provenance. If `current` is
#' not a symlink or resolves outside its source directory, the resolved
#' basename is retained but a warning is emitted at most once per source per
#' session. With `options(warn = 2)`, an outside-source layout aborts instead
#' of returning the path as it did before this warning was introduced.
#'
#' @param source Character scalar naming one refcache source. It must be a
#'   single path segment.
#' @param ... Character scalar path components below the source's `current`
#'   snapshot. Ignored when `path` is supplied.
#' @param path Optional caller-supplied path. Used as given.
#' @param root Optional refcache root. `NULL` reads `REFCACHE_ROOT`.
#' @param must_exist Logical scalar; whether the returned path must exist.
#' @return A character scalar with `snapshot`, `source`, `root`, and
#'   `caller_supplied` attributes. `snapshot` and `root` are `NA_character_`
#'   for a caller-supplied path.
#' @keywords internal
.ref_path <- function(source, ..., path = NULL, root = NULL,
                      must_exist = TRUE) {
  if (!is.character(source) || length(source) != 1L || is.na(source) ||
      !nzchar(source) || source == "." ||
      grepl("[/\\\\]|\\.\\.", source)) {
    stop("`source` must be a single path segment -- no `/`, `\\`, `..` or ",
         "`.`; got ",
         encodeString(paste0(source, collapse = ", "), quote = "\""), ".",
         call. = FALSE)
  }
  source <- unname(source)
  if (!is.logical(must_exist) || length(must_exist) != 1L ||
      is.na(must_exist)) {
    stop("`must_exist` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  previous_resolution <- if (exists(
    source,
    envir = .ref_resolution_state,
    inherits = FALSE
  )) {
    get(source, envir = .ref_resolution_state, inherits = FALSE)
  } else {
    NULL
  }
  warned_layout <- isTRUE(previous_resolution$warned_layout)
  layout_warning <- NULL

  caller_supplied <- !is.null(path)
  if (caller_supplied) {
    if (!is.character(path) || length(path) != 1L || is.na(path) ||
        !nzchar(path)) {
      stop("`path` must be a non-empty character scalar or `NULL`.",
           call. = FALSE)
    }
    if (must_exist && !file.exists(path)) {
      stop("Caller-supplied reference `path` does not exist: ", path, ".",
           call. = FALSE)
    }
    resolved <- unname(as.character(path))
    snapshot <- NA_character_
    resolved_root <- NA_character_
  } else {
    components <- unname(list(...))
    valid_components <- vapply(
      components,
      function(x) is.character(x) && length(x) == 1L && !is.na(x),
      logical(1L)
    )
    if (length(valid_components) && !all(valid_components)) {
      stop("Each path component in `...` must be a character scalar.",
           call. = FALSE)
    }

    if (is.null(root)) {
      root <- Sys.getenv("REFCACHE_ROOT", unset = "")
    }
    if (!is.character(root) || length(root) != 1L || is.na(root)) {
      stop("`root` must be a character scalar or `NULL`.", call. = FALSE)
    }
    if (!nzchar(root)) {
      stop(
        "Reference data \"", source, "\" is not available. ",
        "`REFCACHE_ROOT` is unset. Bind the shared refcache into the ",
        "container and set `REFCACHE_ROOT=/refcache`, or pass `path=` to ",
        "point at a local copy.",
        call. = FALSE
      )
    }

    source_dir <- file.path(root, source)
    if (!dir.exists(source_dir)) {
      stop(
        "Reference data \"", source, "\" is not available. The source ",
        "directory is missing under `REFCACHE_ROOT`: ", source_dir, ". ",
        "Bind the shared refcache into the container and set ",
        "`REFCACHE_ROOT=/refcache`, or pass `path=` to point at a local copy.",
        call. = FALSE
      )
    }

    current <- file.path(source_dir, "current")
    if (!dir.exists(current)) {
      stop(
        "Reference data \"", source, "\" is not available. `current` is ",
        "missing or dangling under the source directory: ", current, ". ",
        "The refcache refresh may have failed. Repair the refcache snapshot, ",
        "or pass `path=` to point at a local copy.",
        call. = FALSE
      )
    }

    current_target <- Sys.readlink(current)
    resolved_target <- normalizePath(
      current,
      winslash = "/",
      mustWork = TRUE
    )
    resolved_source_dir <- normalizePath(
      source_dir,
      winslash = "/",
      mustWork = TRUE
    )
    if (!nzchar(current_target)) {
      layout_warning <- paste0(
        "`current` for reference source ",
        encodeString(source, quote = "\""),
        " is not a symlink, so its snapshot cannot be identified: ",
        current, "."
      )
    } else if (!startsWith(
      resolved_target,
      paste0(resolved_source_dir, "/")
    )) {
      layout_warning <- paste0(
        "Resolved `current` for reference source ",
        encodeString(source, quote = "\""),
        " outside its source directory: ", resolved_target, ". The recorded ",
        "snapshot ", encodeString(basename(resolved_target), quote = "\""),
        " may not identify a refcache snapshot."
      )
    }
    snapshot <- basename(resolved_target)
    resolved_root <- unname(root)
    resolved <- if (length(components)) {
      unname(do.call(file.path, c(list(current), components)))
    } else {
      unname(current)
    }

    if (length(resolved) != 1L) {
      stop("Path components in `...` must resolve to one path.", call. = FALSE)
    }
    if (must_exist && !file.exists(resolved)) {
      stop("`source` reference data path for ",
           encodeString(source, quote = "\""), " does not exist: ", resolved,
           ". Install or mount that snapshot under `REFCACHE_ROOT`.",
           call. = FALSE)
    }
  }

  if (!is.null(layout_warning) && !warned_layout) {
    warned_layout <- TRUE
    assign(
      source,
      list(
        source = source,
        snapshot = snapshot,
        path = resolved,
        warned_layout = warned_layout
      ),
      envir = .ref_resolution_state
    )
    warning(layout_warning, call. = FALSE)
  }

  out <- structure(
    resolved,
    snapshot = snapshot,
    source = source,
    root = resolved_root,
    caller_supplied = caller_supplied
  )
  assign(
    source,
    list(
      source = source,
      snapshot = snapshot,
      path = resolved,
      warned_layout = warned_layout
    ),
    envir = .ref_resolution_state
  )
  out
}

#' Return reference data resolved during this session
#'
#' @return A data frame with character columns `source`, `snapshot`, and
#'   `path`, containing at most one row per source.
#' @keywords internal
.ref_resolutions <- function() {
  keys <- sort(ls(envir = .ref_resolution_state, all.names = TRUE))
  if (!length(keys)) {
    return(data.frame(
      source = character(0),
      snapshot = character(0),
      path = character(0),
      stringsAsFactors = FALSE
    ))
  }

  records <- lapply(
    keys,
    get,
    envir = .ref_resolution_state,
    inherits = FALSE
  )
  data.frame(
    source = unname(vapply(records, `[[`, character(1L), "source")),
    snapshot = unname(vapply(records, `[[`, character(1L), "snapshot")),
    path = unname(vapply(records, `[[`, character(1L), "path")),
    stringsAsFactors = FALSE
  )
}

#' Clear the session reference-data record
#'
#' @return `NULL`, invisibly.
#' @keywords internal
.clear_ref_resolutions <- function() {
  keys <- ls(envir = .ref_resolution_state, all.names = TRUE)
  if (length(keys)) {
    rm(list = keys, envir = .ref_resolution_state)
  }
  invisible(NULL)
}
