#' bulkiRNA: bulk RNA-seq differential expression and gene-set analysis
#'
#' Four layers, one direction of flow:
#'
#' \preformatted{
#'   gsdb_*        gs_test() / gs_score()     gs_result / gs_matrix     gs_plot_*
#'   providers --> compute               -->  data objects         -->  renderers
#' }
#'
#' Compute functions return objects and never plot or write; `gs_write()` /
#' `gs_read()` / `gs_save()` are the only disk access; `gs_plot_*` take an
#' object and return a ggplot.
#'
#' See [gs_result-class] and [gs_matrix-class] for the shared data contracts.
#'
#' @section Reference data:
#' Bundled processed gene-set databases live in `inst/extdata` and resolve via
#' `system.file("extdata", ..., package = "bulkiRNA")`. Raw source files are in
#' the git checkout only and are not shipped.
#'
#' @import ggplot2
#' @keywords internal
"_PACKAGE"
