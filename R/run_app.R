#' Launch the phyloseq Explorer app
#'
#' Opens the interactive Shiny application for exploring any
#' \code{phyloseq} object: composition, alpha/beta diversity with
#' statistics, RDA against environmental variables, taxonomy tables and
#' publication-ready exports. All processing happens locally.
#'
#' @param launch.browser Open the app in the browser? Default \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#' @return Called for its side effect; returns the \code{\link[shiny]{runApp}} result.
#' @examples
#' \dontrun{
#' phyloseqExplorer::run_app()
#' }
#' @export
run_app <- function(launch.browser = TRUE, ...) {
  options(shiny.maxRequestSize = 500 * 1024^2)
  suppressWarnings(try(Sys.setlocale("LC_ALL", "en_US.UTF-8"), silent = TRUE))
  app <- shiny::shinyApp(ui = app_ui, server = app_server)
  shiny::runApp(app, launch.browser = launch.browser, ...)
}
