pane_runApp <- function(app, ...) {
  ## selectively use the RStudio viewer pane (if available)
  viewer <- getOption("viewer")
  if (!is.null(viewer) && is.function(viewer)) {
    runApp(app, launch.browser = viewer, ...)
  } else {
    runApp(app, ...)
  }
}
