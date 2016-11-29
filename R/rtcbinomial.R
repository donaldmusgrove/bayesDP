library(shiny)
library(shinyFiles)
library(shinythemes)

library(survival)
library(ggplot2)
library(MCMCpack)
library(knitr)
library(rmarkdown)

.runApp <- function(app, ...) {
  ## selectively use the RStudio viewer pane (if available)
  viewer <- getOption("viewer")
  if (!is.null(viewer) && is.function(viewer)) {
    runApp(app, launch.browser = viewer, ...)
  } else {
    runApp(app, ...)
  }
}


################################################################################
# UI
################################################################################

options(shiny.launch.browser = .rs.invokeShinyWindowViewer)

ui <- shinyUI(pageWithSidebar(
