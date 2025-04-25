#' module1 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_module1_ui <- function(id) {
  ns <- NS(id)
  tagList(

  )
}

#' module1 Server Functions
#'
#' @noRd
mod_module1_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

  })
}

## To be copied in the UI
# mod_module1_ui("module1_1")

## To be copied in the server
# mod_module1_server("module1_1")
