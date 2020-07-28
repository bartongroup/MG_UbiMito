# scp app.R  gjb-shiny-x@rstudio.compbio.dundee.ac.uk:/homes/gjb-shiny-x/shiny/private/marek_ubimito_all

lib <- "/cluster/gjb_lab/mgierlinski/R_shiny/library/3.6"
if(dir.exists(lib)) .libPaths(lib)

library(shiny)
library(tidyverse)
source("func.R")

dat <- read_rds("data.rds")

###############################################################################

ui <- shinyUI(fluidPage(
  titlePanel("Mitochondrial ubiquitin landscape in neurons"),
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("checks", "Select", choices=c("In MitoCarta" = "in_mito", "In total proteome" = "in_total", "Significant DE" = "sig"), selected=c("in_mito")),
      sliderInput("up_fc", "Upregulated FC limit", min=0, max=5, value=0, step=0.01),
      sliderInput("down_fc", "Downregulated FC limit", min=0, max=5, value=0, step=0.01)
    ),
    mainPanel(
      fluidRow(
        p("Gene list"),
        DT::dataTableOutput("gene_table", width="98%")
        #div(style = 'height: 150px; overflow-y: scroll', tableOutput("gene_table"))
        # br(),
        # p("GO enrichment"),
        # div(style = 'height: 150px; overflow-y: scroll', tableOutput("GOEnrichment")),
        # br(),
        # p("Pathway enrichment"),
        # div(style = 'height: 250px; overflow-y: scroll', tableOutput("ReactomeEnrichment"))
      )
    )
  )
))

###############################################################################

server <- shinyServer(function(input, output, session) {
  
  inputValues <- reactive({
    list(
      checks = input$checks,
      up_fc = input$up_fc,
      down_fc = input$down_fc
    )
  })
  
  getData <- function() {
    vals <- inputValues()
    filter_str <- ifelse(length(vals$checks) > 0, paste(vals$checks, collapse=" & "), "TRUE")
    filter_expr <- rlang::parse_expr(filter_str)
    tab <- dat$kgg %>% 
      filter(!!filter_expr & (log_fc >= vals$up_fc | log_fc <= -vals$down_fc))
  }
  
  output$gene_table <- DT::renderDataTable({
    getData() %>% 
      select(gene_name, description, site_position, log_fc, fdr) %>% 
      DT::datatable(class = 'cell-border strip hover') %>% 
      DT::formatStyle(1:5, fontSize = '80%') %>% 
      DT::formatSignif(c("log_fc", "fdr"), digits=2)
  })
  
})

# Run the application
shinyApp(ui = ui, server = server)

