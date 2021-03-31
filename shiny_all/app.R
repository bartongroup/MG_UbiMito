# scp app.R  gjb-shiny-x@rstudio.compbio.dundee.ac.uk:/homes/gjb-shiny-x/shiny/MitoNUb

lib <- "/cluster/gjb_lab/mgierlinski/R_shiny/library/3.6"
if(dir.exists(lib)) .libPaths(lib)

library(shiny)
library(DT)
library(tidyverse)
source("func.R")

dat <- read_rds("data.rds")
all_genes <- dat$kgg$gene_id %>% unique()

###############################################################################

ui <- shinyUI(fluidPage(
  titlePanel("MitoNUb: mitochondrial ubiquitin landscape in neurons"),
  p("This app allows for quick selection of proteins from the diGly experiment. When Show->Proteins is selected, site position, logFC and FDR come from the peptide with the largest absolute fold change. For each selection of proteins the two tables at the bottom show GO-term and Reactome pathway enrichment. ", em("tot"), " is the total number of proteins with this term/pathway, ", em("sel"),  " - number in selection, ", em("expect"), " - expected count in selection based on random distribution, ", em("enrich"), " - enrichment over random background (observed / expected)."),
  sidebarLayout(
    sidebarPanel(
      radioButtons("show", "Show", choices=c("UB sites", "Proteins"), inline = TRUE),
      checkboxGroupInput("checks", "Select", 
        choices=c("In MitoCarta" = "in_mito", "In total proteome" = "in_total", "Significant DE" = "sig"), selected=c("in_mito")
      ),
      sliderInput("up_fc", "Upregulated FC limit", min=0, max=5, value=0, step=0.01),
      sliderInput("down_fc", "Downregulated FC limit", min=0, max=5, value=0, step=0.01),
      radioButtons("go_selection", "GO database", choices=c("Full", "Slim"), inline=TRUE),
      width = 3
    ),
    mainPanel(
      fluidRow(
        DT::dataTableOutput("gene_table", width="98%"),
        DT::dataTableOutput("go_enrichment", width="98%"),
        DT::dataTableOutput("reactome_enrichment", width="98%")
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
      show = input$show,
      checks = input$checks,
      up_fc = input$up_fc,
      down_fc = input$down_fc,
      go = input$go_selection
    )
  })
  
  getData <- function() {
    vals <- inputValues()
    filter_str <- ifelse(length(vals$checks) > 0, paste(vals$checks, collapse=" & "), "TRUE")
    filter_expr <- rlang::parse_expr(filter_str)
    tab <- dat$kgg %>% 
      filter(!!filter_expr & (log_fc >= vals$up_fc | log_fc <= -vals$down_fc))
    if(vals$show == "UB sites") {
      tab <- tab %>%
        select(gene_name, description, site_position, log_fc, fdr, gene_id, ubi, sub)
    } else {
      tab <- tab %>%
        group_by(gene_name, description) %>%
        summarise(imax = which.max(abs(log_fc)), site_position = site_position[imax], log_fc = log_fc[imax], fdr = fdr[imax], gene_id = gene_id[imax], ubi = ubi[imax], sub = sub[imax]) %>% 
        select(-imax)
    }
    tab
  }
  
  datatable_class <- "compact row-border hover"
  
  output$gene_table <- DT::renderDataTable({
    tab <- getData() %>% select(-gene_id) %>% 
      rename(
        "Gene name" = gene_name,
        "Description" = description,
        "Site position" = site_position,
        "log2 FC" = log_fc,
        "FDR" = fdr,
        "UbiHub" = ubi,
        "Subcompartment" = sub
      )
    tab %>% 
      DT::datatable(class = datatable_class, caption="Selected genes") %>% 
      DT::formatStyle(colnames(tab), fontSize = '80%') %>% 
      DT::formatSignif(c("log2 FC", "FDR"), digits=2)
  })
  
  enrichmentTable <- function(terms, cap) {
    tab <- getData()
    sel <- tab$gene_id %>% unique()
    n <- length(sel)
    fe <- NULL
    if(n > 0) {
      fe <- functionalEnrichment(all_genes, sel, terms, dat$gene2name)
      if(is.null(fe)) return(NULL)
      fe %>%   
        select(-P) %>% 
        rename(
          "Term ID" = term_id,
          "Term name" = term_name,
          "Genes" = ids
        )
    }
    fe %>% 
      DT::datatable(class = datatable_class, caption=cap) %>% 
      DT::formatStyle(colnames(fe), fontSize = '80%')
  }
  
  output$go_enrichment <- DT::renderDataTable({
    vals <- inputValues()
    if(vals$go == "Full") {
      enrichmentTable(dat$bm_go, "GO term enrichment (full)")
    } else {
      enrichmentTable(dat$bm_go_slim, "GO term enrichment (slim)")
    }
  })
  
  output$reactome_enrichment <- DT::renderDataTable({
    enrichmentTable(dat$reactome, "Reactome pathway enrichment")
  })
  
  
})

# Run the application
shinyApp(ui = ui, server = server)

