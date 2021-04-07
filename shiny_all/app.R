# scp app.R  gjb-shiny-x@rstudio.compbio.dundee.ac.uk:/homes/gjb-shiny-x/shiny/MitoNUb

VERSION <- "1.0"
DATE <- "2021-04-07"

lib <- "/cluster/gjb_lab/mgierlinski/R_shiny/library/3.6"
if(dir.exists(lib)) .libPaths(lib)

library(shiny)
library(DT)
library(tidyverse)
source("func.R")

options = list(language = list(emptyTable = 'No data'))

dat <- read_rds("data.rds")
all_genes <- dat$kgg$gene_id %>% unique()

###############################################################################

ui <- shinyUI(fluidPage(
  
  fluidRow(div(
    column(width=6, titlePanel("MitoNUb: mitochondrial ubiquitin landscape in neurons")),
    column(width=6, tags$a(img(src = "full_logo.png", height="60px"), href="https://www.ppu.mrc.ac.uk/"))
  )),
  
  div(paste0("Version: ", VERSION, ", Last updated: ", DATE), style="font-size:8pt"),
  br(),
  
  div(p("This app allows for quick selection of proteins from the diGly analysis of neurons reported in ", a(href = "https://www.biorxiv.org/content/10.1101/2021.04.01.438131v1", "Antico et al. (2021)", .noWS = "outside"), ". When Show->Proteins is selected, site position, logFC and FDR come from the peptide with the largest absolute fold change. For each selection of proteins the two tables at the bottom show GO-term and Reactome pathway enrichment. ", em("tot"), " is the total number of proteins with this term/pathway, ", em("sel"),  " - number in selection, ", em("expect"), " - expected count in selection based on random distribution, ", em("enrich"), " - enrichment over random background (observed / expected)."), style="color:grey"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("show", "Show", choices=c("UB sites", "Proteins"), inline = TRUE),
      checkboxGroupInput("checks", "Select", 
        choices = c(
          "In MitoCarta" = "in_mito",
          "In total proteome" = "in_total",
          "Significant DE" = "sig"
        ), selected=c("in_mito")
      ),
      checkboxGroupInput("ubi", "Select in UbiHub",
        choices = c(
          "E2" = "e2",
          "simple E3" = "e3_simple",
          "complex E3" = "e3_complex",
          "UBIs USP" = "dub_usp",
          "UBIs non-USP" = "dub_nonusp"
        ), selected = NULL
      ),
      sliderInput("up_fc", "Upregulated FC limit", min=0, max=5, value=0, step=0.01),
      sliderInput("down_fc", "Downregulated FC limit", min=0, max=5, value=0, step=0.01),
      radioButtons("go_selection", "GO database", choices=c("Full", "Slim"), inline=TRUE),
      width = 3
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
        tabPanel("Gene selection", DT::dataTableOutput("gene_table", width="98%")),
        tabPanel("GO enrichment", DT::dataTableOutput("go_enrichment", width="98%")),
        tabPanel("Reactome enrichment", DT::dataTableOutput("reactome_enrichment", width="98%"))
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
      ubi = input$ubi,
      up_fc = input$up_fc,
      down_fc = input$down_fc,
      go = input$go_selection
    )
  })
  
  getData <- function() {
    vals <- inputValues()
    checks_str <- ifelse(length(vals$checks) > 0, paste(vals$checks, collapse=" & "), "TRUE")
    ubi_str <-  ifelse(length(vals$ubi) > 0, paste(vals$ubi, collapse=" | "), "TRUE")
    filter_str <- paste("(", checks_str, ") & (", ubi_str, ")")
    filter_expr <- rlang::parse_expr(filter_str)
    
    tab <- dat$kgg %>% 
      filter(!!filter_expr & (log_fc >= vals$up_fc | log_fc <= -vals$down_fc))
    if(nrow(tab) == 0) return(NULL)
    
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
    tab <- getData()
    if(is.null(tab)) {
      shiny::showNotification("Selection resulted in no data", type = "error")
      return(NULL)
    }
    tab <- tab %>% select(-gene_id) %>% 
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
      DT::datatable(class = datatable_class) %>% 
      DT::formatStyle(colnames(tab), fontSize = '80%') %>% 
      DT::formatSignif(c("log2 FC", "FDR"), digits=2)
  })
  
  enrichmentTable <- function(terms, cap) {
    tab <- getData()
    if(is.null(tab)) {
      shiny::showNotification("Selection resulted in no data", type = "error")
      return(NULL)
    }
    sel <- tab$gene_id %>% unique()
    n <- length(sel)
    fe <- NULL
    if(n > 0) {
      fe <- functionalEnrichment(all_genes, sel, terms, dat$gene2name)
      if(is.null(fe)) {
        shiny::showNotification("Selection resulted in no enrichment data", type = "error")
        return(NULL)
      }
      fe <- fe %>%   
        select(-P) %>% 
        rename(
          "Term ID" = term_id,
          "Term name" = term_name,
          "Genes" = ids
        )
    }
    fe %>% 
      DT::datatable(class = datatable_class) %>% 
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

