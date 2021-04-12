# scp app.R  gjb-shiny-x@rstudio.compbio.dundee.ac.uk:/homes/gjb-shiny-x/shiny/MitoNUb

VERSION <- "1.11"
DATE <- "2021-04-12"

lib <- "/cluster/gjb_lab/mgierlinski/R_shiny/library/3.6"
if(dir.exists(lib)) .libPaths(lib)

library(shiny)
library(shinyBS)
library(DT)
library(tidyverse)
source("func.R")

select <- dplyr::select

options = list(language = list(emptyTable = 'No data'), dplyr.summarise.inform = FALSE)

dat <- read_rds("data.rds")

###############################################################################

ui <- shinyUI(fluidPage(
  
  fluidRow(div(
    column(width=5, titlePanel("MitoNUb: mitochondrial ubiquitin landscape in neurons")),
    column(width=7, tags$a(img(src = "full_logo.png", height="60px"), href="https://www.ppu.mrc.ac.uk/"))
  )),
  
  div(paste0("Version: ", VERSION, ", Last updated: ", DATE), style="font-size:8pt"),
  br(),
  
  div(p("This app allows for quick selection of ubiquitilation sites and proteins from the diGly analysis of neurons stimulated with mitochondrial depolarisation, reported in ", a(href = "https://www.biorxiv.org/content/10.1101/2021.04.01.438131v1", "Antico et al. (2021)", .noWS = "outside"), ". For each selection of sites/proteins GO-term and Reactome pathway enrichment are calculated. ", em("tot"), " is the total number of proteins with this term/pathway, ", em("sel"),  " - number in selection, ", em("expect"), " - expected count in selection based on random distribution, ", em("enrich"), " - enrichment over random background (observed / expected)."), style="color:grey"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("select", "Select", choices=c("UB sites", "Total proteome"), inline = TRUE),
      checkboxGroupInput("checks", "Select", 
        choices = c(
          "In MitoCarta" = "in_mito",
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
      sliderInput("down_fc", "Downregulated FC limit", min=-5, max=0, value=0, step=0.01),
      radioButtons("go_selection", "GO database", choices=c("Full", "Slim"), inline=TRUE),
      
      # tooltips
      bsTooltip("select", "Selection between diGlycin ubiquitilation capture and total proteome results", placement = "right", options = list(container = "body")),
      bsTooltip("checks", "Selection of proteins in present MitoCarta; selection of UB sites/proteins significantly differentially expressed (FDR < 0.05) between mitochondrial depolarisation and control", placement = "right", options = list(container = "body")),
      bsTooltip("ubi", "Selection of categories from UbiHub", placement = "right", options = list(container = "body")),
      bsTooltip("up_fc", "Lower limit on positive log2 fold change between mitochondrial depolarisation and control", placement = "right", options = list(container = "body")),
      bsTooltip("down_fc", "Upper limit on negative log2 fold change between mitochondrial depolarisation and control", placement = "right", options = list(container = "body")),
      bsTooltip("go_selection", "Whether to use full GO-term database or a reduced GO-slim database", placement = "right", options = list(container = "body")),
      width = 3
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
        tabPanel("Gene selection", DT::dataTableOutput("gene_table", width="98%")),
        tabPanel("Volcano plot", plotOutput("volcano", width="600px", height="500px")),
        tabPanel("GO enrichment", DT::dataTableOutput("go_enrichment", width="98%")),
        tabPanel("Reactome enrichment", DT::dataTableOutput("reactome_enrichment", width="98%"))
      )
    )
  )
))

###############################################################################

server <- shinyServer(function(input, output, session) {
  
  input_values <- reactive({
    list(
      select = input$select,
      checks = input$checks,
      ubi = input$ubi,
      up_fc = input$up_fc,
      down_fc = -input$down_fc,
      go = input$go_selection
    )
  })
  
  get_data <- function(vals) {
    if(vals$select == "UB sites") {
      tab <- dat$kgg
    } else {
      tab <- dat$tot
    }
    tab
  }
  
  get_data_selection <- function(vals) {
    checks_str <- ifelse(length(vals$checks) > 0, paste(vals$checks, collapse=" & "), "TRUE")
    ubi_str <-  ifelse(length(vals$ubi) > 0, paste(vals$ubi, collapse=" | "), "TRUE")
    filter_str <- paste("(", checks_str, ") & (", ubi_str, ")")
    filter_expr <- rlang::parse_expr(filter_str)
    
    tab <- get_data(vals) %>% 
      filter(!!filter_expr & (log_fc >= vals$up_fc | log_fc <= -vals$down_fc))
    if(nrow(tab) == 0) return(NULL)
    tab
  }
  
  datatable_class <- "compact row-border hover"
  
  output$gene_table <- DT::renderDataTable({
    vals <- input_values()
    tab <- get_data_selection(vals)
    if(is.null(tab)) {
      shiny::showNotification("Selection resulted in no data", type = "error")
      return(NULL)
    }
    ts <- tab %>% 
      select(
        `Gene name` = gene_name,
        Description = description,
        `log2 FC` = log_fc,
        FDR = fdr,
        UbiHub = ubi,
        Subcompartment = sub
      )
    if(vals$select == "UB sites") ts <- add_column(ts, `Site position` = tab$site_position, .after="Gene name")
    ts %>% 
      DT::datatable(class = datatable_class) %>% 
      DT::formatStyle(colnames(ts), fontSize = '80%') %>% 
      DT::formatSignif(c("log2 FC", "FDR"), digits=2)
  })
  
  enrichmentTable <- function(terms, cap) {
    vals <- input_values()
    tab_all <- get_data(vals)
    tab_sel <- get_data_selection(vals)
    
    if(is.null(tab_sel)) {
      shiny::showNotification("Selection resulted in no data", type = "error")
      return(NULL)
    }
    
    all_genes <- tab_all$gene_id %>% unique()
    sel_genes <- tab_sel$gene_id %>% unique()
    
    n <- length(sel_genes)
    fe <- NULL
    if(n > 0) {
      fe <- functionalEnrichment(all_genes, sel_genes, terms, dat$gene2name)
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
    vals <- input_values()
    if(vals$go == "Full") {
      enrichmentTable(dat$bm_go, "GO term enrichment (full)")
    } else {
      enrichmentTable(dat$bm_go_slim, "GO term enrichment (slim)")
    }
  })
  
  output$reactome_enrichment <- DT::renderDataTable({
    enrichmentTable(dat$reactome, "Reactome pathway enrichment")
  })
  
  output$volcano <- renderPlot({
    vals <- input_values()
    d <- get_data_selection(vals)
    if(is.null(d)) {
      sel <- NULL
    } else {
      sel <- d$id %>% unique()
    }
    
    if(vals$select == "UB sites") {
      tab <- dat$kgg
    } else {
      tab <- dat$tot
    }
    pl_volcano(tab, sel)
  })
  
  
})

# Run the application
shinyApp(ui = ui, server = server)

