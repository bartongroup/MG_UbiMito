geneTable <- function(tab) {
 tab %>% 
    select(gene_name, description, site_position, log_fc, fdr) %>% 
    DT::datatable(class = 'cell-border strip hover') %>% 
    DT::formatStyle(0, cursor = 'pointer')
}


# Enrichment function
#
# Performs hypergeometric test.
# Input: genes_sel and genes_all are vectors with gene IDs.
# gene2term: data frame with columns id and term
# gene2name: named vector to change genes ids into gene names
# term_info: data frame with columns term and other columns with name, description and so on.
# These additional columns will be inluded in the output.

functionalEnrichment <- function(genes_all, genes_sel, term_data, gene2name = NULL,
                                 min_count=3, sig_limit=0.05) {
  
  gene2term <- term_data$gene2term
  term_info <- term_data$terms
  
  # select only terms represented in our gene set
  gene2term <- gene2term %>% filter(gene_id %in% genes_all)
  
  # all terms present in the selection
  terms <- gene2term %>% 
    filter(gene_id %in% genes_sel) %>% 
    pull(term_id) %>% 
    unique()
  
  # number of selected genes
  Nsel <- length(genes_sel)
  # size of the universe
  Nuni <- length(genes_all)
  
  # empty line for missing terms
  na_term <- term_info %>% slice(1) %>% mutate_all(~NA)
  
  res <- map_dfr(terms, function(term) {
    info <- term_info %>% filter(term_id == term)
    # returns NAs if no term found
    if(nrow(info) == 0) info <- na_term %>% mutate(term_id = term)
    
    # all genes with the term
    tgenes <- gene2term %>% filter(term_id == term) %>% pull(gene_id)
    # genes from selection with the term
    tgenes_sel <- intersect(tgenes, genes_sel)
    
    nuni <- length(tgenes)
    nsel <- length(tgenes_sel)
    
    expected <- nuni * Nsel / Nuni
    fish <- matrix(c(nsel, nuni - nsel, Nsel - nsel, Nuni + nsel - Nsel - nuni), nrow=2)
    ft <- fisher.test(fish, alternative = "greater")
    p <- as.numeric(ft$p.value)
    
    if(!is.null(gene2name)) tgenes_sel <- gene2name[tgenes_sel] %>% unname()
    
    bind_cols(
      info,
      tibble(
        tot = nuni,
        sel = nsel,
        expect = expected,
        enrich = nsel / expected,
        ids = paste(tgenes_sel, collapse=","),
        P = p
      )
    )
  }) %>% 
    mutate(P = p.adjust(P, method="BH")) %>% 
    filter(sel >= min_count & P <= sig_limit) %>% 
    arrange(desc(enrich)) %>% 
    mutate(enrich = round(enrich, 1), expect = round(expect, 2))
  
  res
}
