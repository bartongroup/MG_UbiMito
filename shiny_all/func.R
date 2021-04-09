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
# term_data - a list with two elements:
#    term2gene: list term_id -> [gene_id1, gene_id2, ...]
#    term_info: data frame with columns term and other columns with name, description and so on.
# gene2name: named vector to change genes ids into gene names
# These additional columns will be inluded in the output.

functionalEnrichment <- function(genes_all, genes_sel, term_data, gene2name = NULL,
                                 min_count=3, sig_limit=0.05) {
  
  term2gene <- term_data$term2gene
  term_info <- term_data$terms
  terms <- names(term2gene)
  
  # number of selected genes
  Nsel <- length(genes_sel)
  # size of the universe
  Nuni <- length(genes_all)
  
  if(Nsel == Nuni) return(NULL)
  
  # empty line for missing terms
  na_term <- term_info %>% slice(1) %>% mutate_all(~NA)
  
  res <- map_dfr(terms, function(term) {
    # all genes with the term 
    tgenes <- intersect(term2gene[[term]], genes_all)
    # genes from selection with the term
    tgenes_sel <- intersect(tgenes, genes_sel)
    
    nuni <- length(tgenes)
    nsel <- length(tgenes_sel)
    
    g <- NULL
    if(nsel >= min_count) {
      info <- term_info[term_info$term_id == term, ]
      # returns NAs if no term found
      if(nrow(info) == 0) info <- na_term %>% mutate(term_id = term)
      expected <- nuni * Nsel / Nuni
      fish <- matrix(c(nsel, nuni - nsel, Nsel - nsel, Nuni + nsel - Nsel - nuni), nrow=2)
      ft <- fisher.test(fish, alternative = "greater")
      p <- as.numeric(ft$p.value)
    
      if(!is.null(gene2name)) tgenes_sel <- sort(unname(gene2name[tgenes_sel]))
    
      g <- bind_cols(
        info,
        tibble(
          tot = nuni,
          sel = nsel,
          expect = expected,
          enrich = nsel / expected,
          ids = paste(tgenes_sel, collapse=" "),
          P = p
        )
      )
    }
    g
  })
  
  if(nrow(res) == 0) return(NULL)
  
  res %>% 
    mutate(P = p.adjust(P, method="BH")) %>% 
    filter(P <= sig_limit) %>% 
    arrange(desc(enrich)) %>% 
    mutate(enrich = round(enrich, 1), expect = round(expect, 2))
}


okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")


pl_volcano <- function(res, sel, point.size=0.7) {
  lmt <- max(res %>% filter(fdr<0.05) %>% pull(p_value))
  ggplot(res, aes(x=log_fc, y=-log10(p_value))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    geom_point(colour = "grey70", size=point.size) +
    geom_point(data = res %>% filter(id %in% sel), colour = "black", size=point.size) +
    theme(legend.position = "none") +
    geom_hline(yintercept = -log10(lmt), colour="red", linetype="dashed") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x=expression(log[2]~FC), y=expression(-log[10]~P))
}
