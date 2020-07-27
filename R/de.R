cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


limma_de <- function(d) {
  tab <- d %>% 
    select(id, ends_with("Ratio")) %>% 
    column_to_rownames("id") %>% 
    as.matrix()
  fnorm <- apply(tab, 2, median)
  tab <- t(t(tab) / fnorm)
  tab <- log2(tab)
  
  meta <- tibble(colname = colnames(tab)) %>% 
    mutate(sample = str_remove(colname, "_.atio")) %>% 
    mutate(group = str_remove(sample, "\\d") %>% as_factor() %>% fct_relevel("UT"))

  design_mat <- model.matrix(~group, data=meta)

  fit <- tab %>% 
    lmFit(design_mat) %>% 
    eBayes()
  
  coef <- colnames(fit$design)[2]
  res <- topTable(fit, coef=coef, adjust="BH", sort.by="none", number=1e9) %>% 
    as_tibble(rownames = "id") %>% 
    select(id, log_fc=logFC, p_value=P.Value, fdr=adj.P.Val, ave_ratio=AveExpr)

  d %>% 
    left_join(res, by="id")
}

diff_expr <- function(prot) {
  prot$kgg <- limma_de(prot$kgg)
  prot$total <- limma_de(prot$total)
  prot
}


plotVolcano <- function(res, fc="logFC", p="PValue", fdr="FDR", group="contrast",
                        alpha=0.05, point.size=0.5, point.alpha=0.5, sel=NULL) {
  r <- res %>%
    mutate(
      x = !!as.name(fc),
      y = !!as.name(p),
      fdr = !!as.name(fdr),
      mrk = 1L
    )
  lmt <- max(r %>% filter(fdr<0.05) %>% pull(y))
  if(!is.null(sel)) r <- r %>% mutate(mrk = if_else(id %in% sel, 2L, mrk))
  ggplot(r, aes(x=x, y=-log10(y), colour=as.factor(mrk))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    geom_point(size=point.size, alpha=point.alpha) +
    scale_colour_manual(values=c("grey70", "black")) +
    #facet_grid(as.formula(glue::glue(". ~ {group}"))) +
    theme(legend.position = "none") +
    geom_hline(yintercept = -log10(lmt), colour="red", linetype="dashed") +
    labs(x=fc, y=glue::glue("-log10({p})"))
}


plot_protein_sites <- function(d, gene, cex=3) {
  nd <- d %>% 
    pivot_longer(cols=UT1_Ratio:AO5_Ratio, names_to="colname") %>%
    select(id, gene_name, colname, value, site_id) %>% 
    mutate(sample = str_remove(colname, "_.atio")) %>% 
    mutate(group = str_remove(sample, "\\d") %>% as_factor() %>% fct_relevel("UT")) %>% 
    group_by(sample) %>% 
    mutate(med = median(value, na.rm=TRUE)) %>% 
    mutate(value_norm = value / med)
    
  nd %>%
    filter(gene_name == gene) %>% 
  ggplot(aes(x=group, y=log2(value))) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_hline(yintercept = 0, colour="grey80") +
    geom_beeswarm(cex=cex) +
    facet_wrap(~site_id, nrow=1) 
}
