find_basal <- function(prot) {
  quant_prot <- prot$total$uniprot %>% unique()
  prot$kgg_basal <- prot$kgg %>% filter(uniprot %in% quant_prot)
  prot
}


merge_prot_mito <- function(dat, mito, ubihub, columns=NULL) {
  dat %>% 
    left_join(mito$carta, by=c("gene_name" = "Symbol")) %>% 
    left_join(ubihub, by="gene_name") %>% 
    select(
      id,
      uniprot,
      gene_name,
      description,
      p_value,
      fdr,
      log_fc,
      MCARTA2_FDR,
      ubi_part,
      all_of(columns)
    ) %>% 
    mutate(in_mito = !is.na(MCARTA2_FDR)) %>% 
    select(-MCARTA2_FDR)
}

merge_all <- function(prot, mito, ubihub) {
  prot$kgg %>% 
    left_join(mito$carta, by=c("gene_name" = "Symbol")) %>% 
    left_join(ubihub, by="gene_name") %>% 
    left_join(prot$total %>% select(gene_name, numpep = `Number of Peptides Quantified`), by="gene_name") %>% 
    select(
      id,
      uniprot,
      gene_name,
      description,
      site_position,
      p_value,
      fdr,
      log_fc,
      MCARTA2_FDR,
      numpep,
      ubi_part
    ) %>% 
    mutate(in_mito = !is.na(MCARTA2_FDR), in_total = !is.na(numpep)) %>% 
    select(-c(MCARTA2_FDR, numpep)) %>% 
    mutate(ubi_part = as_factor(ubi_part)) %>% 
    mutate(change = case_when(fdr < 0.05 & log_fc > 0 ~ "up", fdr < 0.05 & log_fc < 0 ~ "down", TRUE ~ "none") %>% factor(levels=c("none", "down", "up")))
}

plot_mito_change <- function(d, label.size = 8) {
  d <- d %>% 
    filter(in_mito) %>%
    group_by(gene_name) %>%
    mutate(site = rank(site_position)) %>%
    ungroup() %>%
    mutate(gene_name = as_factor(gene_name) %>% fct_rev) %>% 
    mutate(log_fc_sig = if_else(fdr < 0.05, log_fc, as.numeric(NA)))
  ggplot(d, aes(x=site, y=gene_name, fill=log_fc_sig)) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_text(size = label.size)) +
    geom_tile(colour="grey") +
    scale_fill_distiller(type="div", palette="BrBG", limits = c(-1,1)*max(abs(d$log_fc)), na.value="white") +
    labs(x="Site", y=NULL, fill=expression(log[2]~FC)) +
    scale_x_continuous(breaks=1:20, expand=expansion(mult = c(0, 0.05)))
}
