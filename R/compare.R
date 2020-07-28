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

merge_all <- function(prot, mito, ubihub, genes) {
  prot$kgg %>% 
    left_join(mito$carta, by=c("gene_name" = "Symbol")) %>% 
    left_join(ubihub, by="gene_name") %>% 
    left_join(prot$total %>% select(uniprot, numpep = `Number of Peptides Quantified`), by="uniprot") %>% 
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
    mutate(change = case_when(fdr < 0.05 & log_fc > 0 ~ "up", fdr < 0.05 & log_fc < 0 ~ "down", TRUE ~ "none") %>% factor(levels=c("none", "down", "up"))) %>% 
    select(-description) %>% 
    left_join(genes, by="gene_name") 
}

