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