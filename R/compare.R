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
      sub_local = MitoCarta3.0_SubMitoLocalization,
      MitoCarta2.0_FDR,
      ubi_part,
      all_of(columns)
    ) %>% 
    mutate(in_mito = !is.na(MitoCarta2.0_FDR)) %>% 
    select(-MitoCarta2.0_FDR)
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
      sub_local = MitoCarta3.0_SubMitoLocalization,
      MitoCarta2.0_FDR,
      numpep,
      ubi_part
    ) %>% 
    mutate(in_mito = !is.na(MitoCarta2.0_FDR), in_total = !is.na(numpep)) %>% 
    select(-c(MitoCarta2.0_FDR, numpep)) %>% 
    mutate(ubi_part = as_factor(ubi_part)) %>% 
    mutate(change = case_when(fdr < 0.05 & log_fc > 0 ~ "up", fdr < 0.05 & log_fc < 0 ~ "down", TRUE ~ "none") %>% factor(levels=c("none", "down", "up"))) %>% 
    select(-description) %>% 
    left_join(genes, by="gene_name") 
}


merge_ineurons_mito <- function(kgg, ineu) {
  ineu_sel <- ineu %>%
    filter((`Significant WT 6h vs UT` == "+" | `Significant WT 2h vs UT` == "+") & MitoCarta2.0 == "+")
  kgg_sel <- kgg %>% 
    filter(fdr < 0.05 & in_mito)
  
  full_join(ineu_sel, kgg_sel, by=c("gene_name", "site_position"))
}


make_stat_mito <- function(mito, tot, kgg) {
  stat_mito <- mito$carta %>%
    rename(sub_local = MitoCarta3.0_SubMitoLocalization) %>% 
    group_by(sub_local) %>% tally() %>% mutate(stat = "mito")
  stat_kgg <- kgg_mito %>%
    filter(in_mito) %>%
    group_by(sub_local) %>% summarise(n = length(unique(uniprot))) %>% mutate(stat = "kgg")
  stat_kgg_de <- kgg_mito %>%
    filter(in_mito & fdr < 0.05) %>%
    group_by(sub_local) %>% summarise(n = length(unique(uniprot))) %>% mutate(stat = "kgg_de")
  stat_tot <- tot %>%
    filter(in_mito) %>%
    group_by(sub_local) %>% summarise(n = length(unique(uniprot))) %>% mutate(stat = "tot")
  bind_rows(stat_mito, stat_kgg, stat_kgg_de, stat_tot)
}

