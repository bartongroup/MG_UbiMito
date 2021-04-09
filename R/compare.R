find_basal <- function(prot) {
  quant_prot <- prot$total$uniprot %>% unique()
  prot$kgg_basal <- prot$kgg %>% filter(uniprot %in% quant_prot)
  prot
}


merge_prot_mito <- function(dat, mito, ubihub, columns=NULL) {
  dat %>% 
    left_join(mito, by=c("gene_name" = "genes")) %>% 
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


merge_total <- function(prot, mito, ubihub, genes) {
  prot$total %>% 
    left_join(mito, by=c("gene_name" = "genes")) %>% 
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
      ubi_part
    ) %>% 
    mutate(in_mito = !is.na(sub_local)) %>% 
    mutate(change = case_when(fdr < 0.05 & log_fc > 0 ~ "up", fdr < 0.05 & log_fc < 0 ~ "down", TRUE ~ "none") %>% factor(levels=c("none", "down", "up"))) %>% 
    select(-description) %>% 
    left_join(genes, by="gene_name") %>% 
    mutate_at(c("ubi_part", "sub_local", "gene_biotype"), as_factor)
}


merge_kgg <- function(prot, mito, ubihub, genes) {
  prot$kgg_norm %>% 
    left_join(mito, by=c("gene_name" = "genes")) %>% 
    left_join(ubihub, by="gene_name") %>% 
    left_join(prot$total %>%
                select(uniprot, numpep=`Number of Peptides Quantified`),
              by="uniprot"
    ) %>% 
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
      numpep,
      ubi_part
    ) %>% 
    mutate(in_mito = !is.na(sub_local), in_total = !is.na(numpep)) %>% 
    select(-c(numpep)) %>% 
    mutate(change = case_when(fdr < 0.05 & log_fc > 0 ~ "up", fdr < 0.05 & log_fc < 0 ~ "down", TRUE ~ "none") %>% factor(levels=c("none", "down", "up"))) %>% 
    select(-description) %>% 
    left_join(genes, by="gene_name") %>% 
    mutate_at(c("ubi_part", "sub_local", "gene_biotype"), as_factor)
}


merge_ineurons_mito <- function(ineu, mito) {
  ineu %>% 
    left_join(mito, by=c("gene_name" = "genes")) %>%
    rename(sub_local = MitoCarta3.0_SubMitoLocalization) %>% 
    mutate(in_mito = !is.na(sub_local))
}


make_stat_mito <- function(mito, tot, kgg) {
  stat_mito <- mito$carta %>%
    rename(sub_local = MitoCarta3.0_SubMitoLocalization) %>% 
    group_by(sub_local) %>% tally() %>% mutate(stat = "mito")
  stat_kgg <- kgg %>%
    filter(in_mito) %>%
    group_by(sub_local) %>% summarise(n = length(unique(uniprot))) %>% mutate(stat = "kgg")
  stat_kgg_de <- kgg %>%
    filter(in_mito & fdr < 0.05) %>%
    group_by(sub_local) %>% summarise(n = length(unique(uniprot))) %>% mutate(stat = "kgg_de")
  stat_tot <- tot %>%
    filter(in_mito) %>%
    group_by(sub_local) %>% summarise(n = length(unique(uniprot))) %>% mutate(stat = "tot")
  bind_rows(stat_mito, stat_kgg, stat_kgg_de, stat_tot)
}

reduce_sites <- function(x) {
  tibble(full = x) %>% 
    separate(full, c("gene_name", "site_position"), sep=":") %>% 
    mutate(site_num = site_position %>% str_remove(";\\d+$") %>% as.numeric()) %>% 
    group_by(gene_name) %>% 
    arrange(site_num) %>% 
    summarise(sites = paste(site_position, collapse=",") %>% paste0("(", ., ")")) %>% 
    unite(final, c(gene_name, sites), sep=" ") %>% 
    pull(final)
}

kgg_ineurons_overlap <- function(kgg, ineu) {
  kggs <- kgg %>%
    filter(in_mito & fdr < 0.05 & abs(log_fc) > 1.0) %>% 
    unite(gs, c(gene_name, site_position), sep=":") %>% 
    pull(gs) %>% unique()
  ineus <- ineu %>%
    filter(in_mito & (`Significant WT 6h vs UT` == "+" | `Significant WT 2h vs UT` == "+")) %>%
    unite(gs, c(gene_name, site_position), sep=":") %>% 
    pull(gs) %>% unique()
  v <- list(
    kgg = kggs,
    ineu = ineus,
    kgg_only = setdiff(kggs, ineus),
    ineu_only = setdiff(ineus, kggs),
    overlap = intersect(kggs, ineus)
  )
  vp <- map(v, reduce_sites)
  list(
    v = v,
    vp = vp
  )
}
