get_data <- drake_plan(
  mito = read_mitocarta("data/Human.MitoCarta2.0.xls"),
  prot_raw = read_proteomics("data/Mouse Cortical Neurons Proteomics (Protein, Kgg, Phos) - March 2018.xlsx"),
  ubihub = read_ubihub()
)

process_data <- drake_plan(
  prot = prot_raw %>%
    diff_expr() %>% 
    find_basal()
)

get_numbers <- drake_plan(
  kgg_n_prot = prot$kgg$uniprot %>% unique() %>% length(),
  kgg_n_prot_basal = prot$kgg_basal$uniprot %>% unique() %>% length(),
  kgg_n_sites = nrow(prot$kgg),
  total_n_prot = nrow(prot$total),
  n_mito = nrow(mito$carta),
  
  n_kgg_in_mito = kgg_mito %>% filter(in_mito) %>% pull(uniprot) %>% unique() %>% length(),
  n_kgg_ligase_in_mito = kgg_mito %>% filter(in_mito & !is.na(ubi_part)) %>% pull(uniprot) %>% unique() %>% length(),
  
  n_kgg_de_welch = prot$kgg %>% filter(welch_p_value < 0.05) %>% nrow(),
  n_kgg_de_fdr = prot$kgg %>% filter(fdr < 0.05) %>% nrow()
)

compare_data <- drake_plan(
  kgg_mito = merge_prot_mito(prot$kgg_basal, mito, ubihub, "site_position"),
  tot_mito = merge_prot_mito(prot$total, mito, ubihub),
  all_data = merge_all(prot, mito, ubihub)
)

make_figures <- drake_plan(
  kgg_in_mito = kgg_mito %>% filter(in_mito) %>% pull(id),
  fig_kgg_volcano = plotVolcano(prot$kgg, fc="log_fc",  fdr="fdr", p="p_value", sel=kgg_in_mito),
  fig_mito_change = plot_mito_change(all_data)
)

save_tables <- drake_plan(
  save_kgg_de = save_table(prot$kgg, "kgg_de.tsv"),
  save_total_de = save_table(prot$total, "total_de.tsv"),
  save_kgg_mito = save_table(kgg_mito, "kgg_mito.tsv"),
  save_total_mito = save_table(tot_mito, "total_mito.tsv")
)
