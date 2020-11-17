get_biomart <- drake_plan(
  mart = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version="100"),
  bm_genes = bm_fetch_genes_cached(mart, "cache/genes.rds") %>% mutate(gene_name = toupper(gene_name)),
  bm_go = bm_fetch_go_cached(mart, all_data$gene_id, "cache/go_terms.rds"),
  bm_go_slim = bm_fetch_go_cached(mart, all_data$gene_id, "cache/go_terms_slim.rds", slim=TRUE),
  reactome = fetch_reactome_cached(all_data$gene_id, "cache/rectome.rds")
)

get_data <- drake_plan(
  mito = read_mitocarta("3.0"),
  prot_raw = read_proteomics("data/Mouse Cortical Neurons Proteomics (Protein, Kgg, Phos) - March 2018.xlsx"),
  ubihub = read_ubihub(),
  ineurons = read_ineurons("data/mmc2.xlsx", "Kgg Quant NGN2 (AAVS) - WCL"),
  pink_raw = read_pink("data/AO3009_3020 - Protein - For Miratul.xlsx", "Proteins_Final")
)

process_data <- drake_plan(
  prot = prot_raw %>%
    diff_expr() %>% 
    find_basal(),
  pink = pink_raw %>% 
    pink_limma_de()
)

get_numbers <- drake_plan(
  kgg_n_prot = prot$kgg$uniprot %>% unique() %>% length(),
  kgg_n_prot_basal = prot$kgg_basal$uniprot %>% unique() %>% length(),
  kgg_n_sites = nrow(prot$kgg),
  total_n_prot = nrow(prot$total),
  
  n_mito = nrow(mito$carta),
  stat_mito = mito$carta %>% group_by(MitoCarta3.0_SubMitoLocalization) %>% tally() %>% arrange(desc(n)),
  stat_kgg_mito = kgg_mito %>% filter(in_mito) %>% group_by(sub_local) %>% summarise(n = length(unique(uniprot))) %>% arrange(desc(n)),
  
  n_kgg_in_mito = kgg_mito %>% filter(in_mito) %>% pull(uniprot) %>% unique() %>% length(),
  n_kgg_ligase_in_mito = kgg_mito %>% filter(in_mito & !is.na(ubi_part)) %>% pull(uniprot) %>% unique() %>% length(),
  
  n_kgg_de_welch = prot$kgg %>% filter(welch_p_value < 0.05) %>% nrow(),
  n_kgg_de_fdr = prot$kgg %>% filter(fdr < 0.05) %>% nrow()
)

compare_data <- drake_plan(
  kgg_mito = merge_prot_mito(prot$kgg_basal, mito, ubihub, "site_position"),
  tot_mito = merge_prot_mito(prot$total, mito, ubihub),
  all_data = merge_all(prot, mito, ubihub, bm_genes),
  kgg_ineurons = merge_ineurons_mito(kgg_mito, ineurons)
)

make_figures <- drake_plan(
  kgg_in_mito = kgg_mito %>% filter(in_mito) %>% pull(id),
  fig_kgg_volcano = plot_volcano(prot$kgg, fc="log_fc",  fdr="fdr", p="p_value", sel=kgg_in_mito),
  fig_mito_change = plot_mito_change(all_data)
)

save_tables <- drake_plan(
  save_kgg_de = save_table(prot$kgg, "kgg_de.tsv"),
  save_total_de = save_table(prot$total, "total_de.tsv"),
  save_kgg_mito = save_table(kgg_mito, "kgg_mito.tsv"),
  save_total_mito = save_table(tot_mito, "total_mito.tsv")
)

save_for_shiny <- drake_plan(
  shiny_all = shiny_data_all(all_data, bm_go, bm_go_slim, reactome),
  save_shiny_all = saveRDS(shiny_all, file_out("shiny_all/data.rds")),
)