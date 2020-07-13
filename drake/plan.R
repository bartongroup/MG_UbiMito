get_data <- drake_plan(
  mito = read_mitocarta("data/Human.MitoCarta2.0.xls"),
  prot = read_proteomics("data/Mouse Cortical Neurons Proteomics (Protein, Kgg, Phos) - March 2018.xlsx"),
  ubihub = read_ubihub()
)

get_numbers <- drake_plan(
  kgg_n_prot = prot$kgg$`UniProt ID` %>% unique() %>% length(),
  kgg_n_sites = nrow(prot$kgg),
  total_n_prot = nrow(prot$total),
  n_mito = nrow(mito$carta),
  
  n_kgg_in_mito = kgg_mito %>% filter(in_mito) %>% pull(uniprot) %>% unique() %>% length(),
  n_kgg_ligase_in_mito = kgg_mito %>% filter(in_mito & !is.na(ubi_part)) %>% pull(uniprot) %>% unique() %>% length()
)

compare_data <- drake_plan(
  kgg_mito = merge_prot_mito(prot$kgg, mito, ubihub, c("site_position"="Site Position")),
  tot_mito = merge_prot_mito(prot$total, mito, ubihub)  
)