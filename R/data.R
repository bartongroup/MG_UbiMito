read_mitocarta <- function(file) {
  list(
    carta = readxl::read_excel(file, sheet="A Human MitoCarta2.0"),
    ranking = readxl::read_excel(file, sheet="B Human All Genes")
  )
}

read_proteomics <- function(file) {
  list(
    total = readxl::read_excel(file, sheet="Total Protein (Normalized)") %>% 
      select(!starts_with("...")) %>% 
      mutate(gene_name = toupper(`Gene Symbol`)) %>% 
      rename(
        uniprot = `UniProt ID`,
        description = Description,
        welch_p_value = `-Log10(Welch's T-test p-value)`,
        welch_log_fc = `-Log2 ratio (AO/UT)`
      ) %>% 
      mutate(welch_p_value = 10^(-welch_p_value)) %>% 
      mutate(id = uniprot),
    
    kgg = readxl::read_excel(file, sheet="Kgg Proteomic (Normalized)") %>% 
      select(!starts_with("...")) %>% 
      mutate(gene_name = toupper(`Gene symbol`)) %>% 
      rename(
        uniprot = `UniProt ID`,
        site_position = `Site Position`,
        description = `Protein description`,
        welch_p_value = `-Log10(Welch's T-test p-value)`,
        welch_log_fc = `Log2 Ratio (AO/UT)`        
      )%>% 
      mutate(welch_p_value = 10^(-welch_p_value)) %>% 
      unite("id", uniprot, site_position, remove=FALSE),
    
    phos = readxl::read_excel(file, sheet="Phos Proteomics (Normalized)") %>% 
      select(!starts_with("...")) %>%
      mutate(gene_name = toupper(`Gene Symbol`)) %>% 
      rename(
        uniprot = `UniProt ID`,
        site_position = `Site Position`,
        description = prot_description,
        welch_p_value = `-Log10 (Welch's T-test p-value)`,
        welch_log_fc = `-Log2 Ratio (AO/UT)`
      )%>% 
      mutate(welch_p_value = 10^(-welch_p_value)) %>% 
      unite("id", uniprot, site_position, remove=FALSE)
  )
}

read_names <- function(file, name) {
  read_tsv(file, col_names = "gene_name", col_types = "c") %>% 
    mutate(ubi_part = name)
}


read_ubihub <- function() {
  bind_rows(
    read_names("data/E2.txt", "E2"),
    read_names("data/E3_simple.txt", "E3 simple"),
    read_names("data/E3_complex.txt", "E3 complex")
  )
}

save_table <- function(dat, file) {
  if(!dir.exists("tab")) dir.create("tab")
  file <- file.path("tab", file)
  write_tsv(dat, file)
}


shiny_data_all <- function(all_data, bm_go, reactome) {
  gene2name <- set_names(all_data$gene_name, all_data$gene_id)
  list(
    kgg = all_data %>% mutate(sig = fdr < 0.05),
    bm_go = bm_go,
    reactome = reactome,
    gene2name = gene2name
  )
}