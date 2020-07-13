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
        description = Description,
        p_value = `-Log10(Welch's T-test p-value)`,
        log_fc = `-Log2 ratio (AO/UT)`
      ) %>% 
      mutate(p_value = 10^(-p_value)),
    
    kgg = readxl::read_excel(file, sheet="Kgg Proteomic (Normalized)") %>% 
      select(!starts_with("...")) %>% 
      mutate(gene_name = toupper(`Gene symbol`)) %>% 
      rename(
        description = `Protein description`,
        p_value = `-Log10(Welch's T-test p-value)`,
        log_fc = `Log2 Ratio (AO/UT)`        
      )%>% 
      mutate(p_value = 10^(-p_value)),
    
    phos = readxl::read_excel(file, sheet="Phos Proteomics (Normalized)") %>% 
      select(!starts_with("...")) %>%
      mutate(gene_name = toupper(`Gene Symbol`)) %>% 
      rename(
        description = prot_description,
        p_value = `-Log10 (Welch's T-test p-value)`,
        log_fc = `-Log2 Ratio (AO/UT)`
      )%>% 
      mutate(p_value = 10^(-p_value))
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