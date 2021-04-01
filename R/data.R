read_mitocarta <- function(version) {
  file <- glue("data/Human.MitoCarta{version}.xls")
  list(
    carta = readxl::read_excel(file, sheet=glue("A Human MitoCarta{version}")),
    ranking = readxl::read_excel(file, sheet="B Human All Genes")
  )
}

separate_mitocarta_genes <- function(mc) {
  mc %>% 
    unite(genes, c(Symbol, Synonyms), sep="|") %>%
    separate_rows(genes, sep = "\\|") %>% 
    filter(genes != "-") %>% 
    # remove duplicates, select highest score
    group_split(genes) %>% 
    map_dfr(function(w) {
      w %>% 
        arrange(desc(MitoCarta2.0_Score)) %>% 
        head(1)
    })
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
      rename_at(
        vars(contains("_Ratio")),
        ~str_replace(.x, "_Ratio", "_ratio")
      ) %>% 
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


normalise_prot <- function(prt) {
  prt_tot <- prt$tot %>% 
    pivot_longer(cols=UT1_ratio:AO5_ratio, names_to="colname") %>%
    select(uniprot, colname, tot=value)
  prt_kgg <- prt$kgg %>% 
    pivot_longer(cols=UT1_ratio:AO5_ratio, names_to="colname") %>%
    select(id, uniprot, colname, kgg=value) %>% 
    left_join(prt_tot, by=c("uniprot", "colname")) %>% 
    drop_na() %>% 
    mutate(rat = kgg / tot) %>% 
    pivot_wider(id_cols = id, names_from=colname, values_from=rat)
  kgg_norm <- prt$kgg %>% 
    select(-(UT1_ratio:AO5_ratio)) %>% 
    right_join(prt_kgg, by="id")
  prt$kgg_norm <- kgg_norm
  prt
}

read_ineurons <- function(file, sheet) {
  readxl::read_excel(file, sheet=sheet, skip=1) %>% 
    select(!starts_with("...")) %>% 
    mutate(gene_name = toupper(`Gene Symbol`)) %>% 
    rename(
      uniprot = `UniProt ID`,
      site_position = `Site Position`
    )%>% 
    unite("id", uniprot, site_position, remove=FALSE)
}


read_pink <- function(file, sheet) {
  readxl::read_excel(file, sheet=sheet, skip=1) %>% 
    select(!starts_with("...")) %>% 
    mutate(gene_name = toupper(`Gene Symbol`)) %>% 
    rename(
      uniprot = Accession,
      description = Description,
    )%>% 
    mutate(id = uniprot)
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


shiny_data_all <- function(all_data, bm_go, bm_go_slim, reactome) {
  gene2name <- set_names(all_data$gene_name, all_data$gene_id)
  bm_go$terms <- bm_go$terms %>% select(-term_description)
  bm_go_slim$terms <- bm_go_slim$terms %>% select(-term_description)
  list(
    kgg = all_data %>% 
      mutate(
        sig = fdr < 0.05,
        ubi = replace_na(as.character(ubi_part), "-"),
        e2 = ubi == "E2",
        e3_simple = ubi == "E3 simple",
        e3_complex = ubi == "E3 complex",
        sub = replace_na(as.character(sub_local), "-")
      ),
    bm_go = bm_go,
    bm_go_slim = bm_go_slim,
    reactome = reactome,
    gene2name = gene2name
  )
}