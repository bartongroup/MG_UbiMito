# Get all genes

bm_fetch_genes <- function(mart) {
  getBM(attributes = c(
    "gene_biotype",
    "ensembl_gene_id",
    "external_gene_name",
    "description"
  ), mart=mart) %>% 
  dplyr::rename(
    gene_id = ensembl_gene_id,
    gene_name = external_gene_name
  ) %>%
  dplyr::mutate(description = str_remove(description, "\\s\\[.*\\]")) %>% 
  tibble::as_tibble()
}

# Get genes with a cache

bm_fetch_genes_cached <- function(mart, gene_file) {
  if(!file.exists(gene_file)) {
    bm_fetch_genes(mart) %>% write_rds(gene_file)
  }
  read_rds(gene_file)
}



# For a given list of Ensembl gene IDs, fetch all GO-terms

bm_fetch_go_genes <- function(mart, gene_ids, slim=FALSE) {
  id <- ifelse(slim, "goslim_goa_accession", "go_id")
  go <- getBM(
    attributes = c("ensembl_gene_id", id),
    filters = "ensembl_gene_id",
    values = gene_ids,
    mart = mart
  ) %>% 
    dplyr::rename(gene_id = ensembl_gene_id, term_id = !!sym(id)) %>% 
    dplyr::filter(term_id != "") %>% 
    tibble::as_tibble()
  terms <- go$term_id %>% unique()
  map(terms, function(trm) go[go$term_id == trm, ]$gene_id) %>% 
    set_names(terms)
}

# Get all GO-term descriptions

bm_fetch_go_descriptions <- function(mart) {
  # filtering on GO-terms does not work properly, so I have to fetch all terms
  getBM(
    attributes = c("go_id", "name_1006", "definition_1006", "namespace_1003"),
    mart = mart) %>% 
  dplyr::rename(
    term_id = go_id,
    term_name = name_1006,
    term_description = definition_1006,
    term_domain = namespace_1003
  ) %>% 
  dplyr::filter(term_id != "") %>% 
  tibble::as_tibble()
}

# Fetch GO-gene and GO descriptions

bm_fetch_go <- function(mart, gene_ids, slim=FALSE) {
  gene2go <- bm_fetch_go_genes(mart, gene_ids, slim) 
  goterms <- bm_fetch_go_descriptions(mart)
  list(
    term2gene = gene2go,
    terms = goterms
  )
}

# The same, cached

bm_fetch_go_cached <- function(mart, gene_ids, go_file, slim=FALSE) {
  if(!file.exists(go_file)) {
    bm_fetch_go(mart, gene_ids, slim) %>% write_rds(go_file)
  }
  read_rds(go_file)
}




# Reactome raw data

fetch_reactome_raw <- function() {
  url <- "https://reactome.org/download/current/Ensembl2Reactome.txt"
  colms <- c("gene_id", "reactome_id", "URL", "name", "evidence", "species")
  read_tsv(url, col_names = colms, col_types = cols())
}

# Get Reactome data in the same format as GO-data

fetch_reactome <- function(gene_ids) {
  r <- fetch_reactome_raw() %>% filter(gene_id %in% gene_ids)
  g2r <- select(r, gene_id, reactome_id) %>%
    rename(term_id = reactome_id)
  reactometerms <- select(r, reactome_id, name) %>%
    rename(term_id = reactome_id, term_name = name) %>% 
    distinct()
  terms <- g2r$term_id %>% unique()
  gene2reactome <- map(terms, function(trm) g2r[g2r$term_id == trm, ]$gene_id) %>% 
    set_names(terms)
  list(
    term2gene = gene2reactome,
    terms = reactometerms
  )
}

# The same, cached

fetch_reactome_cached <- function(gene_ids, reactome_file) {
  if(!file.exists(reactome_file)) {
    fetch_reactome(gene_ids) %>% write_rds(reactome_file)
  }
  read_rds(reactome_file)
}
