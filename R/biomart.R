biomartGeneDownload <- function(mart) {
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


biomartGOGeneDownload <- function(mart, gene_ids) {
  getBM(
    attributes = c("ensembl_gene_id", "go_id"),
    filters = "ensembl_gene_id",
    values = gene_ids,
    mart = mart
  ) %>% 
  dplyr::rename(gene_id = ensembl_gene_id) %>% 
  dplyr::filter(go_id != "") %>% 
  tibble::as_tibble()
}



biomartGODescriptions <- function(mart) {
  getBM(attributes = c(
    "go_id",
    "name_1006",
    "namespace_1003"
  ), mart=mart) %>% 
  dplyr::rename(
    go_name = name_1006,
    go_domain = namespace_1003
  ) %>% 
  dplyr::filter(go_id != "") %>% 
  tibble::as_tibble()
}

# Reactome

loadReactome <- function() {
  url <- "https://reactome.org/download/current/Ensembl2Reactome.txt"
  colms <- c("gene_id", "reactome_id", "URL", "name", "evidence", "species")
  read_tsv(url, col_names = colms, col_types = cols())
}


# Genes and transcripts with a text file cache

biomartGetGenes <- function(mart, gene_file) {
  if(!file.exists(gene_file)) {
    biomartGeneDownload(mart) %>% write_tsv(gene_file)
  }
  read_tsv(gene_file, col_types=cols())
}

biomartGetTranscripts <- function(mart, transcript_file) {
  if(!file.exists(transcript_file)) {
    biomartTranscriptDownload(mart) %>% write_tsv(transcript_file)
  }
  read_tsv(transcript_file, col_types=cols())
}

biomartGetLengths <- function(mart, gene_length_file, gene_ids) {
  if(!file.exists(gene_length_file)) {
    biomartGeneLength(mart, gene_ids) %>% write_tsv(gene_length_file)
  }
  read_tsv(gene_length_file, col_types=cols())
}

# GO terms - returns a list with gene2term and terms.

biomartGetGODir <- function(mart, dir, gene_ids) {
  gene2go_file <- file.path(dir, "gene2go.txt")
  goterms_file <- file.path(dir, "goterms.txt")
  if(!file.exists(gene2go_file)) {
    biomartGODownload(mart, gene_ids) %>% write_tsv(gene2go_file)
  }
  if(!file.exists(goterms_file)) {
    biomartGODescriptions(mart) %>% write_tsv(goterms_file)
  }
  gene2go <- read_tsv(gene2go_file, col_types=cols()) %>% rename(term_id = go_id)
  goterms <- read_tsv(goterms_file, col_types=cols()) %>% rename(term_id = go_id, term_name = go_name)
  list(
    gene2term = gene2go,
    terms = goterms
  )
}

biomartGetGeneGO <- function(mart, gene_ids) {
  gene2go <- biomartGOGeneDownload(mart, gene_ids) %>% rename(term_id = go_id)
  goterms <- biomartGODescriptions(mart) %>% rename(term_id = go_id, term_name = go_name)
  list(
    gene2term = gene2go,
    terms = goterms
  )
}

biomartGetTranscriptGO <- function(mart, transcript_ids) {
  gene2go <- biomartGOTranscriptDownload(mart, transcript_ids) %>% rename(term_id = go_id)
  goterms <- biomartGODescriptions(mart) %>% rename(term_id = go_id, term_name = go_name)
  list(
    gene2term = gene2go,
    terms = goterms
  )
}

reactomeGetDir <- function(dir, gene_ids) {
  gene2reactome_file <- file.path(dir, "gene2reactome.txt")
  reactometerms_file <- file.path(dir, "reactometerms.txt")
  if(!file.exists(gene2reactome_file) || !file.exists(reactometerms_file)) {
    r <- loadReactome() %>% filter(gene_id %in% gene_ids)
    select(r, gene_id, reactome_id) %>% write_tsv(gene2reactome_file)
    select(r, reactome_id, name) %>% write_tsv(reactometerms_file)
  }
  gene2reactome <- read_tsv(gene2reactome_file, col_types=cols()) %>% rename(term_id = reactome_id)
  reactometerms <- read_tsv(reactometerms_file, col_types=cols()) %>% rename(term_id = reactome_id, term_name = name)
  list(
    gene2term = gene2reactome,
    terms = reactometerms
  )
  
}

reactomeGet <- function(gene_ids) {
  r <- loadReactome() %>% filter(gene_id %in% gene_ids)
  gene2reactome <- select(r, gene_id, reactome_id) %>%
    rename(term_id = reactome_id)
  reactometerms <- select(r, reactome_id, name) %>%
    rename(term_id = reactome_id, term_name = name) %>% 
    distinct()
  list(
    gene2term = gene2reactome,
    terms = reactometerms
  )
}