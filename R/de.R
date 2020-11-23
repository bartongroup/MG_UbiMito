cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


limma_de <- function(d) {
  tab <- d %>% 
    select(id, ends_with("Ratio")) %>% 
    column_to_rownames("id") %>% 
    as.matrix()
  fnorm <- apply(tab, 2, median)
  tab <- t(t(tab) / fnorm)
  tab <- log2(tab)
  
  meta <- tibble(colname = colnames(tab)) %>% 
    mutate(sample = str_remove(colname, "_.atio")) %>% 
    mutate(group = str_remove(sample, "\\d") %>% as_factor() %>% fct_relevel("UT"))

  design_mat <- model.matrix(~group, data=meta)

  fit <- tab %>% 
    lmFit(design_mat) %>% 
    eBayes()
  
  coef <- colnames(fit$design)[2]
  res <- topTable(fit, coef=coef, adjust="BH", sort.by="none", number=1e9) %>% 
    as_tibble(rownames = "id") %>% 
    select(id, log_fc=logFC, p_value=P.Value, fdr=adj.P.Val, ave_ratio=AveExpr)

  d %>% 
    left_join(res, by="id")
}

diff_expr <- function(prot) {
  prot$kgg <- limma_de(prot$kgg)
  prot$total <- limma_de(prot$total)
  prot
}



pink_limma_de <- function(d) {
  tab <- d %>% 
    select(id, ends_with("Ratio_(Total)")) %>% 
    column_to_rownames("id") %>% 
    as.matrix()
  fnorm <- apply(tab, 2, median)
  tab <- t(t(tab) / fnorm)

  meta <- tibble(colname = colnames(tab)) %>% 
    mutate(sample = str_remove(colname, "_Ratio_\\(Total\\)")) %>% 
    separate(sample, c("genotype", "treatment", "replicate"), remove=FALSE) %>% 
    mutate(genotype = genotype %>% as_factor() %>% fct_relevel("WT")) %>% 
    mutate(treatment = treatment %>% as_factor() %>% fct_relevel("UT"))
    

  design_mat <- model.matrix(~genotype + treatment, data=meta)
  
  fit <- tab %>% 
    lmFit(design_mat) %>% 
    eBayes()
  
  coefs <- colnames(fit$design)[-1]
  res <- map_dfr(coefs, function(coef) {
    topTable(fit, coef=coef, adjust="BH", sort.by="none", number=1e9) %>% 
      as_tibble(rownames = "id") %>% 
      select(id, log_fc=logFC, p_value=P.Value, fdr=adj.P.Val, ave_ratio=AveExpr) %>% 
      mutate(contrast = coef)
  })
  
  d %>% 
    left_join(res, by="id")
}


mean_expression <- function(dat, root) {
  dat %>% 
    select(id, ends_with("Ratio")) %>% 
    select(id, starts_with(root)) %>% 
    pivot_longer(-id) %>% 
    group_by(id) %>% 
    summarise(M = mean(value, na.rm=TRUE), SD = sd(value, na.rm=TRUE))
}


