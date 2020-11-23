plot_volcano <- function(res, fc="logFC", p="PValue", fdr="FDR", group="contrast",
                        alpha=0.05, point.size=0.5, point.alpha=0.5, sel=NULL) {
  r <- res %>%
    mutate(
      x = !!as.name(fc),
      y = !!as.name(p),
      fdr = !!as.name(fdr),
      mrk = 1L
    )
  lmt <- max(r %>% filter(fdr<0.05) %>% pull(y))
  if(!is.null(sel)) r <- r %>% mutate(mrk = if_else(id %in% sel, 2L, mrk))
  ggplot(r, aes(x=x, y=-log10(y), colour=as.factor(mrk))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    geom_point(size=point.size, alpha=point.alpha) +
    scale_colour_manual(values=c("grey70", "black")) +
    #facet_grid(as.formula(glue::glue(". ~ {group}"))) +
    theme(legend.position = "none") +
    geom_hline(yintercept = -log10(lmt), colour="red", linetype="dashed") +
    labs(x=fc, y=glue::glue("-log10({p})"))
}

# A 'heatmap' plot of log FC for proteins in MitoCarta
# Each tile is one ub site
plot_mito_change <- function(dat, label.size = 8) {
  d <- dat %>% 
    filter(in_mito) %>%
    group_by(gene_name) %>%
    mutate(pos = str_extract(site_position, "^\\d+") %>% as.integer()) %>% 
    arrange(pos) %>% 
    mutate(site = 1:n()) %>% 
    ungroup() %>%
    mutate(gene_name = as_factor(gene_name) %>% fct_rev) %>% 
    mutate(log_fc_sig = if_else(fdr < 0.05, log_fc, as.numeric(NA)))
  ggplot(d, aes(x=site, y=gene_name, fill=log_fc_sig)) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_text(size = label.size)) +
    geom_tile(colour="grey") +
    scale_fill_distiller(type="div", palette="BrBG", limits = c(-1,1)*max(abs(d$log_fc)), na.value="white") +
    labs(x="Site", y=NULL, fill=expression(log[2]~FC)) +
    scale_x_continuous(breaks=1:20, expand=expansion(mult = c(0, 0.05)))
}


plot_protein_sites <- function(dat, gene, cex=3) {
  nd <- dat %>% 
    pivot_longer(cols=UT1_Ratio:AO5_Ratio, names_to="colname") %>%
    select(id, gene_name, colname, value, site_position, fdr) %>% 
    mutate(sample = str_remove(colname, "_.atio")) %>% 
    mutate(group = str_remove(sample, "\\d") %>% as_factor() %>% fct_relevel("UT")) %>% 
    group_by(sample) %>% 
    mutate(med = median(value, na.rm=TRUE)) %>% 
    mutate(value_norm = value / med)
  
  nd %>%
    filter(gene_name == gene) %>% 
    mutate(pos = str_extract(site_position, "^\\d+") %>% as.integer()) %>% 
    arrange(pos) %>% 
    mutate(site_position = as_factor(site_position))  %>% 
    ggplot(aes(x=group, y=log2(value_norm), colour=fdr<0.05)) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none") +
    geom_hline(yintercept = 0, colour="grey80") +
    geom_beeswarm(cex=cex) +
    facet_wrap(~site_position, nrow=1) +
    scale_colour_manual(values=c("grey", "black")) +
    labs(x=NULL, y=expression(log[2]~Ratio))
}


plot_sub_dist <- function(dat) {
  dat %>% 
    filter(in_mito) %>%
  ggplot(aes(x=sub_local, y=log_fc)) +
    theme_bw() +
    geom_boxplot(outlier.shape = NA) +
    geom_beeswarm() +
    labs("Sub-compartment", y=expression(log[2]~FC))
}


plot_stat_mito <- function(st) {
  st %>%
    filter(!str_detect(sub_local, "\\|")) %>%
    mutate(stat = factor(stat, levels = c("mito", "tot", "kgg", "kgg_de"))) %>%
    mutate(sub_local = factor(sub_local, levels=c("Matrix", "MIM", "MOM", "IMS", "Membrane", "unknown"))) %>% 
  ggplot(aes(x=sub_local, y=n, fill=stat)) +
    geom_bar(stat="identity", position=position_dodge2(padding=0), alpha=1, width=0.8, colour="grey30") +
    theme_bw() +
    scale_fill_manual(values=cbPalette) +
    theme(panel.grid = element_blank()) +
    scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
    labs(x="Sub-compartment", y="Protein count", fill="Selection")
}


plot_mito_fc <- function(dat) {
  dat %>%
    filter(in_mito) %>% 
    mutate(sub_local = factor(sub_local, levels=c("Matrix", "MIM", "MOM", "IMS", "Membrane", "unknown"))) %>% 
  ggplot(aes(x=sub_local, y=log_fc)) +
    theme_bw() +
    geom_boxplot(outlier.shape = NA, colour="grey70") +
    geom_beeswarm(aes(colour = fdr < 0.05)) +
    scale_colour_manual(values=c("grey", "black")) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    labs(x="Mitochondrial subcompartments", y=expression(log[2]~(AO/UT)))
}


plot_ineurons_venn = function(kgg, ineu) {
  kggs <- kgg %>% filter(in_mito & fdr < 0.05 & abs(log_fc) > 1.0) %>% unite(nid, gene_name, site_position) %>% pull(nid) %>% unique()
  ineus <- ineu %>% filter(in_mito & (`Significant WT 6h vs UT` == "+" | `Significant WT 2h vs UT` == "+")) %>% unite(nid, gene_name, site_position) %>% pull(nid) %>% unique()
  v <- list(
    kgg = kggs,
    ineurons = ineus
  )
  list(
    v = v,
    plot = euler(v)
  )
}