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
  d_tot <- d %>% 
    select(gene_name, tot_log_fc) %>% 
    distinct() %>% 
    arrange(tot_log_fc)
  genes <- d_tot$gene_name
  
  d <- d %>% mutate(gene_name = factor(gene_name, levels = genes))
  d_tot <- d_tot %>% mutate(gene_name = factor(gene_name, levels = genes))
  
  lims <- c(-1,1)*max(abs(d$log_fc))
  lims_tot <- c(-1,1)*max(abs(d_tot$tot_log_fc))
  
  g1 <- ggplot(d_tot, aes(x="Total", y=gene_name, fill=tot_log_fc)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = label.size),
      axis.text.y = element_blank(),
      plot.margin = margin(unit(c(5.5, 0, 5.5, 5.5), "pt")),
      legend.position = "left"
    ) +
    scale_y_discrete(position="right") +
    geom_tile(colour="grey") +
    scale_fill_distiller(type="div", palette="RdBu", limits = lims_tot, na.value="white") +
    labs(x=NULL, y=NULL, fill=expression(log[2]~FC[tot]))
  
  g2 <- ggplot(d, aes(x=site, y=gene_name, fill=log_fc_sig)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = label.size),
      axis.text.y = element_text(hjust=0.5)
    ) +
    geom_tile(colour="grey") +
    scale_fill_distiller(type="div", palette="BrBG", limits = c(-1,1)*max(abs(d$log_fc)), na.value="white") +
    labs(x="Site", y=NULL, fill=expression(log[2]~FC)) +
    scale_x_continuous(breaks=1:20, expand=expansion(mult = c(0, 0.05)))
  plot_grid(g1, g2, align="h", rel_widths = c(1,3))
}

plot_mito_change_sep <- function(dat, label.size = 3, site.size=2.5, ncol=10) {
  d <- dat %>% 
    filter(in_mito) %>%
    group_by(gene_name) %>%
    mutate(pos = str_extract(site_position, "^\\d+") %>% as.integer()) %>% 
    ungroup()
   
  
  genes <- d$gene_name %>% unique() %>% sort()
  mx <- max(abs(d$log_fc))
  
  P <- map(genes, function(gene) {
    ds <- d %>%
      filter(gene_name == gene) %>% 
      select(site_position, log_fc, change, pos) %>% 
      distinct() %>% 
      arrange(pos) %>% 
      mutate(idx = row_number())
    ns <- nrow(ds)
    
    dp <- tibble(x=c(0.5,ns+0.5,ns+0.5,0.5, 0.5), y=c(0,0,-1,-1, 0))

          
    ggplot(ds, aes(x=idx, y=abs(log_fc))) +
      theme_nothing() +
      theme(
        #axis.text.x = element_text(size=label.size, angle=90, hjust=1),
        plot.margin = margin(unit(c(5, 2, 1, 2), "pt"))
      ) +
      geom_polygon(data=dp, aes(x,y), fill="grey80") +
      geom_segment(aes(xend=idx, yend=0), colour="grey80") +
      geom_point(aes(colour=change), size=1.3) +
      #geom_hline(yintercept = 0) +
      #geom_hline(yintercept = -1) +
      scale_x_continuous(limits=c(0.5, ns+0.5), expand=c(0,0)) + 
      scale_y_continuous(limits=c(-1, mx*1.05)) +
      scale_colour_manual(values=c("grey", "blue", "red"), drop=FALSE) +
      annotate("text", label=gene, x=(ns+1)/2, y=-0.5, hjust=0.5, size=label.size) +
      geom_text(aes(label=site_position), angle=90, hjust=0, nudge_y = mx/20, size=site.size)
  })
  
  plot_grid(plotlist=P, ncol=ncol)
 
}

plot_protein <- function(dat, gene, cex=3) {
  nd <- dat %>% 
    pivot_longer(cols=UT1_ratio:AO5_ratio, names_to="colname") %>%
    select(id, gene_name, colname, value, fdr) %>% 
    mutate(sample = str_remove(colname, "_ratio")) %>% 
    separate(sample, c("group", "replicate"), 2, remove=FALSE) %>% 
    group_by(sample) %>% 
    mutate(med = median(value, na.rm=TRUE)) %>% 
    mutate(value_norm = value / med)
  
  nd %>%
    filter(gene_name == gene) %>% 
  ggplot(aes(x=group, y=log2(value_norm), shape=fdr<0.05, colour=replicate)) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none") +
    geom_hline(yintercept = 0, colour="grey80") +
    geom_beeswarm(cex=cex) +
    scale_colour_manual(values=cbPalette) +
    scale_shape_manual(values=c(1, 19)) +
    labs(x=NULL, y=expression(log[2]~Ratio))
  
}

plot_protein_sites <- function(dat, gene, cex=3) {
  nd <- dat %>% 
    pivot_longer(cols=UT1_ratio:AO5_ratio, names_to="colname") %>%
    select(id, gene_name, colname, value, site_position, fdr) %>% 
    mutate(sample = str_remove(colname, "_ratio")) %>% 
    separate(sample, c("group", "replicate"), 2, remove=FALSE) %>% 
    mutate(group = as_factor(group) %>% fct_relevel("UT")) %>% 
    group_by(sample) %>% 
    mutate(med = median(value, na.rm=TRUE)) %>% 
    mutate(value_norm = value / med)
  
  nd %>%
    filter(gene_name == gene) %>% 
    mutate(pos = str_extract(site_position, "^\\d+") %>% as.integer()) %>% 
    arrange(pos) %>% 
    mutate(site_position = as_factor(site_position))  %>% 
  ggplot(aes(x=group, y=log2(value_norm), shape=fdr<0.05, colour=replicate)) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none") +
    geom_hline(yintercept = 0, colour="grey80") +
    geom_beeswarm(cex=cex) +
    facet_wrap(~site_position, nrow=1) +
    scale_colour_manual(values=cbPalette) +
    scale_shape_manual(values=c(1, 19)) +
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
  kggs <- kgg %>% filter(in_mito & fdr < 0.05 & abs(log_fc) > 1.0) %>% pull(gene_name) %>% unique()
  ineus <- ineu %>% filter(in_mito & (`Significant WT 6h vs UT` == "+" | `Significant WT 2h vs UT` == "+")) %>% pull(gene_name) %>% unique()
  v <- list(
    kgg = kggs,
    ineurons = ineus
  )
  list(
    v = v,
    plot = euler(v)
  )
}