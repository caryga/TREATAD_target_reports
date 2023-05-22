
library(igraph)
library(tidygraph)
library(ggraph)
library(tidyverse)

synapser::synLogin()

perc.rank <- function(x) trunc(rank(x))/length(x)

theme_set(theme_bw())

# Gather data -------------------------------------------------------------

##
# risk scores

scores <- synapser::synTableQuery('select * from syn25575156')$filepath %>% read_csv() %>%
  mutate(pctl = perc.rank(Overall)*100) %>%
  select(ENSG, GeneName, Overall, rank= Overall_rank, pctl, 
         Omics = OmicsScore, Genetics = GeneticsScore
         # , Literature = LiteratureScore
         # , Neuropath = NeuropathScore,
         # , Druggability = SM_Druggability_bucket
         # , Safety = safety_bucket
         # , Feasibility = feasibility_bucket
         ) %>%
  # mutate(
  #   Druggability = 2*(max(Druggability, na.rm=T) + 1 - Druggability) / max(Druggability, na.rm=T),
  #   # Safety = 2*(max(Safety, na.rm=T) + 1 - Safety) / max(Safety, na.rm=T),
  #   # Feasibility = 2*(max(Feasibility, na.rm=T) + 1 - Feasibility) / max(Feasibility, na.rm=T),
  # ) %>%
  pivot_longer(cols = c(Overall, Omics, Genetics
                        # , Literature
                        # , Neuropath
                        # , Druggability
                        # , Safety, Feasibility
                        )) %>%
  mutate(
    name = fct_relevel(name, c(
      # 'Feasibility','Safety',
      # 'Druggability',
      # 'Neuropath',
      # 'Literature',
      'Omics', 'Genetics','Overall')),
    value = case_when( value == 0 ~ NA_real_, T ~ value)
  ) %>%
  arrange(name)

##
# genetics scores

gen <- synapser::synTableQuery('SELECT * FROM syn26844312')$filepath %>% read_csv() %>%
  filter(!is.na(GeneName), !duplicated(GeneName)) %>%
  mutate(pctl = perc.rank(GeneticsScore)*100) %>%
  mutate(
    meanRank_qtlFDR = if_else( meanRank_qtlFDR == 0, NA_real_, meanRank_qtlFDR)
  )%>%
  select(GeneName, rank=score_rank, pctl,
         meanRank_gwasP, meanRank_qtlFDR,
         coding_variant_summary, noncoding_variant_summary,
         Hsap_pheno_score, Ortholog_pheno_score
  ) %>%
  pivot_longer(-c(GeneName,rank, pctl)) %>%
  mutate( name = fct_relevel(name,
                             c('Ortholog_pheno_score',
                               'Hsap_pheno_score',
                               'noncoding_variant_summary',
                               'coding_variant_summary',
                               'meanRank_qtlFDR',
                               'meanRank_gwasP')),
          value = case_when( value == 0 ~ NA_real_, T ~ value)) %>%
  arrange( name )

##
# omics scores

omics <- read_csv( synapser::synTableQuery('select * from syn22758536')$filepath ) %>%
  select(GName, RNA_TE, RNA_fdr_CorPVal, Pro_TE, Pro_fdr_CorPVal) %>%
  pivot_longer(
    c(RNA_TE, Pro_TE), names_sep = '_', names_to = c('x', 'y'), values_to = 'TE'
  ) %>%
  pivot_longer(
    c(RNA_fdr_CorPVal, Pro_fdr_CorPVal), names_sep = '_', names_to = c('a', 'b'), values_to = 'fdr'
  ) %>%
  distinct() %>% filter(x == a) %>% select(-y, -a, -b) %>%
  mutate(x = str_replace_all(x, 'Pro','Protein'))


## 
# gene-biodom mapping

biodom_genes <- readRDS(synapser::synGet('syn25428992')$path) %>%
  select(Biodomain, GO_ID, hgnc_symbol) %>%
  unnest_longer(hgnc_symbol) %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != '') %>%
  group_by(hgnc_symbol, Biodomain) %>%
  summarise(n_term = length(unique(GO_ID))) %>%
  ungroup() %>%
  left_join(
    .,
    readRDS(synapser::synGet('syn25428992')$path) %>%
      select(Biodomain, GO_ID, hgnc_symbol) %>%
      group_by(Biodomain) %>%
      summarise(bd_terms = length(unique(GO_ID))) %>%
      ungroup(),
    by = 'Biodomain') %>%
  mutate(pct = 100*(n_term / bd_terms) ) %>%
  left_join(
    .,
    readRDS(synapser::synGet('syn25428992')$path) %>%
      select(Biodomain, hgnc_symbol, n_hgncSymbol, GO_ID, GOterm_Name ) %>%
      unnest_longer(hgnc_symbol),
    by = c('hgnc_symbol','Biodomain') ) %>%
  full_join(
    .,
    # domain labels
    read_csv(synapser::synGet('syn26856828')$path, col_types = cols()),
    by = c('Biodomain'='domain')) %>%
  mutate(Biodomain = case_when(Biodomain == 'none' ~ NA_character_, T ~ Biodomain)) %>%
  rename(symbol = hgnc_symbol, n_symbol = n_hgncSymbol)

##
# biodom term nw

biodom = readRDS(synapser::synGet('syn25428992')$path)
domains = biodom$Biodomain %>% unique() %>% sort()
dom.cols = read_csv(synapser::synGet('syn26856828')$path, col_types = cols())

ks = read_csv(synapser::synGet('syn51163744')$path)
ks1 <- ks %>% filter(goid_1 != goid_2, kappa > 0.5)

# # look @ N overlap by kappa score level
# x = map_dfr(
#   seq(.1,1,.1),
#   ~ ks1 %>% arrange(desc(kappa)) %>% filter(kappa <= .x) %>% slice(1:5)
# ) %>%
#   rowwise() %>%
#   mutate(
#     n_genes_1 = biodom$n_symbol[biodom$GO_ID == goid_1] %>% unique(),
#     n_genes_2 = biodom$n_symbol[biodom$GO_ID == goid_2] %>% unique(),
#     n_overlap = length(intersect(
#       biodom %>% filter(GO_ID == goid_1) %>% slice(1) %>%
#         pull(symbol) %>% unlist() %>% unique(),
#       biodom %>% filter(GO_ID == goid_2) %>% slice(1) %>%
#         pull(symbol) %>% unlist() %>% unique()
#     ))
#   )

ks1 <- ks %>%
  left_join(., biodom %>% select(goid_1 = GO_ID, n_gene1 = n_symbol), by = 'goid_1') %>%
  left_join(., biodom %>% select(goid_2 = GO_ID, n_gene2 = n_symbol), by = 'goid_2') %>%
  distinct() %>%
  filter(goid_1 != goid_2,
         kappa > 0.1,
         n_gene1 > 3,
         n_gene2 > 3)

nw <- biodom %>%
  select(n1 = GO_ID, n2 = Biodomain) %>%
  mutate(kappa = 1) %>%
  distinct() %>%
  bind_rows(
    ks1 %>% select(n1 = goid_1, n2 = goid_2, kappa),
    . )

g = igraph::graph_from_data_frame(nw, directed = F)

v = tibble( n = igraph::V(g)$name ) %>%
  left_join(., biodom %>% select(GO_ID, term = GOterm_Name, n_symbol), by = c('n'='GO_ID')) %>%
  left_join(., dom.cols, by = c('n'='domain')) %>%
  mutate(term = if_else(is.na(term), n, term),
         term_size = if_else( is.na(n_symbol), max(n_symbol, na.rm =T), n_symbol ),
         term_size = rank(term_size, ties.method = 'first'),
         term_size = term_size / max(term_size),
         term_size = if_else( n == term, term_size*3, term_size),
         color = case_when( is.na(color) ~ 'grey80', T ~ color)
  ) %>%
  filter(!duplicated(n))

for(i in 2:ncol(v)){
  igraph::vertex_attr(g, names(v)[i], index = igraph::V(g)) <- v %>% pull(i)
}

term_graph = g

##
# SEA-AD snRNA-seq data

seaad <- read_csv( synapser::synGet('syn51197803')$path, col_types = cols() ) %>%
  filter( mean_exp > 0 ) %>%
  mutate(
    group = factor(group, c( 'High',  'Intermediate','Low', 'Not AD', 'Dementia', 'No dementia')),
    broad1 = as.factor(broad) %>% as.numeric(),
    broad2 = case_when(group == 'Dementia' ~ broad1+.25,
                       group == 'No dementia' ~ broad1-.25,
                       group == 'High' ~ broad1+.25,
                       group == 'Intermediate' ~ broad1+.125,
                       group == 'Low' ~ broad1-.125,
                       group == 'Not AD' ~ broad1-.25),
    et = case_when( group %in% c('No dementia','Dementia') ~ 'dem',
                    group %in% c('Not AD','Low','Intermediate','High')~'path' ),
    label = NA_character_
  ) %>%
  arrange(desc(mean_exp), desc(fraction_expressed))

# Generate plots ----------------------------------------------------------

tep.tg = 'FCER1G'

# tg.list = c('APOE','ARHGEF2','C4A','CAPN2','CD44','CNN3','CTSH',
#             'DAG1','DDX1','DHX58','EPHX2','FCER1G','GPNMB','HTRA1','IFIH1',
#             'MDK','MSN','NDUFS2','PLEC','PRDX1','PRDX6','PTN','QPRT',
#             'RABEP1','SDC4','SFRP1','SMOC1','SNX32','STX4','SYK','TICAM1')

for(tep.tg in tg.list){
  
  ##
  # risk scores
  
  score.plot <- ggplot(scores, aes( value, name )) + 
    geom_violin(scale = 'width', aes(fill = name), 
                alpha = .3, color = 'grey50',
                draw_quantiles = c(.5,.9), trim = T)+ 
    geom_point(data = subset(scores, GeneName == tep.tg), 
               shape = 18, size = 4) + 
    viridis::scale_fill_viridis(discrete = T, guide ='none',
                                begin = 0, end = .8,
                                option = 'E', direction = -1) +
    theme(plot.title.position = 'plot')+
    labs(y='', x= 'score', 
         title = paste0('     ',tep.tg, ' TREAT-AD risk scores')
         , subtitle = paste0('               Overall rank #',
                             scores$rank[which(scores$GeneName==tep.tg)] %>% unique(),
                             ' (', scores$pctl[which(scores$GeneName==tep.tg)] %>% unique() %>% signif(3), 
                             ' percentile)'
                             # ,'\nlines: 50th & 90th quantiles'
                             )
         )
  # ggsave(paste0(tep.tg,'_risk_scores.pdf'), device = 'pdf', 
  #        width = 3.5, height = 2, units = 'in')
  # 3.4 in x 2 in
  
  ##
  # genetics
  
  gen.plot <-   ggplot(gen, aes(value, name)) +
    geom_violin(scale = 'width', aes(fill = name), 
                alpha = .3, color = 'grey50',
                draw_quantiles = c(.5,.9), trim = T)+
    geom_point(data = subset(gen, GeneName == tep.tg),
               shape = 18, size = 4)+
    viridis::scale_fill_viridis(discrete = T, guide ='none',
                                begin = 0, end = .8,
                                option = 'E', direction = -1) +
    theme(plot.title.position = 'plot')+
    labs(y='', x= 'score', 
         title = paste0('     ',tep.tg,' genetic risk'),
         subtitle = paste0('               Genetics rank #',
                           gen$rank[which(gen$GeneName==tep.tg)] %>% unique(),
                           ' (', gen$pctl[which(gen$GeneName==tep.tg)] %>% unique() %>% signif(3), 
                           ' percentile)'
                           # ,'\nlines: 50th & 90th quantiles'
                           )
    )
  # ggsave(paste0(tep.tg,'_genetics.pdf'), device = 'pdf', 
  #        width = 3.5, height = 2.7, units = 'in')
  
  ##
  # omics
  
  expr.plot <-     ggplot(omics, aes(x, TE)) + theme_bw() +
    geom_violin(fill = 'grey60', alpha = .3) +
    geom_point(
      data = subset(omics, GName == tep.tg & fdr < 0.05),
      size = 3, color = 'red', alpha = .6
    ) +
    geom_point(
      data = subset(omics, GName == tep.tg & fdr > 0.05),
      size = 3, color = 'grey65', alpha = .6
    ) +
    geom_hline(yintercept = 0, lty = 2) +
    # theme(plot.title=element_text(size=11, hjust = 0),plot.subtitle=element_text(size=9),
    #       plot.margin = margin(t = -2, r=5, b = 2, l=-2, unit = 'pt'),
    #       axis.title.x = element_text(size=9))+
    theme(plot.title.position = 'plot')+
    labs(x = '', 
         y = 'meta-analysis log fold change\nAD vs Control',
         title = paste0('     ',tep.tg, ' differential expression'),
         subtitle = '          red: FDR < 0.05') + 
    coord_flip()
  
  # ggsave(paste0(tep.tg,'_omics.pdf'), device = 'pdf', 
         # width = 3.5, height = 2.5, units = 'in')
  
  risk_plots = cowplot::plot_grid(score.plot, gen.plot, expr.plot, 
                                  ncol = 1, rel_heights = c(.8,1,.8), align = 'v')
  ggsave(plot = risk_plots, filename = paste0(tep.tg,'_risk_plots.pdf'), 
         device = cairo_pdf, width = 3.5, height = 6, units = 'in')
  
  ##
  # bd lolipop
  
  bd.plot <- biodom_genes %>% 
    filter(symbol == tep.tg) %>% 
    select(-GO_ID,-GOterm_Name, -n_symbol) %>%  distinct() %>% 
    full_join(., 
              biodom_genes %>% select(Biodomain, color, label) %>% distinct(), 
              by = c('Biodomain','color', 'label')) %>%
    filter(!is.na(Biodomain)) %>% 
    mutate( pct = if_else(is.na(pct), 0, pct),
            n_term = if_else(is.na(n_term), 0, as.double(n_term)),
            label = fct_reorder(label, pct) ) %>% 
    ggplot(., aes( pct, label )) + 
    geom_segment( aes(yend=label, xend=0), color = 'grey50') + 
    geom_point(alpha = .9, aes(size = n_term, fill = color), shape = 21) + 
    scale_size_continuous('# term', range=c(2,6))+
    scale_fill_identity() + expand_limits(x=0) + 
    theme(legend.position = 'top', plot.title.position = 'plot')+
    labs(y='', x = '% Biodomain GO terms', 
         title = paste0(tep.tg,' Biodomains')) 
  
  ##
  # bd term graph
  
  tg.terms <- biodom_genes %>% 
    filter(symbol == tep.tg) %>% 
    select(GO_ID, Biodomain) %>% 
    unlist() %>% 
    unique()
  
  # remove terms not annotated to term
  sub.g <- igraph::delete.vertices(
    term_graph, igraph::V(term_graph)[ 
      !(igraph::V(term_graph)$name %in% unique(tg.terms)) 
    ]
  )
  
  # remove edges below kappa score threshold
  kappa_thresh = 0.5
  sub.g = igraph::delete.edges(
    sub.g , 
    edges = igraph::E(sub.g)[ igraph::E(sub.g)$kappa < kappa_thresh ] )
  
  simpl <- igraph::simplify(sub.g, edge.attr.comb = c(kappa = 'max'))
  
  bd.nw <- ggraph(simpl, layout = 'fr') +
    geom_edge_link(aes(alpha = kappa), color = 'grey50') +
    geom_node_point( shape = 21, aes(size = n_symbol, fill = color, filter = is.na(.N()$abbr)) )+
    geom_node_point( shape = 21, aes(fill = color, color = color, filter = !is.na(.N()$abbr)), size = 9 )+
    geom_node_label( aes(label = abbr, filter = !is.na(.N()$abbr) ), 
                     size = 2.5, fontface = 'bold', alpha = .7 ) +
    scale_color_identity()+ scale_fill_identity()+ 
    scale_size_continuous('# genes in term', range = c(1,6))+
    scale_alpha_continuous("Cohen's kappa")+
    theme_graph(background = 'white')+
    theme(legend.position = 'bottom', legend.box = 'vertical')
  # +
    # labs(title = paste0(tep.tg, ' Biodomain GO term network'))

  
  bd_plots = cowplot::plot_grid(
    cowplot::plot_grid(NULL, bd.plot, NULL, ncol = 1, rel_heights = c(.1, 1,.1) ),
    bd.nw, nrow = 1, rel_widths = c(.9,1))
  
  ggsave(plot = bd_plots, filename = paste0(tep.tg,'_biodomain_plots.pdf'),
         device = cairo_pdf, width = 7.5, height = 6, units = 'in')
    
  ##
  # cell type graph
  
  # pos = position_jitter(width = 0.3, seed = 2)
  # 
  # seaad2 <- seaad %>% 
  #   filter(fraction_expressed > 0.05, mean_exp > 0) %>% 
  #   group_by(broad, broad1, broad2, et, group, gene) %>% 
  #   summarise(fxn = mean(fraction_expressed), mn = mean(mean_exp))
  
  tmp <- seaad2 %>% 
  filter(gene == tep.tg, 
         et == 'dem' 
         # , fxn > .1
         ) %>% 
    mutate(broad1 = broad2)
  
  cell_types <- ggplot(tmp, aes(broad1, mn)) +
    geom_violin(data = subset(seaad2, et == 'dem' 
                              # & fraction_expressed > 0.1
                              ), 
                aes(fill = broad), draw_quantiles = c(.5,.9), 
                color = 'grey65', alpha = .3, scale = 'width') +
    theme(legend.position = 'bottom', legend.box = 'vert')+
    viridis::scale_fill_viridis(discrete = T, guide ='none',
                                begin = 0, end = .8,option = 'D') +
    scale_x_continuous(breaks = 1:8,
                       labels = seaad2$broad %>% unique() %>% sort()) +
    labs(x = '', y = 'mean expression, ln(UP10K+1)' )+
    theme(axis.title.x = element_blank())+
    scale_color_discrete(name = 'cell subset')+
    scale_size_continuous('fraction of cells expressing', range = c(.5,3), limits = c(0,1) )+
    geom_line(data = tmp, #subset(z1, is.na(label) ),
              aes(group = broad),
              color = 'grey35') +
    geom_point(data = tmp, #subset(z1, is.na(label) ),
               aes(
                 group = broad,
                 color = group,
                   size = fxn)
               ,shape = 16, alpha = .8)+
    labs(title = paste0(tep.tg, ' SEA-AD cell-type expression'))
  
  ggsave(plot = cell_types, filename = paste0(tep.tg,'_celltypes_plots.pdf'), 
         device = cairo_pdf, width = 7.5, height = 4, units = 'in')
  
  # # Plot dimensions: 5.84 in x 7.10 in
  # p = cowplot::plot_grid(
  #   cowplot::plot_grid( 
  #     cowplot::plot_grid(score.plot, expr.plot, nrow=2, align = 'hv', rel_heights = c(1,.8)), 
  #                       bd.plot, nrow = 1, rel_widths = c(.875,1)),
  #   allen.plot, nrow = 2, rel_heights = c(1,.5), rel_widths = c(1,.7))
  
  # View(tad.scores %>% filter(GeneName == tep.tg))
  
  t = tad.scores %>% filter(GeneName == tep.tg) %>% 
    select(Overall, Overall_rank, Overall_pctl, GeneticsScore, OmicsScore, LiteratureScore) %>% 
    mutate(across(everything(), ~ signif(.x, digits = 4))) %>% 
    gridExtra::tableGrob(rows = NULL, 
                         theme = gridExtra::ttheme_default(base_size = 9))
  
  ggsave(
    filename = paste0(here::here(), '/results/TEP_report_plots/', tep.tg,'.png'),
    plot = cowplot::plot_grid(gridExtra::grid.arrange(t), p, ncol = 1, rel_heights = c(.1,.9)),
    width = 6.5, height = 8
  )
  
  }


