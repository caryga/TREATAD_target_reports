library(tidyverse)

synapser::synLogin()

perc.rank <- function(x) trunc(rank(x))/length(x)

theme_set(theme_bw())

# Gather data -------------------------------------------------------------

##
# risk scores

scores <- synapser::synTableQuery('select * from syn25575156')$filepath %>% 
  read_csv() %>%
  mutate(pctl = perc.rank(Overall)*100) %>%
  select(ENSG, GeneName, Overall, rank= Overall_rank, pctl, 
         Omics = OmicsScore, Genetics = GeneticsScore
         ) %>%
  pivot_longer(cols = c(Overall, Omics, Genetics)) %>%
  mutate(
    name = fct_relevel(name, c('Omics', 'Genetics','Overall')),
    value = case_when( value == 0 ~ NA_real_, T ~ value)
  ) %>%
  arrange(name)

##
# genetics scores

gen <- synapser::synTableQuery('SELECT * FROM syn26844312')$filepath %>% 
  read_csv() %>%
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
# SEA-AD snRNA-seq data

seaad <- read_csv( synapser::synGet('syn51197803')$path, col_types = cols() ) %>%
  filter( 
    mean_exp > 0, 
    group %in% c('Dementia', 'No dementia')) %>%
  mutate(
    group = factor(group, c('Dementia', 'No dementia')),
    broad1 = as.factor(broad) %>% as.numeric(),
    broad2 = case_when(group == 'Dementia' ~ broad1+.25,
                       group == 'No dementia' ~ broad1-.25),
    et = case_when( group %in% c('No dementia','Dementia') ~ 'dem'),
    label = NA_character_
  ) %>%
  arrange(desc(mean_exp), desc(fraction_expressed))

# Generate plots ----------------------------------------------------------

tep.tg = 'FCER1G'

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
                           )
       )

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
                         )
  )

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
  theme(plot.title.position = 'plot')+
  labs(x = '', 
       y = 'meta-analysis log fold change\nAD vs Control',
       title = paste0('     ',tep.tg, ' differential expression'),
       subtitle = '          red: FDR < 0.05') + 
  coord_flip()



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

ggsave(plot = bd.plot, filename = paste0(tep.tg,'_biodomain_plots.pdf'),
       device = cairo_pdf, width = 7.5, height = 6, units = 'in')
  
##
# cell type graph

seaad2 <- seaad %>%
  filter(fraction_expressed > 0.05, mean_exp > 0) %>%
  group_by(broad, broad1, broad2, et, group, gene) %>%
  summarise(fxn = mean(fraction_expressed), mn = mean(mean_exp))

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


