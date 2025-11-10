# --- 1. LOAD SETUP, HELPERS, AND QC ---
source(here::here("src","setup_main.r"),echo = F,chdir = T)
source(here::here("src","helpers_analysis.r"),echo = F,chdir = T)

# --- 1b. QC AND LOAD DATA ---
data_list <- load_csat_data()
ATAR_IDR <- data_list$ATAR_IDR
#ATAR_IDR = read_tsv("~/Downloads/ATAR-UMAP-DATA.tsv")
# .RData objects (hs_aacount, df_info, features_to_use, etc.) 
# are loaded into GlobalEnv by the function

# --- 2. CORRELATION WITH CSAT ---
# The saturation concentration values are taken from the micDROP-red system (obtained from Figure S6)
.info$log("Running analysis: C-SAT Correlation...")
CSAT = tribble( ~PROTEIN, ~csat_conc,
                "Ddx4",752,
                "DYRK3", 2275,
                "FUS", 780,
                "hnRNPA1", 271,
                "RBM-14", 148,
                "TAF-15", 306,
                "TDP-43", 442,
                "UBQ2", 81,
                "HSPB8", 423,
                "Erα", 2242) %>% 
  mutate(csat_conc_log10 = log10(csat_conc))

csat_features = left_join(CSAT,ATAR_IDR) %>% 
  dplyr::select(colnames(df_info),where(is.numeric),
                csat_conc,csat_conc_log10,
                -tot, -X, -fr_X, -fc_X, -peptide_len, -sum_aa, -LEN)
df_feature_csat <- csat_features %>% 
  filter(!is.na(csat_conc)) %>%
  dplyr::select(-PROTEIN,-AC,-IDR_id,-START,-END,-from_atar)

features = df_feature_csat %>% dplyr::select(-starts_with('csat')) %>% colnames

features_pearson = map(features, ~pearson.toplot(df_feature_csat$csat_conc_log10,df_feature_csat[[.x]])) %>% bind_rows()
features_pearson$name = features
features_pearson$is_umap = features_pearson$name %in% features_to_use
features_pearson$signif = symnum(features_pearson$pv, corr = FALSE, cutpoints = c(0,  .001,.01,.05, .1, 1), symbols = c("(***)","(**)","(*)","(.)"," "))
features_pearson$p.adj = p.adjust(features_pearson$pv, method="BH")


features_spearman = map(features, ~spearman.toplot(df_feature_csat$csat_conc_log10,df_feature_csat[[.x]])) %>% bind_rows()
features_spearman$name = features
features_spearman$is_umap = features_spearman$name %in% features_to_use
features_spearman$signif =symnum(features_spearman$pv, corr = FALSE, cutpoints = c(0,  .001,.01,.05, .1, 1), symbols = c("(***)","(**)","(*)","(.)"," "))
features_spearman$p.adj =  p.adjust(features_spearman$pv, method="BH")

plot_csat_pearson = features_pearson %>% 
  filter( abs(r) > 0.3 ) %>% 
  ggplot(aes(y=reorder(name,-r), x=r, label=reorder(name,-r))) + 
  geom_col(orientation='y', width=0.7, linewidth=0.5, aes(col=is_umap)) + 
  ggfittext::geom_bar_text() + 
  geom_text(mapping = aes(label=sprintf(" %.2f ",r)),col='gray',hjust='outward', size=3,check_overlap = T) + 
  geom_text(mapping = aes(label=sprintf(" %.2f",p.adj)),x=-1,hjust=0, size=2.5) +   
  geom_text(mapping = aes(label=sprintf(" %.2f %s",pv,signif)),x=-0.9,col='red',hjust='left', size=3) +   
  ggeasy::easy_remove_axes('y') +
  scale_color_manual(values=c("TRUE"='dodgerblue',"FALSE"='transparent')) +
  xlab("Pearson correlation (Csat vs. feature)") + xlim(-1,1)

plot_csat_spearman = features_spearman %>% 
  filter( abs(r) > 0.3 ) %>%
  ggplot(aes(y=reorder(name,-r), x=r, label=reorder(name,-r))) + 
  geom_col(orientation='y', width=0.7, linewidth=0.5, aes(col=is_umap)) + 
  ggfittext::geom_bar_text() + 
  geom_text(mapping = aes(label=sprintf(" %.2f ",r)),col='gray',hjust='outward', size=3,check_overlap = T) + 
  geom_text(mapping = aes(label=sprintf(" %.2f",p.adj)),x=-1,hjust=0, size=2.5) +   
  geom_text(mapping = aes(label=sprintf(" %.2f %s",pv,signif)),x=-0.9,col='red',hjust='left', size=3) +   
  ggeasy::easy_remove_axes('y') +
  scale_color_manual(values=c("TRUE"='dodgerblue',"FALSE"='transparent')) +
  xlab("Spearman correlation (Csat vs. feature)") + xlim(-1,1)

plot_csat_cor = patchwork::wrap_plots(plot_csat_spearman,plot_csat_pearson)

ggsave(plot_csat_cor, path=here::here("plots"),
       filename = 'csat_correlation.pdf', scale=1.2,
       device = 'pdf', height=12, width=16, bg='white')
ggsave(plot_csat_cor, path=here::here("plots"),
       filename = 'csat_correlation.png', scale=1.2,
       device = 'png', height=12, width=16, bg='white')

features_cor = bind_rows(features_pearson,features_spearman) %>% 
  pivot_wider(id_cols=name, names_from = cor_type, values_from=c(r,pv,p,signif) ) %>% 
  mutate(is_signif = (pv_pearson<=0.05 | pv_spearman <=0.05),
         toshow = sprintf("r=%.2f p=%.2f %s\ns=%.2f p=%.2f %s",r_pearson,pv_pearson,signif_pearson,r_spearman,pv_spearman,signif_spearman))

features_plot = pivot_longer(df_feature_csat, cols = -csat_conc_log10) %>% 
  left_join(features_cor) %>%
  arrange(csat_conc_log10, desc(r_pearson)) %>% 
  ggplot(aes(x=value,y=csat_conc_log10,col=is_signif)) + 
  geom_point() + 
  geom_line(aes(group=name),linewidth=0.1) + 
  geom_smooth(method='lm') +
  geom_text(mapping=aes(label=toshow),
            x=-Inf,y=Inf, hjust='inward',vjust='inward',size=2.5, check_overlap = T) +
  scale_color_manual(values=c("TRUE"='red',"FALSE"="black")) +
  scale_y_log10() +
  facet_wrap(~reorder(name,abs(r_pearson)), scales = 'free_x',nrow = 8,strip.position = 'right') + 
  theme(strip.text.y = element_text(size = 7),panel.spacing = unit(0.1, "cm"),
        axis.text.x = element_text(size=5),axis.text.y = element_text(size=5))


ggsave(features_plot, path=here::here("plots"),
       filename = 'csat_vs_features.pdf', scale=1.2,
       device = 'pdf', height=12, width=12, bg='white')
ggsave(features_plot, path=here::here("plots"),
       filename = 'csat_vs_features.png', scale=1.2,
       device = 'png', height=12, width=12, bg='white')


# save UMAP features + CSAT 
right_join(CSAT,df_features) %>% 
  filter(from_atar) %>%
  write_tsv(file = here::here('data','ATAR-UMAP-DATA.tsv'))


cor_Vsticky = features_cor %>% filter(name=='voronoi_stickiness')


plot_vsticky = csat_features %>%
  ggplot(aes(y=csat_conc_log10, x=voronoi_stickiness)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label=PROTEIN)) +
  geom_smooth(method='lm') + 
  geom_text(data=cor_Vsticky,aes(label=toshow),x=Inf,y=Inf,hjust='inward',vjust='inward') +
  xlab('Voronoi Stickiness (sum)') + ylab('IDR saturating concentration (log10)') +
  theme(axis.title.y = element_text(size=8,family = 'Helvetica'),
        axis.title.x = element_text(size=8,family = 'Helvetica'), 
        axis.text.y = element_text(size=8,family = 'Helvetica'),
        axis.text.x = element_text(size=8,family = 'Helvetica')) 


ggsave(plot_vsticky, path=here::here("plots"),
       filename = 'csat_voronoi_stickiness.png', scale=1,
       device = 'png', height=5, width=5, bg='white')

ggsave(plot_vsticky, path=here::here("plots"),
       filename = 'csat_voronoi_stickiness.pdf', scale=1,
       device = 'pdf', height=5, width=5, bg='white')

.succ$log("C-SAT analysis complete. Plots saved to plots/ directory.")


