# LOAD DATASETS ################################################################
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source(here::here("src","setup_idr_analysis.r"),echo = F,chdir = T)
library(Biostrings)
library(hablar)

# BUILD HUMAN IDR DATASET ######################################################
# Filter consensus disordered regions (at least 50% agreement of predictors)
hs_mobidb = get_human_mobidb() 
hs_mobidb_merged = hs_mobidb |> dplyr::rename(S=START,E=END) |> merge_mobidb(gap_min=3L)
hs_diso =  dplyr::filter(hs_mobidb, is_uniref & feature == 'disorder' & source %in% 'th_50') |>
           dplyr::rename(IDR_id=feature_id) |> 
           arrange(IDR_id)

cluster <- new_cluster(14)
cluster_library(cluster, "dplyr")
hs_idr = add_idr_sequence(hs_diso, cl=cluster) |>  arrange(IDR_id)
hs_peptides = get_peptstats(hs_idr,cl=cluster) |>  arrange(IDR_id)
rm(cluster)

hs_idr_info = hs_idr |> dplyr::select(acc,IDR_id,feature_len,feature_seq)
# Compute amino acid counts in human IDRs
hs_aacount = get_aa_count(hs_idr,'feature_seq')
# Compute molecular features based on human IDR amino acid counts
hs_aafreq = get_aa_freq(hs_aacount)
hs_charge  = get_aa_charge(hs_aacount)
hs_scores  = get_aa_scores(hs_aacount)
hs_topfreq = get_aa_topfreq(hs_aacount)
hs_foldchange = get_aa_foldchange(hs_aacount)
hs_topfc = get_aa_topfc(hs_foldchange)

# Combine molecular features of human IDRs
hs_features = bind_cols(hs_aacount, hs_aafreq, hs_charge,hs_scores,
                        hs_peptides, hs_foldchange,hs_topfreq, hs_topfc)

# Phase-separation proteins with uniprot infos
hs_ps = get_human_phase_separation(mobidb=hs_mobidb) |>
        # merging phase-separating regions by positions
        merge_ps_region() |> 
        dplyr::select(-GN,-LEN,-PN,-KW,-PE,-Reviewed) |>
        group_by(acc) |>
        add_count(acc,name='PS_n')
  

# Dataset for human IDRs with molecular features + phase-separating regions
hs_idr_ps = left_join(hs_idr,hs_features,by=c("IDR_id")) |>
            ungroup() |>
            left_join(hs_ps, by = 'acc') # |>
            #rowwise() |>
            #mutate(has_ps = has_overlap(START,END,PS_START,PS_END))
                   #ps_overlap = which_overlap(START,END,PS_START,PS_END))
# Subset human IDR based on length (at least 35 residues)
# Rename and reorganize the columns
HS_IDR = inner_join(hs_uni,hs_idr_ps,by=c('AC'='acc')) %>%
  filter(feature_len>35) %>%
  mutate(PS_n = replace_na(PS_n,0), has_PS = PS_n > 0 ) %>%
  dplyr::select(-DB,-id_cdna,-OS,-OX,-length,-PE,-SV) %>%
  dplyr::rename(IDR_len=feature_len,IDR_seq=feature_seq,
                IDR_frac=content_fraction,IDR_count=content_count) %>%
  relocate(ncbi_taxid,AC,ID,GN,ensp,NAME, is_uniref, uniprot_seq,
           uniprot_len,IDR_frac,IDR_count, tot,
           IDR_id,IDR_len,START,END,source,feature,evidence,IDR_seq,
           colnames(hs_charge),starts_with('peptide_'),
           sum_aa, all_of(AA3),"X", all_of(names(AA_PROP)),
           starts_with("fr_"),starts_with("fc_"), 
           starts_with("frtop_"), starts_with("fctop_"),
           aggrescan:wimleywhite,
           starts_with("mean_"),
           region,has_PS,starts_with("PS_")
           ) %>%
  ungroup() %>% distinct()

summary(HS_IDR)
head(HS_IDR)
write_tsv(HS_IDR,file=here::here('data','HUMAN_MOBIDB_FEATURES.tsv'))

# BUILD ATAR IDR DATASET #######################################################
# LAST UPDATE OF INPUT DATA *FEB 2023*
ATAR_CANDIDATES = rio::import(here::here('data',"ATAR_candidates_table.xlsx"))
df_atar = left_join(ATAR_CANDIDATES, df_hs_seq, by=c("UNIPROT"='uniprot_id')) |>
          dplyr::select(PROTEIN,UNIPROT,uniprot_seq,atar_sequence=SEQUENCE,
                        START="Start Position",END="End Position", LEN = SEQ_LENGTH) |>
          mutate( IDR_id  = paste0(UNIPROT,"_",START,"..",END), atar_len = nchar(atar_sequence) ) |>
          arrange(IDR_id)

cluster <- new_cluster(10)
atar_idr = add_idr_sequence(df_atar,cluster) |> arrange(IDR_id)
atar_peptides = get_peptstats(atar_idr,col_sequence = 'atar_sequence',cl=cluster) |> arrange(IDR_id)
rm(cluster)

# Compute amino acid counts in candidate IDRs
atar_aacount = get_aa_count(atar_idr,'atar_sequence')
# Compute molecular features based on candidate IDR amino acid counts
atar_aafreq = get_aa_freq(atar_aacount,col_aa = get.AAA())
atar_charge  = get_aa_charge(atar_aacount)
atar_scores  = get_aa_scores(atar_aacount)
atar_topfreq = get_aa_topfreq(atar_aacount)

hs_naa = sum(rowSums(hs_aacount[,c(get.AAA(),"X")]))
hs_idr_aafreq = hs_aacount[,c(get.AAA(),"X")] |> colSums() |> magrittr::divide_by(hs_naa)

atar_foldchange = get_aa_foldchange(atar_aacount,ref_aa_freq=hs_idr_aafreq)
atar_topfc = get_aa_topfc(atar_foldchange,col_aa = paste0("fc_",get.AAA()))

# Combine molecular features of candidate IDRs
atar_features = bind_cols(atar_aacount, atar_aafreq, atar_charge, atar_scores,
                        atar_peptides, atar_foldchange, atar_topfreq, atar_topfc)

## QUALITY CONTROL OF COMPUTED FEATURES ########################################
atar_seq = atar_idr$atar_sequence %>% str_split("")

# checking stickiness score
sticky = get.stickiness()
sum_stickiness = map_dbl(atar_seq, ~ sum(sticky[.x],na.rm=T) )
avg_stickiness = map_dbl(atar_seq, ~ mean(sticky[.x],na.rm=T) )

all.equal(sum_stickiness,atar_features$stickiness)
all.equal(avg_stickiness,atar_features$mean_stickiness)

# Atar's candidate IDRs with molecular features + phase-separating regions
ATAR_IDR = left_join(atar_idr,atar_features,by="IDR_id") |> 
  group_by(UNIPROT,LEN) |>
  #left_join(hs_ps, by = c('UNIPROT'='acc')) |>
  left_join(hs_diso[,c('acc','content_fraction','content_count','length')], by=c('UNIPROT'='acc')) |>
  mutate(from_atar=TRUE) |>
  distinct() |>
  dplyr::rename(AC=UNIPROT, IDR_frac = content_fraction, IDR_count = content_count,
                IDR_len = atar_len) |>
  relocate(AC,PROTEIN, LEN,
           IDR_id,START,END, IDR_len, atar_sequence, 
           colnames(hs_charge),
           starts_with('pep_'),starts_with("fr_"),starts_with("fc_"), aggrescan:wimleywhite,
           starts_with("mean_"))

#write_tsv(ATAR_IDR, file = here::here('data','ATAR-IDR-FEATURES.tsv'))
save.image(file = here::here('data', 'IDR-features-data.rdata'))
#load(here::here('data', 'IDR-features-data.rdata'))
#### IDR FEATURES ####

##### FEATURES COMPARISON #####

# Check correlogram of numeric features to remove redundancy
# (high absoluter correlation == redundant features)

col_mobidb = c('IDR_len','IDR_count','IDR_frac','uniprot_len')
numeric_features = c( colnames(hs_aacount), colnames(hs_aafreq),
                      colnames(hs_charge), colnames(hs_peptides),
                      colnames(hs_scores), colnames(hs_foldchange),
                      col_mobidb)

df_all_features = bind_rows(HS_IDR,ATAR_IDR) %>% 
  dplyr::select(PROTEIN,AC,IDR_id,START,END, all_of(numeric_features), has_PS, from_atar) %>%
  mutate(across(all_of(c('IDR_len',"IDR_count",'uniprot_len')), .fns=~log10(.x), .names="{col}.log10"), 
         from_atar = replace_na(from_atar,FALSE)) %>%
  dplyr::select(-all_of(c('IDR_len',"IDR_count",'uniprot_len'))) %>%
  distinct() 

p_all = make_features_correlation(df_all_features)
ggsave(p_all,path=here::here("plots"),
       filename ='correlation-human-idr-all-features.png',
       height=12,width=12,bg='white')
ggsave(p_all,path=here::here("plots"),
       filename ='correlation-human-idr-all-features.pdf',
       height=12,width=12,bg='white')

##### FEATURES SELECTION ######
# AA Frequencies + by chemical group
# AA Foldchanges
# AA scores (stickiness, roseman aggrescan)
# Peptide stats (netcharge, PI, IDR_frac)

features_to_use = c(paste0("fr_",c(get.AAA(),"X")),
                    colnames(hs_foldchange), 
                    c("netcharge_residue","charge_asymetry","peptide_PI"),
                    c("mean_stickiness","mean_roseman","mean_aggrescan"),
                    "IDR_len","IDR_frac","IDR_count")

df_features = bind_rows(HS_IDR,ATAR_IDR) %>% 
          dplyr::select(PROTEIN,AC,IDR_id,START,END, all_of(features_to_use) ,from_atar) %>%
          mutate(across(all_of(c('IDR_len',"IDR_count")), .fns=~log10(.x), .names="{col}.log10"), 
          from_atar = replace_na(from_atar,FALSE)) %>%
          dplyr::select(-all_of(c('IDR_len',"IDR_count"))) %>%
          distinct() 

df_num = df_features %>%
         dplyr::select(where(~ is.numeric(.x))) %>%
         dplyr::select(-START,-END)

df_info = df_features %>% dplyr::select( -colnames(df_num) )

df_scaled = bind_cols(df_info,as_tibble(scale(df_num))) %>% 
                        dplyr::rename_with(.cols=colnames(df_num),.fn = xxS, sx='scaled',s='.')
# Check correlogram of selected features
p_used = make_features_correlation(df_features)
ggsave(p_used,path=here::here("plots"),
       filename ='correlation-human-idr-selected-features.png',
       height=12,width=12,bg='white')
ggsave(p_used,path=here::here("plots"),
       filename ='correlation-human-idr-selected-features.pdf',
       height=12,width=12,bg='white')

# IDR UMAP #####################################################################
make_umap(df_scaled, K = 30, seed = 142, is_scaled = T)
make_umap(df_features, K = 30, seed = 142, is_scaled = F)

k_neighbors = c(10,20,30,40,50,100)
UMAP_SCALED = lapply(k_neighbors, function(x){ make_umap(df_scaled, K = x, seed = 142, is_scaled = T)} ) |>
              patchwork::wrap_plots(nrow = 2)

# save umap in PNG/PDF
plot_name = 'umap-idr-human-scaled'
ggsave(UMAP_SCALED, path=here::here("plots"),
       filename = paste0(plot_name,'.png'),scale=1.3,
       device = 'png', height=12, width=12, bg='white')
ggsave(UMAP_SCALED, path=here::here("plots"),
       filename = paste0(plot_name,'.pdf'), scale=1.3,
       device = 'pdf', height=12, width=12, bg='white')

UMAP_UNSCALED = lapply(k_neighbors, function(x){ make_umap(df_features, K = x, seed = 142, is_scale = F)} ) |>
  patchwork::wrap_plots(nrow = 2)

plot_name = 'umap-idr-human-unscaled'
ggsave(UMAP_UNSCALED, path=here::here("plots"),
       filename = paste0(plot_name,'.png'),scale=1.3,
       device = 'png', height=12, width=12, bg='white')
ggsave(UMAP_UNSCALED, path=here::here("plots"),
       filename = paste0(plot_name,'.pdf'), scale=1.3,
       device = 'pdf', height=12, width=12, bg='white')


# CORRELATION WITH CSAT ########################################################
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


features_spearman = map(features, ~spearman.toplot(df_feature_csat$csat_conc_log10,df_feature_csat[[.x]])) %>% bind_rows()
features_spearman$name = features
features_spearman$is_umap = features_spearman$name %in% features_to_use
features_spearman$signif =symnum(features_spearman$pv, corr = FALSE, cutpoints = c(0,  .001,.01,.05, .1, 1), symbols = c("(***)","(**)","(*)","(.)"," "))


plot_csat_pearson = features_pearson %>% 
  filter( abs(r) > 0.3 ) %>% 
  ggplot(aes(y=reorder(name,-r), x=r, label=reorder(name,-r))) + 
    geom_col(orientation='y', width=0.7, linewidth=0.5, aes(col=is_umap)) + 
    ggfittext::geom_bar_text() + 
    geom_text(mapping = aes(label=sprintf(" %.2f ",r)),col='gray',hjust='outward', size=3) + 
    geom_text(mapping = aes(label=sprintf(" %.2f %s",pv,signif)),x=-Inf,col='red',hjust='inward', size=3) +   
    ggeasy::easy_remove_axes('y') +
    scale_color_manual(values=c("TRUE"='dodgerblue',"FALSE"='transparent')) +
    xlab("Pearson correlation (Csat vs. feature)") + xlim(-1,1)

plot_csat_spearman = features_spearman %>% 
  filter( abs(r) > 0.3 ) %>%
  ggplot(aes(y=reorder(name,-r), x=r, label=reorder(name,-r))) + 
  geom_col(orientation='y', width=0.7, linewidth=0.5, aes(col=is_umap)) + 
  ggfittext::geom_bar_text() + 
  geom_text(mapping = aes(label=sprintf(" %.2f ",r)),col='gray',hjust='outward', size=3) + 
  geom_text(mapping = aes(label=sprintf(" %.2f %s",pv,signif)),x=-Inf,col='red',hjust='inward', size=3) +   
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

features_cor = bind_rows(features_pearson,features_spearman) %>% mutate(cor_type=subname(method,"'")) %>%
               pivot_wider(id_cols=name, names_from = cor_type, values_from=c(r,pv,p,signif) ) %>% 
               mutate(is_signif = (pv_Pearson<=0.05 | pv_Spearman <=0.05),
                      toshow = sprintf("r=%.2f p=%.2f %s\ns=%.2f p=%.2f %s",r_Pearson,pv_Pearson,signif_Pearson,r_Spearman,pv_Spearman,signif_Spearman))
features_plot = pivot_longer(df_feature_csat, cols = -csat_conc_log10) %>% 
                left_join(features_cor) %>%
  arrange(csat_conc_log10, desc(r_Pearson)) %>% 
  ggplot(aes(x=value,y=csat_conc_log10,col=is_signif)) + 
  geom_point() + 
  geom_line(aes(group=name),linewidth=0.1) + 
  geom_smooth(method='lm') +
  geom_text(mapping=aes(label=toshow),
            x=-Inf,y=Inf, hjust='inward',vjust='inward',size=2.5, check_overlap = T) +
  scale_color_manual(values=c("TRUE"='red',"FALSE"="black")) +
  scale_y_log10() +
  facet_wrap(~reorder(name,abs(r_Pearson)), scales = 'free_x',nrow = 8,strip.position = 'right') + 
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
