# LOAD DATASETS ################################################################
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source(here::here("setup_idr_analysis.r"),echo = F,chdir = T)
library(Biostrings)
library(hablar)

# BUILD HUMAN IDR DATASET ######################################################
# Filter consensus disordered regions (at least 50% agreement of predictors)
hs_mobidb = get_human_mobidb()
hs_diso =  dplyr::filter(hs_mobidb, is_uniref & feature == 'disorder' & source %in% 'th_50') |>
           dplyr::rename(IDR_id=feature_id)

cluster <- new_cluster(14)
cluster_library(cluster, "dplyr")
hs_idr = add_idr_sequence(hs_diso, cl=cluster)
hs_peptides = get_peptstats(hs_idr,cl=cluster)
rm(cluster)

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
  dplyr::select(-DB,-id_cdna,-OS,-OX,-length,-PE,-SV) %>%
  dplyr::rename(IDR_len=feature_len,IDR_seq=feature_seq,
                IDR_frac=content_fraction,IDR_count=content_count) %>%
  relocate(ncbi_taxid,AC,ID,GN,ensp,NAME,
           uniprot_len,IDR_frac,IDR_count,
           IDR_id,IDR_len,START,END,source,feature,IDR_seq,
           positive,negative,netcharge,
           starts_with('pep_'),starts_with("fr_"),starts_with("fc_"), 
           aggrescan:wimleywhite,
           starts_with("mean_")) %>%
  ungroup() %>% distinct()

summary(HS_IDR)
#write_tsv(HS_IDR,file=here::here('data','HUMAN_MOBIDB_FEATURES.tsv'))

# BUILD ATAR IDR DATASET #######################################################
# LAST UPDATE OF INPUT DATA *FEB 2023*
ATAR_CANDIDATES = rio::import(here::here('data',"ATAR_candidates_table.xlsx"))
df_atar = left_join(ATAR_CANDIDATES, df_hs_seq, by=c("UNIPROT"='uniprot_id')) |>
          dplyr::select(PROTEIN,UNIPROT,uniprot_seq,atar_sequence=SEQUENCE,
                        START="Start Position",END="End Position", LEN = SEQ_LENGTH) |>
          mutate( IDR_id  = paste0(UNIPROT,"_",START,"..",END), atar_len = nchar(atar_sequence) )

cluster <- new_cluster(10)
atar_idr = add_idr_sequence(df_atar,cluster)
atar_peptides = get_peptstats(atar_idr,cl=cluster)
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

# Atar's candidate IDRs with molecular features + phase-separating regions
ATAR_IDR = left_join(atar_idr,atar_features,by="IDR_id") |> 
  group_by(UNIPROT,LEN) |>
  #left_join(hs_ps, by = c('UNIPROT'='acc')) |>
  left_join(hs_diso[,c('acc','content_fraction','content_count','length')], by=c('UNIPROT'='acc')) |>
  mutate(from_atar=TRUE) |>
  distinct() |>
  dplyr::rename(AC=UNIPROT, IDR_frac = content_fraction, IDR_len = content_count) |>
  relocate(AC,PROTEIN, LEN,
           IDR_id,START,END, atar_len, atar_sequence, 
           positive,negative,netcharge,
           starts_with('pep_'),starts_with("fr_"),starts_with("fc_"), aggrescan:wimleywhite,
           starts_with("mean_"))

#write_tsv(ATAR_IDR, file = here::here('data','ATAR-IDR-FEATURES.tsv'))
#save.image(file = here::here('data', 'IDR-features-data.rdata'))


# IDR UMAP #####################################################################
# IDR FEATURES TO USE
# AA Frequencies + by chemical group
# AA Foldchanges
# AA scores (stickiness, roseman aggrescan)
# Peptide stats (average MW, Netcharge, PI, IDR_frac)

features_to_use = c(paste0("fr_",c(get.AAA(),"X")),
                    colnames(hs_foldchange), 
                    c("mean_stickiness","mean_roseman","mean_aggrescan"),
                    c("peptide_mw_avg","peptide_netcharge","peptide_PI"),
                    "IDR_len","IDR_frac")

df_features = bind_rows(HS_IDR,ATAR_IDR) %>% 
          dplyr::select(PROTEIN,AC,IDR_id,START,END,features_to_use,from_atar) %>%
          mutate(IDR_len.log10 = log10(IDR_len), from_atar = replace_na(from_atar,FALSE)) %>%
          dplyr::select(-IDR_len) %>%
          distinct() 

### Select features for umap
# all numeric
df_features_num = df_features %>%
              dplyr::select(where(~ is.numeric(.x))) %>%
              dplyr::select(-START,-END)

# Check correlogram of numeric features to remove redundancy
# (high absoluter correlation == redundant features)
cor_features = ggcorrplot::ggcorrplot(cor(df_features_num),outline.color = 'transparent', tl.cex = 8, tl.srt = 70)
ggsave(cor_features,filename = here::here('plots','correlation-human-idr-features.png'),height=12,width=12,bg='white')

# Get non-features columns (protein id, IDR boundaries...)
df_info = df_features %>% dplyr::select( -colnames(df_features_num) )

K=20
library(umap)
set.seed(142)
umap.config = umap.defaults
umap.config$n_neighbors = K
# Compute umap based on scaled features
idr_features_map <- df_features_num %>% umap::umap(seed = 142, config=umap.config)
#df_features_num %>% t() %>% scale() %>% t() %>% umap::umap(seed = 142, config=umap.config)

# mark outliers on umap coordinates (in 1/99% percentiles or outside the interval defined below for X1/X2)
df_umap = idr_features_map$layout %>% magrittr::set_colnames(c('X1','X2')) %>%
  bind_cols(idr_features_map$data %>% as_tibble() %>% dplyr::rename_with(.fn=xxS, s=".", sx='scaled')) %>%
  bind_cols(df_features_num) %>%
  bind_cols(df_info) %>%
  mutate(outliers_x1 = !between(percent_rank(X1),0.001,0.999),
         outliers_x2 = !between(percent_rank(X2),0.001,0.999),
         outliers = !between(X1,-8,8) | !between(X2,-6,6))
summary(df_umap[,c('X1','X2')])
#library(ggtrace)
library(ggalt)
library(ggiraph)
library(ggforce)
umap_data = df_umap %>% 
            filter( !outliers )  %>%
            mutate(idr_lab=paste0(PROTEIN," ",START,"-",END))

summary(umap_data[,c('X1','X2')])
# Plot the umap
UMAP = ggplot(data=umap_data ,aes(x = X1,y = X2)) +
  # all idrs
  geom_point(size=1.5,shape=16,color='gray',alpha=1) +
  # circle PS idr
  #geom_point(data=subset(umap_data,PS_overlap ), size=3, shape=21, color='black',stroke=0.5) +

  # highlight atar idr
  geom_point(data=subset(umap_data,from_atar), aes(color=PROTEIN), shape=16, size=4,alpha=0.9) +
  # geom_text_repel(data=atar_idr,aes(label=idr_lab,color=PROTEIN),
  #                 size=4,fontface='bold', force=5,  max.overlaps=15,force_pull = -1,
  #                 seed=291222) +

  # highlight idr in atar's protein (not overlapping)
  #geom_point(data=atar_neighbor, aes(color=PROTEIN), fill='white',shape=21, stroke=1, size=4, alpha=0.7) +
  ggrepel::geom_text_repel(data=subset(umap_data,from_atar),aes(label=idr_lab,color=PROTEIN),
                 size=3,fontface='italic', force=5, max.overlaps=15, force_pull =0,
                seed=291222) +

  # graphical parameters
  labs(x = "UMAP1", y = "UMAP2", subtitle="")+
  see::scale_color_metro(palette = 'rainbow', discrete = T) +
  theme(legend.position="bottom") + #theme_blackboard()+
  ggeasy::easy_text_size(c("axis.title","axis.text.x", "axis.text.y"), size = 20)+
  ggeasy::easy_remove_legend() +
  ggeasy::easy_remove_x_axis('text') + ggeasy::easy_remove_y_axis('text') +
  theme(aspect.ratio = 1)

# Make interactive points
plot(UMAP)

# save umap in PNG/PDF
ggsave(UMAP,filename = here::here('plots',sprintf('umap-idr-human-k%s.png',K)), height=12,width=12,bg='white')
ggsave(UMAP, filename = here::here('plots',sprintf('umap-idr-human-k%s.pdf',K)), height=12,width=12,bg='white')

#subset(umap_data,atar_proteins) %>% arrange(IDR_id) %>%
#  dplyr::relocate(PROTEIN,acc,IDR_id,S,E,S_atar,E_atar,atar_proteins,atar_overlap,PS_overlap) %>% View()

