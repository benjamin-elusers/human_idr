# LOAD DATASETS ################################################################
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source(here::here("src","setup_idr_analysis.r"),echo = F,chdir = T)
library(Biostrings)
library(hablar)

# BUILD HUMAN IDR DATASET ######################################################
# Filter consensus disordered regions (at least 50% agreement of predictors)
hs_mobidb = get_human_mobidb()
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
  relocate(ncbi_taxid,AC,ID,GN,ensp,NAME,
           uniprot_len,IDR_frac,IDR_count,
           IDR_id,IDR_len,START,END,source,feature,IDR_seq,
           colnames(hs_charge),
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
          mutate( IDR_id  = paste0(UNIPROT,"_",START,"..",END), atar_len = nchar(atar_sequence) ) |>
          arrange(IDR_id)

cluster <- new_cluster(10)
atar_idr = add_idr_sequence(df_atar,cluster) |> arrange(IDR_id)
atar_peptides = get_peptstats(atar_idr,cl=cluster) |> arrange(IDR_id)
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
  dplyr::rename(AC=UNIPROT, IDR_frac = content_fraction, IDR_count = content_count,
                IDR_len = atar_len) |>
  relocate(AC,PROTEIN, LEN,
           IDR_id,START,END, IDR_len, atar_sequence, 
           colnames(hs_charge),
           starts_with('pep_'),starts_with("fr_"),starts_with("fc_"), aggrescan:wimleywhite,
           starts_with("mean_"))

#write_tsv(ATAR_IDR, file = here::here('data','ATAR-IDR-FEATURES.tsv'))
#save.image(file = here::here('data', 'IDR-features-data.rdata'))


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
summary(df_features)
# Check correlogram of selected features

p_used = make_features_correlation(df_features)
ggsave(p_used,path=here::here("plots"),
       filename ='correlation-human-idr-selected-features.png',
       height=12,width=12,bg='white')
ggsave(p_used,path=here::here("plots"),
       filename ='correlation-human-idr-selected-features.pdf',
       height=12,width=12,bg='white')

# IDR UMAP #####################################################################

k_neighbors = c(10,20,30,40,50,100)
UMAP_SCALED = lapply(k_neighbors, function(x){ make_umap(df_scaled, K = x, seed = 142, scale = T)} ) |>
              patchwork::wrap_plots(nrow = 2)

# save umap in PNG/PDF
plot_name = 'umap-idr-human-scaled'
ggsave(UMAP_SCALED, path=here::here("plots"),
       filename = paste0(plot_name,'.png'),scale=1.3,
       device = 'png', height=12, width=12, bg='white')
ggsave(UMAP_SCALED, path=here::here("plots"),
       filename = paste0(plot_name,'.pdf'), scale=1.3,
       device = 'pdf', height=12, width=12, bg='white')

UMAP_UNSCALED = lapply(k_neighbors, function(x){ make_umap(df_features, K = x, seed = 142, scale = F)} ) |>
  patchwork::wrap_plots(nrow = 2)

plot_name = 'umap-idr-human-unscaled'
ggsave(UMAP_UNSCALED, path=here::here("plots"),
       filename = paste0(plot_name,'.png'),scale=1.3,
       device = 'png', height=12, width=12, bg='white')
ggsave(UMAP_UNSCALED, path=here::here("plots"),
       filename = paste0(plot_name,'.pdf'), scale=1.3,
       device = 'pdf', height=12, width=12, bg='white')

