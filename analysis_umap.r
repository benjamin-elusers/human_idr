# --- 1. LOAD SETUP, HELPERS, AND QC ---
source(here::here("src","setup_main.r"),echo = F,chdir = T)
source(here::here("src","helpers_analysis.r"),echo = F,chdir = T)

# --- 1b. QC AND LOAD DATA ---
data_list <- load_umap_data()
HS_IDR <- data_list$HS_IDR
ATAR_IDR <- data_list$ATAR_IDR
# .RData objects (hs_aacount, etc.) are loaded into GlobalEnv by the function

# --- 2. FEATURES COMPARISON (All) ---
.info$log("Running analysis: Features Comparison...")
col_mobidb = c('IDR_len','IDR_count','IDR_frac','uniprot_len')
numeric_features = c( colnames(hs_aacount), colnames(hs_aafreq),
                      colnames(hs_charge), colnames(hs_peptides),
                      colnames(hs_scores), colnames(hs_foldchange),
                      col_mobidb)

df_all_features = bind_rows(HS_IDR,ATAR_IDR) %>% 
  dplyr::select(PROTEIN,AC,IDR_id,START,END, any_of(numeric_features),has_PS, from_atar) %>%
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

# --- 3. FEATURES SELECTION ---
.info$log("Running analysis: Feature Selection & Correlation...")

features_to_use = c(paste0("fr_",c(get.AAA(),"X")),
                    colnames(hs_foldchange), 
                    c("netcharge_residue","charge_asymmetry","peptide_PI"),
                    c("stickiness","mean_stickiness","mean_roseman","mean_aggrescan"),
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

# --- 4. IDR UMAP ---
.info$log("Running analysis: UMAP...")
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

# --- 4. SAVE OUTPUT FOR C-SAT ANALYSIS ---
.info$log("Saving UMAP analysis output for C-SAT script...")
save(df_features, df_info, df_scaled, features_to_use, 
     file=here::here('data', 'analysis_umap_output.RData'))

.succ$log("UMAP analysis complete. Plots saved to plots/ directory.")


