# LOAD DATASETS ################################################################
# This script builds the primary datasets from source and saves them as TSV files.
# It should be run before the analysis scripts.

# Use here::here to ensure file paths are relative to the project root
# This sources all packages, helpers, and base data
source(here::here("src","setup_main.r"),echo = F,chdir = T)

# Helper function to write TSV with version/date header
write_tsv_with_header <- function(data, file_path, db_version_info) {
  # Create header with database info and timestamp
  header <- paste(
    "# DATASET VERSION",
    paste0("# Source data: ", db_version_info),
    paste0("# Last generated: ", Sys.time()),
    sep = "\n"
  )
  
  # Write the header to a new file (overwrites existing)
  readr::write_lines(header, file_path)
  
  # Append the data to the file
  readr::write_tsv(data, file_path, append = TRUE, col_names = TRUE)
  
  .info$log(sprintf("Successfully saved dataset to %s", file_path))
}

# --- Create data/ directory if it doesn't exist ---
data_dir <- here::here('data')
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
  .info$log(sprintf("Created directory: %s", data_dir))
}

# BUILD HUMAN IDR DATASET ######################################################
# Filter consensus disordered regions (at least 50% agreement of predictors)
.info$log("Building Human IDR Dataset...")
hs_mobidb = get_human_mobidb() 
hs_mobidb_merged = hs_mobidb |> dplyr::rename(S=START,E=END) |> merge_mobidb(gap_min=3L)
hs_diso =  dplyr::filter(hs_mobidb, is_uniref & feature == 'disorder' & source %in% 'th_50') |>
  dplyr::rename(IDR_id=feature_id) |> 
  arrange(IDR_id)

cluster <- multidplyr::new_cluster(14)
multidplyr::cluster_library(cluster, "dplyr")
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
.info$log("Writing HUMAN_MOBIDB_FEATURES.tsv...")
write_tsv_with_header(HS_IDR,
                      file=here::here('data','HUMAN_MOBIDB_FEATURES.tsv'),
                      db_version_info = "mobiDB (retrieved via setup_main.r)")

# BUILD ATAR IDR DATASET #######################################################
# LAST UPDATE OF INPUT DATA *FEB 2023*
.info$log("Building ATAR IDR Dataset...")
ATAR_CANDIDATES = rio::import(here::here('data',"ATAR_candidates_table.xlsx"))
df_atar = left_join(ATAR_CANDIDATES, df_hs_seq, by=c("UNIPROT"='uniprot_id')) |>
  dplyr::select(PROTEIN,UNIPROT,uniprot_seq,atar_sequence=SEQUENCE,
                START="Start Position",END="End Position", LEN = SEQ_LENGTH) |>
  mutate( IDR_id  = paste0(UNIPROT,"_",START,"..",END), atar_len = nchar(atar_sequence) ) |>
  arrange(IDR_id)

cluster <- multidplyr::new_cluster(10)
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

.info$log("Writing ATAR_IDR_FEATURES.tsv...")
write_tsv_with_header(ATAR_IDR,
                      file = here::here('data','ATAR_IDR_FEATURES.tsv'),
                      db_version_info = "ATAR_candidates_table.xlsx (Feb 2023)")


.info$log("Saving RData image for analysis scripts...")
# This .RData file contains calculated objects needed for the analysis scripts
save.image(file = here::here('data', 'IDR-FEATURES-DATA.rdata'))

.info$log("Data generation complete.")

