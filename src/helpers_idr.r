# Helper Functions: IDR Data Retrieval

# Get human data on intrinsically disordered regions from mobidb (matching reference uniprot)
# This function depends on `df_hs_seq` and `hs_uniref` which are loaded in `setup_main.r`
get_human_mobidb = function(){
  dataset_mobidb_human = here::here('data','mobidb-human-features.rds')
  hs_mobidb = preload(dataset_mobidb_human,
                      load.mobidb(9606), # should be 9606
                      'retrieve human mobidb features...') |>
    dplyr::rename(START=S, END=E) |>
    mutate(feature_id = paste0(acc,"_",START,"..",END),
           is_uniref = (acc %in% hs_uniref)) |>
    left_join(df_hs_seq,by=c('acc'='uniprot_id')) |>
    filter(is_uniref & length == uniprot_len & !duplicated(feature_id) ) |>
    distinct()
}

merge_mobidb = function(mobidb, gap_min=1L){
  library(GenomicRanges)
  if(missing(mobidb)){ stop("requires mobidb data...") }
  
  merged = mobidb %>%
    dplyr::select(seqnames=acc,start=S,end=E) %>%
    as("GRanges") %>%
    reduce(min.gapwidth=gap_min) %>% # maximum gap for merging two intervals
    as_tibble() %>%
    dplyr::rename(acc=seqnames,S=start,E=end,feature_len=width) %>%
    mutate(acc=as.character(acc)) %>%
    dplyr::select(-strand) %>%
    arrange(acc)
  return(merged)
}

add_idr_sequence = function(.data,cl){
  # This function expects a mobidb formatted dataset
  .info$log("Add sequences of mobidb intrinsically disordered regions...")
  tictoc::tic("Add sequences of mobidb intrinsically disordered regions...")
  df_diso = group_by(.data, IDR_id) |>
    multidplyr::partition(cl) |>
    mutate( feature_seq = stringr::str_sub(uniprot_seq, start=START,end=END) ) |>
    ungroup() |>
    dplyr::collect()
  tictoc::toc()
  return(df_diso)
}
