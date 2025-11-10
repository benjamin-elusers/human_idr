# Helper Functions: Phase Separation Data

# Get data from source database for LLPS related regions in human
load_phasep_human = function(){
  .info$log("loading data from phasepDB on human phase-separating regions")
  #url_phasepdb_data = "http://db.phasep.pro/static/db/database/phaseodbv2_1_llps.xlsx"
  # 09.11.2025 DATA ARE STORED IN A ZIP FOLDER. URL CHANGED TO: "https://db.phasep.pro/data/phasepdbv2_1_data.zip" 
  # THE DATASET WAS DOWNLOADED AND EXTRACTED TO THE DATA FOLDER
  url_phasepdb_data = "data/phasepdbv2_1_llps.xlsx"
  
  phasepdb_llps = rio::import(url_phasepdb_data, na = c("","_")) |>
    mutate(region = str_replace_all(region,"–","-") ) |>
    # Keep only human
    filter(str_detect(organism,'Homo sapiens') & !is.na(region)) |>
    # split regions to rows
    separate_rows(region,sep='[,\\+]') |>
    # some regions are not given as positions
    mutate( ambiguous_region = str_detect(region,negate = T, pattern = "^\\d+-\\d+$")) |>
    filter( !ambiguous_region ) |>
    # Split region boundaries between start and end positions
    separate(region,remove = F, into=c('PS_START','PS_END'),sep='-') |>
    rowwise() |>
    # Some regions may be indicated by a domain name or a loose indication (e.g. C-terminal)
    mutate( PS_START = parse_integer(PS_START), PS_END =parse_integer(PS_END)) |>
    type_convert(col_types=cols(PMID = 'c')) |>
    mutate(PS_db='phasepdb',
           PS_id = ifelse( is_number(PS_START) & is_number(PS_END),
                           paste0(uniprot_entry,"_",PS_START,"..",PS_END),
                           paste0(uniprot_entry,"_",str_replace_all(region," ","-")))
    ) |>
    # remove phase separating regions with undefined or loosely defined boundaries
    filter( !is.na(PS_START) & !is.na(PS_END) | !is.na(region))  |>
    arrange(PS_id)
  return(phasepdb_llps)
}

load_phasepro_human = function(all_col=F){
  url_phaseprodb_data = "https://phasepro.elte.hu/download_full.json"
  col_id = c("taxon","region_ref","id","phase_id")
  col_ref = c("pmids","annotator","date","under_annot")
  col_ctrl = c("rna_req","ptm_affect","discrete_oligo","splice","membrane_clust")
  col_desc = c("description","segment","interaction","functional_class","disease")
  col_exp = c("determinants","partners","experiment_llps","experiment_state")
  
  phasepro = RJSONIO::fromJSON(url_phaseprodb_data) |>
    bind_rows() |>
    separate_rows(boundaries, sep="; ",convert = T) |>
    mutate( boundaries = str_extract_all(boundaries,"[0-9]+\\-[0-9]+")) |>
    separate(boundaries, into=c('PS_START','PS_END')) |>
    type_convert() |>
    mutate(ncbi_taxid=9606,evidence="curated",feature='phase_separation',source="phasepro",
           protein_len = nchar(sequence),
           feature_seq = str_sub(sequence,PS_START,PS_END),
           feature_len = nchar(feature_seq)) |>
    relocate(organism,accession,gene,name,common_name,
             PS_START,PS_END,segment,forms,organelles,
             domain_dep,ptm_dep,rna_dep,partner_dep,in_vivo,in_vitro) |>
    dplyr::rename(protein_seq = sequence) |>
    dplyr::select(-any_of(c(col_id,col_ref,col_ctrl,col_desc,col_exp,'boundaries'))) |>
    filter(organism == "Homo sapiens") |> 
    group_by(accession) |>
    mutate(content_count=sum(feature_len),
           content_fraction=content_count/protein_len,
           PS_id = paste0(accession,"_",PS_START,"..",PS_END),
           region=sprintf("%s-%s",PS_START,PS_END))
  return(phasepro)
}


# Get uniprot data about proteins with phase-separation regions
merge_ps_region = function(df_ps_region,maxgap=1){
  .info$log("Merge overlapping and nearly-contiguous phase-separating regions...")
  # make a dataset with only phase separation regions
  PS_all = df_ps_region |>
    arrange(acc) |>
    mutate(PS_len = PS_END-PS_START+1, PS_full = LEN == PS_len) |>
    distinct() |>
    add_count(acc,name='PS_n') |>
    ungroup() |>
    mutate(row=row_number())
  library(GenomicRanges)
  g0 = PS_all |> filter(PS_full | PS_n==1 )
  g1 = PS_all |> filter(!(PS_full | PS_n==1 ))
  
  g1.1 = g1 |>
    dplyr::select(db=PS_db,seqnames=acc,start=PS_START,end=PS_END) |>
    as("GRanges") %>%
    reduce(min.gapwidth=maxgap) %>% # maximum gap for merging two intervals
    as_tibble() |>
    dplyr::select(-strand) |>
    dplyr::rename(acc=seqnames, PS_START=start,PS_END=end,PS_len=width) |>
    mutate(PS_db ='merged',PS_id = paste0(acc,"_",PS_START,"..",PS_END)) |>
    mutate(region = paste0(PS_START,"-",PS_END))
  
  info_prot = df_ps_region |> dplyr::select(acc,GN,PN,LEN,KW,PE,Reviewed) |> distinct()
  
  # Merged phase-separation regions
  ps_merged = bind_rows(g1.1,g0) |>
    dplyr::select(-c(LEN,PS_full,PS_n,row) ) |>
    left_join(info_prot,by=c('acc','GN','PN','KW','PE','Reviewed')) |>
    mutate(PS_full = LEN == PS_len) |>
    ungroup() %>% distinct() |>
    arrange(acc,PS_id) |>
    group_by(acc) |>
    fill(c('GN','PN','KW','PE','Reviewed','LEN'),.direction = "updown") |>
    add_count(acc,name='PS_n') |>
    mutate(PS_full = LEN == PS_len, PS_n = PS_n - ((PS_n>1)*PS_full)) |>
    dplyr::select(-PS_db) |> distinct()
  
  n_ps = n_distinct(ps_merged$PS_id)
  n_psprot = n_distinct(ps_merged$acc)
  .info$log(sprintf("Phase-separating MERGED regions: %s (in %s proteins)",n_ps,n_psprot))
  
  return(ps_merged)
}

which_overlap = function(s1,e1,s2,e2){
  case_when(
    s2 < s1 & between_(e2,s1,e1) ~ "left",
    s2 < e1 & (s1 < s2 & e2 > e1) ~ "right",
    s1 < s2 & e2 < e1 ~ "around",
    s1 > s2 & e1 < e2 ~ "inside",
    s1 == s2 & e1 == e2 ~ "superposed",
    is.na(s2) | is.na(e2) ~ NA,
    .default = 'none'
  )
}

has_overlap = function(s1,e1,s2,e2){
  which_overlap(s1,e1,s2,e2) != "none"
}

get_human_phase_separation = function(mobidb=get_human_mobidb()){
  cat("loading data from mobidb on human phase-separating regions (annotated in phaseproDB)\n")
  ### phasepro -------------------------------------------------------------------
  
  # From mobidb with feature phase_separation and source phasepro
  phasepro_db  = load_phasepro_human() |>
    mutate(is_uniref = (accession %in% hs_uniref)) |>
    dplyr::select(ncbi_taxid,acc=accession,evidence,feature,source,
                  PS_id,PS_START,PS_END,content_fraction, content_count, 
                  feature_len, is_uniref, 
                  uniprot_seq = protein_seq, uniprot_len = protein_len) 
  
  phasepro_mobidb = mobidb |> 
    dplyr::select(-length) |>
    filter(feature == 'phase_separation' & source=='phasepro') |>
    dplyr::rename(PS_id = feature_id, PS_START=START, PS_END=END)
  
  # Turn phasepro as molecular features (columns)
  phasepro_wide = full_join(phasepro_mobidb,phasepro_db) |>
    mutate( region =  paste0(PS_START,"-",PS_END))|>
    dplyr::select(acc,PS_db=source,PS_id,PS_START,PS_END,region) |> 
    distinct() 
  ### phasepdb -------------------------------------------------------------------
  
  # Turn phasepdb as molecular features (columns)
  phasep_wide = load_phasep_human() |>
    dplyr::select(acc=uniprot_entry,PS_db,PS_id,PS_START,PS_END,region)
  
  # Combine human phase-separation regions
  PS_raw = full_join(phasepro_wide,phasep_wide) |> distinct()
  
  # Get uniprot annotations for phase-separating proteins
  dataset_uniprot_phase = here::here('data','uniprot-human-phase_separating_proteins.rds')
  UNI_PS = preload(dataset_uniprot_phase,
                   get_uniprot_ids(PS_raw$acc) |> distinct(),
                   doing = 'loading uniprot annotations for PS proteins')
  
  hs_PS = full_join(PS_raw,UNI_PS, by=c('acc'='AC')) |>
    arrange(acc,PS_id)
  
  n_ps = n_distinct(hs_PS$PS_id)
  n_psprot = n_distinct(hs_PS$acc)
  .info$log(sprintf("Phase-separating regions: %s (in %s proteins)",n_ps,n_psprot))
  
  return(hs_PS)
}
