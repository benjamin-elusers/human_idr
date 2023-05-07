# Source Yeastomics functions (recommended) ----------------------------------------------
url_yeastomics = "https://raw.githubusercontent.com/benjamin-elusers/yeastomics/main/src/"
#source(paste0(url_yeastomics,"__setup_yeastomics__.r"))
# NOTE: Missing packages would be installed on first run.
# Those two scripts should contain most functions used
source(paste0(url_yeastomics,"utils.r"))
source(paste0(url_yeastomics,"function_datapub.r"))

# defaul options for R environment
options(dplyr.summarise.inform=FALSE,
        dplyr.width=Inf,
        max.print=2e4,
        timeout = max(600, getOption("timeout")))

# Loading custom functions -----------------------------------------------------
normalize_sequence = function(BS){
  #data("BLOSUM62")
  BS_norm = BS |>
            Biostrings::chartr(old="U",new="C") |>
            Biostrings::chartr(old="O",new="K") |>
            Biostrings::chartr(old="J",new="L")
  return(BS_norm)
}

get_aa_scales=function(AA=1){
  scores = data.frame(row.names = get.AA1(),
    aggrescan=get.aggrescan(),
    camsol=get.camsol(),
    foldamyloid=get.foldamyloid(),
    kytedoolittle=get.kytedoolittle(),
    pawar=get.pawar_ph7(),
    roseman=get.roseman(),
    stickiness=get.stickiness(),
    voronoi_sickiness=get.voronoi_stickiness(),
    wimleywhite=get.wimleywhite()
  )
  if(AA==3){ rownames(scores)=get.AAA() }
  return(scores)
}


count_aa = function(string,verbose=F){
  if(verbose){
    .info$log("Count amino acids in sequence...")
    tic("Count amino acids in sequence...")
  }
 aminoacids = c(Biostrings::AA_STANDARD,"X")

  m_aacount = Biostrings::AAStringSet(x=string) |>
              letterFrequency(letters = aminoacids) 
  toc()
  return(m_aacount)
}

calculate_aascore = function(aa,verbose=F,scores=get_aa_scales(AA=3)){
  if(verbose){
    .info$log("Calculate sum of amino acid scores from residues count...")
    tic("Calculate sum of amino acid scores from residues count...")
  }
  aas = intersect(rownames(scores),colnames(aa))
  sum_scores = map(scores[aas,],
              ~magrittr::multiply_by(aa,.x) |>  rowSums()) |> bind_rows()
  
  toc()
  df_scores = bind_cols(sum_aa=rowSums(aa[,aas]),sum_scores) |>
              mutate( across(all_of(names(scores)), ~ .x / sum_aa, .names = 'mean_{col}') )
  return(df_scores)
}

# Get data from source database for LLPS related regions in human
load_phasep_human = function(){
  .info$log("loading data from phasepDB on human phase-separating regions")
  url_phasepdb_data = "http://db.phasep.pro/static/db/database/phaseodbv2_1_llps.xlsx"
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

# Loading needed packages ------------------------------------------------------
library(tidyverse)
library(multidplyr) # Using this library to speed up using parallel computing
library(log)
.info  = infoLog()
.error =  errorLog()
.warn  = warningLog()
.succ  = successLog()
.dbg = Logger$new("DEBUG")$
  date()$
  time()$
  hook(crayon::bgWhite)

pkg <- c('here','tidyverse','dplyr','furrr','progressr','pbmcapply','Biostrings')
xfun::pkg_load2(pkg)
NCPUS <- parallel::detectCores() - 2

# Loading precomputed data -----------------------------------------------------
## amino acids -----------------------------------------------------------------
AA1=get.AA1()
AA3=get.AAA()
AA_PROP = get.aa.poperties()
AA_SCALES = get_aa_scales()

## human datasets --------------------------------------------------------------
# human = find.uniprot_refprot(c('9606','HUMAN','homo sapiens'))  %>%
#   arrange(desc(keyword_matched))

### uniprot --------------------------------------------------------------------
# Get human uniprot reference proteome dataset and sequences
# human = find.uniprot_refprot(c('HUMAN','homo sapiens'))
dataset_uniprot_human =  here::here('data','uniprot-human-proteome_data.rds')
hs_uni = preload(dataset_uniprot_human,
                 get_uniprot_reference(9606),
                 "load human uniprot reference proteome data")

sequence_uniprot_human = here::here('data','uniprot-human-proteome_sequences.rds')
hs_aa = preload(sequence_uniprot_human,
                  get.uniprot.proteome(taxid = 9606, DNA = F),
                  "load human uniprot reference proteome sequences")
hs_uniref = names(hs_aa)

df_hs_seq = normalize_sequence(as.vector(hs_aa)) |>
            enframe(name = 'uniprot_id',value='uniprot_seq') |>
            mutate(uniprot_len = nchar(uniprot_seq))


### mobidb --------------------------------------------------------------------
# Get human data on intrinsically disordered regions from mobidb (matching reference uniprot)
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

# Shorthands for building IDR dataset -------------------------------------------
# They should be used in series to generate a complete dataset

add_idr_sequence = function(.data,cl){
  # This function expects a mobidb formatted dataset
  .info$log("Add sequences of mobidb intrinsically disordered regions...")
  tic("Add sequences of mobidb intrinsically disordered regions...")
  df_diso = group_by(.data, IDR_id) |>
    partition(cl) |>
    mutate( feature_seq = stringr::str_sub(uniprot_seq, start=START,end=END) ) |>
    ungroup() |>
    collect()
  toc()
  return(df_diso)
}

get_peptstats = function(.data=.,cl){
  .info$log("Add pepstats in disordered regions...")
  tic("Add pepstats in disordered regions...")

  df_peptide = group_by(.data,IDR_id) |>
    partition(cl) |>
    summarise(
      peptide_len = Peptides::lengthpep(feature_seq),
      peptide_mw  = Peptides::mw(feature_seq),
      peptide_mw_avg = peptide_mw / peptide_len,
      peptide_netcharge = Peptides::charge(feature_seq),
      peptide_PI = Peptides::pI(feature_seq)
    ) |>
    ungroup() |>
    collect()
  toc()
  return( df_peptide )
}

get_aa_count =function(.data,col_sequence="feature_seq"){
  .info$log("Count amino-acids (single + chemical properties) in sequence...")
  tic("Count amino-acids (single + chemical properties) in sequence...")
  aa_properties = get.aa.poperties()

  df_aacount = count_aa(.data[[col_sequence]]) |> as_tibble()
  df_aaprop = map(aa_properties, function(x){ dplyr::select(df_aacount,any_of(x)) |> rowSums() }) |> bind_cols()
  df_aaacount = df_aacount |> dplyr::rename(invert(get.AA3()))
  toc()
  df_count = bind_cols(df_aaacount, df_aaprop)
  return(df_count)
}

get_aa_freq = function(aa_count,col_aa=c(get.AAA(),"X")){
  aas = aa_count |> dplyr::select(all_of(col_aa))
  df_freq = 100* aa_count / rowSums(aas)
  colnames(df_freq) = paste0("fr_",colnames(df_freq))
  return(df_freq)
}

get_aa_scores =function(aa_count,col_aa_count=c(get.AAA())){
  .info$log("Compute amino-acids scores based on amino-acids counts...")
  tic("Compute amino-acids scores  based on amino-acids counts...")
  df_aascore = calculate_aascore(aa_count |> dplyr::select(all_of(col_aa_count))) 
  toc()
  return( df_aascore )
}

get_aa_topfreq = function(aa_count,col_aa_count=c(get.AAA(),"X"),topn=4){
  .info$log("Top-ranked amino-acids by frequencies...")
  tic("Top-ranked amino-acids by frequencies...")

  # Add top4 most frequent amino acids in IDR (with their frequency in IDR)
  df_topfreq = aa_count |>
    dplyr::select(all_of(col_aa_count)) |>
    mutate(id=row_number()) |>
    pivot_longer(cols = -id, names_to = 'aa', values_to='count') |>
    group_by(id) |>
    mutate( rk = rank(-count, ties = 'first'), tot=sum(count), freq = round(100*count/tot,2),
            name.val = sprintf("%s:%2.1f:%s",aa,freq,count)) |>
    slice_max(order_by = count, n = topn, with_ties = F) |>
    group_by(id) %>%
    reframe(frtop_aa = paste0(aa, collapse=":"),
              frtop_fr=paste0(freq, collapse=":"),
              frtop_count = paste0(count, collapse = ":"),
              tot) |>
    distinct() |>
    ungroup() |>
    dplyr::select(-id)
  toc()
  return( df_topfreq )
}

get_aa_topfc= function(aa_fc,col_aa=paste0("fc_",c(get.AAA(),"X")),topn=4){
  .info$log("Top-ranked amino-acids by fold-changes....")
  tic("Top-ranked amino-acids by fold-changes...")

  # Add top4 most frequent amino acids in IDR (with their frequency in IDR)
  df_topfc = aa_fc |>
    dplyr::select(all_of(col_aa)) |>
    mutate(id=row_number()) |>
    pivot_longer(cols = -id, names_to = 'aa', values_to='fc') |>
    group_by(id) |>
    mutate( rk = rank(-fc, ties = 'first'), fc=round(fc,1), lfc = round(log2(fc),1) ,
            name.val = sprintf("%s:%2.1f:%2.1f",aa,fc,lfc)) |>
    slice_max(order_by = fc, n = topn, with_ties = F) |>
    group_by(id) %>%
    reframe(fctop_aa = paste0(str_replace_all(aa,"fc_",""), collapse=":"),
            fctop_fc=paste0(fc, collapse=":"),
            fctop_lfc = paste0(lfc, collapse = ":")) |>
    distinct() |>
    ungroup() |>
    dplyr::select(-id)
  toc()
  return( df_topfc )
}

get_aa_charge = function(aa_count,col_pos=c("LYS","ARG","HIS"), col_neg= c("ASP","GLU") ){
  .info$log("Count of charged residues...")
  tic("Count of charged residues...")
  df_charged = tibble( positive=  aa_count[,col_pos] |> rowSums(),
            negative= aa_count[,col_neg] |> rowSums(),
            netcharge=(positive-negative) )
  toc()
  return( df_charged )
}

get_aa_foldchange = function(aa_count,ref_aa_freq,col_aa=c(get.AAA(),"X")){
  .info$log("Count of fold-change of residues within disordered regions...")
  tic("Count of fold-change of residues within disordered regions...")

  aas = aa_count |> dplyr::select(all_of(col_aa))
  naa=  rowSums(aas)
  # calculate the frequencies of residues in disordered regions
  aas_freq = tibble(aas,tot=naa) |> rowwise() |> mutate(across(all_of(col_aa), ~.x/tot))
  if(missing(ref_aa_freq)){
    # calculate the frequencies of residues across disordered proteome
    tot_aa = sum(naa)
    ref_aa_freq = aas |> colSums() |> magrittr::divide_by(tot_aa)
  }

  df_foldchange = sweep( aas_freq[,col_aa],2,ref_aa_freq,"/")
  toc()
  colnames(df_foldchange) = paste0("fc_",colnames(df_foldchange))

  return(df_foldchange)
}

make_features_correlation= function(df_data, path_plot = here::here('plots')){
  
  df_num = df_data %>%
    dplyr::select(where(~ is.numeric(.x))) %>%
    dplyr::select(-START,-END)
  
  cor_features = ggcorrplot::ggcorrplot(cor(df_num),
                                        outline.color = 'transparent', 
                                        tl.cex = 8, tl.srt = 70)
  ggsave(cor_features,path=path_plot,
         filename ='correlation-human-idr-features.png',
         height=12,width=12,bg='white')
  
}

make_umap = function(df_data, K=10, seed=142, scale=F, path_plot = here::here('plots')){
  
  library(umap)
  library(ggalt)
  library(ggiraph)
  library(ggforce)
  set.seed(seed)
  umap.config = umap.defaults
  umap.config$n_neighbors = K
  
  tag_neighbor = sprintf("-K%s",K)
  
  df_num = df_data %>%
    dplyr::select(where(~ is.numeric(.x))) %>%
    dplyr::select(-START,-END)
  
  df_info = df_data %>% dplyr::select( -colnames(df_num) )
  # Compute umap based on scaled features
  cat(sprintf("Compute umap with K=%s neighbors using ",K))
  
  if(scale){
    cat(": scaled features...\n")
    tag_scale = "-scaled"
    features_map <- df_num %>% scale() %>% umap::umap(seed = seed, config=umap.config) 
  }else{
    cat(": raw features (unscaled)...\n")
    tag_scale = "-raw"
    features_map <- df_num %>% umap::umap(seed = seed, config=umap.config)
  }
  
  # Mark outliers on umap coordinates (in 1/99% percentiles or outside the interval defined below for X1/X2)
  df_umap = features_map$layout %>% magrittr::set_colnames(c('X1','X2')) %>%
    bind_cols(features_map$data %>% as_tibble() %>% rename_with(.fn = xxS, sx=tag_scale,s='')) %>%
    bind_cols(df_num) %>%
    bind_cols(df_info) %>%
    mutate(outliers_x1 = !between(percent_rank(X1),0.01,0.99),
           outliers_x2 = !between(percent_rank(X2),0.01,0.99)) %>%
    mutate(pt_lab=paste0(PROTEIN," ",START,"-",END)) |>
    filter(between(X1,-10,10), between(X2,-10,10)) 
  
  x1.q = quantile(df_umap$X1,seq(0,1,len=101))
  x2.q = quantile(df_umap$X2,seq(0,1,len=101))
  print(summary(df_umap[,c('X1','X2')]))
  
  n_idr= n_distinct(df_info$IDR_id)
  n_idr_umap= n_distinct(df_umap$IDR_id)
  
  n_prot= n_distinct(df_info$AC)
  n_prot_umap= n_distinct(df_umap$AC)
  sample_size_text = sprintf("IDR=%s (%s) | PROT=%s (%s) ",n_idr_umap,n_idr,n_prot_umap,n_prot)
  
  umap_atar = subset(df_umap,from_atar)
  # Plot the umap
  UMAP = ggplot(data=df_umap ,aes(x = X1,y = X2)) +
    # all idrs
    geom_point(size=1.5,shape=16,color='gray',alpha=1) +
    geom_text(data=NULL,label=sample_size_text, x=Inf, y=Inf, hjust='right',vjust='top',check_overlap = T) +
    # highlight atar idr
    geom_point(data=umap_atar, aes(color=PROTEIN), shape=16, size=4,alpha=0.9) +
    # annotate atar idr 
    ggrepel::geom_text_repel(data=umap_atar,aes(label=pt_lab,color=PROTEIN), size=3, fontface='italic',
                             max.overlaps=15, force=5, force_pull=0, seed=291222) +
    # graphical parameters
    labs(x = "UMAP1", y = "UMAP2", subtitle="")+
    see::scale_color_metro(palette = 'rainbow', discrete = T) +
    theme(legend.position="bottom") +
    ggeasy::easy_text_size(c("axis.title","axis.text.x", "axis.text.y"), size = 20)+
    ggeasy::easy_remove_legend() +
    ggeasy::easy_remove_x_axis('ticks') + ggeasy::easy_remove_y_axis('ticks') +
    theme(aspect.ratio = 1)
  
  # Make interactive points
  plot(UMAP)
  
  # save umap in PNG/PDF
  plot_name = sprintf('umap-idr-human%s%s',tag_neighbor,tag_scale)
  ggsave(UMAP, path = path_plot,
         filename = paste0(plot_name,'.png'),
         device = 'png', height=12, width=12, bg='white')
  ggsave(UMAP, path = path_plot, 
         filename = paste0(plot_name,'.pdf'), 
         device = 'pdf', height=12, width=12, bg='white')
  return(UMAP)
}

