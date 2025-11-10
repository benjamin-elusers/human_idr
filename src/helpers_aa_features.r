# Helper Functions: Amino Acid Features

# --- 1. AA Constants ---
AA1=get.AA1()
AA3=get.AAA()
AA_PROP = get.aa.poperties()
AA_SCALES = get_aa_scales()

# --- 2. AA Helper Functions ---
normalize_sequence = function(BS){
  #data("BLOSUM62")
  BS_norm = BS |>
    Biostrings::chartr(old="U",new="C") |>
    Biostrings::chartr(old="O",new="K") |>
    Biostrings::chartr(old="J",new="L")
  return(BS_norm)
}

count_aa = function(string,verbose=F){
  if(verbose){
    .info$log("Count amino acids in sequence...")
    tic("Count amino acids in sequence...")
  }
  aminoacids = c(Biostrings::AA_STANDARD,"X")
  
  m_aacount = Biostrings::AAStringSet(x=string) |>
    Biostrings::letterFrequency(letters = aminoacids) 
  toc()
  return(m_aacount)
}

calculate_aascore = function(aa,verbose=F,scores=get_aa_scales(AA=3)){
  if(verbose){
    .info$log("Calculate sum of amino acid scores from residues count...")
    tic("Calculate sum of amino acid scores from residues count...")
  }
  aas = intersect(colnames(aa),rownames(scores))
  sum_scores = map(scores[aas,],
                   ~magrittr::multiply_by(t(aa),.x) |>  colSums()) %>% bind_rows()
  
  toc()
  df_scores = bind_cols(sum_aa=rowSums(aa[,aas]),sum_scores) |>
    mutate( across(all_of(names(scores)), ~ .x / sum_aa, .names = 'mean_{col}') )
  return(df_scores)
}

# --- 3. Shorthands for building IDR dataset ---
# These are the main functions used by build_datasets.r

get_peptstats = function(.data=.,col_sequence="feature_seq",cl){
  .info$log("Add pepstats in disordered regions...")
  tictoc::tic("Add pepstats in disordered regions...")
  
  df_peptide = group_by(.data,IDR_id) |>
    multidplyr::partition(cl) |>
    summarise(
      peptide_len = Peptides::lengthpep(!!sym(col_sequence)),
      peptide_mw  = Peptides::mw(!!sym(col_sequence)),
      peptide_mw_avg = peptide_mw / peptide_len,
      peptide_netcharge = Peptides::charge(!!sym(col_sequence)),
      peptide_PI = Peptides::pI(!!sym(col_sequence))
    ) |>
    ungroup() |>
    dplyr::collect()
  tictoc::toc()
  return( df_peptide )
}

get_aa_count =function(.data,col_sequence="feature_seq"){
  .info$log("Count amino-acids (single + chemical properties) in sequence...")
  tictoc::tic("Count amino-acids (single + chemical properties) in sequence...")
  aa_properties = get.aa.poperties()
  
  df_aacount = count_aa(.data[[col_sequence]]) |> as_tibble()
  df_aaprop = map(aa_properties, function(x){ dplyr::select(df_aacount,any_of(x)) |> rowSums() }) |> bind_cols()
  df_aaacount = df_aacount |> dplyr::rename(invert(get.AA3()))
  tictoc::toc()
  df_count = bind_cols(df_aaacount, df_aaprop)
  return(df_count)
}

get_aa_freq = function(aa_count,col_aa=c(get.AAA(),"X"),to100=T){
  aas = aa_count |> dplyr::select(all_of(col_aa))
  df_freq = df_freq = aa_count / rowSums(aas)  
  if(to100){ df_freq = 100 * df_freq  }
  colnames(df_freq) = paste0("fr_",colnames(df_freq))
  return(df_freq)
}

get_aa_scores =function(aa_count,col_aa_count=c(get.AAA())){
  .info$log("Compute amino-acids scores based on amino-acids counts...")
  tictoc::tic("Compute amino-acids scores  based on amino-acids counts...")
  df_aascore = calculate_aascore(aa_count |> dplyr::select(all_of(col_aa_count))) 
  tictoc::toc()
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
  tictoc::tic("Count of charged residues...")
  col_aa = c(get.AAA(),"X")
  POSITIVE =  aa_count[, col_pos] 
  NEGATIVE =  aa_count[, col_neg] 
  
  df_charged = tibble( sum_aa = aa_count[, col_aa] |> rowSums(),
                       positive= POSITIVE |> rowSums(),
                       negative= NEGATIVE |> rowSums()) |>
    rowwise() |>
    mutate(
      fr_positive = positive/sum_aa,
      fr_negative = negative/sum_aa,
      netcharge = positive - negative,
      charged = positive + negative,
      fr_charged = fr_positive+fr_negative,
      netcharge_residue =  abs(fr_positive-fr_negative),
      charge_asymmetry = charge.asym(fr_positive,fr_negative)) |>
    dplyr::select(-sum_aa)
  tictoc::toc()
  return( df_charged )
}

get_aa_foldchange = function(aa_count,ref_aa_freq,col_aa=c(get.AAA(),"X")){
  .info$log("Count of fold-change of residues within disordered regions...")
  tictoc::tic("Count of fold-change of residues within disordered regions...")
  
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
  tictoc::toc()
  colnames(df_foldchange) = paste0("fc_",colnames(df_foldchange))
  
  return(df_foldchange)
}
