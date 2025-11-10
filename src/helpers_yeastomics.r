# Helper Functions: Yeastomics

# Strings ----------------------------------------------------------------------
str2chr  <- function(x){ return( unlist(strsplit(x,split='')) ) } # split string by character
concat   <- function(...){ return( paste(...,collapse='') ) }     # concatenate characters to string

towords  <- function(x){ return( unlist(strsplit(x,split="\\s")) ) } # split string by any white space
xxS      <- function(x,sx,s='.'){ paste0(x,s,sx) } # Add suffix to a string
Pxx      <- function(x,px,s='.'){ paste0(px,s,x) } # Add prefix to a string
starting <- function(str,pre){ str[ startsWith(str,pre) ] }
ending   <- function(str,suf){ str[ endsWith(str,suf) ] }

trim.lead  <- function (x){ sub("^\\s+", "", x) } # returns string w/o leading whitespace
trim.trail <- function (x){ sub("\\s+$", "", x) }# returns string w/o trailing whitespace
trim       <- function (x){ gsub("^\\s+|\\s+$", "", x) } # returns string w/o leading or trailing whitespace
trim.NaN   <- function(d){ d[apply(d,2,function(x) any(is.nan(x))),] } # returns a vector without NaN

paste.even <- function(x,s='-'){ paste(evens(x),odds(x),sep=s) }
paste.odd  <- function(x,s='-'){ c(paste.even(but.last(x),s), last(x)) }

rm.dup.str = function(str,sep='.'){ # removes duplicated word in string
  words = strsplit(str,paste0("\\",sep))
  words.unique = lapply(words,unique)
  new.str = unlist(lapply(words.unique,paste,collapse=sep))
  return(new.str)
}

repchar <- function(char,times){ # replicates a character multiple times as a string
  return(paste( rep(char, times), collapse='' ))
}

subname=function(name,sep="\\.",lc=F){ # extracts substring until first separator
  b4sep = sprintf("^([^%s]+)%s",sep,sep)
  part1 = sub(b4sep, "\\1", x=name)
  if(lc){ tolower(part1) }
  return(part1)
}

get.longest = function(S, s='\\.'){
  # get the longest string in a list of splitted string
  library(stringr)
  L = str_split(string = S, pattern = s)
  long=sapply(L,function(x){ nc=nchar(x); which.max(nc)})
  sapply(1:length(L),function(i){ L[[i]][long[i]] })
}

strfind = function(strings, patterns, index=F){ # search multiple patterns in character vectors
  sapply(patterns,  function(p){ grep(x = strings, pattern = p, value = !index) })
}

# Testing/Subsetting -----------------------------------------------------------
is.whole  <- function(x){ all(floor(x) == x) }      # checks if a value has decimal part (not necessarily integer e.g. 1.0 is whole)
is.string <- function(x){  length(x) == 1 & is.character(x) } # checks if a value is a single string of letters
#is.binary <- function(x){ all( (1*x) %in% 0:1) }    # checks if a value is 0/1 (binary/logical)

is_num_bin <- function(x){ is.numeric(x) & length(unique(na.omit(x))) == 2 } # checks if a numeric vector has only 2 values
is_fac_bin <- function(x){ is.factor(x) & nlevels(na.omit(x)) == 2 }         # checks if a factor vector has only 2 levels
# checks if a value is binary (2 unique outcomes -- by default numeric values must be either 0 and 1)
is_binary  <- function(x,xvals=c(0,1)){ is_num_bin(x) & all(range(na.omit(x)) %in% xvals[1:2]) | is_fac_bin(x) | is.logical(x) }
is_number <- function(x){ grepl("^[0-9]+$",x) } # checks if it contains only number
is_frequency <- function(x){ all(is.numeric(x) & x>=0 & x<=100) } # checks if values are between 0 and 100

is.even   <- function(x){ as.integer(x) %% 2 == 0 } # checks if a value is even
is.odd    <- function(x){ as.integer(x) %% 2 != 0 } # checks if a value is odd
is.dup    <- function(x){ x %in% x[duplicated(x)] } # detects duplicates

is.in = function(el,set,case.sensitive=T,withNames=T){
  found = el %in% set
  if(!case.sensitive){ found = toupper(el) %in% toupper(set) }
  if(withNames) return(setNames(found,el))
  return(found)
}
invert    <- function(x){ setNames(names(x),make.unique(x)) }

load_datafile = function(datafile) {
  ext <- tools::file_ext(tolower(datafile))
  if (ext == "rds") {
    return(readRDS(datafile))
  } else if (ext %in% c("rda", "rdata", "data")) {
    load(datafile, envir = .GlobalEnv) # Loads into the global environment
    return(invisible(NULL)) # Return NULL since the data is in the environment
  } else {
    stop("Unsupported file format: ", ext)
  }
}

preload = function(saved.file, loading.call, doing = "Creating data...") {
  library(tictoc)
  cat(doing, "\n")
  
  # If the file doesn't exist, create and save it
  if (!file.exists(saved.file)) {
    tic(doing)
    res <- eval(substitute(loading.call))
    
    # Ensure the directory exists
    save_dir <- dirname(saved.file)
    if (!dir.exists(save_dir)) {
      message("Directory does not exist. Creating path: ", save_dir)
      dir.create(save_dir, recursive = TRUE)
    }
    
    # Save the result in RDS format
    saveRDS(res, saved.file)
    toc()
    return(res)
  } else {
    # Load the existing file
    return(load_datafile(saved.file))
  }
}
# Amino acids ------------------------------------------------------------------
# SEQ must be a character vector
is.pos  <- function(SEQ){ SEQ == "K" | SEQ == "R" }
is.neg  <- function(SEQ){ SEQ == "D" | SEQ == "E" }
pos     <- function(SEQ){ sum(is.pos(SEQ)) }
neg     <- function(SEQ){ sum(is.neg(SEQ)) }
fpos    <- function(SEQ){ mean(is.pos(SEQ)) }
fneg    <- function(SEQ){ mean(is.neg(SEQ)) }
charged <- function(SEQ){ return( pos(SEQ) + neg(SEQ) ) }
fcr     <- function(SEQ){ return( fpos(SEQ) + fneg(SEQ)  ) }
npcr    <- function(SEQ){ return( abs( fpos(SEQ) - fneg(SEQ) )  ) }
charge.asym <- function(p,n){ #p=positives n=negatives (sums)
  if( (p+n) == 0){ return(0) }
  return( (p-n)**2 / (p+n) )
}
get.AA1 = function(){ unlist(strsplit("ACDEFGHIKLMNPQRSTVWY","")) }
get.AAA = function(){ c("ALA","CYS","ASP","GLU","PHE",
                        "GLY","HIS","ILE","LYS","LEU",
                        "MET","ASN","PRO","GLN","ARG",
                        "SER","THR","VAL","TRP","TYR") }
get.AA3 = function(){ c('A'="ALA",'C'="CYS",'D'="ASP",'E'="GLU",'F'="PHE",
                        'G'="GLY",'H'="HIS",'I'="ILE",'K'="LYS",'L'="LEU",
                        'M'="MET",'N'="ASN",'P'="PRO",'Q'="GLN",'R'="ARG",
                        'S'="SER",'T'="THR",'V'="VAL",'W'="TRP",'Y'="TYR") }

get.AA.df = function(){
  aaa = get.AA3()
  return( data.frame( a=names(aaa), aaa ) )
}
get.aggrescan = function(){
  setNames(
    object=c(-0.036, 0.604, -1.836, -1.412, 1.754, -0.535, -1.033, 1.822, -0.931, 1.38,
             0.91, -1.302, -0.334, -1.231, -1.24, -0.294, -0.159, 1.594, 1.037, 1.159),
    nm=get.AA1()
  )
}

get.camsol = function(){
  setNames(
    object=c(-0.4533509211, -3.039164569,   3.806164055,  4.112172118, -3.922859742,
             0.809691392,  -0.2864452149, -2.986464178,  3.64837007,  -2.150991033,
             -1.649718164,   1.353162125,   1.405662227,  0.2038292233, 4.093823461,
             0.5372203047, -0.9718594221, -3.593600781, -3.759140236, -2.931244491),
    nm = get.AA1()
  )
}

get.foldamyloid = function(){
  setNames(
    object = c(0.086, 0.568, -0.776, -0.632, 0.958, -1.088, 0.025, 1.217, -0.565, 1.015, 0.725,
               -0.713, -2.303, -0.271, 0.032, -0.73, -0.349, 0.92, 1.027, 0.851),
    nm = get.AA1()
  )
}

get.kytedoolittle = function(){
  setNames(
    object=c(1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, -1.6, -3.5,
             -4.5, -0.8, -0.7, 4.2, -0.9, -1.3),
    nm=get.AA1()
  )
}

get.pawar_ph7 = function(){
  setNames(
    object=c(-3.31, 1.61, -9.42, -10.38, 2.8, -3.96, -4.31, 0.93, -9.55, -0.25, -1.06, -6.02,
             -11.96, -6, -11.93, -5.08, -2.12, 0.49, 2.92, 1.03),
    nm=get.AA1()
  )
}

get.roseman = function(){
  setNames(
    object=c(0.39, 0.25, -3.81, -2.91, 2.27, 0, -0.64, 1.82, -2.77, 1.82,  0.96, -1.91, 0.99,
             -1.3, -3.95, -1.24, -1, 1.3, 2.13, 1.47),
    nm=get.AA1()
  )
}

get.stickiness = function(){
  setNames(
    object=c(0.0062, 1.0372, -0.7485, -0.7893, 1.2727, -0.1771, 0.1204, 1.1109, -1.1806,
             0.9138, 1.0124, -0.2693, -0.1799, -0.4114, -0.0876, 0.1376, 0.1031, 0.7599,
             0.7925, 0.8806),
    nm=get.AA1()
  )
}

get.voronoi_stickiness = function(){
  setNames(
    object=c(-0.040638565,  0.739811788, -0.281053408, -0.392609711, 0.855804199,
             -0.015564908,  0.189364177,  0.593731643, -0.546173776, 0.560032835,
             0.569250951, -0.163722093, -0.089014045, -0.142636678, 0.027334912,
             -0.101035135, -0.024028749,  0.372857295,  0.913761361, 0.79208657),
    nm=get.AA1()
  )
}

get.wimleywhite = function(){
  setNames(
    object=c(4.08, 4.49, 3.02, 2.23, 5.38, 4.24, 4.08, 4.52, 3.77, 4.81, 4.48, 3.83, 3.8,
             3.67, 3.91, 4.12, 4.11, 4.18, 6.1, 5.19),
    nm=get.AA1()
  )
}


get_aa_scales=function(AA=1){
  scores = data.frame(row.names = get.AA1(),
                      aggrescan= get.aggrescan(),
                      camsol=get.camsol(),
                      foldamyloid=get.foldamyloid(),
                      kytedoolittle=get.kytedoolittle(),
                      pawar=get.pawar_ph7(),
                      roseman=get.roseman(),
                      stickiness=get.stickiness(),
                      voronoi_stickiness=get.voronoi_stickiness(),
                      wimleywhite=get.wimleywhite()
  )
  if(AA==3){ rownames(scores)=get.AAA() }
  return(scores)
}

get.aa.poperties = function(){
  aa_prop = list(
    tiny      = c('A','C','G','S','T'),
    small     = c('A','B','C','D','G','N','P','S','T','V'),
    aliphatic = c('I','L','V'),
    aromatic  = c('F','H','W','Y'),
    nonpolar  = c('A','C','F','G','I','L','M','P','V','W','Y'),
    polar     = c('D','E','H','K','N','Q','R','S','T','Z'),
    charged   = c('B','D','E','H','K','R','Z'),
    basic     = c('H','K','R'),
    acidic    = c('B','D','E','Z'),
    alcohol   = c('S','T'),
    turnlike  = c('A','C','D','E','G','H','K','N','Q','R','S','T')
  )
  
  aa_grouped = sapply(aa_prop,paste0,collapse="")
  names(aa_prop) = paste0(names(aa_prop),"_",aa_grouped)
  return(aa_prop)
}