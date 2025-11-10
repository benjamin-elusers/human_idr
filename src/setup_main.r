# Main Setup Script
# This script loads all packages, sets global options, and sources all helper scripts.
library(here)
library(tidyverse)
# --- 1. LOAD YEASTOMICS FUNCTIONS ---
#url_yeastomics = "https://raw.githubusercontent.com/benjamin-elusers/yeastomics/main/src/"
# NOTE: Missing packages would be installed on first run.
#source(paste0(url_yeastomics,"utils.r"))
#source(paste0(url_yeastomics,"function_datapub.r"))
source(here::here("src", "helpers_yeastomics.r"))

# --- 2. SET GLOBAL R OPTIONS ---
options(dplyr.summarise.inform=FALSE,
        dplyr.width=Inf,
        max.print=2e4,
        timeout = max(600, getOption("timeout")))

# --- 3. SETUP LOGGING ---
library(log)
.info  = infoLog()
.error =  errorLog()
.warn  = warningLog()
.succ  = successLog()
.dbg = Logger$new("DEBUG")$
  date()$
  time()$
  hook(crayon::bgWhite)

# --- 4. LOAD ALL REQUIRED PACKAGES ---
pkg <- c(
  'here', 'tidyverse', 'dplyr', 'furrr', 'progressr', 'pbmcapply', 'multidplyr', 'Biostrings', 
  'GenomicRanges', 'RJSONIO', 'Peptides', 'ggcorrplot', 'umap', 'ggalt', 'ggplot2',
  'ggiraph', 'ggforce', 'ggrepel', 'see', 'ggpubr', 'hrbrthemes', 
  'ggeasy', 'patchwork', 'rio', 'magrittr', 'hablar', 'multidplyr', 'readr'
)
xfun::pkg_load2(pkg)
NCPUS <- parallel::detectCores() - 2

# --- 5. SET GLOBAL GGPLOT THEME ---
AXIS_TITLE_SIZE=14
AXIS_TEXT_SIZE=10
TEXT_SIZE=4

th_txt_size = hrbrthemes::theme_ipsum(
  axis_title_just = 'm', axis = 'xy', axis_col = 'black',
  grid = F, base_family = 'Helvetica',
  axis_text_size = AXIS_TEXT_SIZE,
  axis_title_size = AXIS_TITLE_SIZE)

default_look = th_txt_size +
  ggpubr::grids(axis = 'xy') +
  ggplot2::theme(aspect.ratio = 1,
        plot.background = ggplot2::element_rect(fill=NA,linewidth = NA,color=NA),
        panel.grid.major = ggplot2::element_line(linewidth=0.5,linetype='solid'),
        panel.grid.minor = ggplot2::element_line(linewidth=0.5,linetype='dashed'),
        plot.margin = ggplot2::margin(0,3,0,3,unit = 'pt')
  ) +
  ggeasy::easy_remove_legend()

ggplot2::theme_set(default_look)

# --- 6. SOURCE HELPER FUNCTION SCRIPTS ---
.info$log("Sourcing helper functions...")
tryCatch({
  source(here::here("src", "helpers_aa_features.r"))
  source(here::here("src", "helpers_idr.r"))
  source(here::here("src", "helpers_phasesep.r"))
}, error = function(e) {
  .error$log("Error sourcing helper files. Make sure helpers_aa_features.r, helpers_idr.r, and helpers_phasesep.r are in the src/ directory.")
  stop(e)
})

# --- 7. LOAD BASE UNIPROT/PROTEOME DATA ---
# This data is required by the build_datasets.r script
.info$log("Loading base Uniprot data...")
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
  dplyr::mutate(uniprot_len = nchar(uniprot_seq))

.succ$log("Main setup complete.")

# --- 8. QC AND FILE CHECKING FUNCTION ---
# This function is called by analysis scripts
check_data_files <- function(files_needed, build_script_message) {
  files_missing <- files_needed[!file.exists(files_needed)]
  if (length(files_missing) > 0) {
    .error$log("Missing required data files:")
    for (f in files_missing) { .error$log(f) }
    .error$log(build_script_message)
    stop("Missing required data files.", call. = FALSE)
  }
  .succ$log("QC Check: All required data files found.")
  return(TRUE)
}

# --- 9b. DATA LOADING FUNCTIONS ---

load_umap_data <- function() {
  .info$log("Running QC check for UMAP analysis...")
  files_for_umap <- c(
    here::here('data', 'HUMAN_MOBIDB_FEATURES.tsv'),
    here::here('data', 'ATAR_IDR_FEATURES.tsv'),
    here::here('data', 'IDR-FEATURES-DATA.rdata')
  )
  check_data_files(
    files_needed = files_for_umap,
    build_script_message = "Please run src/build_datasets.r to generate missing data files."
  )
  
  .info$log("Loading UMAP data...")
  
  # Load .RData first, which populates the global environment
  .info$log("Loading .RData image for feature definitions...")
  load(here::here('data', 'IDR-FEATURES-DATA.rdata'), envir = .GlobalEnv)
  
  # Load .tsv files and return them
  .info$log("Loading TSV files...")
  hs_idr <- readr::read_tsv(
    here::here('data', 'HUMAN_MOBIDB_FEATURES.tsv'),
    comment = "#"
  )
  .info$log("Loaded HUMAN_MOBIDB_FEATURES.tsv")
  
  atar_idr <- readr::read_tsv(
    here::here('data', 'ATAR_IDR_FEATURES.tsv'),
    comment = "#"
  )
  .info$log("Loaded ATAR_IDR_FEATURES.tsv")
  
  .succ$log("UMAP data loading complete.")
  return(list(HS_IDR = hs_idr, ATAR_IDR = atar_idr))
}

load_csat_data <- function() {
  .info$log("Running QC check for C-SAT analysis...")
  files_for_csat <- c(
    here::here('data', 'ATAR_IDR_FEATURES.tsv'),
    here::here('data', 'IDR-FEATURES-DATA.rdata'),
    here::here('data', 'analysis_umap_output.RData') # Output from analysis_umap.r
  )
  check_data_files(
    files_needed = files_for_csat,
    build_script_message = "Please run src/build_datasets.r and then analysis_umap.r first."
  )
  
  .info$log("Loading C-SAT data...")
  
  # Load .RData files, which populate the global environment
  .info$log("Loading .RData image for feature definitions...")
  load(here::here('data', 'IDR-FEATURES-DATA.rdata'), envir = .GlobalEnv)
  load(here::here('data', 'analysis_umap_output.RData'), envir = .GlobalEnv)
  
  # Load .tsv files and return them
  .info$log("Loading TSV files...")
  atar_idr <- readr::read_tsv(
    here::here('data', 'ATAR_IDR_FEATURES.tsv'),
    comment = "#"
  )
  .info$log("Loaded ATAR_IDR_FEATURES.tsv")
  
  .succ$log("C-SAT data loading complete.")
  return(list(ATAR_IDR = atar_idr))
}


# --- 10. CREATE PLOTS DIRECTORY ---
plots_dir <- here::here('plots')
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
  .info$log(sprintf("Created directory: %s", plots_dir))
}


