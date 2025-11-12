# =============================================================================
# Generate HTAN Cell-Level and Sample-Level Metadata
# =============================================================================
#
# This script generates comprehensive metadata files for Human Tumor Atlas Network
# (HTAN) single-cell data. It processes h5ad files created by ~/git_control/cellNexus/dev/HTA_parse_one_count_to_sce.R
# and HTAN metadata files to create two output files:
#
# 1. Cell-level metadata (hta_cell_metadata.parquet):
#    - Extracts cell information from all saved h5ad files
#    - Creates three key columns:
#      * cell_id: Unique cell identifier (format: barcode___Biospecimen)
#      * sample_id: Biospecimen identifier (or Biospecimen___channel for multi-channel)
#      * file_id_cellNexus_single_cell: MD5 hash of counts_path + ".h5ad"
#    - Memory efficient: Only reads colData from h5ad files, not expression data
#
# 2. Sample-level metadata (hta_sample_metadata.parquet):
#    - Combines file metadata, samples metadata, and donors metadata
#    - Handles multi-channel data (one Biospecimen can have multiple channels)
#    - Creates sample_id that matches cell-level metadata
#    - Includes all annotations from HTAN samples and donors metadata
#
# Key features:
# - Processes all h5ad files in the save directory
# - Handles multi-channel samples (combines channel information)
# - Creates file_id using MD5 hash for consistent file identification
# - Joins sample and donor-level annotations
# - Saves as parquet format for efficient storage and querying
#
# Output files:
#   - hta_cell_metadata.parquet: Cell-level metadata with cell_id, sample_id, file_id
#   - hta_sample_metadata.parquet: Sample-level metadata with all HTAN annotations
library(targets)
library(dplyr)
library(duckdb)
# =============================================================================
# Create cell-level metadata
# =============================================================================

store_file_hta_cell_metadata <- "/vast/scratch/users/shen.m/htan/hta_lung_scrnaseq_cell_metadata_target_store"

tar_script({
  library(ensembldb) 
  library(dplyr)
  library(purrr)
  library(stringr)
  library(fs)
  library(readr)
  library(Matrix)
  library(R.utils)
  library(tidyr)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(zellkonverter)
  library(Seurat)
  library(tibble)
  library(targets)
  library(crew)
  library(crew.cluster)
  
  tar_option_set(
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue",
    workspace_on_error = TRUE,
    cue = tar_cue(mode = "thorough"),
    controller = crew_controller_slurm(
      name = "elastic",
      workers = 300,
      tasks_max = 20,
      seconds_idle = 30,
      crashes_error = 10,
      options_cluster = crew_options_slurm(
        memory_gigabytes_required = c(20, 40, 60, 100, 150), 
        cpus_per_task = c(2),
        time_minutes = c(60*24),
        verbose = T
      )
    )
  )
  
# Extract cell metadata from each h5ad file
  get_h5ad_cell_metadata <- function(h5ad_file) {
    # Read h5ad file (only colData to be memory efficient)
    sce <- readH5AD(
      file = h5ad_file,
      reader = "R",
      use_hdf5 = TRUE,
      obs = TRUE,  # Read cell metadata
      var = FALSE, # Don't read gene metadata
      raw = FALSE, # Don't read raw data
      layers = FALSE # Don't read layers
    )
    
    # Extract metadata
    col_data <- colData(sce) |>
      as.data.frame(row.names = NULL) |>
      tibble::rownames_to_column("cell_id") |>
      as_tibble() |> 
      mutate(
        # THERE SHOULD BE ANOTHE COLUMN IN THE METADATA TO SPECIFIY HTA DATA
        file_id_cellNexus_single_cell = paste0(sample_id, ".h5ad")
      ) |>
      dplyr::select(cell_id, sample_id, file_id_cellNexus_single_cell)
    
    col_data
  }
  
  list(
    tar_target(
      h5ad_files,
      # Get all saved h5ad files
      list.files("/vast/scratch/users/shen.m/htan/09-11-2025/counts/", pattern = "\\.h5ad$", full.names = TRUE)
      ),
    tar_target(
      cell_data_list,
      get_h5ad_cell_metadata(h5ad_files),
      pattern = map(h5ad_files),
      iteration = "list"
    )
  )
}, script = paste0(store_file_hta_cell_metadata, "_target_script.R"), ask = FALSE)


job::job({
  
  tar_make(
    script = paste0(store_file_hta_cell_metadata, "_target_script.R"), 
    store = store_file_hta_cell_metadata, 
    reporter = "summary" #, callr_function = NULL
  )
  
})

cell_metadata <- tar_read(cell_data_list,  store = store_file_hta_cell_metadata) |> bind_rows()

# Save cell metadata
save_directory <- "/vast/scratch/users/shen.m/htan/"
arrow::write_parquet(
  cell_metadata,
  file.path(save_directory, "hta_cell_metadata.parquet")
)

cell_metadata <-  tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/htan/hta_cell_metadata.parquet')")
)



# =============================================================================
# Create Sample, Donor metadata
# =============================================================================
file_metadata <- read.csv("/home/users/allstaff/shen.m/projects/HTAN/files_metadata_2025_10_21.tsv",
                          sep = "\t", na.strings = c("NA",""), header = TRUE) |>
  

  mutate(
    Filename_basename = basename(Filename),
    channel_number = str_extract(Filename_basename, "channel[0-9]+"),
    sample_id = if_else(
      is.na(channel_number),
      Biospecimen,
      paste(Biospecimen, channel_number, sep = "___") # Because there are cases One Biospecimen match more than one counts (Channel)
    )) |>
  
  # Because only counts from these files are generated 
  filter(str_detect(Filename_basename, "molecule_counts.mtx$|_matrix.mtx.gz$") ) |> 
  mutate(
    file_id_cellNexus_single_cell = paste0(sample_id, ".h5ad")
  ) |> 
  dplyr::select(sample_id, file_id_cellNexus_single_cell, Biospecimen, channel_number, Assay, Organ
                #, Treatment
                ) |> 
  as_tibble()


sample_metadata <-  read.csv(
  "/home/users/allstaff/shen.m/projects/HTAN/samples_metadata_2025_10_21.tsv",
  sep = "\t",
  na.strings = c("NA", ""),
  header = TRUE
)
 
donors_metadata <- read.csv(
  "/home/users/allstaff/shen.m/projects/HTAN/donors_metadata_2025_10_21.tsv",
  sep = "\t",
  na.strings = c("NA", ""),
  header = TRUE
)


# Joining -----------------------------------------------------------------
sample_metadata <- file_metadata |> left_join(
  sample_metadata,
  by = c("Biospecimen" = "HTAN.Biospecimen.ID"),
  copy = TRUE
)

sample_metadata <- sample_metadata |> 
  left_join(
    donors_metadata, 
    by = c("Participant.ID" = "HTAN.Participant.ID", "Atlas.Name", "Publications", "Synapse.ID",
           "Atlas.ID", "Level")
  ) |> 
  dplyr::rename(
    donor_id = Participant.ID,
    diagnosis_age = Age.at.Diagnosis..years.,
    title = Publications,
    center = Atlas.Name,
    center_id = Atlas.ID,
    sex = Gender,
    tissue = Organ,
    url = Protocol.Link,
    primary_diagnosis = Primary.Diagnosis,
    assay = Assay
    
  ) |> 
  mutate(self_reported_ethnicity = ifelse(!is.na(Race) & !is.na(Ethnicity) , paste(Race, Ethnicity, sep = "___"), NA),
         sex = tolower(sex),
         age_days = diagnosis_age * 365,
         tissue = tolower(tissue),
         atlas_id = "hta/09-11-2025" # For directory hierarchy
         )



# Subset necessary columns ------------------------------------------------
cols <- c(
  "sample_id",
  "donor_id",
  "self_reported_ethnicity",
  "age_days",
  "tissue",
  "sex",
  "assay",
  "file_id_cellNexus_single_cell",
  "atlas_id", 
  #"Biospecimen",
 # "HTAN.Parent.ID",
 # "Source.HTAN.Biospecimen.ID",
  "diagnosis_age",
  "primary_diagnosis",
  "url",
  "center",
  "center_id",
  "title",
  "Last.Known.Disease.Status",
  "Progression.or.Recurrence.Type",
  "Vital.Status",
  "Cause.of.Death",
  "Year.Of.Birth",
  "Year.of.Death",
  "Year.Of.Diagnosis",
  "Tumor.Grade",
  "Progression.or.Recurrence",
  "Tissue.or.Organ.of.Origin"
)

sample_metadata <- sample_metadata |> dplyr::select(any_of(cols))
#sample_metadata |> glimpse()
# Save Sample metadata
arrow::write_parquet(
  sample_metadata,
  file.path(save_directory, "hta_sample_metadata.parquet")
)


sample_metadata <-  tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/htan/hta_sample_metadata.parquet')")
) 

sample_metadata |> distinct(sample_id) |> dplyr::count()



# Testing with cellNexus --------------------------------------------------
library(cellNexus)
library(dplyr)
library(testthat)
library(duckdb)
test_that("get_metadata and get_single_cell_experiment return expected SCE for test sample", {
  
  cache <- "/vast/scratch/users/shen.m/htan/"
  save_directory <- tempdir()   # or your defined directory
  test_sample_id <- c("HTA1_203_332101___channel3" ,
                      "dc1a2e1504a4b71427b682a6300d02d3___1.h5ad") # One of cellNexus file
  
  cell_metadata <-  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/htan/hta_cell_metadata.parquet')")
  )
  
  sample_metadata <-  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/scratch/users/shen.m/htan/hta_sample_metadata.parquet')")
  ) 
  
  # Build test data
  test_data <- cell_metadata |>
    dplyr::filter(sample_id %in% test_sample_id) |>
    dplyr::left_join(
      sample_metadata,
      by = c("sample_id", "file_id_cellNexus_single_cell"),
      copy = TRUE
    )
  
  # Save parquet test metadata
  test_data |>
    dplyr::collect() |>
    arrow::write_parquet(
      file.path(save_directory, "test_metadata.parquet")
    )
  
  # Read metadata and construct SCE
  sce <- get_metadata(
    cache_directory = cache,
    local_metadata = file.path(save_directory, "test_metadata.parquet")
  ) |>
    dplyr::filter(sample_id %in% test_sample_id ) |>
    get_single_cell_experiment(cache_directory = cache)
  
  # Basic structural checks
  expect_s4_class(sce, "SingleCellExperiment")
  expect_true(ncol(sce) > 0)
  expect_true(nrow(sce) > 0)
  
  # Check assay column exists
  assay_counts <- sce |> dplyr::count(assay)
  expect_true("assay" %in% colnames(colData(sce)))
  expect_s3_class(assay_counts, "data.frame")
  
})
