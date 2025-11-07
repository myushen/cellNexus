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

# =============================================================================
# Create cell-level metadata
# =============================================================================

save_directory = "/vast/scratch/users/shen.m/htan/"
# Get all saved h5ad files
h5ad_files <- list.files(save_directory, pattern = "\\.h5ad$", full.names = TRUE)

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
    as_tibble() |> 
    mutate(
      cell_id = rownames(colData(sce)),
      sample_id = Biospecimen,
      # THERE SHOULD BE ANOTHE COLUMN IN THE METADATA TO SPECIFIY HTA DATA
      file_id_cellNexus_single_cell = paste0(sapply(counts_path, digest::digest), ".h5ad")
    ) |>
    select(cell_id, sample_id, file_id_cellNexus_single_cell)
  
  col_data
}

cell_metadata <- map_dfr(h5ad_files, get_h5ad_cell_metadata, .progress = T)
cell_metadata

# Save cell metadata
arrow::write_parquet(
  cell_metadata,
  file.path(save_directory, "hta_cell_metadata.parquet")
)


# =============================================================================
# Create Sample, Donor metadata
# =============================================================================
file_metadata <- read.csv("/home/users/allstaff/shen.m/projects/HTAN/files_metadata_2025_10_21.tsv",
                          sep = "\t", na.strings = c("NA",""), header = TRUE) |>
  

  mutate(
    Filename_basename = basename(Filename),
    channel_number = str_extract(Filename_basename, "(?<=channel)[0-9]+"),
    sample_id = if_else(
      is.na(channel_number),
      Biospecimen,
      paste(Biospecimen, channel_number, sep = "___") # Because there are cases One Biospecimen match more than one counts (Channel)
    )) |>
  filter(str_detect(Filename_basename, "molecule_counts.mtx$|_matrix.mtx.gz$") ) |> 
  mutate(
    file_id_cellNexus_single_cell = paste0(sapply(Filename_basename, digest::digest), ".h5ad")
  ) |> 
  select(sample_id, file_id_cellNexus_single_cell, Biospecimen, channel_number ) |> 
  as_tibble()


samples_metadata <-  read.csv(
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
  samples_metadata,
  by = c("Biospecimen" = "HTAN.Biospecimen.ID"),
  copy = TRUE
)

sample_metadata <- sample_metadata |> 
  left_join(
    donors_metadata, 
    by = c("Participant.ID" = "HTAN.Participant.ID", "Atlas.Name", "Publications", "Synapse.ID",
           "Atlas.ID", "Level")
  )

# Save Sample metadata
arrow::write_parquet(
  sample_metadata,
  file.path(save_directory, "hta_sample_metadata.parquet")
)

