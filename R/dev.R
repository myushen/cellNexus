# Utility scripts for development purposes, that are not exported to users

#' Upload a file to the Nectar object store
#' @param source A character scalar indicating the local path to the file to
#'   upload
#' @param container A character scalar indicating the name of the container to
#'   upload to
#' @param name An optional character scalar indicating the name the file should
#'   have after being uploaded. Defaults to being the basename of the source
#'   file.
#' @param credential_id The OpenStack application credential ID as a character
#'   scalar. This is optional because you can alternatively source a
#'   `-openrc.sh` file instead of providing it here.
#' @param credential_id The OpenStack application credential secret as a
#'   character scalar
#' @return `NULL`, invisibly
#' @keywords internal
#' @examples
#' \dontrun{
#' upload_swift(
#'   "/vast/projects/cellxgene_curated/metadata_parquet_0.2",
#'   "cellNexus-metadata",
#'   credential_id = "ABCDEFGHIJK",
#'   credential_secret = "ABCD1234EFGH-5678IJK"
#' )
#' }
upload_swift <- function(
  source,
  container,
  name = basename(source),
  credential_id = NULL,
  credential_secret = NULL
) {
  swift_exe <- Sys.which("swift")
  if (nchar(swift_exe) == 0) {
    stop("'swift' executable not found on PATH. Did you source your openrc.sh?")
  }

  if (!is.null(credential_id) && !is.null(credential_secret)) {
    auth <- c(
      "--os-auth-type",                     "v3applicationcredential",
      "--os-application-credential-id",     credential_id,
      "--os-application-credential-secret", credential_secret
    )
  } else {
    auth <- character()
  }

  args <- c(
    "--os-auth-url", "https://keystone.rc.nectar.org.au:5000/v3/",
    "--os-project-id", "06d6e008e3e642da99d806ba3ea629c5",
    auth,
    "upload",
    container,
    source,
    "--object-name", name
  )

  ret <- system2(swift_exe, args = args)
  if (ret != 0L) stop("swift upload failed with exit code: ", ret)

  invisible(NULL)
}

#' Register a new atlas data version in the remote changelog
#'
#' Downloads the current `atlas_versions.parquet` registry from Nectar, appends
#' a new entry, and re-uploads it. Call this once per data release, after the
#' new metadata parquets and counts files have been uploaded.
#'
#' @param atlas_id A character scalar. The new atlas identifier, e.g.
#'   `"cellxgene_2024/0.1.0"`. The directory portion (before the slash)
#'   encodes the source collection and year; the version portion follows
#'   semantic versioning: major i.e. counts class change, minor = file IDs or
#'   schema change, patch = bug fix.
#' @param census_version A character scalar. The CellxGene Census snapshot this
#'   atlas was built from, e.g. `"01-07-2024"`.
#' @param container A character scalar indicating the name of the container to
#'   upload to
#' @param change_type A character scalar. One of `"initial"`, `"patch"`,
#'   `"minor"`, or `"major"`.
#' @param description A character scalar. Free-text summary of what changed in
#'   this release.
#' @inheritDotParams upload_swift credential_id credential_secret
#' @return `NULL`, invisibly
#' @keywords internal
#' @noRd
#' @examples
#' \dontrun{
#' register_atlas_version(
#   atlas_id = "cellxgene_2024/0.1.0",
#   census_version = "01-07-2024",
#   change_type = "initial",
#   description = "Initial release linked to CellxGene Census 01-07-2024.",
#   credential_id = "ABCDEFGHIJK",
#   credential_secret = "ABCD1234EFGH-5678IJK"
# )
#' }
register_atlas_version <- function(
  atlas_id,
  census_version,
  container,
  change_type,
  description,
  credential_id = NULL,
  credential_secret = NULL
) {
  registry_name <- "atlas_versions.parquet"
  local_path <- file.path(tempdir(), registry_name)

  registry_url <- paste0(
    "https://object-store.rc.nectar.org.au/v1/",
    "AUTH_06d6e008e3e642da99d806ba3ea629c5/",
    container, "/", registry_name
  )

  # Download the existing registry; silently skip if it does not exist yet
  tryCatch(
    sync_remote_file(registry_url, local_path, overwrite = TRUE),
    error = function(e) NULL
  )

  new_row <- tibble::tibble(
    atlas_id = atlas_id,
    census_version = census_version,
    change_type = change_type,
    description = description,
    modified_at = as.character(Sys.Date())
  )
  existing <- if (file.exists(local_path)) {
    arrow::read_parquet(local_path)
  } else {
    new_row[0L, ]
  }
  
  updated <- dplyr::bind_rows(existing, new_row)
  arrow::write_parquet(updated, local_path)

  upload_swift(
    local_path,
    container         = "cellNexus-metadata",
    name              = registry_name,
    credential_id     = credential_id,
    credential_secret = credential_secret
  )
  
  invisible(NULL)
}

#' Update the unharmonised parquet files
#' @param unharmonised_parquet_dir The path to a directory containing parquet
#'   files, one for each dataset, e.g.
#'   /vast/projects/cellxgene_curated/metadata_non_harmonised_parquet_0.2
#' @inheritDotParams upload_swift
#' @inherit upload_swift return
#' @keywords internal
#' @noRd
#' @examples
#' \dontrun{
#' update_unharmonised(
#'   "/vast/projects/cellxgene_curated/metadata_non_harmonised_parquet_0.2",
#'   credential_id = "ABCDEFGHIJK",
#'   credential_secret = "ABCD1234EFGH-5678IJK"
#' )
#' }
update_unharmonised <- function(unharmonised_parquet_dir, ...) {
  # name="/" forces it have no prefix, ie be at the top level in the bucket
  upload_swift(
    unharmonised_parquet_dir,
    container = "unharmonised_metadata",
    name = "/",
    ...
  )
}

#' Converts a series of HDF5Array-serialized SingleCellExperiments to AnnData
#' @param input_directory A character scalar. The path to a directory containing one or more
#'  directories created by [HDF5Array::saveHDF5SummarizedExperiment()].
#' @param output_directory A character scalar. The path to a directory in which to save the
#'  created anndata files.
#' @keywords internal
#' @noRd
#' @return A character vector of the newly-created anndata files
#' @examples
#' \dontrun{
#' hdf5_to_anndata(
#'   "/vast/projects/cellxgene_curated/splitted_DB2_data_0.2.1",
#'   "/vast/projects/cellxgene_curated/splitted_DB2_anndata_0.2.1"
#' )
#' hdf5_to_anndata(
#'   "/vast/projects/cellxgene_curated/splitted_DB2_data_scaled_0.2.1",
#'   "/vast/projects/cellxgene_curated/splitted_DB2_anndata_scaled_0.2.1"
#' )
#' }
hdf5_to_anndata <- function(input_directory, output_directory) {
  dir.create(output_directory, showWarnings = FALSE)
  # This is a quick utility script to convert the SCE files into AnnData format for use in Pythonlist.files("/vast/projects/RCP/human_cell_atlas/splitted_DB2_data", full.names = FALSE) |>  purrr::walk(function(dir){
  basilisk::basiliskRun(fun = function(sce) {
    list.dirs(input_directory)[-1] |>
      purrr::map_chr(function(sce_dir) {
        cli::cli_alert_info("Processing {sce_dir}.")
        prefix <- basename(sce_dir)
        out_path <- glue::glue("{prefix}.h5ad") |>
          file.path(output_directory, name = _)

        if (file.exists(out_path)) {
          cli::cli_alert_info("{out_path} already exists. Skipping")
        } else {
          sce <- HDF5Array::loadHDF5SummarizedExperiment(sce_dir)
          single_column <- length(colnames(sce)) == 1
          if (single_column) {
            # Hack, so that single-column SCEs will convert
            # correctly
            cli::cli_alert_info(
              "{sce_dir} has only 1 column. Duplicating column."
            )
            sce <- cbind(sce, sce)
            single_column <- TRUE
          }
          ad <- zellkonverter::SCE2AnnData(sce)
          if (single_column) {
            # Remove the duplicate column
            sce$X <- sce$X[1]
          }
          # TODO: customize chunking here, when anndata supports it
          # (see https://github.com/scverse/anndata/issues/961)
          ad$write_h5ad(out_path)
        }
        out_path
      }, .progress = "Converting files")
  }, env = zellkonverter::zellkonverterAnnDataEnv())
}

#' Makes a "downsampled" metadata file that only contains the minimal data
#' needed to run the vignette and unit tests
#' @keywords internal
#' @noRd
#' @param output Character scalar. Path to the output file.
#' @return NULL
downsample_metadata <- function(output = "sample_metadata.2.1.0.parquet") {
  metadata <- get_metadata() |>
    join_census_table()

  # Make a table of rows per dataset
  dataset_sizes <- metadata |>
    dplyr::group_by(.data$file_id_cellNexus_single_cell) |>
    summarise(n = dplyr::n()) |>
    dplyr::collect()

  # Find a minimal set of file_id_dbs we need
  minimal_file_ids <- rlang::exprs(
    # Used by the vignette
    .data$self_reported_ethnicity == "African" &
      stringr::str_like(.data$assay, "%10x%") &
      .data$tissue == "lung parenchyma" &
      stringr::str_like(.data$cell_type, "%CD4%"),
    .data$cell_type_unified_ensemble == "nk",
    .data$cell_type_unified_ensemble == "cd14 mono",
    .data$tissue == "kidney blood vessel",
    .data$file_id_cellNexus_single_cell == "e52795dec7b626b6276b867d55328d9f___1.h5ad",
    # Used by tests
    .data$file_id_cellNexus_single_cell == "4164d0eb972ad5e12719b6858c9559ea___1.h5ad",
    .data$file_id_cellNexus_single_cell == "7ddd6775d704d6826539abaee8d22f65___1.h5ad",
    .data$file_id_cellNexus_single_cell == "6f0d7a93cff864a1cf029a28e92057e6___1.h5ad",
    .data$file_id_cellNexus_single_cell == "51250a91a93bfe3c6ba22c9f7d8e04ee___1.h5ad",
    .data$file_id_cellNexus_single_cell == "490bfce7a70edde5d4fc6352374ead0c___1.h5ad",
    .data$file_id_cellNexus_single_cell == "25807d1aa8bb7fc63d587a64361fd2db___1.h5ad",
    .data$file_id_cellNexus_single_cell == "0ad7715d6d9052559d660aaa61bd69d2___1.h5ad",
    .data$file_id_cellNexus_single_cell == "d2ab9251c9670f0de3657130aa79d7e6___1.h5ad",
    .data$file_id_cellNexus_single_cell == "88a378f0cd39ed31bd4b83ee4d2fba5a___1.h5ad",
    .data$file_id_cellNexus_single_cell == "5257022969fc37563e82b608d3908565___1.h5ad",
    .data$file_id_cellNexus_single_cell == "1d90787ac66a5ed18222e2a9851eabbc___1.h5ad",
    .data$file_id_cellNexus_single_cell == "4414dffc701125c467adad7977adcf21___1.h5ad",
  ) |>
    purrr::map(function(filter) {
      all_ids <- metadata |>
        dplyr::filter(!!filter) |>
        dplyr::group_by(.data$file_id_cellNexus_single_cell) |>
        dplyr::pull(.data$file_id_cellNexus_single_cell) |>
        unique()

      dataset_sizes |>
        dplyr::filter(.data$file_id_cellNexus_single_cell %in% all_ids) |>
        dplyr::slice_min(n = 50, order_by = .data$n) |>
        dplyr::pull(.data$file_id_cellNexus_single_cell)
    }) |>
    purrr::reduce(union)

  metadata |>
    dplyr::filter(.data$file_id_cellNexus_single_cell %in% minimal_file_ids) |>
    dplyr::arrange(.data$file_id_cellNexus_single_cell, .data$sample_id) |>
    dplyr::collect() |>
    arrow::write_parquet(output)

  NULL
}
