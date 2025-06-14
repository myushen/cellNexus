# Functions that relate to the harmonised metadata database

#' @include utils.R
NULL

#' Environment that we use to cache the DuckDB connections
#' @noRd
cache <- rlang::env(
  metadata_table = rlang::env()
)

#' Returns the URLs for all metadata files 
#' @param databases A character vector specifying the names of the metadata files. 
#'   Download the specific metadata by defining the metadata version. The default is 
#'   metadata.1.0.12.parquet
#' @export
#' @return A character vector of URLs to parquet files to download
#' @examples
#' get_metadata_url("metadata.1.0.12.parquet")
#' @references Mangiola, S., M. Milton, N. Ranathunga, C. S. N. Li-Wai-Suen, 
#'   A. Odainic, E. Yang, W. Hutchison et al. "A multi-organ map of the human 
#'   immune system across age, sex and ethnicity." bioRxiv (2023): 2023-06.
#'   doi:10.1101/2023.06.08.542671.
#' @source [Mangiola et al.,2023](https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3)
get_metadata_url <- function(databases = c("metadata.1.0.12.parquet")) {
  clear_old_metadata(updated_data = databases)
  glue::glue(
    "https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/cellNexus-metadata/{databases}")
}

#' URL pointing to the sample metadata file, which is smaller and for test,
#' demonstration, and vignette purposes only
#' @export
#' @return A character scalar consisting of the URL
#' @examples
#' get_metadata(cloud_metadata = SAMPLE_DATABASE_URL, cache_directory = tempdir())
#' @references Mangiola, S., M. Milton, N. Ranathunga, C. S. N. Li-Wai-Suen, 
#'   A. Odainic, E. Yang, W. Hutchison et al. "A multi-organ map of the human 
#'   immune system across age, sex and ethnicity." bioRxiv (2023): 2023-06.
#'   doi:10.1101/2023.06.08.542671.
#' @source [Mangiola et al.,2023](https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3)
SAMPLE_DATABASE_URL <- single_line_str(
  "https://object-store.rc.nectar.org.au/v1/
    AUTH_06d6e008e3e642da99d806ba3ea629c5/metadata/
    sample_metadata.0.2.3.parquet"
)

#' Gets the Curated Atlas metadata as a data frame.
#'
#' Downloads a parquet database of the Human Cell Atlas metadata to a local
#' cache, and then opens it as a data frame. It can then be filtered and passed
#' into [get_single_cell_experiment()] to obtain a
#' [`SingleCellExperiment::SingleCellExperiment-class`]
#'
#' @param cloud_metadata Optional character vector of any length. HTTP URL/URLs pointing
#'   to the name and location of parquet database/databases. By default, it points to 
#'   cellNexus ARDC Nectar Research Cloud. Assign `NULL` to query local_metadata only 
#'   if exists. 
#' @param local_metadata Optional character vector of any length representing the local
#'   path of parquet database(s).
#' @param cache_directory Optional character vector of length 1. A file path on
#'   your local system to a directory (not a file) that will be used to store
#'   `metadata.parquet`
#' @param use_cache Optional logical scalar. If `TRUE` (the default), and this
#'   function has been called before with the same parameters, then a cached
#'   reference to the table will be returned. If `FALSE`, a new connection will
#'   be created no matter what.
#' @return A lazy data.frame subclass containing the metadata. You can interact
#'   with this object using most standard dplyr functions. For string matching,
#'   it is recommended that you use `stringr::str_like` to filter character
#'   columns, as `stringr::str_match` will not work.
#' @export
#' @examples
#' library(dplyr)
#' filtered_metadata <- get_metadata() |>
#'     filter(
#'         ethnicity == "African" &
#'             assay %LIKE% "%10x%" &
#'             tissue == "lung parenchyma" &
#'             cell_type %LIKE% "%CD4%"
#'     )
#'
#' @importFrom DBI dbConnect
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl
#' @importFrom httr progress
#' @importFrom cli cli_alert_info hash_sha256
#' @importFrom glue glue
#' @importFrom purrr walk
#'
#' @details
#'
#' The metadata was collected from the Bioconductor package `cellxgenedp`.
#' `vignette("using_cellxgenedp", package="cellxgenedp")` provides an overview of the columns in the
#' metadata. The data for which the column `organism_name` included "Homo
#' sapiens" was collected collected from `cellxgenedp`.
#'
#' The columns `dataset_id` and `file_id_cellNexus_single_cell` link the datasets explorable through
#' `cellNexus` and `cellxgenedp`to the CELLxGENE portal.
#'
#' Our representation, harmonises the metadata at dataset, sample and cell
#' levels, in a unique coherent database table.
#'
#' Dataset-specific columns (definitions available at cellxgene.cziscience.com):
#' `cell_count`, `collection_id`, `created_at.x`, `created_at.y`,
#' `dataset_deployments`, `dataset_id`, `file_id_cellNexus_single_cell`, `filename`, `filetype`,
#' `is_primary_data.y`, `is_valid`, `linked_genesets`, `mean_genes_per_cell`,
#' `name`, `published`, `published_at`, `revised_at`, `revision`, `s3_uri`,
#' `schema_version`, `tombstone`, `updated_at.x`, `updated_at.y`,
#' `user_submitted`, `x_normalization`
#'
#' Sample-specific columns (definitions available at cellxgene.cziscience.com):
#' `sample_id`, `.sample_name`, `age_days`, `assay`, `assay_ontology_term_id`,
#' `development_stage`, `development_stage_ontology_term_id`, `ethnicity`,
#' `ethnicity_ontology_term_id`, `experiment___`, `organism`,
#' `organism_ontology_term_id`, `sample_placeholder`, `sex`,
#' `sex_ontology_term_id`, `tissue`, `tissue_harmonised`,
#' `tissue_ontology_term_id`, `disease`, `disease_ontology_term_id`,
#' `is_primary_data.x`
#'
#' Cell-specific columns (definitions available at cellxgene.cziscience.com):
#' `cell_id`, `cell_type`, `cell_type_ontology_term_idm`, `cell_type_harmonised`,
#' `confidence_class`, `cell_annotation_azimuth_l2`,
#' `cell_annotation_blueprint_singler`
#'
#' Through harmonisation and curation we introduced custom columns not present
#' in the original CELLxGENE metadata:
#'
#' - `tissue_harmonised`: a coarser tissue name for better filtering
#' - `age_days`: the number of days corresponding to the age
#' - `cell_type_harmonised`: the consensus call identity (for immune cells)
#'   using the original and three novel annotations using Seurat Azimuth and 
#'   SingleR
#' - `confidence_class`: an ordinal class of how confident
#'   `cell_type_harmonised` is. 1 is complete consensus, 2 is 3 out of four and
#'   so on.
#' - `cell_annotation_azimuth_l2`: Azimuth cell annotation
#' - `cell_annotation_blueprint_singler`: SingleR cell annotation using 
#'   Blueprint reference
#' - `cell_annotation_blueprint_monaco`: SingleR cell annotation using Monaco 
#'   reference
#' - `sample_id_db`: Sample subdivision for internal use
#' - `file_id_db`: File subdivision for internal use
#' - `sample_id`: Sample ID
#' - `.sample_name`: How samples were defined
#'
#' **Possible cache path issues**
#'
#' If your default R cache path includes non-standard characters (e.g. dash
#' because of your user or organisation name), the following error can occur.
#'
#' ```
#' Error in `db_query_fields.DBIConnection()`: ! Can't query fields. Caused by
#' error: ! Parser Error: syntax error at or near "/" LINE 2: FROM
#' /Users/bob/Library/Caches...
#' ```
#'
#' The solution is to choose a different cache, for example
#' ```R
#' get_metadata(cache_directory = path.expand('~'))
#' ```
#' 
#' @inheritDotParams read_parquet
#' 
#' @references Mangiola, S., M. Milton, N. Ranathunga, C. S. N. Li-Wai-Suen, 
#'   A. Odainic, E. Yang, W. Hutchison et al. "A multi-organ map of the human 
#'   immune system across age, sex and ethnicity." bioRxiv (2023): 2023-06.
#'   doi:10.1101/2023.06.08.542671.
#' @source [Mangiola et al.,2023](https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3)
get_metadata <- function(
    cloud_metadata = get_metadata_url(),
    local_metadata = NULL,
    cache_directory = get_default_cache_dir(),
    use_cache = TRUE,
    ...
) {
  # Synchronize remote files
  walk(cloud_metadata, function(url) {
    # Calculate the file path from the URL
    path <- file.path(cache_directory, url |> basename())
    if (!file.exists(path)) {
      report_file_sizes(url)
      sync_remote_file(url,
                       path,
                       progress(type = "down", con = stderr()))
    }
  })
  
  if (is.null(cloud_metadata)) all_parquet <- c(local_metadata)
  if (!is.null(cloud_metadata)) all_parquet <- c(file.path(cache_directory, 
                                                           cloud_metadata |> basename()), 
                                                 local_metadata)
  
  # We try to avoid re-reading a set of parquet files 
  # that is identical to a previous set by hashing the file list
  hash <- all_parquet |> paste0(collapse="") |>
    hash_sha256()
  cached_connection <- cache$metadata_table[[hash]]
  
  if (!is.null(cached_connection) && isTRUE(use_cache)) {
    # If the file list is identical, just re-use the database table
    cached_connection
  }
  else {
    table <- duckdb() |>
      dbConnect(drv = _, read_only = TRUE) |>
      read_parquet(path = all_parquet, ...)
    cache$metadata_table[[hash]] <- table
    table
  }
}

