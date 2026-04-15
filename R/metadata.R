# Functions that relate to the harmonised metadata database

#' Environment that we use to cache the DuckDB connections
#' @keywords internal
#' @noRd
cache <- rlang::env(
  metadata_table = rlang::env()
)

#' Returns the URLs for all metadata files
#' @param databases A character vector specifying the names of the metadata files.
#'   Download the specific metadata by defining the metadata version.
#' @export
#' @return A character vector of URLs to parquet files to download
#' @examples
#' get_metadata_url("cellnexus_metadata.2.2.0.parquet")
#' @references Mangiola, S., M. Milton, N. Ranathunga, C. S. N. Li-Wai-Suen,
#'   A. Odainic, E. Yang, W. Hutchison et al. "A multi-organ map of the human
#'   immune system across age, sex and ethnicity." bioRxiv (2023): 2023-06.
#'   doi:10.1101/2023.06.08.542671.
#' @source [Mangiola et al.,2023](https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3)
get_metadata_url <- function(databases = c(
                               "cellnexus_metadata.2.2.0.parquet",
                               "census_cell_metadata.2.2.0.parquet"
                             )) {
  glue::glue(
    "https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/cellNexus-metadata/{databases}"
  )
}

#' URL pointing to the sample metadata file, which is smaller and for test,
#' demonstration, and vignette purposes only
#' @export
#' @return Character scalar consisting of the URL/URLs
#' @examples
#' get_metadata(cloud_metadata = SAMPLE_DATABASE_URL["cellnexus"], cache_directory = tempdir())
#' @references Mangiola, S., M. Milton, N. Ranathunga, C. S. N. Li-Wai-Suen,
#'   A. Odainic, E. Yang, W. Hutchison et al. "A multi-organ map of the human
#'   immune system across age, sex and ethnicity." bioRxiv (2023): 2023-06.
#'   doi:10.1101/2023.06.08.542671.
#' @source [Mangiola et al.,2023](https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3)
SAMPLE_DATABASE_URL <- c(
  cellnexus = paste0(
    "https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/",
    "cellNexus-metadata/cellnexus_sample_metadata.2.2.0.parquet"
  ),
  census = paste0(
    "https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/",
    "cellNexus-metadata/census_sample_metadata.2.2.0.parquet"
  )
)

#' Gets the CellNexus metadata as a data frame.
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
#'   your local system to a directory (not a file) that will be used to store metadata.
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
#' # For fast build purpose only, you do not need to specify anything in cloud_metadata.
#' filtered_metadata <- get_metadata(cloud_metadata = SAMPLE_DATABASE_URL) |>
#'   filter(
#'     self_reported_ethnicity == "African" &
#'       assay %LIKE% "%10x%" &
#'       tissue == "lung parenchyma" &
#'       cell_type %LIKE% "%CD4%"
#'   )
#'
#' @importFrom DBI dbConnect
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl
#' @importFrom httr progress
#' @importFrom cli cli_alert_info hash_sha256
#' @importFrom glue glue glue_sql
#' @importFrom purrr walk
#' @importFrom dbplyr sql
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
#' @references Mangiola, S., M. Milton, N. Ranathunga, C. S. N. Li-Wai-Suen,
#'   A. Odainic, E. Yang, W. Hutchison et al. "A multi-organ map of the human
#'   immune system across age, sex and ethnicity." bioRxiv (2023): 2023-06.
#'   doi:10.1101/2023.06.08.542671.
#' @source [Mangiola et al.,2023](https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3)
get_metadata <- function(cloud_metadata = get_metadata_url("cellnexus_metadata.2.2.0.parquet"),
                         local_metadata = NULL,
                         cache_directory = get_default_cache_dir(),
                         use_cache = TRUE) {
  # Synchronize remote files using parallel downloads
  if (!is.null(cloud_metadata) && length(cloud_metadata) > 0) {
    output_paths <- file.path(cache_directory, basename(cloud_metadata))
    # Filter to only files that need downloading
    to_download <- !file.exists(output_paths)
    if (any(to_download)) {
      report_file_sizes(cloud_metadata[to_download])
      sync_remote_files(cloud_metadata[to_download],
        output_paths[to_download],
        progress = TRUE
      )
    }
  }

  if (is.null(cloud_metadata)) all_parquet <- c(local_metadata)
  if (!is.null(cloud_metadata)) {
    all_parquet <- c(
      file.path(
        cache_directory,
        cloud_metadata |>
          basename()
      ),
      local_metadata
    )
  }

  # We try to avoid re-reading a set of parquet files
  # that is identical to a previous set by hashing the file list
  hash <- all_parquet |>
    paste0(collapse = "") |>
    hash_sha256()
  cached_connection <- cache$metadata_table[[hash]]

  if (!is.null(cached_connection) && isTRUE(use_cache)) {
    # If the file list is identical, just re-use the database table
    cached_connection
  } else {
    table <- duckdb() |>
      dbConnect(drv = _, read_only = TRUE) |>
      duckdb_read_parquet(path = all_parquet)
    cache$metadata_table[[hash]] <- table
    table
  }
}

#' Retrieve cellNexus cell communication ligand–receptor strength as a data frame.
#'
#' Downloads a parquet database of the cell communication strength to a local
#' cache, and then opens it as a data frame. It can then be filtered.
#'
#' @param cloud_metadata Character vector of any length. HTTP URL/URLs pointing
#'   to the name and location of parquet database/databases. By default, it points to
#'   cell communication metadata in cellNexus ARDC Nectar Research Cloud.
#'   Assign `NULL` to query local_metadata only if exists.
#' @param local_metadata Optional character vector of any length representing the local
#'   path of parquet database(s).
#' @param cache_directory Optional character vector of length 1. A file path on
#'   your local system to a directory (not a file) that will be used to store
#'   `metadata.parquet`
#' @param use_cache Optional logical scalar. If `TRUE` (the default), and this
#'   function has been called before with the same parameters, then a cached
#'   reference to the table will be returned. If `FALSE`, a new connect138/4.7ion will
#'   be created no matter what.
#' @return A lazy data.frame subclass containing the metadata. You can interact
#'   with this object using most standard dplyr functions. For string matching,
#'   it is recommended that you use `stringr::str_like` to filter character
#'   columns, as `stringr::str_match` will not work.
#' @details
#' The returned table integrates three levels of cell communication
#'   inference from `CellChat`, for each sample:
#'   (i) ligand–receptor–level communication (\code{lr_prob}, \code{lr_pval}),
#'   (ii) pathway-level aggregated signaling (\code{pathway_prob}, \code{pathway_pval}),
#'   (iii) cell-pair–level summaries of communication breadth
#'   (\code{interaction_count} - number of significant LR interactions) and
#'   intensity (\code{interaction_weight} - overall communication strength).
#'
#'   Together, these metrics allow simultaneous assessment of signaling specificity,
#'   pathway dominance, and global communication structure between cell populations.
#' @examples
#' # For fast build purpose only, you do not need to specify anything in the function.
#' communication_meta <- get_cell_communication_strength(cloud_metadata = SAMPLE_DATABASE_URL)
#' @export
get_cell_communication_strength <- function(
  cloud_metadata = get_metadata_url("cellNexus_lr_signaling_pathway_strength.parquet"),
  local_metadata = NULL,
  cache_directory = get_default_cache_dir(),
  use_cache = TRUE
) {
  get_metadata(cloud_metadata, local_metadata, cache_directory, use_cache)
}

#' Join metacell metadata to an existing data frame
#'
#' Downloads and joins the metacell metadata with cellNexus metadata.
#' This function creates indexed tables for efficient joining and returns a data frame.
#'
#' @param tbl A `tbl_sql` object (from get_metadata) or a database connection
#' @param cache_directory A character string specifying the local cache
#'   directory where remote parquet files will be stored. Defaults to
#'   [get_default_cache_dir()].
#' @param join_keys A character vector of column names used for the join.
#'   Defaults to `c("sample_id", "dataset_id", "cell_id")`.
#' @return A lazy SQL table with metacell metadata joined to the cellNexus metadata.
#' @importFrom dplyr left_join
#' @keywords internal
#' @noRd
join_metacell_table <- function(tbl,
                                cache_directory = get_default_cache_dir(),
                                join_keys = c("sample_id", "dataset_id", "cell_id")) {
  # Temporary only because internal
  cloud_metadata <- SAMPLE_DATABASE_URL
  # Synchronize remote files
  walk(cloud_metadata, function(url) {
    # Calculate the file path from the URL
    path <- file.path(cache_directory, url |>
      basename())
    if (!file.exists(path)) {
      report_file_sizes(url)
      sync_remote_file(
        url,
        path,
        progress(type = "down", con = stderr())
      )
    }
  })
  parquet_path <- file.path(cache_directory, cloud_metadata |>
    basename())
  # Fetch current connection
  conn <- dbplyr::remote_con(tbl)
  # Register the metacell parquet as a lazy table
  metacell_tbl <- duckdb_read_parquet(conn, parquet_path)
  # Join to the incoming tbl_lazy
  tbl |>
    left_join(metacell_tbl, by = join_keys)
}

#' Join Census metadata to an existing data frame
#'
#' Downloads and joins the Census metadata with cellNexus metadata.
#' This function creates indexed tables for efficient joining and returns a data frame.
#'
#' @param tbl A `tbl_sql` object (from get_metadata) or a database connection.
#' @param cloud_metadata  HTTP URL/URLs pointing to the census parquet database name and location.
#' @param cache_directory A character string specifying the local cache
#'   directory where remote parquet files will be stored. Defaults to
#'   [get_default_cache_dir()].
#' @param join_keys A character vector of column names used for the join.
#'   Defaults to `c("sample_id", "dataset_id", "observation_joinid")`.
#' @examples
#' library(dplyr)
#' get_metadata(cloud_metadata = SAMPLE_DATABASE_URL) |> head(2) |>
#'   # You do not need to specify anything in cloud_metadata
#'   join_census_table(
#'     cloud_metadata = SAMPLE_DATABASE_URL,
#'     cache_directory = tempdir()
#'   )
#' @return A lazy SQL table with Census metadata joined to the cellNexus metadata.
#' @importFrom dplyr left_join
#' @export
join_census_table <- function(tbl,
                              cloud_metadata = get_metadata_url("census_cell_metadata.2.2.0.parquet"),
                              cache_directory = get_default_cache_dir(),
                              join_keys = c("sample_id", "dataset_id", "observation_joinid")) {
  # Synchronize remote files
  walk(cloud_metadata, function(url) {
    # Calculate the file path from the URL
    path <- file.path(cache_directory, url |>
      basename())
    if (!file.exists(path)) {
      report_file_sizes(url)
      sync_remote_file(
        url,
        path,
        progress(type = "down", con = stderr())
      )
    }
  })
  parquet_path <- file.path(cache_directory, cloud_metadata |>
    basename())
  # Fetch current connection
  conn <- dbplyr::remote_con(tbl)
  # Register the census parquet as a lazy table
  census_tbl <- duckdb_read_parquet(conn, parquet_path)
  # Join to the incoming tbl_lazy
  tbl |>
    left_join(census_tbl, by = join_keys)
}

#' Returns the atlas version changelog as a tibble
#'
#' Downloads the `atlas_versions.parquet` registry from the cellNexus metadata
#' store, caches it locally, and returns it as an in-memory tibble. Each row
#' describes one atlas data release and its relationship to a CellxGene Census
#' snapshot.
#'
#' The `atlas_id` column in this table corresponds directly to the `atlas_id`
#' column returned by [get_metadata()], so you can join them to find which
#' Census snapshot any cell in your query came from.
#'
#' @param cache Optional character scalar. A local directory used to
#'   cache the downloaded parquet file. Defaults to a temporary directory to 
#'   separate from the main cache directory.
#' @return A tibble with columns:
#'   \describe{
#'     \item{atlas_id}{Atlas version identifier, e.g. `"cellxgene_2024/0.1.0"`.
#'       Matches the `atlas_id` column in [get_metadata()].}
#'     \item{census_version}{The CellxGene Census snapshot this atlas was built
#'       from, e.g. `"01-07-2024"`.}
#'     \item{change_type}{One of `"initial"`, `"patch"`, `"minor"`, or
#'       `"major"`. See `ATLAS_VERSIONS.md` for the conventions standards.}
#'     \item{description}{Summary text of what changed in this release.}
#'     \item{modified_at}{Modification date as a character scalar (`"YYYY-MM-DD"`).
#'       By default use `Sys.Date()`}
#'   }
#' @seealso
#' \itemize{
#'   \item CellxGene Census data releases (LTS):
#'   \url{https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_docsite_data_release_info.html}
#' }
#' @export
#' @examples
#' get_atlas_versions()
#' @references Mangiola, S., M. Milton, N. Ranathunga, C. S. N. Li-Wai-Suen,
#'   A. Odainic, E. Yang, W. Hutchison et al. "A multi-organ map of the human
#'   immune system across age, sex and ethnicity." bioRxiv (2023): 2023-06.
#'   doi:10.1101/2023.06.08.542671.
#' @source [Mangiola et al.,2023](https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3)
get_atlas_versions <- function(cache = tempdir()) {
  registry_url <- get_metadata_url("atlas_versions.parquet")
  local_path <- file.path(cache, "atlas_versions.parquet")
  sync_remote_file(registry_url,
    local_path,
    overwrite = TRUE
  )
  arrow::read_parquet(local_path)
}
