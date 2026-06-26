# Functions that relate to unharmonised metadata

#' Base URL for all the unharmonised data
#' @keywords internal
#' @noRd
UNHARMONISED_URL <- "https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/unharmonised_metadata"

#' Returns unharmonised metadata for selected datasets.
#'
#' Various metadata fields are *not* common between datasets, so it does not
#' make sense for these to live in the main metadata table. This function is a
#' utility that allows easy fetching of this data if necessary.
#'
#' @param dataset_id A character vector, where each entry is a dataset ID
#'   obtained from the `$file_id_cellNexus_single_cell` column of the table returned from
#'   [get_metadata()]
#' @param cells An optional character vector of cell IDs. If provided, only
#'   metadata for those cells will be returned.
#' @param conn An optional DuckDB connection object. If provided, it will re-use
#'   the existing connection instead of opening a new one.
#' @param remote_url Optional character vector of length 1. An HTTP URL pointing
#'   to the root URL under which all the unharmonised dataset files are located.
#' @param cache_directory Optional character vector of length 1. A file path on
#'   your local system to a directory (not a file) that will be used to store
#'   the unharmonised metadata files.
#' @importFrom purrr map set_names
#' @importFrom glue glue
#' @importFrom DBI dbConnect
#' @importFrom duckdb duckdb
#' @importFrom dplyr tbl filter
#' @importFrom rlang .data
#' @keywords internal
#' @return A named list, where each name is a dataset file ID, and each value is
#'   a "lazy data frame", ie a `tbl`.
#' @examples
#' \dontrun{
#' dataset <- "838ea006-2369-4e2c-b426-b2a744a2b02b"
#' harmonised_meta <- get_metadata() |>
#'   dplyr::filter(file_id_cellNexus_single_cell == dataset) |>
#'   dplyr::collect()
#' unharmonised_meta <- get_unharmonised_dataset(dataset)
#' unharmonised_tbl <- dplyr::collect(unharmonised_meta[[dataset]])
#' dplyr::left_join(harmonised_meta, unharmonised_tbl,
#'   by = c("file_id_cellNexus_single_cell", "cell_id")
#' )
#' }
#' @references Shen, M., Y. Gao, N. Liu, D. Bhuva, M. Milton, J. Henao,
#'   J. Andrews, E. Yang, C. Zhan, N. Liu, S. Si, J. W. Hutchison,
#'   M. H. Shakeel, M. Morgan, A. T. Papenfuss, J. Iskander, J. M. Polo,
#'   and S. Mangiola. "cellNexus: Quality control, annotation, aggregation
#'   and analytical layers for the Human Cell Atlas data." bioRxiv (2026).
#'   doi:10.64898/2026.04.14.718336.
#' @source [Shen et al.,2026](https://www.biorxiv.org/content/10.64898/2026.04.14.718336v3)
get_unharmonised_dataset <- function(
  dataset_id,
  cells = NULL,
  conn = duckdb() |>
    dbConnect(drv = _, read_only = TRUE),
  remote_url = UNHARMONISED_URL,
  cache_directory = get_default_cache_dir()
) {
  unharmonised_root <- file.path(
    cache_directory,
    "unharmonised"
  )
  file_name <- glue::glue("{dataset_id}.parquet")
  local_path <- file.path(unharmonised_root, file_name)
  glue("{remote_url}/{file_name}") |>
    sync_remote_file(
      local_path,
      progress(type = "down", con = stderr())
    )

  duckdb_read_parquet(conn, local_path) |>
    filter(.data$cell_id %in% cells)
}
