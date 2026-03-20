# Utility scripts that are used internally by the package at runtime

#' Gets the file size of a number of remote files
#' @param urls A character vector containing URLs
#' @return The file size of each of the files pointed to by the provided URL,
#' in gigabytes, as double vector
#' @importFrom curl multi_run new_pool curl_fetch_memory parse_headers_list
#' @importFrom purrr map_dbl
#' @keywords internal
#' @noRd
url_file_size <- function(urls) {
    if (length(urls) == 0) return(numeric(0))
    
    # Use curl for parallel HEAD requests
    env <- new.env(parent = emptyenv())
    env$sizes <- rep(NA_real_, length(urls))
    pool <- new_pool()
    
    for (i in seq_along(urls)) {
        curl::curl_fetch_multi(
            urls[i],
            done = function(res) {
                idx <- which(urls == res$url)
                if (length(idx) > 0) {
                    headers <- parse_headers_list(res$headers)
                    content_length <- headers[["content-length"]]
                    if (!is.null(content_length)) {
                        env$sizes[idx[1]] <- as.numeric(content_length) / 10^9
                    }
                }
            },
            fail = function(msg) { },
            pool = pool,
            handle = curl::new_handle(nobody = TRUE)
        )
    }
    
    multi_run(pool = pool)
    env$sizes[is.na(env$sizes)] <- 0
    env$sizes
}

#' Prints a message indicating the size of a download
#' @param urls A character vector containing URLs
#' @importFrom cli cli_alert_info
#' @keywords internal
#' @noRd
#' @return `NULL`, invisibly
report_file_sizes <- function(urls) {
  total_size <- url_file_size(urls) |>
    sum() |>
    round(digits = 2)

  "Downloading {length(urls)} file{?s}, totalling {total_size} GB" |>
    cli_alert_info()
  
  invisible(NULL)
}

#' Returns the default cache directory with a version number
#' @export
#' @return A length one character vector.
#' @importFrom tools R_user_dir
#' @importFrom utils packageName
#' @examples
#' get_metadata(cloud_metadata = SAMPLE_DATABASE_URL, cache_directory = get_default_cache_dir())
#' @references Mangiola, S., M. Milton, N. Ranathunga, C. S. N. Li-Wai-Suen, 
#'   A. Odainic, E. Yang, W. Hutchison et al. "A multi-organ map of the human 
#'   immune system across age, sex and ethnicity." bioRxiv (2023): 2023-06.
#'   doi:10.1101/2023.06.08.542671.
#' @source [Mangiola et al.,2023](https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3)
get_default_cache_dir <- function() {
  packageName() |>
    R_user_dir(
      "cache"
    ) |>
    normalizePath() |>
    suppressWarnings()
}

#' Clear the default cache directory
#' @return A length one character vector.
#' @keywords internal
#' @noRd
clear_cache <- function() {
  get_default_cache_dir() |> unlink(TRUE, TRUE)
}

#' Clear the outdated metadata in the default cache directory.
#' @param updated_data A character vector of new metadata file name
#' @return `NULL`, invisibly
#' @keywords internal
#' @noRd
keep_updated_metadata <- function(updated_data) {
  cache_directory <- get_default_cache_dir()
  files_in_cache <- list.files(cache_directory)
  pattern <- "\\.parquet$"
  parquet_files <- grep(pattern, files_in_cache, value = TRUE)
  files_to_delete <- setdiff(parquet_files, updated_data)
  unlink(file.path(cache_directory, files_to_delete))
} 

#' Synchronises a single remote file with a local path
#' @importFrom httr write_disk GET stop_for_status
#' @importFrom cli cli_abort cli_alert_info
#' @return `NULL`, invisibly
#' @keywords internal
#' @noRd
sync_remote_file <- function(full_url, output_file, ...) {
  if (!file.exists(output_file)) {
    output_dir <- dirname(output_file)
    dir.create(output_dir,
               recursive = TRUE,
               showWarnings = FALSE
    )
    cli_alert_info("Downloading {full_url} to {output_file}")
    
    tryCatch(
      GET(full_url, write_disk(output_file), ...) |> stop_for_status(),
      error = function(e) {
        # Clean up if we had an error
        file.remove(output_file)
        cli_abort("File {full_url} could not be downloaded. {e}")
      }
    )
  }
  invisible(NULL)
}

#' Synchronises multiple remote files with local paths in parallel
#' 
#' Uses curl::multi_download for concurrent async I/O downloads.
#' Falls back to sequential downloads if parallel downloads are disabled.
#' 
#' @param urls A character vector of URLs to download
#' @param output_files A character vector of local file paths (same length as urls)
#' @param progress Whether to show a progress bar (default TRUE)
#' @importFrom curl multi_download
#' @importFrom cli cli_alert_info cli_alert_warning cli_abort
#' @return The output_files vector, invisibly
#' @keywords internal
#' @noRd
sync_remote_files <- function(urls, output_files, progress = TRUE) {
    if (length(urls) == 0) return(invisible(character(0)))
    if (length(urls) != length(output_files)) {
        cli_abort("urls and output_files must have the same length")
    }
    
    # Filter to only files that don't exist
    to_download <- !file.exists(output_files)
    urls_to_download <- urls[to_download]
    files_to_download <- output_files[to_download]
    
    if (length(urls_to_download) == 0) {
        return(invisible(output_files))
    }
    
    # Create directories for all output files
    unique_dirs <- unique(dirname(files_to_download))
    for (dir in unique_dirs) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    use_parallel <- getOption("cellNexus.parallel_downloads", TRUE)
    
    if (use_parallel && length(urls_to_download) > 1) {
        cli_alert_info("Downloading {length(urls_to_download)} file{?s} in parallel...")
        
        # Use curl::multi_download for parallel async I/O
        # Note: multiplex=TRUE enables HTTP/2 multiplexing for concurrent streams
        results <- multi_download(
            urls = urls_to_download,
            destfiles = files_to_download,
            progress = progress,
            multiplex = TRUE
        )
        
        # Check for failures
        failed <- results$success == FALSE | results$status_code >= 400
        if (any(failed, na.rm = TRUE)) {
            failed_urls <- urls_to_download[failed]
            failed_files <- files_to_download[failed]
            # Clean up failed downloads
            for (f in failed_files) {
                if (file.exists(f)) file.remove(f)
            }
            cli_alert_warning("{sum(failed)} file{?s} failed to download")
            if (sum(failed) == length(urls_to_download)) {
                cli_abort("All downloads failed. Check your network connection.")
            }
        }
    } else {
        # Sequential fallback
        for (i in seq_along(urls_to_download)) {
            sync_remote_file(urls_to_download[i], files_to_download[i])
        }
    }
    
    invisible(output_files)
}

#' Returns a tibble from a parquet file path
#' Since dbplyr 2.4.0, raw file paths aren't handled very well
#' See: <https://github.com/duckdb/duckdb-r/issues/38>
#' Hence the need for this method
#' @param conn A DuckDB connection.
#' @param path Path(s) to parquet file(s).
#' @importFrom dplyr tbl
#' @importFrom dbplyr sql
#' @importFrom glue glue_sql
#' @return An SQL data frame
#' @keywords internal
#' @noRd
duckdb_read_parquet <- function(conn, path) {
  from_clause <- glue_sql("FROM read_parquet([{`path`*}], union_by_name=true)", .con = conn) |> sql()
  tbl(conn, from_clause)
}

#' Deletes specific counts and metadata from cache
#' @importFrom purrr map
#' @importFrom dplyr filter distinct pull collect
#' @return `NULL`, invisibly
#' @keywords internal
#' @noRd
delete_counts <- function(data, 
                          assay = c("original","cpm"), 
                          cache_directory = get_default_cache_dir()){
  data <- collect(data)
  ids <- data |> distinct(file_id_db) |> pull(file_id_db)
  counts_path <- file.path(cache_directory, assay, ids)
  # counts
  map(counts_path, ~ .x |> unlink(recursive = TRUE))
  
  # metadata
  filename <- get_metadata(cache_directory = cache_directory, use_cache = FALSE) |>
    filter(file_id_db %in% ids) |>
    distinct(meta_filename) |>
    pull(meta_filename)
  arrow::read_parquet(filename) |>
    filter(!file_id_db %in% ids) |>
    arrow::write_parquet(filename)
}

#' Write table to Parquet file using DuckDB
#' This function takes an SQL table, renders it into an SQL query, and writes the output directly
#' to a Parquet file without loading into disk.
#' @param .tbl_sql A data frame loaded by DuckDB.
#' @param path A string specifying the file path where the Parquet file will be saved.
#' @param con A DuckDB connection object that allows executing SQL queries on the DuckDB database.
#' @return An integer specifying the row number of the output data frame.
#' @importFrom DBI dbConnect dbExecute
#' @keywords internal
#' @noRd
duckdb_write_parquet <- function(.tbl_sql, 
                                 path, 
                                 con = dbConnect(duckdb::duckdb(),  dbdir = ":memory:")) {
  sql_tbl <- 
    .tbl_sql |>
    dbplyr::sql_render()
  
  sql_call <- glue::glue("COPY ({sql_tbl}) TO '{path}' (FORMAT 'parquet')")
  
  res <- dbExecute(con, sql_call)
  res
}

#' Check whether a column is all NA in a dataframe, drop the column and warn users
#' @return A data frame where all values are not NA
#' @keywords internal
#' @noRd
clean_and_report_NA_columns <- function(df) {
  na_column_names <- df |>
    dplyr::select(dplyr::where(~ all(is.na(.)))) |>
    names()
  if (length(na_column_names) > 0) {
    cli::cli_alert_info(
      "Dropping {.field {na_column_names}} because they contain only NA values."
    )
    df <- df |>
      dplyr::select(-dplyr::all_of(na_column_names))
  }
  df
}

#' Save a SingleCellExperiment as an AnnData file
#'
#' Reports and removes `colData` columns that are entirely `NA` before writing to `.h5ad`.
#' If the object contains only one column, the assay is duplicated to avoid
#' single-column export issues.
#'
#' @param sce A `SingleCellExperiment` object.
#' @param path Output file path for the `.h5ad` file.
#' @param ... Additional arguments passed to [zellkonverter::writeH5AD()].
#'
#' @return Called for its side effect of writing an `.h5ad` file.
#' @inheritDotParams zellkonverter::writeH5AD
#' @keywords internal
#' @noRd
save_sce_as_h5ad <- function(sce, path, ...) {
  # Remove columns in colData that are all NA and warn the user
  cleaned_coldata <- SummarizedExperiment::colData(sce) |>
    as.data.frame() |>
    clean_and_report_NA_columns()
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(cleaned_coldata)

  if (ncol(SummarizedExperiment::assay(sce)) == 1) {
    sce <- sce |> duplicate_single_column_assay()
  }

  # anndataR does not support writing DelayedArray yet. Issue: https://github.com/scverse/anndataR/pull/387
  sce |> zellkonverter::writeH5AD(path, compression = "gzip", ...)
}

#' Duplicate Single-Column Assay in SingleCellExperiment Object
#'
#' This function handles SingleCellExperiment (SCE) objects where a specified assay 
#' contains only one column. It duplicates the single-column assay to avoid potential 
#' errors during saving or downstream analysis that require at least two columns. 
#' The duplicated column is marked with a prefix `DUMMY___` to distinguish it. 
#' Corresponding entries in the column metadata (`colData`) are also duplicated.
#'
#' @param sce A `SingleCellExperiment` object.
#' @importFrom SummarizedExperiment assay assays colData
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom rlang set_names
#' @keywords internal
#' @noRd
duplicate_single_column_assay <- function(sce) {

  assay_name <- (sce |> assays() |> names())[[1L]]

  if (ncol(assay(sce)) == 1) {

    # Duplicate the assay to prevent saving errors due to single-column matrices
    my_assay <- cbind(assay(sce), assay(sce))
    # Rename the second column to distinguish it
    colnames(my_assay)[2] <- paste0("DUMMY", "___", colnames(my_assay)[2])

    cd <- colData(sce)
    cd <- cd |> rbind(cd)
    rownames(cd)[2] <- paste0("DUMMY", "___", rownames(cd)[2])

    sce <- SingleCellExperiment(assay = list(my_assay) |> set_names(assay_name), colData = cd)
    sce
  } 
  sce
}

#' Synchronize metadata assay files with remote repository
#' 
#' @description Downloads and caches assay files from a remote repository based on metadata specifications.
#'
#' @param data A data frame or tbl_sql containing metadata with required columns `cell_id`, `atlas_id`, and a grouping column
#' @param assays Character vector specifying which assays to sync
#' @param repository URL of the remote repository containing the assay files
#' @param cell_aggregation Character string specifying the cell aggregation strategy 
#' @param grouping_column Column name in data used to group files for synchronization
#' @param cache_directory Local directory path where files will be cached. Uses default cache if not specified
#'
#' @return `NULL`, invisibly. Progress messages are displayed for the downloads.
#' 
#' @importFrom dplyr pull transmute distinct
#' @importFrom checkmate assert check_subset check_true
#' @importFrom cli cli_alert_info
#' @importFrom purrr pmap
#' @importFrom httr parse_url
#' @importFrom rlang .data
#' @keywords internal
#' @noRd
sync_metadata_assay_files <- function(data,
                                      assays = "counts",
                                      repository = COUNTS_URL,
                                      cell_aggregation = "",
                                      grouping_column,
                                      cache_directory = get_default_cache_dir()
) {
  assert(check_subset(c("cell_id", "atlas_id", grouping_column), colnames(data)))
  atlas_name <- data |> pull(atlas_id) |> unique()
  
  subdirs <- assay_map[assays]
  
  if (!is.null(repository)) {
    cli_alert_info("Synchronising files")
    parsed_repo <- parse_url(repository)
    assert(check_true(parsed_repo$scheme %in% c("http", "https")))
    
    files_to_read <-
      data |>
      transmute(
        files = .data[[grouping_column]],
        atlas_name = atlas_id,
        cache_dir = cache_directory
      ) |>
      distinct() |>
      # Convert to tibble here to avoid breaking the memory after loading all metadata
      tibble::as_tibble() |>
      pmap(function(files, atlas_name, cache_dir) {
        sync_assay_files(
          files = files,
          atlas_name = atlas_name,
          cache_dir = cache_dir,
          url = parsed_repo,
          cell_aggregation = cell_aggregation,
          subdirs = subdirs
        )
      })
  }
}

#' Keep high-quality cells based on QC columns
#'
#' @param data A data frame or tibble containing single-cell metadata.
#' @param empty_droplet_col A string specifying the column name 
#'   that indicates empty droplets (default: `"empty_droplet"`). 
#'   Expected logical vector
#' @param alive_col A string specifying the column name 
#'   that indicates whether cells are alive (default: `"alive"`). 
#'   Expected logical vector
#' @param doublet_col A string specifying the column name 
#'   that indicates doublets (default: `"scDblFinder.class"`). 
#'   Expected character vector: `"doublet"` and/or `"singlet"` and/or `"unknown"`.
#'
#' @return A filtered data frame containing only cells that pass all QC checks.
#' @examples
#' get_metadata(cloud_metadata = SAMPLE_DATABASE_URL, cache_directory = tempdir()) |> 
#'   head(2) |>
#'   keep_quality_cells()
#' @export
#' @importFrom rlang .data
#' @importFrom dplyr filter
#' @source [Mangiola et al.,2023](https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3)
keep_quality_cells <- function(data,
                               empty_droplet_col = "empty_droplet",
                               alive_col = "alive",
                               doublet_col = "scDblFinder.class") {
  data |>
    filter(
      .data[[empty_droplet_col]] == FALSE,
      .data[[alive_col]] == TRUE,
      .data[[doublet_col]] != "doublet"
    )
}

