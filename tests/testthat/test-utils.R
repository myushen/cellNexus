test_that("url_file_size() returns the correct sizes", {
    c(
        "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz",
        "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz"
    ) |>
        url_file_size() |>
        expect_equal(c(
            0.973,
            0.944
        ), tolerance = 0.001)
})

test_that("get_default_cache_dir() returns a character path", {
  d <- get_default_cache_dir()
  expect_type(d, "character")
  expect_length(d, 1L)
  expect_true(nzchar(d))
  expect_true(grepl("cellNexus", d, fixed = TRUE))
})

test_that("keep_updated_metadata() removes stale parquet files from the cache directory", {
  cache_dir <- get_default_cache_dir()
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  current_file <- file.path(cache_dir, "sample_metadata.2.0.0.parquet")
  stale_file   <- file.path(cache_dir, "metadata.old.parquet")
  file.create(current_file)
  file.create(stale_file)
  on.exit(unlink(c(current_file, stale_file)), add = TRUE)
  
  expect_invisible(cellNexus:::keep_updated_metadata("sample_metadata.2.0.0.parquet"))
  expect_true( file.exists(current_file))
  expect_false(file.exists(stale_file))
})
