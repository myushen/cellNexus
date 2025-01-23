library(dplyr)

test_that("get_default_cache_dir() returns the correct directory on Linux", {
    grepl("linux", version$platform, fixed = TRUE) |>
        skip_if_not()

    "~/.cache/R/cellNexus" |>
        normalizePath() |>
        expect_equal(
            get_default_cache_dir(),
        )
})

test_that("sync_assay_files() syncs appropriate files", {
    temp <- tempfile()
    test_file <- "4164d0eb972ad5e12719b6858c9559ea___1.h5ad"
    
    atlas_name <- "cellxgene/19-12-2024"

    sync_assay_files(
      atlas_name = atlas_name,
      cell_aggregation = "",
      cache_dir = temp,
      files = test_file,
      subdirs = "cpm"
    )
    
    temp_subdir <- file.path(temp, atlas_name, "cpm")
    
    test_file %in% list.files(temp_subdir) |>
        expect(failure_message = "The file was not downloaded")
})

test_that("get_SingleCellExperiment() syncs appropriate files", {
    temp <- tempfile()
    test_file <- "4164d0eb972ad5e12719b6858c9559ea___1.h5ad"

    meta <- get_metadata() |> head(2)

    # The remote dataset should have many genes
    sce <- get_SingleCellExperiment(meta, cache_directory = temp)
    sce |>
        row.names() |>
        length() |>
        expect_gt(1)
})

test_that(
    "The assays argument to get_SingleCellExperiment controls the number
  of returned assays",
    {
        # We need this for the assays() function
        library(SummarizedExperiment)

        meta <- get_metadata() |> head(2)
        
        atlas_name <- meta |> distinct(atlas_id) |> pull()
        
        assay_names <- c("counts", "cpm")
        
        temp <- tempfile()

        # If we request both assays, we get both assays
        get_SingleCellExperiment(meta, cache_directory = temp,
                                 assays = assay_names) |>
            assays() |>
            names() |>
            expect_setequal(assay_names)

        # If we request one assay, we get one assays
        get_SingleCellExperiment(meta, cache_directory = temp,
                                 assays = "counts") |>
            assays() |>
            names() |>
            expect_setequal("counts")
        get_SingleCellExperiment(meta, cache_directory = temp,
                                 assays = "cpm") |>
            assays() |>
            names() |>
            expect_setequal("cpm")
    }
)

test_that("The features argument to get_SingleCellExperiment subsets genes", {
    meta <- get_metadata() |> head(2)

    # The un-subset dataset should have many genes
    sce_full <- get_SingleCellExperiment(meta) |>
        row.names() |>
        length()
    expect_gt(sce_full, 1)

    # The subset dataset should only have one gene
    sce_subset <- get_SingleCellExperiment(meta, features = "ENSG00000254123") |>
        row.names() |>
        length()
    expect_equal(sce_subset, 1)

    expect_gt(sce_full, sce_subset)
})

test_that("get_seurat() returns the appropriate data in Seurat format", {
    meta <- get_metadata() |> head(2)

    sce <- get_SingleCellExperiment(meta, features = "ENSG00000254123")
    seurat <- get_seurat(meta, features = "ENSG00000254123")

    # The output should be a Seurat object
    expect_s4_class(seurat, "Seurat")
    # Both methods should have appropriately subset genes
    expect_equal(
        rownames(sce),
        rownames(seurat)
    )
})

test_that("get_SingleCellExperiment() assigns the right cell ID to each cell", {
    atlas_id = "cellxgene/19-12-2024"
    id = "a65bcc2d-4243-44c1-a262-ab7dcddfcf86"
    file_id <- "7ddd6775d704d6826539abaee8d22f65___1.h5ad"
    
    # Force the file to be cached
    get_metadata() |>
        filter(dataset_id == id) |>
        get_SingleCellExperiment()
    
    # Load the SCE from cache directly
    assay_1 = cellNexus:::get_default_cache_dir() |>
        file.path(atlas_id, "counts", file_id ) |>
        zellkonverter::readH5AD(reader = "R", use_hdf5 = TRUE) |>
        assay("counts") |>
        as.matrix()
    
    # Make a SCE that has the right column names, but reversed
    assay_2 = 
        assay_1 |>
        colnames() |>
        tibble::tibble(
            file_id_cellNexus_single_cell = file_id,
            dataset_id = id,
            cell_id = _,
            atlas_id = atlas_id
        ) |>
        arrange(-row_number()) |>
        get_SingleCellExperiment(assays = "counts") |>
        assay("counts") |>
        as.matrix()
        
    colnames(assay_2) = sub("_1", "", x=colnames(assay_2))

    expect_equal(
        assay_1,
        assay_2[, colnames(assay_1)]
    )
})

test_that("get_metadata() is cached", {
    table = get_metadata()
    table_2 = get_metadata()
    
    identical(table, table_2) |> expect_true()
})

test_that("database_url() expect character ", {
  get_metadata_url() |>
    expect_s3_class("character")
})

test_that("get_metadata() expect a unique cell_type `abnormal cell` is present", {
  n_cell <- get_metadata() |> filter(cell_type == 'abnormal cell') |> as_tibble() |> nrow()
  expect_true(n_cell > 0)
})

test_that("get_metadata() expect to combine local and cloud metadata", {
  data("sample_sce_obj")
  
  cache = tempfile()
  if (!dir.exists(cache)) dir.create(cache, recursive = TRUE)
  metadata(sample_sce_obj)$data |> arrow::write_parquet(glue("{cache}/my_metadata.parquet"))
  
  atlas_name <- metadata(sample_sce_obj)$data |> pull(atlas_id) |> unique()
  assay <- "counts"
  file_id_from_cloud <- "68aea852584d77f78e5c91204558316d___1.h5ad"
  file_path <- file.path(cache, atlas_name, assay)
  
  n_cell = get_metadata(local_metadata = glue("{cache}/my_metadata.parquet")) |>
    filter(file_id_cellNexus_single_cell %in% c("sample_sce_obj.h5ad",
                                                file_id_from_cloud)) |>
    dplyr::count() |> pull() |> as.integer()
  expect_gt(n_cell , 1)
})

test_that("get_metadata() expect to query local metadata only", {
  data("sample_sce_obj")
  
  cache = tempfile()
  if (!dir.exists(cache)) dir.create(cache, recursive = TRUE)
  metadata(sample_sce_obj)$data |> arrow::write_parquet(glue("{cache}/my_metadata.parquet"))
  
  n_cell = get_metadata(local_metadata = glue("{cache}/my_metadata.parquet"),
                        cloud_metadata = NULL) |>
    dplyr::count() |> pull() |> as.integer()
  expect_equal(n_cell, 2)
})

test_that("get_single_cell_experiment() expect to combine local and cloud counts", {
  data("sample_sce_obj")
  
  cache = tempfile()
  if (!dir.exists(cache)) dir.create(cache, recursive = TRUE)
  metadata(sample_sce_obj)$data |> arrow::write_parquet(glue("{cache}/my_metadata.parquet"))
  
  atlas_name <- metadata(sample_sce_obj)$data |> pull(atlas_id) |> unique()
  assay <- "counts"
  file_id_from_cloud <- "68aea852584d77f78e5c91204558316d___1.h5ad"
  file_path <- file.path(cache, atlas_name, assay)
  
  if (!dir.exists(file_path)) dir.create(file_path, recursive = TRUE)
  sample_sce_obj |> writeH5AD(glue("{cache}/{atlas_name}/{assay}/sample_sce_obj.h5ad"), 
                              compression = "gzip")
  
  sce = get_metadata(local_metadata = glue("{cache}/my_metadata.parquet")) |>
    filter(file_id_cellNexus_single_cell %in% c("sample_sce_obj.h5ad",
                                                file_id_from_cloud)) |>
    get_single_cell_experiment(cache = cache)
    
  expect_contains(colnames(sce), "Liver_cDNA_CCTATTAGTTTGTGAC-1_2" )
})

test_that("get_single_cell_experiment() expect to get local counts only", {
  data("sample_sce_obj")
  
  cache = tempfile()
  if (!dir.exists(cache)) dir.create(cache, recursive = TRUE)
  metadata(sample_sce_obj)$data |> arrow::write_parquet(glue("{cache}/my_metadata.parquet"))
  
  atlas_name <- metadata(sample_sce_obj)$data |> pull(atlas_id) |> unique()
  assay <- "counts"
  file_path <- file.path(cache, atlas_name, assay)
  
  if (!dir.exists(file_path)) dir.create(file_path, recursive = TRUE)
  sample_sce_obj |> writeH5AD(glue("{cache}/{atlas_name}/{assay}/sample_sce_obj.h5ad"), 
                              compression = "gzip")
  
  n_cell_in_sce = get_metadata(local_metadata = glue("{cache}/my_metadata.parquet"),
                     cloud_metadata = NULL) |> 
    get_single_cell_experiment(cache = cache) |>
    colnames() |> length()
  
  expect_equal(n_cell_in_sce, 2)
})

# unharmonised_data is not implemented yet
# test_that("get_unharmonised_dataset works with one ID", {
#     dataset_id = "838ea006-2369-4e2c-b426-b2a744a2b02b"
#     unharmonised_meta = get_unharmonised_dataset(dataset_id)
# 
#     expect_s3_class(unharmonised_meta, "tbl")
# })

# test_that("get_unharmonised_metadata() returns the appropriate data", {
#     harmonised <- get_metadata() |> dplyr::filter(tissue == "kidney blood vessel")
#     unharmonised <- get_unharmonised_metadata(harmonised)
#     
#     unharmonised |> is.data.frame() |> expect_true()
#     expect_setequal(colnames(unharmonised), c("file_id", "unharmonised"))
#     
#     # The number of cells in both harmonised and unharmonised should be the same
#     expect_equal(
#         dplyr::collect(harmonised) |> nrow(),
#         unharmonised$unharmonised |> purrr::map_int(function(df) dplyr::tally(df) |> dplyr::pull(n)) |> sum()
#     )
#     
#     # The number of datasets in both harmonised and unharmonised should be the same
#     expect_equal(
#         harmonised |> dplyr::group_by(file_id) |> dplyr::n_groups(),
#         nrow(unharmonised)
#     )
# })


# test_that("import_one_sce() loads metadata from a SingleCellExperiment object into a parquet file and generates pseudobulk", {
#   # Test both functionalities together because if import them independently,
#   # the sample data will be loaded into the cache, which causes the second import to fail the unique file check
#   data(sample_sce_obj)
#   temp <- tempfile()
#   dataset_id <- "GSE122999"
#   import_one_sce(sce_obj = sample_sce_obj,
#                          cache_dir = temp,
#                  pseudobulk = TRUE)
#   
#   dataset_id %in% (get_metadata(cache_directory = temp) |> 
#                     dplyr::distinct(dataset_id) |> 
#                     dplyr::pull()) |>
#     expect(failure_message = "The correct metadata was not created")
#   
#   sme <- get_metadata(cache_directory = temp) |> filter(file_id == "id1") |>
#     get_pseudobulk(cache_directory = file.path(temp, "pseudobulk"))
#   sme |>
#     row.names() |>
#     length() |>
#     expect_gt(1)
# })
# 
# test_that("get_pseudobulk() syncs appropriate files", {
#   temp <- tempfile()
#   id <- "0273924c-0387-4f44-98c5-2292dbaab11e"
#   meta <- get_metadata(cache_directory = temp) |> filter(file_id == id)
#   
#   # The remote dataset should have many genes
#   sme <- get_pseudobulk(meta, cache_directory = temp)
#   sme |>
#     row.names() |>
#     length() |>
#     expect_gt(1)
# })
# 
# test_that("get_pseudobulk() syncs appropriate fixed file", {
#   temp <- tempfile()
#   ids <- c(
#     "b50b15f1-bf19-4775-ab89-02512ec941a6",
#     "bffedc04-5ba1-46d4-885c-989a294bedd4",
#     "cc3ff54f-7587-49ea-b197-1515b6d98c4c",
#     "0af763e1-0e2f-4de6-9563-5abb0ad2b01e",
#     "51f114ae-232a-4550-a910-934e175db814",
#     "327927c7-c365-423c-9ebc-07acb09a0c1a",
#     "3ae36927-c188-4511-88cc-572ee1edf906",
#     "6ed2cdc2-dda8-4908-ad6c-cead9afee85e",
#     "56e0359f-ee8d-4ba5-a51d-159a183643e5",
#     "5c64f247-5b7c-4842-b290-65c722a65952"
#   )
#   meta <- get_metadata(cache_directory = temp) |> dplyr::filter(file_id %in% ids)
#   
#   # The remote dataset should have many genes
#   sme <- get_pseudobulk(meta, cache_directory = temp)
#   sme |>
#     row.names() |>
#     length() |>
#     expect_gt(1)
# })
# 
# test_that("get_single_cell_experiment() syncs prostate atlas", {
#   temp <- tempfile()
#   # A sample from prostate atlas
#   sample <- "GSM4089151"
#   meta <- get_metadata(cache_directory = temp) |> filter(sample_ == sample)
#   sce <- meta |> get_single_cell_experiment(cache_directory = temp)
#   sce |>
#     row.names() |>
#     length() |>
#     expect_gt(1)
# })
