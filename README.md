cellNexus
================
Mangiola et al.

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#maturing)
<!-- badges: end -->

# Introduction

`cellNexus` builds upon and extends the functionality of the previously
released `CuratedAtlasQueryR`, providing a unified query and access
interface to the harmonised, curated, and reannotated CELLxGENE human
cell atlas. It enables reproducible and programmatic exploration of
large-scale single-cell data resources, supporting retrieval at the
cell, sample, and dataset levels through flexible filtering by tissue,
cell type, experimental condition, or other metadata features. The
retrieved data are returned for downstream analysis.

`cellNexus` integrates over 40 million human cells processed with
standardised quality control, consistent normalisation, and unified
abundance representations—including single-cell, counts-per-million,
pseudobulk, and metacell layers. This harmonised design facilitates
efficient cross-dataset analyses and downstream integration.

Data are hosted on the ARDC Nectar Research Cloud, and most `cellNexus`
functions interact with Nectar via web requests, so a network connection
is required for most functionality.

`cellNexus` builds on top of package CuratedAtlasQueryR. While both rely
on pre-computed expression layers, they differ in how these layers are
generated. cellNexus implements a more standardised workflow, including
explicit removal of empty droplets and dead cells, followed by
harmonised quality control, normalisation, and multi-layer data
generation. Through this process, it produces updated datasets that
remain aligned with the evolving CELLxGENE releases.

<img src="man/figures/logo.png" width="120x" height="139px" />

<img src="man/figures/svcf_logo.jpeg" width="155x" height="58px" /><img src="man/figures/czi_logo.png" width="129px" height="58px" /><img src="man/figures/bioconductor_logo.jpg" width="202px" height="58px" /><img src="man/figures/vca_logo.png" width="219px" height="58px" /><img src="man/figures/nectar_logo.png" width="180px" height="58px" />

# Query interface

## Installation

``` r
devtools::install_github("MangiolaLaboratory/cellNexus")
```

## Load the package

``` r
library(cellNexus)
```

## Load additional packages

``` r
suppressPackageStartupMessages({
  library(ggplot2)
})
```

## Load and explore the metadata

### Load the metadata

``` r
metadata <- get_metadata()
metadata
```

    #> # Source:   SQL [?? x 89]
    #> # Database: DuckDB 1.4.3 [unknown@Linux 5.14.0-362.24.1.el9_3.x86_64:R 4.5.2/:memory:]
    #>    cell_id dataset_id       observation_joinid sample_id cell_type cell_type_ontology_t…¹ sample_ assay assay_ontology_term_id cell_count
    #>      <dbl> <chr>            <chr>              <chr>     <chr>     <chr>                  <chr>   <chr> <chr>                       <int>
    #>  1      81 cda2c8cd-be1c-4… *NUPW@J{c2         034f0fb1… monocyte  CL:0000576             034f0f… 10x … EFO:0011025                255901
    #>  2      82 cda2c8cd-be1c-4… KIV>qGFIS?         034f0fb1… monocyte  CL:0000576             034f0f… 10x … EFO:0011025                255901
    #>  3      83 cda2c8cd-be1c-4… p5e=WoIq0d         034f0fb1… monocyte  CL:0000576             034f0f… 10x … EFO:0011025                255901
    #>  4      84 cda2c8cd-be1c-4… I6>u{Gb-J_         034f0fb1… monocyte  CL:0000576             034f0f… 10x … EFO:0011025                255901
    #>  5      85 cda2c8cd-be1c-4… lx`7Bo-&7n         034f0fb1… monocyte  CL:0000576             034f0f… 10x … EFO:0011025                255901
    #>  6      86 cda2c8cd-be1c-4… 6mRCZW}rOM         034f0fb1… monocyte  CL:0000576             034f0f… 10x … EFO:0011025                255901
    #>  7      87 cda2c8cd-be1c-4… -NL-OH3!IA         034f0fb1… monocyte  CL:0000576             034f0f… 10x … EFO:0011025                255901
    #>  8      88 cda2c8cd-be1c-4… zHCZWNmUHu         034f0fb1… monocyte  CL:0000576             034f0f… 10x … EFO:0011025                255901
    #>  9      89 cda2c8cd-be1c-4… *_#lQ<oUnT         034f0fb1… monocyte  CL:0000576             034f0f… 10x … EFO:0011025                255901
    #> 10      99 cda2c8cd-be1c-4… IdHwp1GBZm         03ddfd57… monocyte  CL:0000576             03ddfd… 10x … EFO:0009899                255901
    #> # ℹ more rows
    #> # ℹ abbreviated name: ¹​cell_type_ontology_term_id
    #> # ℹ 79 more variables: citation <chr>, collection_id <chr>, dataset_version_id <chr>, default_embedding <chr>, development_stage <chr>,
    #> #   development_stage_ontology_term_id <chr>, disease <chr>, disease_ontology_term_id <chr>, donor_id <chr>, experiment___ <chr>,
    #> #   explorer_url <chr>, feature_count <int>, filesize <dbl>, filetype <chr>, is_primary_data <chr>, mean_genes_per_cell <dbl>,
    #> #   organism <chr>, organism_ontology_term_id <chr>, primary_cell_count <chr>, published_at <chr>, raw_data_location <chr>,
    #> #   revised_at <chr>, run_from_cell_id <chr>, sample_heuristic <chr>, schema_version <chr>, self_reported_ethnicity <chr>, …

Metadata is saved to `get_default_cache_dir()` unless a custom path is
provided via the cache_directory argument. The `metadata` variable can
then be re-used for all subsequent queries.

### Explore the tissue

``` r
metadata |>
  dplyr::distinct(tissue, cell_type_unified_ensemble)
#> # Source:   SQL [?? x 2]
#> # Database: DuckDB 1.4.3 [unknown@Linux 5.14.0-362.24.1.el9_3.x86_64:R 4.5.2/:memory:]
#>    tissue              cell_type_unified_ensemble
#>    <chr>               <chr>                     
#>  1 thymus              cd14 mono                 
#>  2 breast              nk                        
#>  3 renal pelvis        epithelial                
#>  4 kidney              epithelial                
#>  5 kidney blood vessel epithelial                
#>  6 lung parenchyma     cd4 th2 em                
#>  7 respiratory airway  cd4 th2 em                
#>  8 lung                cd4 th1 em                
#>  9 lung                cd4 fh em                 
#> 10 nose                cd4 fh em                 
#> # ℹ more rows
```

## Quality control

cellNexus metadata applies standardised quality control to filter out
empty droplets, dead or damaged cells, doublets, and samples with low
gene counts.

``` r
metadata <- metadata |>
  dplyr::filter(feature_count >= 5000) |>
  keep_quality_cells()
```

## Download single-cell RNA sequencing counts

### Query raw counts

``` r
single_cell_counts <-
  metadata |>
  dplyr::filter(
    self_reported_ethnicity == "African" &
      assay |>
        stringr::str_like("%10x%") &
      tissue == "lung parenchyma" &
      cell_type |>
        stringr::str_like("%CD4%")
  ) |>
  head() |>
  get_single_cell_experiment()

single_cell_counts
```

    #> class: SingleCellExperiment 
    #> dim: 56239 6 
    #> metadata(0):
    #> assays(1): counts
    #> rownames(56239): ENSG00000121410 ENSG00000268895 ... ENSG00000135605 ENSG00000109501
    #> rowData names(0):
    #> colnames(6): LAP92_CATTCTAGTGCGGATA-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_1
    #>   LAP92_CTCATGCCACCTGATA-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_1 ...
    #>   GCTCCTAAGGGTATCG_F02607___9f222629-9e39-47d0-b83f-e08d610c7479_1
    #>   AACACGTCACGCATCG_F01853___9f222629-9e39-47d0-b83f-e08d610c7479_2
    #> colData names(98): dataset_id observation_joinid ... dir_prefix original_cell_
    #> reducedDimNames(0):
    #> mainExpName: NULL
    #> altExpNames(0):

### Query counts scaled per million

``` r
single_cell_cpm <-
  metadata |>
  dplyr::filter(
    self_reported_ethnicity == "African" &
      assay |>
        stringr::str_like("%10x%") &
      tissue == "lung parenchyma" &
      cell_type |>
        stringr::str_like("%CD4%")
  ) |>
  head() |>
  get_single_cell_experiment(assays = "cpm")

single_cell_cpm
```

    #> class: SingleCellExperiment 
    #> dim: 1 6 
    #> metadata(0):
    #> assays(1): cpm
    #> rownames(1): ENSG00000134644
    #> rowData names(0):
    #> colnames(6): LAP92_CATTCTAGTGCGGATA-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_1
    #>   LAP92_CTCATGCCACCTGATA-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_1 ...
    #>   GCTCCTAAGGGTATCG_F02607___9f222629-9e39-47d0-b83f-e08d610c7479_1
    #>   AACACGTCACGCATCG_F01853___9f222629-9e39-47d0-b83f-e08d610c7479_2
    #> colData names(98): dataset_id observation_joinid ... dir_prefix original_cell_
    #> reducedDimNames(0):
    #> mainExpName: NULL
    #> altExpNames(0):

### Query pseudobulk

``` r
pseudobulk_counts <-
  metadata |>
  dplyr::filter(
    self_reported_ethnicity == "African" &
      assay |>
        stringr::str_like("%10x%") &
      tissue == "lung parenchyma" &
      cell_type |>
        stringr::str_like("%CD4%")
  ) |>
  head() |>
  get_pseudobulk()

pseudobulk_counts
```

    #> class: SingleCellExperiment 
    #> dim: 56239 3 
    #> metadata(0):
    #> assays(1): counts
    #> rownames(56239): ENSG00000000003 ENSG00000000005 ... ENSG00000290292 ENSG00000291237
    #> rowData names(0):
    #> colnames(3): a2459ad4272363e6eb775e8e99607c3e___cd4 th1 em 9c8fa5a8d2ae37179b579a0217670512___LAP92_1_duong___cd4 th2 em
    #>   e4d7f8162faf68a85f61bdbd81dae627___cd4 th2 em
    #> colData names(59): dataset_id sample_id ... dir_prefix sample_identifier
    #> reducedDimNames(0):
    #> mainExpName: NULL
    #> altExpNames(0):

### Extract only a subset of genes

This is helpful if just few genes are of interest (e.g ENSG00000134644
(PUM1)), as they can be compared across samples. cellNexus uses ENSEMBL
gene ID(s).

``` r
single_cell_cpm <-
  metadata |>
  dplyr::filter(
    self_reported_ethnicity == "African" &
      assay |>
        stringr::str_like("%10x%") &
      tissue == "lung parenchyma" &
      cell_type |>
        stringr::str_like("%CD4%")
  ) |>
  head() |>
  get_single_cell_experiment(assays = "cpm", features = "ENSG00000134644")

single_cell_counts
```

    #> class: SingleCellExperiment 
    #> dim: 1 6 
    #> metadata(0):
    #> assays(1): cpm
    #> rownames(1): ENSG00000134644
    #> rowData names(0):
    #> colnames(6): LAP92_CATTCTAGTGCGGATA-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_1
    #>   LAP92_CTCATGCCACCTGATA-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_1 ...
    #>   GCTCCTAAGGGTATCG_F02607___9f222629-9e39-47d0-b83f-e08d610c7479_1
    #>   AACACGTCACGCATCG_F01853___9f222629-9e39-47d0-b83f-e08d610c7479_2
    #> colData names(98): dataset_id observation_joinid ... dir_prefix original_cell_
    #> reducedDimNames(0):
    #> mainExpName: NULL
    #> altExpNames(0):

### Extract the counts as a Seurat object

This convert the H5 SingleCellExperiment to Seurat so it might take long
time and occupy a lot of memory depending on how many cells you are
requesting.

``` r
seurat_counts <-
  metadata |>
  dplyr::filter(
    self_reported_ethnicity == "African" &
      assay |>
        stringr::str_like("%10x%") &
      tissue == "lung parenchyma" &
      cell_type |>
        stringr::str_like("%CD4%")
  ) |>
  head() |>
  get_seurat()

seurat_counts
```

    #> An object of class Seurat 
    #> 56239 features across 6 samples within 1 assay 
    #> Active assay: originalexp (56239 features, 0 variable features)
    #>  2 layers present: counts, data

By default, data is downloaded to `get_default_cache_dir()` output. If
memory is a concern, users can specify a custom cache directory to
metadata and counts functions:

## Load metadata from the custom cache directory

``` r
metadata <- get_metadata(cache_directory = "/MY/CUSTOM/PATH")
```

## Query raw counts from the custom cache directory

``` r
single_cell_counts <-
  metadata |>
  dplyr::filter(
    self_reported_ethnicity == "African" &
      assay |>
        stringr::str_like("%10x%") &
      tissue == "lung parenchyma" &
      cell_type |>
        stringr::str_like("%CD4%")
  ) |>
  get_single_cell_experiment(cache_directory = "/MY/CUSTOM/PATH")

single_cell_counts
```

Same strategy can be applied for functions `get_pseuodbulk()`,
`get_metacell()`, `get_seurat()` by passing your custom directory
character to “cache_directory” parameter.

## Save your `SingleCellExperiment`

The returned `SingleCellExperiment` can be saved with three modalities,
as `.rds` or as `HDF5` or as `H5AD`.

### Saving as RDS (fast saving, slow reading)

Saving as `.rds` has the advantage of being fast, and the `.rds` file
occupies very little disk space as it only stores the links to the files
in your cache.

However it has the disadvantage that for big `SingleCellExperiment`
objects, which merge a lot of HDF5 from your
`get_single_cell_experiment`, the display and manipulation is going to
be slow. In addition, an `.rds` saved in this way is not portable: you
will not be able to share it with other users.

``` r
single_cell_counts |>
  saveRDS("single_cell_counts.rds")
```

### Saving as HDF5 (slow saving, fast reading)

Saving as `.hdf5` executes any computation on the `SingleCellExperiment`
and writes it to disk as a monolithic `HDF5`. Once this is done,
operations on the `SingleCellExperiment` will be comparatively very
fast. The resulting `.hdf5` file will also be totally portable and
sharable.

However this `.hdf5` has the disadvantage of being larger than the
corresponding `.rds` as it includes a copy of the count information, and
the saving process is going to be slow for large objects.

``` r
# ! IMPORTANT if you save 200K+ cells
HDF5Array::setAutoBlockSize(size = 1e+09)

single_cell_counts |>
  HDF5Array::saveHDF5SummarizedExperiment(
    "single_cell_counts",
    replace = TRUE,
    as.sparse = TRUE,
    verbose = TRUE
  )
```

### Saving as H5AD (slow saving, fast reading)

Saving as `.h5ad` executes any computation on the `SingleCellExperiment`
and writes it to disk as a monolithic `H5AD`. The `H5AD` format is the
HDF5 disk representation of the AnnData object and is well-supported in
Python.

However this `.h5ad` saving strategy has a bottleneck of handling
columns with only NA values of a `SingleCellExperiment` metadata.

``` r
# ! IMPORTANT if you save 200K+ cells
HDF5Array::setAutoBlockSize(size = 1e+09)

single_cell_counts |>
  anndataR::write_h5ad("single_cell_counts.h5ad",
    compression = "gzip",
    verbose = TRUE
  )
```

## Visualise gene transcription

We can gather all CD14 monocytes cells and plot the distribution of
ENSG00000085265 (FCN1) across all tissues

``` r
# Plots with styling
counts <- metadata |>
  # Filter and subset
  dplyr::filter(cell_type_unified_ensemble == "cd14 mono") |>
  # Get counts per million for FCN1 gene
  get_single_cell_experiment(assays = "cpm", features = "ENSG00000085265") |>
  suppressMessages() |>
  # Add feature to table
  tidySingleCellExperiment::join_features("ENSG00000085265", shape = "wide") |>
  # Rank x axis
  tibble::as_tibble() |>
  # Rename to gene symbol
  dplyr::rename(FCN1 = ENSG00000085265)

# Plot by disease
counts |>
  dplyr::with_groups(disease, ~ .x |>
    dplyr::mutate(median_count = median(`FCN1`, rm.na = TRUE))) |>
  # Plot
  ggplot(aes(forcats::fct_reorder(disease, median_count, .desc = TRUE), `FCN1`, color = dataset_id)) +
  geom_jitter(shape = ".") +

  # Style
  guides(color = "none") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
  xlab("Disease") +
  ggtitle("FCN1 in CD14 monocytes by disease. Coloured by datasets")
```

![](man/figures/FCN1_disease_plot.png)<!-- -->

``` r
# Plot by tissue
counts |>
  dplyr::with_groups(tissue, ~ .x |>
    dplyr::mutate(median_count = median(`FCN1`, rm.na = TRUE))) |>
  # Plot
  ggplot(aes(
    forcats::fct_reorder(tissue,
      median_count,
      .desc = TRUE
    ),
    `FCN1`,
    color = dataset_id
  )) +
  geom_jitter(shape = ".") +

  # Style
  guides(color = "none") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
  xlab("Tissue") +
  ggtitle("FCN1 in CD14 monocytes by tissue. Colored by datasets") +
  theme(legend.position = "none", axis.text.x = element_text(size = 6.5))
```

![](man/figures/FCN1_tissue_plot.png)<!-- -->

## Integrate cloud and local metadata

`cellNexus` not only enables users to query our metadata but also allows
integration with your local metadata. Additionally, users can integrate
with your metadata stored in the cloud.

To enable this feature, users must include
`file_id_cellNexus_single_cell` and `atlas_id` (e.g cellxgene/dd-mm-yy)
columns in the metadata. See metadata structure in cellNexus::pbmc3k_sce

``` r
# Set up local cache and paths
local_cache <- tempdir()
layer <- "counts"
meta_path <- file.path(local_cache, "pbmc3k_metadata.parquet")
data(pbmc3k_sce)

# Extract and prepare metadata
pbmc3k_metadata <- pbmc3k_sce |>
  S4Vectors::metadata() |>
  purrr::pluck("data") |>
  dplyr::mutate(
    counts_directory = file.path(tempdir(), atlas_id, layer),
    sce_path = file.path(counts_directory, file_id_cellNexus_single_cell)
  )

# Get unique paths
counts_directory <- pbmc3k_metadata |>
  dplyr::pull(counts_directory) |>
  unique()

sce_path <- pbmc3k_metadata |>
  dplyr::pull(sce_path) |>
  unique()

# Create directory structure
dir.create(counts_directory, recursive = TRUE, showWarnings = FALSE)

# Save data to disk
pbmc3k_sce |>
  S4Vectors::metadata() |>
  purrr::pluck("data") |>
  arrow::write_parquet(meta_path)

# Save SCE object
pbmc3k_sce |>
  anndataR::write_h5ad(sce_path, compression = "gzip", mode = "w")
```

``` r
# A cellNexus file
file_id_from_cloud <- "e52795dec7b626b6276b867d55328d9f___1.h5ad"
file_id_local <- basename(sce_path)

get_metadata(
  cloud_metadata = METADATA_URL,
  local_metadata = meta_path,
  cache_directory = local_cache
) |>
  # For illustration purpose, only filter a selected cloud metadata and the saved metadata
  dplyr::filter(file_id_cellNexus_single_cell %in% c(file_id_from_cloud, file_id_local)) |>
  dplyr::select(cell_id, sample_id, dataset_id, cell_type_unified_ensemble, atlas_id, file_id_cellNexus_single_cell) |>
  get_single_cell_experiment(cache_directory = local_cache)
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ! The number of cells in the SingleCellExperiment will be less than the number of cells you have selected from the metadata. Are cell IDs duplicated? Or, do cell IDs correspond to the counts file?
#> ! cellNexus says: Not all genes completely overlap across the provided objects. Counts are generated by genes intersection.
#> ℹ Compiling Experiment.
#> class: SingleCellExperiment 
#> dim: 12795 500 
#> metadata(1): data
#> assays(1): counts
#> rownames(12795): ENSG00000228463 ENSG00000228327 ... ENSG00000273748 ENSG00000278384
#> rowData names(0):
#> colnames(500): AAACATACAACCAC_1 AAACATTGAGCTAC_1 ... AGGAGTCTTGTCAG_1 AGGATAGACATTTC_1
#> colData names(6): sample_id dataset_id ... file_id_cellNexus_single_cell original_cell_
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

# Cell metadata

Dataset-specific columns (definitions available at
cellxgene.cziscience.com)

`cell_count`, `collection_id`, `filetype`, `is_primary_data`,
`mean_genes_per_cell`, `published_at`, `revised_at`, `schema_version`,
`tombstone`, `x_normalization`, `explorer_url`, `dataset_id`,
`dataset_version_id`

Sample-specific columns (definitions available at
cellxgene.cziscience.com)

`sample_id`, `sample_`, `age_days`, `assay`, `assay_ontology_term_id`,
`development_stage`, `development_stage_ontology_term_id`,
`self_reported_ethnicity`, `self_reported_ethnicity_ontology_term_id`,
`experiment___`, `organism`, `organism_ontology_term_id`,
`sample_placeholder`, `sex`, `sex_ontology_term_id`, `tissue`,
`tissue_type`, `tissue_ontology_term_id`, `tissue_groups`, `disease`,
`disease_ontology_term_id`, `is_primary_data`, `donor_id`, `is_immune`

Cell-specific columns (definitions available at
cellxgene.cziscience.com)

`cell_id`, `cell_type`, `cell_type_ontology_term_id`,
`cell_annotation_azimuth_l2`, `cell_annotation_blueprint_singler`,
`observation_joinid`, `empty_droplet`, `alive`, `scDblFinder.class`

Through harmonisation and curation we introduced custom column, not
present in the original CELLxGENE metadata

- `age_days`: donors’ age in days
- `cell_type_unified_ensemble`: the consensus call identity (for immune
  cells) using the original and three novel annotations using Seurat
  Azimuth and SingleR
- `cell_annotation_azimuth_l2`: Azimuth cell annotation
- `cell_annotation_blueprint_singler`: SingleR cell annotation using
  Blueprint reference
- `cell_annotation_blueprint_monaco`: SingleR cell annotation using
  Monaco reference
- `sample_heuristic`: sample subdivision for internal use
- `file_id_cellNexus_single_cell`: file subdivision for internal use
- `file_id_cellNexus_pseudobulk`: file subdivision for internal use
- `sample_id`: sample ID
- `nCount_RNA`: total number of RNA detected in a cell per sample
- `nFeature_expressed_in_sample`: total number of genes expressed in a
  cell per sample

# RNA abundance

The `counts` assay includes RNA abundance in the positive real scale
(not transformed with non-linear functions, e.g. log sqrt). Originally
CELLxGENE include a mix of scales and transformations specified in the
`x_normalization` column.

The `cpm` assay includes counts per million.

# Other representations

The `rank` assay is the representation of each cell’s gene expression
profile where genes are ranked by expression intensity.

The `pseudobulk` assay includes aggregated RNA abundance for sample and
cell type combination.

The metacell (e.g `metacell_2`, `metacell_4` etc) assays represent
hierarchical partitions of cells into metacell groups.

# Session Info

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Red Hat Enterprise Linux 9.3 (Plow)
#> 
#> Matrix products: default
#> BLAS:   /stornext/System/data/software/rhel/9/base/tools/R/4.5.2/lib64/R/lib/libRblas.so 
#> LAPACK: /stornext/System/data/software/rhel/9/base/tools/R/4.5.2/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.1
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Australia/Melbourne
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] ggplot2_4.0.1     BiocStyle_2.38.0  anndataR_1.0.0    shiny_1.12.1      testthat_3.3.1    cellNexus_0.99.14 lintr_3.3.0-1    
#> [8] dplyr_1.1.4      
#> 
#> loaded via a namespace (and not attached):
#>   [1] fs_1.6.6                    matrixStats_1.5.0           spatstat.sparse_3.1-0       xopen_1.0.1                
#>   [5] fontawesome_0.5.3           devtools_2.4.6              httr_1.4.7                  RColorBrewer_1.1-3         
#>   [9] tools_4.5.2                 sctransform_0.4.3           backports_1.5.0             utf8_1.2.6                 
#>  [13] R6_2.6.1                    HDF5Array_1.38.0            lazyeval_0.2.2              uwot_0.2.4                 
#>  [17] rhdf5filters_1.22.0         withr_3.0.2                 sp_2.2-0                    prettyunits_1.2.0          
#>  [21] gridExtra_2.3               progressr_0.18.0            cli_3.6.5                   Biobase_2.70.0             
#>  [25] spatstat.explore_3.6-0      fastDummies_1.7.5           sass_0.4.10                 Seurat_5.4.0               
#>  [29] arrow_22.0.0                S7_0.2.1                    spatstat.data_3.1-9         ggridges_0.5.7             
#>  [33] pbapply_1.7-4               commonmark_2.0.0            R.utils_2.13.0              parallelly_1.46.0          
#>  [37] sessioninfo_1.2.3           styler_1.11.0               rstudioapi_0.17.1           generics_0.1.4             
#>  [41] ica_1.0-3                   spatstat.random_3.4-3       Matrix_1.7-4                waldo_0.6.2                
#>  [45] S4Vectors_0.48.0            rclipboard_0.2.1            abind_1.4-8                 R.methodsS3_1.8.2          
#>  [49] lifecycle_1.0.4             yaml_2.3.12                 SummarizedExperiment_1.40.0 rhdf5_2.54.1               
#>  [53] SparseArray_1.10.6          Rtsne_0.17                  grid_4.5.2                  blob_1.2.4                 
#>  [57] promises_1.5.0              crayon_1.5.3                dir.expiry_1.18.0           miniUI_0.1.2               
#>  [61] lattice_0.22-7              beachmat_2.26.0             cowplot_1.2.0               pillar_1.11.1              
#>  [65] knitr_1.50                  GenomicRanges_1.62.1        future.apply_1.20.1         codetools_0.2-20           
#>  [69] glue_1.8.0                  spatstat.univar_3.1-5       rex_1.2.1                   data.table_1.17.8          
#>  [73] remotes_2.5.0               vctrs_0.6.5                 png_0.1-8                   spam_2.11-1                
#>  [77] gtable_0.3.6                rcmdcheck_1.4.0             assertthat_0.2.1            cachem_1.1.0               
#>  [81] xfun_0.55                   S4Arrays_1.10.1             mime_0.13                   Seqinfo_1.0.0              
#>  [85] rsconnect_1.7.0             survival_3.8-3              SingleCellExperiment_1.32.0 ellipsis_0.3.2             
#>  [89] fitdistrplus_1.2-4          ROCR_1.0-11                 nlme_3.1-168                usethis_3.2.1              
#>  [93] bit64_4.6.0-1               filelock_1.0.3              RcppAnnoy_0.0.22            GenomeInfoDb_1.46.2        
#>  [97] rprojroot_2.1.1             R.cache_0.17.0              bslib_0.9.0                 irlba_2.3.5.1              
#> [101] KernSmooth_2.23-26          otel_0.2.0                  BiocGenerics_0.56.0         DBI_1.2.3                  
#> [105] zellkonverter_1.20.1        duckdb_1.4.3                tidyselect_1.2.1            processx_3.8.6             
#> [109] bit_4.6.0                   compiler_4.5.2              curl_7.0.0                  h5mread_1.2.1              
#> [113] xml2_1.5.1                  desc_1.4.3                  DelayedArray_0.36.0         plotly_4.11.0              
#> [117] bookdown_0.46               checkmate_2.3.3             scales_1.4.0                lmtest_0.9-40              
#> [121] callr_3.7.6                 stringr_1.6.0               digest_0.6.39               goftest_1.2-3              
#> [125] spatstat.utils_3.2-0        rmarkdown_2.30              basilisk_1.22.0             XVector_0.50.0             
#> [129] htmltools_0.5.9             pkgconfig_2.0.3             MatrixGenerics_1.22.0       dbplyr_2.5.1               
#> [133] fastmap_1.2.0               rlang_1.1.6                 htmlwidgets_1.6.4           UCSC.utils_1.6.0           
#> [137] farver_2.1.2                jquerylib_0.1.4             zoo_1.8-14                  jsonlite_2.0.0             
#> [141] BiocParallel_1.44.0         R.oo_1.27.1                 magrittr_2.0.4              scuttle_1.20.0             
#> [145] dotCall64_1.2               patchwork_1.3.2             Rhdf5lib_1.32.0             Rcpp_1.1.0                 
#> [149] reticulate_1.44.1           stringi_1.8.7               brio_1.1.5                  MASS_7.3-65                
#> [153] plyr_1.8.9                  pkgbuild_1.4.8              parallel_4.5.2              listenv_0.10.0             
#> [157] ggrepel_0.9.6               deldir_2.0-4                splines_4.5.2               tensor_1.5.1               
#> [161] ps_1.9.1                    igraph_2.2.1                spatstat.geom_3.6-1         RcppHNSW_0.6.0             
#> [165] reshape2_1.4.5              stats4_4.5.2                pkgload_1.4.1               evaluate_1.0.5             
#> [169] SeuratObject_5.2.0          BiocManager_1.30.27         httpuv_1.6.16               RANN_2.6.2                 
#> [173] tidyr_1.3.1                 purrr_1.2.0                 polyclip_1.10-7             future_1.68.0              
#> [177] scattermore_1.2             xtable_1.8-4                RSpectra_0.16-2             roxygen2_7.3.3             
#> [181] later_1.4.4                 viridisLite_0.4.2           tibble_3.3.0                memoise_2.0.1              
#> [185] IRanges_2.44.0              cluster_2.1.8.1             shinyWidgets_0.9.0          globals_0.18.0             
#> [189] xmlparsedata_1.0.5
```
