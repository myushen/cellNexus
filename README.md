cellNexus
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

`cellNexus` is a query interface that allow the programmatic exploration
and retrieval of the harmonised, curated and reannotated CELLxGENE
single-cell human cell atlas. Data can be retrieved at cell, sample, or
dataset levels based on filtering criteria.

Harmonised data is stored in the ARDC Nectar Research Cloud, and most
`cellNexus` functions interact with Nectar via web requests, so a
network connection is required for most functionality.

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

Load additional packages

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

    #> # Source:   SQL [?? x 94]
    #> # Database: DuckDB v1.2.2 [shen.m@Darwin 23.3.0:R 4.5.0/:memory:]
    #>    cell_id     dataset_id observation_joinid sample_id cell_type cell_type_ontology_t…¹
    #>    <chr>       <chr>      <chr>              <chr>     <chr>     <chr>                 
    #>  1 10X383_6:C… 0ee5ae70-… *|Vt;4ATOI         e732d3bf… leukocyte CL:0000738            
    #>  2 10X383_6:G… 0ee5ae70-… #z5ft!kwbn         e732d3bf… leukocyte CL:0000738            
    #>  3 10X383_5:C… 0ee5ae70-… kd(R>I7G#-         e732d3bf… leukocyte CL:0000738            
    #>  4 10X383_5:G… 0ee5ae70-… mSylksfIn3         e732d3bf… leukocyte CL:0000738            
    #>  5 10X218_3:C… 0ee5ae70-… XEl?wkHAcg         979206e1… leukocyte CL:0000738            
    #>  6 10X383_5:C… 0ee5ae70-… T_4$MfQ!Ss         e732d3bf… leukocyte CL:0000738            
    #>  7 10X218_3:T… 0ee5ae70-… &yS%7ee<`?         979206e1… leukocyte CL:0000738            
    #>  8 10X218_3:C… 0ee5ae70-… 0*ArkI0ju;         979206e1… leukocyte CL:0000738            
    #>  9 10X218_3:A… 0ee5ae70-… r<pyq>qS=W         979206e1… leukocyte CL:0000738            
    #> 10 10X383_6:C… 0ee5ae70-… )tfjoT&}j{         e732d3bf… leukocyte CL:0000738            
    #> # ℹ more rows
    #> # ℹ abbreviated name: ¹​cell_type_ontology_term_id
    #> # ℹ 88 more variables: sample_ <chr>, assay <chr>, assay_ontology_term_id <chr>,
    #> #   cell_count <chr>, citation <chr>, collection_id <chr>, dataset_version_id <chr>,
    #> #   default_embedding <chr>, development_stage <chr>,
    #> #   development_stage_ontology_term_id <chr>, disease <chr>,
    #> #   disease_ontology_term_id <chr>, donor_id <chr>, experiment___ <chr>, …

Metadata is saved to `get_default_cache_dir()` unless a custom path is
provided via the cache_directory argument. The `metadata` variable can
then be re-used for all subsequent queries.

### Explore the tissue

``` r
metadata |>
    dplyr::distinct(tissue, cell_type_unified_ensemble) 
#> # Source:   SQL [?? x 2]
#> # Database: DuckDB v1.2.2 [shen.m@Darwin 23.3.0:R 4.5.0/:memory:]
#>    tissue                  cell_type_unified_ensemble
#>    <chr>                   <chr>                     
#>  1 lower lobe of left lung cd8 tem                   
#>  2 lower lobe of left lung t cd4                     
#>  3 trachea                 muscle                    
#>  4 trachea                 macrophage                
#>  5 trachea                 pdc                       
#>  6 trachea                 cd16 mono                 
#>  7 trachea                 tgd                       
#>  8 lymph node              cd8 tcm                   
#>  9 lymph node              treg                      
#> 10 lung                    macrophage                
#> # ℹ more rows
```

## Quality Control

cellNexus metadata applies standardised quality control to filter out
empty droplets, dead or damaged cells, doublets, and samples with low
gene counts.

``` r
metadata = metadata |> 
  dplyr::filter(empty_droplet == FALSE,
         alive == TRUE,
         scDblFinder.class != "doublet",
         feature_count >= 5000)
```

## Download single-cell RNA sequencing counts

### Query raw counts

``` r
single_cell_counts = 
    metadata |>
    dplyr::filter(
        self_reported_ethnicity == "African" &
        assay |> stringr::str_like("%10x%") &
        tissue == "lung parenchyma" &
        cell_type |> stringr::str_like("%CD4%")
    ) |>
    get_single_cell_experiment()
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Experiment.

single_cell_counts
#> class: SingleCellExperiment 
#> dim: 56239 1780 
#> metadata(0):
#> assays(1): counts
#> rownames(56239): ENSG00000121410 ENSG00000268895 ... ENSG00000135605
#>   ENSG00000109501
#> rowData names(0):
#> colnames(1780):
#>   LAP87_ATCGTAGAGACTAGAT-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_1
#>   LAP87_ATCTTCATCCTGTTAT-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_1 ...
#>   LAP87_GCCGTGACAGGAGGAG-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_76
#>   LAP87_TGCATGACAACCGACC-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_76
#> colData names(95): dataset_id observation_joinid ... dir_prefix
#>   original_cell_
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

### Query counts scaled per million

``` r
single_cell_counts = 
    metadata |>
    dplyr::filter(
        self_reported_ethnicity == "African" &
        assay |> stringr::str_like("%10x%") &
        tissue == "lung parenchyma" &
        cell_type |> stringr::str_like("%CD4%")
    ) |>
    get_single_cell_experiment(assays = "cpm")
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Experiment.

single_cell_counts
#> class: SingleCellExperiment 
#> dim: 56239 1780 
#> metadata(0):
#> assays(1): cpm
#> rownames(56239): ENSG00000121410 ENSG00000268895 ... ENSG00000135605
#>   ENSG00000109501
#> rowData names(0):
#> colnames(1780):
#>   LAP87_ATCGTAGAGACTAGAT-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_1
#>   LAP87_ATCTTCATCCTGTTAT-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_1 ...
#>   LAP87_GCCGTGACAGGAGGAG-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_76
#>   LAP87_TGCATGACAACCGACC-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_76
#> colData names(95): dataset_id observation_joinid ... dir_prefix
#>   original_cell_
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

### Query pseudobulk

``` r
pseudobulk_counts = 
   metadata |>
    dplyr::filter(
        self_reported_ethnicity == "African" &
        assay |> stringr::str_like("%10x%") &
        tissue == "lung parenchyma" &
        cell_type |> stringr::str_like("%CD4%")
    ) |>
    get_pseudobulk()
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Experiment.

pseudobulk_counts
#> class: SingleCellExperiment 
#> dim: 56239 98 
#> metadata(0):
#> assays(1): counts
#> rownames(56239): ENSG00000000003 ENSG00000000005 ... ENSG00000290292
#>   ENSG00000291237
#> rowData names(0):
#> colnames(98): bfe624d44f7e5868cc22e11ad0f13866___t cd4
#>   bfe624d44f7e5868cc22e11ad0f13866___cd4 tcm ...
#>   4f067f7e5f960bc72b0710684a521e84____SC86___cd4 th1/th17 em
#>   e4d7f8162faf68a85f61bdbd81dae627___other
#> colData names(56): dataset_id sample_id ... dir_prefix sample_identifier
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

### Query metacell

The metadata includes a series of metacell aggregation levels, beginning
with 2, 4, 8, and so on. For example, the value of metacell_2 represents
a grouping of cells that can be split into two distinct metacells.

``` r
metacell_counts = 
   metadata |>
    dplyr::filter(!is.na(metacell_2)) |>
    dplyr::filter(
        self_reported_ethnicity == "African" &
        assay |> stringr::str_like("%10x%") &
        tissue == "lung parenchyma" &
        cell_type |> stringr::str_like("%CD4%")
    ) |>
    get_metacell(cell_aggregation = "metacell_2")
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Experiment.

metacell_counts
#> class: SingleCellExperiment 
#> dim: 56239 914 
#> metadata(0):
#> assays(1): counts
#> rownames(56239): ENSG00000121410 ENSG00000268895 ... ENSG00000135605
#>   ENSG00000109501
#> rowData names(0):
#> colnames(914): 9c8fa5a8d2ae37179b579a0217670512___LAP87_1_duong___3
#>   9c8fa5a8d2ae37179b579a0217670512___LAP87_1_duong___1 ...
#>   9c8fa5a8d2ae37179b579a0217670512___LAP87_1_duong___2
#>   9c8fa5a8d2ae37179b579a0217670512___LAP87_1_duong___1
#> colData names(39): metacell_2 dataset_id ... dir_prefix metacell_identifier
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

### Extract only a subset of genes

This is helpful if just few genes are of interest (e.g ENSG00000134644
(PUM1)), as they can be compared across samples. cellNexus uses ENSEMBL
gene ID(s).

``` r
single_cell_counts = 
    metadata |>
    dplyr::filter(
        self_reported_ethnicity == "African" &
        assay |> stringr::str_like("%10x%") &
        tissue == "lung parenchyma" &
        cell_type |> stringr::str_like("%CD4%")
    ) |>
    get_single_cell_experiment(assays = "cpm", features = "ENSG00000134644")
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Experiment.

single_cell_counts
#> class: SingleCellExperiment 
#> dim: 1 1780 
#> metadata(0):
#> assays(1): cpm
#> rownames(1): ENSG00000134644
#> rowData names(0):
#> colnames(1780):
#>   LAP87_ATCGTAGAGACTAGAT-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_1
#>   LAP87_ATCTTCATCCTGTTAT-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_1 ...
#>   LAP87_GCCGTGACAGGAGGAG-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_76
#>   LAP87_TGCATGACAACCGACC-1_duong___9f222629-9e39-47d0-b83f-e08d610c7479_76
#> colData names(95): dataset_id observation_joinid ... dir_prefix
#>   original_cell_
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

### Extract the counts as a Seurat object

This convert the H5 SingleCellExperiment to Seurat so it might take long
time and occupy a lot of memory depending on how many cells you are
requesting.

``` r
seurat_counts = 
    metadata |>
    dplyr::filter(
        self_reported_ethnicity == "African" &
        assay |> stringr::str_like("%10x%") &
        tissue == "lung parenchyma" &
        cell_type |> stringr::str_like("%CD4%")
    ) |>
    get_seurat()
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Experiment.

seurat_counts
#> An object of class Seurat 
#> 56239 features across 1780 samples within 1 assay 
#> Active assay: originalexp (56239 features, 0 variable features)
#>  2 layers present: counts, data
```

By default, data is downloaded to `get_default_cache_dir()` output. If
memory is a concern, users can specify a custom cache directory to
metadata and counts functions:

## Load metadata from the custom cache directory

``` r
metadata <- get_metadata(cache_directory = "/MY/CUSTOM/PATH")
```

## Query raw counts from the custom cache directory

``` r
single_cell_counts = 
    metadata |>
    dplyr::filter(
        self_reported_ethnicity == "African" &
        assay |> stringr::str_like("%10x%") &
        tissue == "lung parenchyma" &
        cell_type |> stringr::str_like("%CD4%")
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
single_cell_counts |> saveRDS("single_cell_counts.rds")
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

single_cell_counts |> zellkonverter::writeH5AD("single_cell_counts.h5ad", 
                                               compression = "gzip",
                                               verbose = TRUE)
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
  dplyr::with_groups(disease, ~ .x |> dplyr::mutate(median_count = median(`FCN1`, rm.na=TRUE))) |> 
  
  # Plot
  ggplot(aes(forcats::fct_reorder(disease, median_count,.desc = TRUE), `FCN1`,color = dataset_id)) +
  geom_jitter(shape=".") +
    
  # Style
  guides(color="none") +
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
  dplyr::with_groups(tissue, ~ .x |> dplyr::mutate(median_count = median(`FCN1`, rm.na=TRUE))) |> 
  
  # Plot
  ggplot(aes(forcats::fct_reorder(tissue, median_count,.desc = TRUE), `FCN1`,color = dataset_id)) +
  geom_jitter(shape=".") +
    
  # Style
  guides(color="none") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) + 
  xlab("Tissue") + 
  ggtitle("FCN1 in CD14 monocytes by tissue. Colored by datasets") + 
  theme(legend.position = "none", axis.text.x = element_text(size = 6.5))
```

![](man/figures/FCN1_tissue_plot.png)<!-- -->

## Integrate Cloud Metadata with Local Metadata

CellNexus not only enables users to query our metadata but also allows
integration with your local metadata. Additionally, users can integrate
with your metadata stored in the cloud.

To enable this feature, users must include
`file_id_cellNexus_single_cell` (e.g sample_sce_obj.h5ad) and `atlas_id`
(e.g cellxgene/dd-mm-yy) columns in the metadata.

``` r
# Define local cache and metadata path
local_cache <- tempdir()
meta_path <- paste0(local_cache,"/my_metadata.parquet")

# Extract the date from example sce metadata 
date <- cellNexus::sample_sce_obj |> S4Vectors::metadata() |> 
  purrr::pluck("data") |> dplyr::pull(atlas_id) |> unique()

# Create a directory for storing counts
counts_path = file.path(local_cache, date, "counts")
dir.create(counts_path, recursive = TRUE, showWarnings = FALSE)

# Define the SCE file path for saving
sce_path = file.path(counts_path, "sample_sce_obj.h5ad")

# Save metadata to disk 
cellNexus::sample_sce_obj |> S4Vectors::metadata() |> purrr::pluck("data") |>
  arrow::write_parquet(meta_path)

# Save SCE to disk
cellNexus::sample_sce_obj |> zellkonverter::writeH5AD(file = sce_path,
                                                      compression = "gzip")
#> ℹ Using the 'counts' assay as the X matrix
```

``` r
# A file from cloud 
file_id_from_cloud <- "68aea852584d77f78e5c91204558316d___1.h5ad"

get_metadata(local_metadata = meta_path,
             cache_directory = local_cache) |>
  
  # For illustration purpose, only filter a selected cloud metadata and the saved metadata 
  dplyr::filter(file_id_cellNexus_single_cell %in% c("sample_sce_obj.h5ad", file_id_from_cloud)) |>
  get_single_cell_experiment(cache_directory = local_cache)
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ! cellNexus says: Not all genes completely overlap across the provided objects.Counts are generated by genes intersection.
#> ℹ Compiling Experiment.
#> class: SingleCellExperiment 
#> dim: 155 4 
#> metadata(1): data
#> assays(1): counts
#> rownames(155): ENSG00000187634 ENSG00000188976 ... ENSG00000187144
#>   ENSG00000157191
#> rowData names(0):
#> colnames(4):
#>   GCTGCGAGTGGGTCAA-1-1-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-1-0-0-0-0___218acb0f-9f2f-4f76-b90b-15a4b7c7f629_1
#>   ACATCAGGTCTGCCAG-1-1-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-1-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0___218acb0f-9f2f-4f76-b90b-15a4b7c7f629_1
#>   Liver_cDNA_CCTATTAGTTTGTGAC-1_2 Liver_cDNA_CCTATTAGTTTGTGTG-1_2
#> colData names(95): dataset_id observation_joinid ... dir_prefix
#>   original_cell_
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
#> R version 4.5.0 (2025-04-11)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Sonoma 14.3
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: Australia/Melbourne
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] ggplot2_3.5.2    cellNexus_1.6.13
#> 
#> loaded via a namespace (and not attached):
#>   [1] RcppAnnoy_0.0.22            splines_4.5.0              
#>   [3] later_1.4.2                 filelock_1.0.3             
#>   [5] tibble_3.3.0                polyclip_1.10-7            
#>   [7] basilisk.utils_1.20.0       preprocessCore_1.70.0      
#>   [9] fastDummies_1.7.5           lifecycle_1.0.4            
#>  [11] rprojroot_2.0.4             globals_0.18.0             
#>  [13] lattice_0.22-6              MASS_7.3-65                
#>  [15] backports_1.5.0             magrittr_2.0.3             
#>  [17] plotly_4.10.4               sass_0.4.10                
#>  [19] rmarkdown_2.29              jquerylib_0.1.4            
#>  [21] yaml_2.3.10                 httpuv_1.6.16              
#>  [23] Seurat_5.3.0                sctransform_0.4.2          
#>  [25] spam_2.11-1                 sp_2.2-0                   
#>  [27] spatstat.sparse_3.1-0       reticulate_1.42.0          
#>  [29] cowplot_1.1.3               pbapply_1.7-2              
#>  [31] DBI_1.2.3                   RColorBrewer_1.1-3         
#>  [33] abind_1.4-8                 Rtsne_0.17                 
#>  [35] GenomicRanges_1.60.0        purrr_1.0.4                
#>  [37] BiocGenerics_0.54.0         GenomeInfoDbData_1.2.14    
#>  [39] IRanges_2.42.0              S4Vectors_0.46.0           
#>  [41] ggrepel_0.9.6               irlba_2.3.5.1              
#>  [43] listenv_0.9.1               spatstat.utils_3.1-4       
#>  [45] goftest_1.2-3               RSpectra_0.16-2            
#>  [47] spatstat.random_3.4-1       fitdistrplus_1.2-2         
#>  [49] parallelly_1.45.0           codetools_0.2-20           
#>  [51] DelayedArray_0.34.1         tidyselect_1.2.1           
#>  [53] UCSC.utils_1.4.0            farver_2.1.2               
#>  [55] shinyWidgets_0.9.0          matrixStats_1.5.0          
#>  [57] stats4_4.5.0                spatstat.explore_3.4-3     
#>  [59] duckdb_1.2.2                jsonlite_2.0.0             
#>  [61] progressr_0.15.1            ggridges_0.5.6             
#>  [63] survival_3.8-3              tools_4.5.0                
#>  [65] ica_1.0-3                   Rcpp_1.0.14                
#>  [67] glue_1.8.0                  gridExtra_2.3              
#>  [69] SparseArray_1.8.0           xfun_0.52                  
#>  [71] MatrixGenerics_1.20.0       GenomeInfoDb_1.44.0        
#>  [73] HDF5Array_1.36.0            dplyr_1.1.4                
#>  [75] withr_3.0.2                 fastmap_1.2.0              
#>  [77] basilisk_1.20.0             rhdf5filters_1.20.0        
#>  [79] ttservice_0.5.3             digest_0.6.37              
#>  [81] R6_2.6.1                    mime_0.13                  
#>  [83] colorspace_2.1-1            scattermore_1.2            
#>  [85] tensor_1.5                  spatstat.data_3.1-6        
#>  [87] h5mread_1.0.1               utf8_1.2.6                 
#>  [89] tidyr_1.3.1                 generics_0.1.4             
#>  [91] data.table_1.17.4           httr_1.4.7                 
#>  [93] htmlwidgets_1.6.4           S4Arrays_1.8.1             
#>  [95] uwot_0.2.3                  pkgconfig_2.0.3            
#>  [97] gtable_0.3.6                rsconnect_1.4.2            
#>  [99] blob_1.2.4                  lmtest_0.9-40              
#> [101] SingleCellExperiment_1.30.1 XVector_0.48.0             
#> [103] htmltools_0.5.8.1           dotCall64_1.2              
#> [105] SeuratObject_5.1.0          scales_1.4.0               
#> [107] Biobase_2.68.0              png_0.1-8                  
#> [109] spatstat.univar_3.1-3       knitr_1.50                 
#> [111] rstudioapi_0.17.1           tzdb_0.5.0                 
#> [113] reshape2_1.4.4              checkmate_2.3.2            
#> [115] nlme_3.1-168                curl_6.3.0                 
#> [117] rhdf5_2.52.1                cachem_1.1.0               
#> [119] zoo_1.8-14                  stringr_1.5.1              
#> [121] KernSmooth_2.23-26          parallel_4.5.0             
#> [123] miniUI_0.1.2                arrow_20.0.0.2             
#> [125] zellkonverter_1.18.0        pillar_1.10.2              
#> [127] grid_4.5.0                  vctrs_0.6.5                
#> [129] RANN_2.6.2                  promises_1.3.3             
#> [131] dbplyr_2.5.0                xtable_1.8-4               
#> [133] cluster_2.1.8.1             evaluate_1.0.3             
#> [135] readr_2.1.5                 cli_3.6.5                  
#> [137] compiler_4.5.0              rlang_1.1.6                
#> [139] crayon_1.5.3                future.apply_1.20.0        
#> [141] tidybulk_1.21.1             plyr_1.8.9                 
#> [143] stringi_1.8.7               viridisLite_0.4.2          
#> [145] deldir_2.0-4                assertthat_0.2.1           
#> [147] lazyeval_0.2.2              spatstat.geom_3.4-1        
#> [149] Matrix_1.7-3                dir.expiry_1.16.0          
#> [151] RcppHNSW_0.6.0              hms_1.1.3                  
#> [153] patchwork_1.3.0             bit64_4.6.0-1              
#> [155] future_1.58.0               Rhdf5lib_1.30.0            
#> [157] shiny_1.10.0                SummarizedExperiment_1.38.1
#> [159] ROCR_1.0-11                 igraph_2.1.4               
#> [161] bslib_0.9.0                 bit_4.6.0
```
