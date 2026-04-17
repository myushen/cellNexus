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
normalised expression, and pseudobulk layers. This harmonised design
facilitates efficient cross-dataset analyses and downstream integration.

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

By default, `get_metadata()` loads harmonised annotations. Users can
retrieve original Census annotations by the function
`join_census_table()`.

``` r
metadata <- get_metadata() |>
  join_census_table()
metadata
```

    #> # Source:   SQL [?? x 73]
    #> # Database: DuckDB 1.4.3 [unknown@Linux 5.14.0-362.24.1.el9_3.x86_64:R 4.5.3/:memory:]
    #>    cell_id observation_joinid dataset_id               sample_id sample_ cell_count citation collection_id dataset_version_id
    #>      <dbl> <chr>              <chr>                    <chr>     <chr>        <int> <chr>    <chr>         <chr>             
    #>  1      15 TjgA2vJ1;{         842c6f5d-4a94-4eef-8510… 1119f482… 1119f4…     714331 Publica… 4195ab4c-20b… b8b5be07-061b-439…
    #>  2      19 lNmuO5xs~3         842c6f5d-4a94-4eef-8510… 1119f482… 1119f4…     714331 Publica… 4195ab4c-20b… b8b5be07-061b-439…
    #>  3      18 h22!#$}SJ*         842c6f5d-4a94-4eef-8510… 1119f482… 1119f4…     714331 Publica… 4195ab4c-20b… b8b5be07-061b-439…
    #>  4      17 QRMCN*8*|#         842c6f5d-4a94-4eef-8510… 1119f482… 1119f4…     714331 Publica… 4195ab4c-20b… b8b5be07-061b-439…
    #>  5      14 qxl7HJjL$L         842c6f5d-4a94-4eef-8510… 1119f482… 1119f4…     714331 Publica… 4195ab4c-20b… b8b5be07-061b-439…
    #>  6      20 6Eu5c&aEH;         842c6f5d-4a94-4eef-8510… 1119f482… 1119f4…     714331 Publica… 4195ab4c-20b… b8b5be07-061b-439…
    #>  7      16 j}0<Y>a#X~         842c6f5d-4a94-4eef-8510… 1119f482… 1119f4…     714331 Publica… 4195ab4c-20b… b8b5be07-061b-439…
    #>  8       5 N>_|{;6_6N         842c6f5d-4a94-4eef-8510… 1f755b9b… 1f755b…     714331 Publica… 4195ab4c-20b… b8b5be07-061b-439…
    #>  9       4 CVfj>iWZ&P         842c6f5d-4a94-4eef-8510… 1f755b9b… 1f755b…     714331 Publica… 4195ab4c-20b… b8b5be07-061b-439…
    #> 10       3 5jp7Uu{#@#         842c6f5d-4a94-4eef-8510… 1f755b9b… 1f755b…     714331 Publica… 4195ab4c-20b… b8b5be07-061b-439…
    #> # ℹ more rows
    #> # ℹ 64 more variables: default_embedding <chr>, experiment___ <chr>, explorer_url <chr>, feature_count <int>,
    #> #   filesize <dbl>, filetype <chr>, mean_genes_per_cell <dbl>, primary_cell_count <chr>, published_at <chr>,
    #> #   raw_data_location <chr>, revised_at <chr>, run_from_cell_id <chr>, sample_heuristic <chr>, schema_version <chr>,
    #> #   suspension_type <chr>, tissue_type <chr>, title <chr>, tombstone <lgl>, url <chr>, x_approximate_distribution <chr>,
    #> #   age_days <int>, tissue_groups <chr>, nFeature_expressed_in_sample <int>, nCount_RNA <dbl>, empty_droplet <lgl>,
    #> #   cell_type_unified_ensemble <chr>, is_immune <lgl>, subsets_Mito_percent <int>, subsets_Ribo_percent <int>, …

Metadata is saved to `get_default_cache_dir()` unless a custom path is
provided via the cache_directory argument. The `metadata` variable can
then be re-used for all subsequent queries.

### Explore the tissue

``` r
metadata |>
  dplyr::distinct(tissue, cell_type_unified_ensemble)
#> # Source:   SQL [?? x 2]
#> # Database: DuckDB 1.4.3 [unknown@Linux 5.14.0-362.24.1.el9_3.x86_64:R 4.5.3/:memory:]
#>    tissue           cell_type_unified_ensemble
#>    <chr>            <chr>                     
#>  1 cortex of kidney cd16 mono                 
#>  2 cortex of kidney cd14 mono                 
#>  3 bone marrow      nk                        
#>  4 blood            cytotoxic                 
#>  5 blood            cd4 tem                   
#>  6 blood            granulocyte               
#>  7 blood            cd14 mono                 
#>  8 breast           monocytic                 
#>  9 blood            Unknown                   
#> 10 forebrain        nk                        
#> # ℹ more rows
```

## Quality control

cellNexus metadata applies standardised quality control to filter out
empty droplets, dead or damaged cells, doublets, and samples with low
gene counts.

``` r
metadata <- metadata |>
  keep_quality_cells()

metadata <- metadata |>
  dplyr::filter(feature_count >= 5000)
```

## Download single-cell RNA sequencing counts

### Query raw counts

``` r
single_cell_counts <-
  metadata |>
  dplyr::filter(
    self_reported_ethnicity == "African American" &
      assay == "10x 3' v3" &
      tissue == "breast" &
      cell_type == "T cell"
  ) |>
  get_single_cell_experiment()
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Experiment.

single_cell_counts
#> class: SingleCellExperiment 
#> dim: 33145 3 
#> metadata(0):
#> assays(1): counts
#> rownames(33145): ENSG00000243485 ENSG00000237613 ... ENSG00000277475 ENSG00000268674
#> rowData names(0):
#> colnames(3): 1_1 2_1 3_1
#> colData names(73): observation_joinid dataset_id ... tissue_ontology_term_id original_cell_
```

### Query counts scaled per million

``` r
single_cell_cpm <-
  metadata |>
  dplyr::filter(
    self_reported_ethnicity == "African American" &
      assay == "10x 3' v3" &
      tissue == "breast" &
      cell_type == "T cell"
  ) |>
  get_single_cell_experiment(assays = "cpm")
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Experiment.

single_cell_cpm
#> class: SingleCellExperiment 
#> dim: 33145 3 
#> metadata(0):
#> assays(1): cpm
#> rownames(33145): ENSG00000243485 ENSG00000237613 ... ENSG00000277475 ENSG00000268674
#> rowData names(0):
#> colnames(3): 1_1 2_1 3_1
#> colData names(73): observation_joinid dataset_id ... tissue_ontology_term_id original_cell_
```

### Query SCT normalised counts

``` r
single_cell_sct <-
  metadata |>
  dplyr::filter(
    self_reported_ethnicity == "African American" &
      assay == "10x 3' v3" &
      tissue == "breast" &
      cell_type == "T cell"
  ) |>
  get_single_cell_experiment(assays = "sct")
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Experiment.

single_cell_sct
#> class: SingleCellExperiment 
#> dim: 33145 3 
#> metadata(0):
#> assays(1): sct
#> rownames(33145): ENSG00000243485 ENSG00000237613 ... ENSG00000277475 ENSG00000268674
#> rowData names(0):
#> colnames(3): 1_1 2_1 3_1
#> colData names(73): observation_joinid dataset_id ... tissue_ontology_term_id original_cell_
```

### Query pseudobulk

``` r
pseudobulk_counts <-
  metadata |>
  dplyr::filter(
      assay  == "10x 5' v1" &
      tissue == "lung" &
      cell_type == "classical monocyte"
  ) |>
  get_pseudobulk()
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Experiment.

pseudobulk_counts
#> class: SingleCellExperiment 
#> dim: 19979 1 
#> metadata(0):
#> assays(1): counts
#> rownames(19979): ENSG00000243485 ENSG00000238009 ... ENSG00000275063 ENSG00000271254
#> rowData names(0):
#> colnames(1): 2e8c9911c9bfbffc07288adef93d3cf2___cd14 mono
#> colData names(56): dataset_id sample_id ... tissue_ontology_term_id sample_identifier
```

## Download cell communication metadata

Cell communication metadata was generated based on post-QC cells per
sample using `CellChat v2` method. It uses our harmonised cell type
annotation (cell_type_unified_ensemble) to infer the communication. It
captures inferred communication at both the ligand–receptor pair level
and the signalling pathway level.

interaction_count: The number of inferred interactions between each pair
of cell groups.  
interaction_weight: The aggregated communication strength between each
pair of cell groups.

For definitions of additional annotations, please refer to the CellChat
v2 documentation: <https://github.com/jinworks/CellChat>.

``` r
get_cell_communication_strength()
```

    #> # Source:   SQL [?? x 16]
    #> # Database: DuckDB 1.4.3 [unknown@Linux 5.14.0-362.24.1.el9_3.x86_64:R 4.5.3/:memory:]
    #>   source    target ligand receptor   lr_prob lr_pval interaction_name    interaction_name_2  pathway_name annotation evidence
    #>   <chr>     <chr>  <chr>  <chr>        <dbl>   <dbl> <chr>               <chr>               <chr>        <chr>      <chr>   
    #> 1 b         b      TGFB1  TGFbR1_R2 0.000116    1    TGFB1_TGFBR1_TGFBR2 TGFB1 - (TGFBR1+TG… TGFb         Secreted … KEGG: h…
    #> 2 b memory  b      TGFB1  TGFbR1_R2 0.000865    1    TGFB1_TGFBR1_TGFBR2 TGFB1 - (TGFBR1+TG… TGFb         Secreted … KEGG: h…
    #> 3 b naive   b      TGFB1  TGFbR1_R2 0.000696    0.99 TGFB1_TGFBR1_TGFBR2 TGFB1 - (TGFBR1+TG… TGFb         Secreted … KEGG: h…
    #> 4 cd14 mono b      TGFB1  TGFbR1_R2 0.00240     0.81 TGFB1_TGFBR1_TGFBR2 TGFB1 - (TGFBR1+TG… TGFb         Secreted … KEGG: h…
    #> 5 cd4 naive b      TGFB1  TGFbR1_R2 0.000957    1    TGFB1_TGFBR1_TGFBR2 TGFB1 - (TGFBR1+TG… TGFb         Secreted … KEGG: h…
    #> 6 cd4 tem   b      TGFB1  TGFbR1_R2 0.00242     0.76 TGFB1_TGFBR1_TGFBR2 TGFB1 - (TGFBR1+TG… TGFb         Secreted … KEGG: h…
    #> # ℹ 5 more variables: pathway_prob <dbl>, pathway_pval <dbl>, sample_id <chr>, interaction_count <dbl>,
    #> #   interaction_weight <dbl>

### Extract only a subset of genes

This is helpful if just few genes are of interest (e.g ENSG00000134644
(PUM1)), as they can be compared across samples. cellNexus uses ENSEMBL
gene ID(s).

``` r
single_cell_cpm <-
  metadata |>
  dplyr::filter(
    self_reported_ethnicity == "African American" &
      assay == "10x 3' v3" &
      tissue == "breast" &
      cell_type == "T cell"
  ) |>
  get_single_cell_experiment(assays = "cpm", features = "ENSG00000134644")
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Experiment.

single_cell_cpm
#> class: SingleCellExperiment 
#> dim: 1 3 
#> metadata(0):
#> assays(1): cpm
#> rownames(1): ENSG00000134644
#> rowData names(0):
#> colnames(3): 1_1 2_1 3_1
#> colData names(73): observation_joinid dataset_id ... tissue_ontology_term_id original_cell_
```

### Extract the counts as a Seurat object

This convert the H5 SingleCellExperiment to Seurat so it might take long
time and occupy a lot of memory depending on how many cells you are
requesting.

``` r
seurat_counts <-
  metadata |>
  dplyr::filter(
    self_reported_ethnicity == "African American" &
      assay == "10x 3' v3" &
      tissue == "breast" &
      cell_type == "T cell"
  ) |>
  get_seurat()
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Experiment.

seurat_counts
#> An object of class Seurat 
#> 33145 features across 3 samples within 1 assay 
#> Active assay: originalexp (33145 features, 0 variable features)
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
single_cell_counts <-
  metadata |>
  dplyr::filter(
    self_reported_ethnicity == "African American" &
      assay == "10x 3' v3" &
      tissue == "breast" &
      cell_type == "T cell"
  ) |>
  get_single_cell_experiment(cache_directory = "/MY/CUSTOM/PATH")

single_cell_counts
```

Same strategy can be applied for functions `get_pseuodbulk()` and
`get_seurat()` by passing your custom directory character to
“cache_directory” parameter.

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
#> ℹ Compiling Experiment.
#> class: SingleCellExperiment 
#> dim: 13132 500 
#> metadata(1): data
#> assays(1): counts
#> rownames(13132): ENSG00000228463 ENSG00000228327 ... ENSG00000273748 ENSG00000278384
#> rowData names(0):
#> colnames(500): AAACATACAACCAC_1 AAACATTGAGCTAC_1 ... AGGAGTCTTGTCAG_1 AGGATAGACATTTC_1
#> colData names(6): sample_id dataset_id ... file_id_cellNexus_single_cell original_cell_
```

# Cell metadata

The complete metadata dictionary (dataset-, sample-, and cell-level
columns, plus cellNexus-specific harmonised fields) is available on the
documentation site: [cellNexus documentation](https://cellnexus.org/).

# RNA abundance

The `counts` assay includes RNA abundance in the positive real scale
(not transformed with non-linear functions, e.g. log sqrt). Originally
CELLxGENE include a mix of scales and transformations specified in the
`x_normalization` column.

The `cpm` assay includes counts per million.

The `sct` assay includes normalised counts by `sctranform`.

# Other representations

The `rank` assay is the representation of each cell’s gene expression
profile where genes are ranked by expression intensity.

The `pseudobulk` assay includes aggregated RNA abundance for sample and
cell type combination.

The detailed documentation for RNA abundance is available on the
documentation site: [cellNexus documentation](https://cellnexus.org/).

# Session Info

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Red Hat Enterprise Linux 9.3 (Plow)
#> 
#> Matrix products: default
#> BLAS:   /stornext/System/data/software/rhel/9/base/tools/R/4.5.3/lib64/R/lib/libRblas.so 
#> LAPACK: /stornext/System/data/software/rhel/9/base/tools/R/4.5.3/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.1
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
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] cellNexus_0.99.20           arrow_22.0.0                shiny_1.12.1                anndataR_1.0.0             
#>  [5] testthat_3.3.1              tidyr_1.3.1                 ttservice_0.5.3             SingleCellExperiment_1.32.0
#>  [9] SummarizedExperiment_1.40.0 Biobase_2.70.0              GenomicRanges_1.62.1        Seqinfo_1.0.0              
#> [13] ggplot2_4.0.1               BiocStyle_2.38.0            HDF5Array_1.38.0            h5mread_1.2.1              
#> [17] rhdf5_2.54.1                DelayedArray_0.36.0         SparseArray_1.10.6          S4Arrays_1.10.1            
#> [21] IRanges_2.44.0              abind_1.4-8                 S4Vectors_0.48.0            MatrixGenerics_1.22.0      
#> [25] matrixStats_1.5.0           BiocGenerics_0.56.0         generics_0.1.4              Matrix_1.7-4               
#> [29] BiocParallel_1.44.0         crew.cluster_0.4.0          targets_1.11.4              dplyr_1.1.4                
#> 
#> loaded via a namespace (and not attached):
#>   [1] igraph_2.2.1           graph_1.88.1           ica_1.0-3              plotly_4.11.0          Formula_1.2-5         
#>   [6] scater_1.38.0          devtools_2.4.6         BiocBaseUtils_1.12.0   tidyselect_1.2.1       bit_4.6.0             
#>  [11] doParallel_1.0.17      clue_0.3-66            lattice_0.22-9         rjson_0.2.21           blob_1.2.4            
#>  [16] stringr_1.6.0          rngtools_1.5.2         rclipboard_0.2.1       hunspell_3.0.6         parallel_4.5.3        
#>  [21] png_0.1-8              cli_3.6.5              registry_0.5-1         goftest_1.2-3          textshaping_1.0.4     
#>  [26] purrr_1.2.0            BiocNeighbors_2.4.0    ggnetwork_0.5.14       RUnit_0.4.33.1         stringdist_0.9.15     
#>  [31] uwot_0.2.4             curl_7.0.0             mime_0.13              evaluate_1.0.5         ComplexHeatmap_2.26.0 
#>  [36] stringi_1.8.7          backports_1.5.0        desc_1.4.3             XML_3.99-0.20          httpuv_1.6.16         
#>  [41] magrittr_2.0.4         rappdirs_0.3.3         splines_4.5.3          rcmdcheck_1.4.0        CellChat_2.2.0.9001   
#>  [46] nanonext_1.7.2         biocViews_1.78.0       sctransform_0.4.3      ggbeeswarm_0.7.3       xopen_1.0.1           
#>  [51] sessioninfo_1.2.3      DBI_1.2.3              jquerylib_0.1.4        withr_3.0.2            systemfonts_1.3.1     
#>  [56] rprojroot_2.1.1        lmtest_0.9-40          RBGL_1.86.0            brio_1.1.5             spelling_2.3.2        
#>  [61] BiocManager_1.30.27    duckdb_1.4.3           htmlwidgets_1.6.4      fs_1.6.6               ggrepel_0.9.6         
#>  [66] statnet.common_4.12.0  diffobj_0.3.6          reticulate_1.44.1      zoo_1.8-14             XVector_0.50.0        
#>  [71] knitr_1.50             network_1.19.0         UCSC.utils_1.6.0       secretbase_1.0.5       foreach_1.5.2         
#>  [76] fansi_1.0.7            patchwork_1.3.2        grid_4.5.3             data.table_1.17.8      R.oo_1.27.1           
#>  [81] RSpectra_0.16-2        irlba_2.3.5.1          fastDummies_1.7.5      commonmark_2.0.0       ellipsis_0.3.2        
#>  [86] base64url_1.4          lazyeval_0.2.2         yaml_2.3.12            survival_3.8-6         scattermore_1.2       
#>  [91] crayon_1.5.3           mirai_2.5.3            RcppAnnoy_0.0.22       RColorBrewer_1.1-3     progressr_0.18.0      
#>  [96] later_1.4.4            ggridges_0.5.7         codetools_0.2-20       GlobalOptions_0.1.3    Seurat_5.4.0          
#> [101] Rtsne_0.17             shape_1.4.6.1          filelock_1.0.3         pkgconfig_2.0.3        xml2_1.5.1            
#> [106] spatstat.univar_3.1-5  ggpubr_0.6.2           spatstat.sparse_3.1-0  viridisLite_0.4.2      gridBase_0.4-7        
#> [111] xtable_1.8-4           car_3.1-3              plyr_1.8.9             httr_1.4.7             tools_4.5.3           
#> [116] globals_0.18.0         SeuratObject_5.2.0     pkgbuild_1.4.8         beeswarm_0.4.0         broom_1.0.11          
#> [121] checkmate_2.3.3        nlme_3.1-168           dbplyr_2.5.1           assertthat_0.2.1       digest_0.6.39         
#> [126] bookdown_0.46          dir.expiry_1.18.0      farver_2.1.2           reshape2_1.4.5         viridis_0.6.5         
#> [131] glue_1.8.0             cachem_1.1.0           BiocFileCache_3.0.0    polyclip_1.10-7        ggalluvial_0.12.5     
#> [136] parallelly_1.46.0      pkgload_1.4.1          R.cache_0.17.0         RcppHNSW_0.6.0         ScaledMatrix_1.18.0   
#> [141] carData_3.0-5          pbapply_1.7-4          httr2_1.2.2            spam_2.11-1            waldo_0.6.2           
#> [146] utf8_1.2.6             basilisk_1.22.0        fontawesome_0.5.3      ggsignif_0.6.4         crew_1.3.0            
#> [151] gridExtra_2.3          R.utils_2.13.0         rhdf5filters_1.22.0    RCurl_1.98-1.17        memoise_2.0.1         
#> [156] rmarkdown_2.30         scales_1.4.0           R.methodsS3_1.8.2      future_1.68.0          svglite_2.2.2         
#> [161] RANN_2.6.2             spatstat.data_3.1-9    rstudioapi_0.17.1      cluster_2.1.8.2        zellkonverter_1.20.1  
#> [166] spatstat.utils_3.2-0   BiocCheck_1.46.3       fitdistrplus_1.2-4     cowplot_1.2.0          colorspace_2.1-2      
#> [171] FNN_1.1.4.1            rlang_1.1.6            GenomeInfoDb_1.46.2    styler_1.11.0          shinyWidgets_0.9.0    
#> [176] dotCall64_1.2          collapse_2.1.6         circlize_0.4.17        scuttle_1.20.0         xfun_0.55             
#> [181] coda_0.19-4.1          sna_2.8                remotes_2.5.0          iterators_1.0.14       tibble_3.3.0          
#> [186] Rhdf5lib_1.32.0        bitops_1.0-9           ps_1.9.1               promises_1.5.0         praise_1.0.0          
#> [191] RSQLite_2.4.5          compiler_4.5.3         prettyunits_1.2.0      beachmat_2.26.0        listenv_0.10.0        
#> [196] Rcpp_1.1.0             roxygen2_7.3.3         BiocSingular_1.26.1    tensor_1.5.1           usethis_3.2.1         
#> [201] MASS_7.3-65            spatstat.random_3.4-3  R6_2.6.1               fastmap_1.2.0          rstatix_0.7.3         
#> [206] vipor_0.4.7            ROCR_1.0-11            rsvd_1.0.5             gtable_0.3.6           KernSmooth_2.23-26    
#> [211] miniUI_0.1.2           deldir_2.0-4           htmltools_0.5.9        bit64_4.6.0-1          spatstat.explore_3.6-0
#> [216] lifecycle_1.0.4        S7_0.2.1               processx_3.8.6         callr_3.7.6            sass_0.4.10           
#> [221] vctrs_0.6.5            rsconnect_1.7.0        spatstat.geom_3.6-1    NMF_0.28               sp_2.2-0              
#> [226] future.apply_1.20.1    bslib_0.9.0            pillar_1.11.1          collections_0.3.9      otel_0.2.0            
#> [231] jsonlite_2.0.0         GetoptLong_1.1.0
```
