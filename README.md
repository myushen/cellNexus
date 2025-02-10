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

    #> # Source:   SQL [?? x 61]
    #> # Database: DuckDB v1.1.0 [unknown@Linux 5.14.0-362.24.1.el9_3.x86_64:R 4.4.1/:memory:]
    #>    cell_id        dataset_id observation_joinid sample_id cell_type cell_type_ontology_t…¹ sample_ assay assay_ontology_term_id
    #>    <chr>          <chr>      <chr>              <chr>     <chr>     <chr>                  <chr>   <chr> <chr>                 
    #>  1 10X229_4:GGAA… 182f6a56-… bPG(iFL7v%         1ff18682… neuron    CL:0000540             1ff186… 10x … EFO:0009922           
    #>  2 10X229_3:ATGT… 182f6a56-… B`gPCpeFM&         1ff18682… neuron    CL:0000540             1ff186… 10x … EFO:0009922           
    #>  3 10X229_4:CATG… 182f6a56-… ?nZ@+pg6>J         1ff18682… neuron    CL:0000540             1ff186… 10x … EFO:0009922           
    #>  4 10X229_3:CATC… 182f6a56-… P^ds`8fcb2         1ff18682… neuron    CL:0000540             1ff186… 10x … EFO:0009922           
    #>  5 10X229_3:TCAT… 182f6a56-… T$lzLt*&U&         1ff18682… neuron    CL:0000540             1ff186… 10x … EFO:0009922           
    #>  6 10X229_3:ATTA… 182f6a56-… O3nX^8L>sD         1ff18682… neuron    CL:0000540             1ff186… 10x … EFO:0009922           
    #>  7 10X229_3:ACAC… 182f6a56-… l)-Q8?=OrJ         1ff18682… neuron    CL:0000540             1ff186… 10x … EFO:0009922           
    #>  8 10X229_3:TGTA… 182f6a56-… W)c5vgJa#v         1ff18682… neuron    CL:0000540             1ff186… 10x … EFO:0009922           
    #>  9 10X229_3:CAGA… 182f6a56-… vq&v@>l`se         1ff18682… neuron    CL:0000540             1ff186… 10x … EFO:0009922           
    #> 10 10X229_4:TCCT… 182f6a56-… mpJX1@{T%B         1ff18682… neuron    CL:0000540             1ff186… 10x … EFO:0009922           
    #> # ℹ more rows
    #> # ℹ abbreviated name: ¹​cell_type_ontology_term_id
    #> # ℹ 52 more variables: cell_count <chr>, citation <chr>, collection_id <chr>, dataset_version_id <chr>,
    #> #   development_stage <chr>, development_stage_ontology_term_id <chr>, disease <chr>, disease_ontology_term_id <chr>,
    #> #   donor_id <chr>, experiment___ <chr>, explorer_url <chr>, feature_count <int>, filesize <dbl>, filetype <chr>,
    #> #   is_primary_data <chr>, mean_genes_per_cell <dbl>, organism <chr>, organism_ontology_term_id <chr>,
    #> #   primary_cell_count <chr>, published_at <date>, raw_data_location <chr>, revised_at <date>, sample_heuristic <chr>, …

The `metadata` variable can then be re-used for all subsequent queries.

### Explore the tissue

``` r
metadata |>
    dplyr::distinct(tissue, cell_type_unified_ensemble) 
#> # Source:   SQL [?? x 2]
#> # Database: DuckDB v1.1.0 [unknown@Linux 5.14.0-362.24.1.el9_3.x86_64:R 4.4.1/:memory:]
#>    tissue              cell_type_unified_ensemble
#>    <chr>               <chr>                     
#>  1 blood               nk                        
#>  2 lung                cd4 tcm                   
#>  3 lung                endothelial               
#>  4 lung                cd4 th17 em               
#>  5 lung                cd4 fh em                 
#>  6 lung                immune                    
#>  7 parietal lobe       glial                     
#>  8 blood               cd4 th2 em                
#>  9 thoracic lymph node tgd                       
#> 10 thoracic lymph node macrophage                
#> # ℹ more rows
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
#> # A SingleCellExperiment-tibble abstraction: 1,800 × 63
#> # [90mFeatures=56239 | Cells=1800 | Assays=counts[0m
#>    .cell          dataset_id observation_joinid sample_id cell_type cell_type_ontology_t…¹ sample_ assay assay_ontology_term_id
#>    <chr>          <chr>      <chr>              <chr>     <chr>     <chr>                  <chr>   <chr> <chr>                 
#>  1 LAP87_GCTGGGT… 9f222629-… YK~|ZXL}LN         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x … EFO:0009922           
#>  2 LAP87_GACCCAG… 9f222629-… em#eUayTqZ         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x … EFO:0009922           
#>  3 LAP87_GCTGGGT… 9f222629-… 8(o)3D`Q%_         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x … EFO:0009922           
#>  4 TCATTTGCAAGCT… 9f222629-… {_UNJ6zNR<         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#>  5 CATCAAGGTCGAA… 9f222629-… z6RMY?)7(D         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#>  6 AGACGTTAGTAGA… 9f222629-… xF*TjD66E+         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#>  7 CCTTTCTTCCCAA… 9f222629-… w@YkRixyk@         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#>  8 GAAGCAGGTCGAC… 9f222629-… !Jd)V-~~>v         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#>  9 AAGCCGCTCTCGA… 9f222629-… zZEWSDl)l<         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#> 10 TACTCGCCAAAGG… 9f222629-… b`)P9|FF1D         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#> # ℹ 1,790 more rows
#> # ℹ abbreviated name: ¹​cell_type_ontology_term_id
#> # ℹ 54 more variables: cell_count <chr>, citation <chr>, collection_id <chr>, dataset_version_id <chr>,
#> #   development_stage <chr>, development_stage_ontology_term_id <chr>, disease <chr>, disease_ontology_term_id <chr>,
#> #   donor_id <chr>, experiment___ <chr>, explorer_url <chr>, feature_count <int>, filesize <dbl>, filetype <chr>,
#> #   is_primary_data <chr>, mean_genes_per_cell <dbl>, organism <chr>, organism_ontology_term_id <chr>,
#> #   primary_cell_count <chr>, published_at <date>, raw_data_location <chr>, revised_at <date>, sample_heuristic <chr>, …
```

### Query counts scaled per million

This is helpful if just few genes are of interest, as they can be
compared across samples.

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
#> # A SingleCellExperiment-tibble abstraction: 1,800 × 63
#> # [90mFeatures=56239 | Cells=1800 | Assays=cpm[0m
#>    .cell          dataset_id observation_joinid sample_id cell_type cell_type_ontology_t…¹ sample_ assay assay_ontology_term_id
#>    <chr>          <chr>      <chr>              <chr>     <chr>     <chr>                  <chr>   <chr> <chr>                 
#>  1 LAP87_GCTGGGT… 9f222629-… YK~|ZXL}LN         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x … EFO:0009922           
#>  2 LAP87_GACCCAG… 9f222629-… em#eUayTqZ         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x … EFO:0009922           
#>  3 LAP87_GCTGGGT… 9f222629-… 8(o)3D`Q%_         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x … EFO:0009922           
#>  4 TCATTTGCAAGCT… 9f222629-… {_UNJ6zNR<         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#>  5 CATCAAGGTCGAA… 9f222629-… z6RMY?)7(D         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#>  6 AGACGTTAGTAGA… 9f222629-… xF*TjD66E+         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#>  7 CCTTTCTTCCCAA… 9f222629-… w@YkRixyk@         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#>  8 GAAGCAGGTCGAC… 9f222629-… !Jd)V-~~>v         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#>  9 AAGCCGCTCTCGA… 9f222629-… zZEWSDl)l<         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#> 10 TACTCGCCAAAGG… 9f222629-… b`)P9|FF1D         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#> # ℹ 1,790 more rows
#> # ℹ abbreviated name: ¹​cell_type_ontology_term_id
#> # ℹ 54 more variables: cell_count <chr>, citation <chr>, collection_id <chr>, dataset_version_id <chr>,
#> #   development_stage <chr>, development_stage_ontology_term_id <chr>, disease <chr>, disease_ontology_term_id <chr>,
#> #   donor_id <chr>, experiment___ <chr>, explorer_url <chr>, feature_count <int>, filesize <dbl>, filetype <chr>,
#> #   is_primary_data <chr>, mean_genes_per_cell <dbl>, organism <chr>, organism_ontology_term_id <chr>,
#> #   primary_cell_count <chr>, published_at <date>, raw_data_location <chr>, revised_at <date>, sample_heuristic <chr>, …
```

### Extract only a subset of genes

``` r
single_cell_counts = 
    metadata |>
    dplyr::filter(
        self_reported_ethnicity == "African" &
        assay |> stringr::str_like("%10x%") &
        tissue == "lung parenchyma" &
        cell_type |> stringr::str_like("%CD4%")
    ) |>
    get_single_cell_experiment(assays = "cpm", features = "PUM1")
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Experiment.

single_cell_counts
#> # A SingleCellExperiment-tibble abstraction: 1,800 × 63
#> # [90mFeatures=0 | Cells=1800 | Assays=cpm[0m
#>    .cell          dataset_id observation_joinid sample_id cell_type cell_type_ontology_t…¹ sample_ assay assay_ontology_term_id
#>    <chr>          <chr>      <chr>              <chr>     <chr>     <chr>                  <chr>   <chr> <chr>                 
#>  1 LAP87_GCTGGGT… 9f222629-… YK~|ZXL}LN         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x … EFO:0009922           
#>  2 LAP87_GACCCAG… 9f222629-… em#eUayTqZ         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x … EFO:0009922           
#>  3 LAP87_GCTGGGT… 9f222629-… 8(o)3D`Q%_         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x … EFO:0009922           
#>  4 TCATTTGCAAGCT… 9f222629-… {_UNJ6zNR<         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#>  5 CATCAAGGTCGAA… 9f222629-… z6RMY?)7(D         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#>  6 AGACGTTAGTAGA… 9f222629-… xF*TjD66E+         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#>  7 CCTTTCTTCCCAA… 9f222629-… w@YkRixyk@         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#>  8 GAAGCAGGTCGAC… 9f222629-… !Jd)V-~~>v         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#>  9 AAGCCGCTCTCGA… 9f222629-… zZEWSDl)l<         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#> 10 TACTCGCCAAAGG… 9f222629-… b`)P9|FF1D         d0a88566… CD4-posi… CL:0000624             d0a885… 10x … EFO:0011025           
#> # ℹ 1,790 more rows
#> # ℹ abbreviated name: ¹​cell_type_ontology_term_id
#> # ℹ 54 more variables: cell_count <chr>, citation <chr>, collection_id <chr>, dataset_version_id <chr>,
#> #   development_stage <chr>, development_stage_ontology_term_id <chr>, disease <chr>, disease_ontology_term_id <chr>,
#> #   donor_id <chr>, experiment___ <chr>, explorer_url <chr>, feature_count <int>, filesize <dbl>, filetype <chr>,
#> #   is_primary_data <chr>, mean_genes_per_cell <dbl>, organism <chr>, organism_ontology_term_id <chr>,
#> #   primary_cell_count <chr>, published_at <date>, raw_data_location <chr>, revised_at <date>, sample_heuristic <chr>, …
```

### Extract the counts as a Seurat object

This convert the H5 SingleCellExperiment to Seurat so it might take long
time and occupy a lot of memory depending on how many cells you are
requesting.

``` r
single_cell_counts_seurat = 
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

single_cell_counts_seurat
#> An object of class Seurat 
#> 56239 features across 1800 samples within 1 assay 
#> Active assay: originalexp (56239 features, 0 variable features)
#>  2 layers present: counts, data
```

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
HLA-A across all tissues

``` r

# Plots with styling
counts <- metadata |>
  # Filter and subset
  dplyr::filter(cell_type_unified_ensemble == "cd14 mono") |>
  dplyr::filter(file_id_cellNexus_single_cell != "c5a05f23f9784a3be3bfa651198a48eb") |> 
  
  # Get counts per million for HCA-A gene
  get_single_cell_experiment(assays = "cpm", features = "HLA-A") |> 
  suppressMessages() |>
  
  # Add feature to table
  tidySingleCellExperiment::join_features("HLA-A", shape = "wide") |> 
    
  # Rank x axis
  tibble::as_tibble()

# Plot by disease
counts |>
  dplyr::with_groups(disease, ~ .x |> dplyr::mutate(median_count = median(`HLA.A`, rm.na=TRUE))) |> 
  
  # Plot
  ggplot(aes(forcats::fct_reorder(disease, median_count,.desc = TRUE), `HLA.A`,color = file_id)) +
  geom_jitter(shape=".") +
    
  # Style
  guides(color="none") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) + 
  xlab("Disease") + 
  ggtitle("HLA-A in CD14 monocytes by disease") 
```

![](man/figures/HLA_A_disease_plot.png)<!-- -->

``` r
# Plot by tissue
counts |> 
  dplyr::with_groups(tissue_harmonised, ~ .x |> dplyr::mutate(median_count = median(`HLA.A`, rm.na=TRUE))) |> 
  
  # Plot
  ggplot(aes(forcats::fct_reorder(tissue_harmonised, median_count,.desc = TRUE), `HLA.A`,color = file_id_cellNexus_single_cell)) +
  geom_jitter(shape=".") +
    
  # Style
  guides(color="none") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) + 
  xlab("Tissue") + 
  ggtitle("HLA-A in CD14 monocytes by tissue") + 
  theme(legend.position = "none")
```

![](man/figures/HLA_A_tissue_plot.png)<!-- -->

## Integrate Cloud Metadata with Local Metadata

CellNexus not only enables users to query our metadata but also allows
integration with your local metadata. Additionally, users can integrate
with your metadata stored in the cloud.

To enable this feature, users must include
`file_id_cellNexus_single_cell` (e.g my_sce.h5ad) and `atlas_id` (e.g
cellxgene/dd-mm-yy) columns in the metadata.

``` r
# Example of creating local metadata 
local_cache = tempdir()
meta_path = paste0(local_cache,"/my_metadata.parquet")
sce_path = paste0(local_cache, "/my_sce.h5ad")
cellNexus::sample_sce_obj |> S4Vectors::metadata() |> purrr::pluck("data") |> arrow::write_parquet(meta_path)
cellNexus::sample_sce_obj |> zellkonverter::writeH5AD(file = sce_path,
                                                      compression = "gzip")
#> ℹ Using the 'counts' assay as the X matrix
```

``` r
file_id_from_cloud <- "68aea852584d77f78e5c91204558316d___1.h5ad"
get_metadata(local_metadata = meta_path) |>
    dplyr::filter(file_id_cellNexus_single_cell %in% c("my_sce.h5ad",
                                                file_id_from_cloud)) |>
    get_single_cell_experiment()
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Experiment.
#> # A SingleCellExperiment-tibble abstraction: 17,791 × 63
#> # [90mFeatures=30867 | Cells=17791 | Assays=counts[0m
#>    .cell          dataset_id observation_joinid sample_id cell_type cell_type_ontology_t…¹ sample_ assay assay_ontology_term_id
#>    <chr>          <chr>      <chr>              <chr>     <chr>     <chr>                  <chr>   <chr> <chr>                 
#>  1 ATGGGAGGTGACT… 218acb0f-… AJ2j@WTK(8         b0588c32… CD8-posi… CL:0000625             b0588c… 10x … EFO:0009899           
#>  2 TTCTTAGTCAACA… 218acb0f-… IzL^|bq19d         ebc92036… CD8-posi… CL:0000625             ebc920… 10x … EFO:0009899           
#>  3 GAACGGATCAGGC… 218acb0f-… segTgRGKmZ         bc8e1a19… CD8-posi… CL:0000625             bc8e1a… 10x … EFO:0009899           
#>  4 GACAGAGCAAGTA… 218acb0f-… XF)oKe`?(T         c08b4008… CD8-posi… CL:0000625             c08b40… 10x … EFO:0009899           
#>  5 ATCATGGAGGAGT… 218acb0f-… 1#;Lh!~ntg         fefe4f5c… CD8-posi… CL:0000625             fefe4f… 10x … EFO:0009899           
#>  6 CTCTGGTCAAATT… 218acb0f-… ao^Q4xhsY|         49b29ac6… CD8-posi… CL:0000625             49b29a… 10x … EFO:0009899           
#>  7 CGGACACGTGATG… 218acb0f-… xwKwskwc{@         a304b977… CD8-posi… CL:0000625             a304b9… 10x … EFO:0009899           
#>  8 CATATTCTCACGA… 218acb0f-… H@+Lyu0Zu5         49b29ac6… CD8-posi… CL:0000625             49b29a… 10x … EFO:0009899           
#>  9 GTGTTAGTCGTTA… 218acb0f-… !t=WUM2sL<         891879fa… CD8-posi… CL:0000625             891879… 10x … EFO:0009899           
#> 10 AAGGTTCAGTAGC… 218acb0f-… <?rN+@p=y1         df563ecd… CD8-posi… CL:0000625             df563e… 10x … EFO:0009899           
#> # ℹ 17,781 more rows
#> # ℹ abbreviated name: ¹​cell_type_ontology_term_id
#> # ℹ 54 more variables: cell_count <chr>, citation <chr>, collection_id <chr>, dataset_version_id <chr>,
#> #   development_stage <chr>, development_stage_ontology_term_id <chr>, disease <chr>, disease_ontology_term_id <chr>,
#> #   donor_id <chr>, experiment___ <chr>, explorer_url <chr>, feature_count <int>, filesize <dbl>, filetype <chr>,
#> #   is_primary_data <chr>, mean_genes_per_cell <dbl>, organism <chr>, organism_ontology_term_id <chr>,
#> #   primary_cell_count <chr>, published_at <date>, raw_data_location <chr>, revised_at <date>, sample_heuristic <chr>, …
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
`observation_joinid`, `empty_droplet`

Through harmonisation and curation we introduced custom column, not
present in the original CELLxGENE metadata

- `tissue_harmonised`: a coarser tissue name for better filtering
- `age_days`: the number of days corresponding to the age
- `cell_type_unified_ensemble`: the consensus call identity (for immune
  cells) using the original and three novel annotations using Seurat
  Azimuth and SingleR
- `cell_annotation_azimuth_l2`: Azimuth cell annotation
- `cell_annotation_blueprint_singler`: SingleR cell annotation using
  Blueprint reference
- `cell_annotation_blueprint_monaco`: SingleR cell annotation using
  Monaco reference
- `sample_heuristic`: Sample subdivision for internal use
- `file_id_cellNexus_single_cell`: File subdivision for internal use
- `file_id_cellNexus_pseudobulk`: File subdivision for internal use
- `sample_id`: Sample ID
- `.sample_name`: How samples were defined

# RNA abundance

The `counts` assay includes RNA abundance in the positive real scale
(not transformed with non-linear functions, e.g. log sqrt). Originally
CELLxGENE include a mix of scales and transformations specified in the
`x_normalization` column.

The `cpm` assay includes counts per million.

The `rank` assay includes cells rank based on their total transcript
count.

# Session Info

``` r
sessionInfo()
#> R version 4.4.1 (2024-06-14)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Red Hat Enterprise Linux 9.5 (Plow)
#> 
#> Matrix products: default
#> BLAS:   /stornext/System/data/software/rhel/9/base/tools/R/4.4.1/lib64/R/lib/libRblas.so 
#> LAPACK: /stornext/System/data/software/rhel/9/base/tools/R/4.4.1/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.0
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
#> [1] cellNexus_1.5.10 testthat_3.2.1.1 assertthat_0.2.1 dplyr_1.1.4      ggplot2_3.5.1   
#> 
#> loaded via a namespace (and not attached):
#>   [1] fs_1.6.4                        matrixStats_1.3.0               spatstat.sparse_3.0-3          
#>   [4] devtools_2.4.5                  httr_1.4.7                      RColorBrewer_1.1-3             
#>   [7] profvis_0.3.8                   tools_4.4.1                     sctransform_0.4.1              
#>  [10] backports_1.5.0                 utf8_1.2.4                      R6_2.5.1                       
#>  [13] HDF5Array_1.32.1                lazyeval_0.2.2                  uwot_0.2.2                     
#>  [16] rhdf5filters_1.16.0             urlchecker_1.0.1                withr_3.0.1                    
#>  [19] sp_2.1-4                        gridExtra_2.3                   preprocessCore_1.66.0          
#>  [22] progressr_0.14.0                cli_3.6.3                       Biobase_2.64.0                 
#>  [25] spatstat.explore_3.2-7          fastDummies_1.7.3               sass_0.4.9                     
#>  [28] Seurat_5.1.0                    arrow_17.0.0                    spatstat.data_3.0-4            
#>  [31] readr_2.1.5                     ggridges_0.5.6                  pbapply_1.7-2                  
#>  [34] commonmark_1.9.1                parallelly_1.37.1               sessioninfo_1.2.2              
#>  [37] rstudioapi_0.16.0               generics_0.1.3                  ica_1.0-3                      
#>  [40] spatstat.random_3.2-3           Matrix_1.7-0                    fansi_1.0.6                    
#>  [43] S4Vectors_0.42.1                abind_1.4-5                     lifecycle_1.0.4                
#>  [46] yaml_2.3.10                     SummarizedExperiment_1.34.0     rhdf5_2.48.0                   
#>  [49] SparseArray_1.4.3               Rtsne_0.17                      grid_4.4.1                     
#>  [52] blob_1.2.4                      promises_1.3.0                  crayon_1.5.2                   
#>  [55] dir.expiry_1.12.0               miniUI_0.1.1.1                  lattice_0.22-6                 
#>  [58] cowplot_1.1.3                   pillar_1.9.0                    knitr_1.48                     
#>  [61] GenomicRanges_1.56.0            future.apply_1.11.2             codetools_0.2-20               
#>  [64] leiden_0.4.3.1                  glue_1.8.0                      data.table_1.16.2              
#>  [67] remotes_2.5.0                   tidySingleCellExperiment_1.15.4 vctrs_0.6.5                    
#>  [70] png_0.1-8                       spam_2.10-0                     gtable_0.3.5                   
#>  [73] cachem_1.0.8                    xfun_0.48                       S4Arrays_1.4.0                 
#>  [76] mime_0.12                       survival_3.7-0                  SingleCellExperiment_1.26.0    
#>  [79] ellipsis_0.3.2                  fitdistrplus_1.1-11             ROCR_1.0-11                    
#>  [82] nlme_3.1-166                    usethis_2.2.3                   bit64_4.0.5                    
#>  [85] filelock_1.0.3                  RcppAnnoy_0.0.22                GenomeInfoDb_1.40.0            
#>  [88] rprojroot_2.0.4                 bslib_0.7.0                     irlba_2.3.5.1                  
#>  [91] KernSmooth_2.23-24              colorspace_2.1-1                BiocGenerics_0.50.0            
#>  [94] DBI_1.2.2                       zellkonverter_1.12.1            duckdb_1.1.0                   
#>  [97] tidyselect_1.2.1                bit_4.0.5                       compiler_4.4.1                 
#> [100] curl_5.2.1                      basilisk.utils_1.14.1           xml2_1.3.6                     
#> [103] desc_1.4.3                      DelayedArray_0.30.1             plotly_4.10.4                  
#> [106] stringfish_0.16.0               checkmate_2.3.1                 scales_1.3.0                   
#> [109] lmtest_0.9-40                   stringr_1.5.1                   digest_0.6.35                  
#> [112] goftest_1.2-3                   spatstat.utils_3.0-5            rmarkdown_2.26                 
#> [115] basilisk_1.14.3                 XVector_0.44.0                  htmltools_0.5.8.1              
#> [118] pkgconfig_2.0.3                 MatrixGenerics_1.16.0           highr_0.11                     
#> [121] dbplyr_2.5.0                    fastmap_1.1.1                   rlang_1.1.4                    
#> [124] htmlwidgets_1.6.4               UCSC.utils_1.0.0                shiny_1.8.1.1                  
#> [127] jquerylib_0.1.4                 zoo_1.8-12                      jsonlite_1.8.8                 
#> [130] magrittr_2.0.3                  GenomeInfoDbData_1.2.12         dotCall64_1.1-1                
#> [133] patchwork_1.2.0                 Rhdf5lib_1.26.0                 munsell_0.5.1                  
#> [136] Rcpp_1.0.13                     reticulate_1.36.1               stringi_1.8.4                  
#> [139] brio_1.1.5                      zlibbioc_1.50.0                 MASS_7.3-61                    
#> [142] plyr_1.8.9                      pkgbuild_1.4.4                  parallel_4.4.1                 
#> [145] listenv_0.9.1                   ggrepel_0.9.5                   deldir_2.0-4                   
#> [148] splines_4.4.1                   tensor_1.5                      hms_1.1.3                      
#> [151] igraph_2.1.1                    spatstat.geom_3.2-9             RcppHNSW_0.6.0                 
#> [154] reshape2_1.4.4                  stats4_4.4.1                    pkgload_1.3.4                  
#> [157] tidybulk_1.17.3                 evaluate_1.0.1                  ttservice_0.4.0                
#> [160] SeuratObject_5.0.2              RcppParallel_5.1.7              tzdb_0.4.0                     
#> [163] httpuv_1.6.15                   RANN_2.6.2                      tidyr_1.3.1                    
#> [166] purrr_1.0.2                     polyclip_1.10-7                 future_1.33.2                  
#> [169] scattermore_1.2                 xtable_1.8-4                    RSpectra_0.16-1                
#> [172] roxygen2_7.3.2                  later_1.3.2                     viridisLite_0.4.2              
#> [175] tibble_3.2.1                    memoise_2.0.1                   IRanges_2.38.0                 
#> [178] cluster_2.1.6                   globals_0.16.3
```
