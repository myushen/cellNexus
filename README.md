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

    #> # Source:   SQL [?? x 81]
    #> # Database: DuckDB v1.1.0 [unknown@Linux 5.14.0-362.24.1.el9_3.x86_64:R 4.4.1/:memory:]
    #>    cell_id                     dataset_id observation_joinid sample_id cell_type cell_type_ontology_t…¹ sample_ assay
    #>    <chr>                       <chr>      <chr>              <chr>     <chr>     <chr>                  <chr>   <chr>
    #>  1 ACTGAGTGTCACACGC-1-1-0-0-0… 218acb0f-… 9Wnf#P{Z$o         6c5d8c2d… classica… CL:0000860             6c5d8c… 10x …
    #>  2 ACACCAATCGGTGTTA-1-1-0-0-0… 218acb0f-… bo`(S3yH2r         6e3d4243… classica… CL:0000860             6e3d42… 10x …
    #>  3 TACGGATTCTTGTCAT-1-1-0-0-0… 218acb0f-… *K*WvP4DrD         c1bba717… conventi… CL:0000990             c1bba7… 10x …
    #>  4 CGGACGTTCACCCGAG-1-1-0-0-0… 218acb0f-… mY6vX<Sny~         bc368cd3… classica… CL:0000860             bc368c… 10x …
    #>  5 GGGCATCAGGTGTGGT-1-1-0-0-0… 218acb0f-… zsf!9I9U%i         447ff240… classica… CL:0000860             447ff2… 10x …
    #>  6 GTAGGCCCAAGAAAGG-1-1-0-0-0… 218acb0f-… ;fGeB3DSt?         25f5de36… conventi… CL:0000990             25f5de… 10x …
    #>  7 ATCCGAATCCTAGAAC-1-1-0-0-0… 218acb0f-… H)Uc-KXi2I         8f0e44ba… B cell    CL:0000236             8f0e44… 10x …
    #>  8 CGCGGTAAGCCACGCT-1-1-0-0-0… 218acb0f-… 8+|oWr+3G}         da8caef9… CD8-posi… CL:0000625             da8cae… 10x …
    #>  9 TGTATTCAGGCATGTG-1-1-0-0-0… 218acb0f-… iQFv3r%pLn         2041da57… CD8-posi… CL:0000625             2041da… 10x …
    #> 10 GTCTCGTGTCGAGTTT-1-1-0-0-0… 218acb0f-… wb2SsYxb}<         d9578857… classica… CL:0000860             d95788… 10x …
    #> # ℹ more rows
    #> # ℹ abbreviated name: ¹​cell_type_ontology_term_id
    #> # ℹ 73 more variables: assay_ontology_term_id <chr>, cell_count <chr>, citation <chr>, collection_id <chr>,
    #> #   dataset_version_id <chr>, default_embedding <chr>, development_stage <chr>,
    #> #   development_stage_ontology_term_id <chr>, disease <chr>, disease_ontology_term_id <chr>, donor_id <chr>,
    #> #   experiment___ <chr>, explorer_url <chr>, feature_count <int>, filesize <dbl>, filetype <chr>,
    #> #   is_primary_data <chr>, mean_genes_per_cell <dbl>, organism <chr>, organism_ontology_term_id <chr>, …

The `metadata` variable can then be re-used for all subsequent queries.

### Explore the tissue

``` r
metadata |>
    dplyr::distinct(tissue, cell_type_unified_ensemble) 
#> # Source:   SQL [?? x 2]
#> # Database: DuckDB v1.1.0 [unknown@Linux 5.14.0-362.24.1.el9_3.x86_64:R 4.4.1/:memory:]
#>    tissue cell_type_unified_ensemble
#>    <chr>  <chr>                     
#>  1 liver  macrophage                
#>  2 liver  tgd                       
#>  3 liver  cd16 mono                 
#>  4 lung   macrophage                
#>  5 lung   cd16 mono                 
#>  6 lung   ilc                       
#>  7 lung   pdc                       
#>  8 lung   tgd                       
#>  9 spleen cdc                       
#> 10 spleen cd4 th1/th17 em           
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
#> # A SingleCellExperiment-tibble abstraction: 1,800 × 83
#> # [90mFeatures=56239 | Cells=1800 | Assays=counts[0m
#>    .cell                       dataset_id observation_joinid sample_id cell_type cell_type_ontology_t…¹ sample_ assay
#>    <chr>                       <chr>      <chr>              <chr>     <chr>     <chr>                  <chr>   <chr>
#>  1 LAP87_ATCGTAGAGACTAGAT-1_d… 9f222629-… T-<6W?u>IS         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  2 LAP87_ATCTTCATCCTGTTAT-1_d… 9f222629-… (2|&?rOm9L         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  3 LAP87_GCTGGGTGTCATAGTC-1_d… 9f222629-… YK~|ZXL}LN         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  4 LAP87_GCTGGGTCAGCGTTTA-1_d… 9f222629-… 8(o)3D`Q%_         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  5 LAP87_TGCTTGCTCACCCTGT-1_d… 9f222629-… N=rI@P7%M#         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  6 LAP87_GACCCAGTCAAACGTC-1_d… 9f222629-… em#eUayTqZ         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  7 LAP87_TAGCACAGTTACGTAC-1_d… 9f222629-… N;1@x5XKt^         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  8 LAP87_AGCATCAAGTCGAAAT-1_d… 9f222629-… l1Cad2IZ8!         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  9 LAP87_CGTAATGGTTGGGTAG-1_d… 9f222629-… 26;xzkCsNW         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#> 10 LAP87_TACCGGGAGAAGCTCG-1_d… 9f222629-… 7XAbqq9MWA         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#> # ℹ 1,790 more rows
#> # ℹ abbreviated name: ¹​cell_type_ontology_term_id
#> # ℹ 75 more variables: assay_ontology_term_id <chr>, cell_count <chr>, citation <chr>, collection_id <chr>,
#> #   dataset_version_id <chr>, default_embedding <chr>, development_stage <chr>,
#> #   development_stage_ontology_term_id <chr>, disease <chr>, disease_ontology_term_id <chr>, donor_id <chr>,
#> #   experiment___ <chr>, explorer_url <chr>, feature_count <int>, filesize <dbl>, filetype <chr>,
#> #   is_primary_data <chr>, mean_genes_per_cell <dbl>, organism <chr>, organism_ontology_term_id <chr>, …
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
#> # A SingleCellExperiment-tibble abstraction: 1,800 × 83
#> # [90mFeatures=56239 | Cells=1800 | Assays=cpm[0m
#>    .cell                       dataset_id observation_joinid sample_id cell_type cell_type_ontology_t…¹ sample_ assay
#>    <chr>                       <chr>      <chr>              <chr>     <chr>     <chr>                  <chr>   <chr>
#>  1 LAP87_ATCGTAGAGACTAGAT-1_d… 9f222629-… T-<6W?u>IS         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  2 LAP87_ATCTTCATCCTGTTAT-1_d… 9f222629-… (2|&?rOm9L         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  3 LAP87_GCTGGGTGTCATAGTC-1_d… 9f222629-… YK~|ZXL}LN         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  4 LAP87_GCTGGGTCAGCGTTTA-1_d… 9f222629-… 8(o)3D`Q%_         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  5 LAP87_TGCTTGCTCACCCTGT-1_d… 9f222629-… N=rI@P7%M#         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  6 LAP87_GACCCAGTCAAACGTC-1_d… 9f222629-… em#eUayTqZ         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  7 LAP87_TAGCACAGTTACGTAC-1_d… 9f222629-… N;1@x5XKt^         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  8 LAP87_AGCATCAAGTCGAAAT-1_d… 9f222629-… l1Cad2IZ8!         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  9 LAP87_CGTAATGGTTGGGTAG-1_d… 9f222629-… 26;xzkCsNW         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#> 10 LAP87_TACCGGGAGAAGCTCG-1_d… 9f222629-… 7XAbqq9MWA         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#> # ℹ 1,790 more rows
#> # ℹ abbreviated name: ¹​cell_type_ontology_term_id
#> # ℹ 75 more variables: assay_ontology_term_id <chr>, cell_count <chr>, citation <chr>, collection_id <chr>,
#> #   dataset_version_id <chr>, default_embedding <chr>, development_stage <chr>,
#> #   development_stage_ontology_term_id <chr>, disease <chr>, disease_ontology_term_id <chr>, donor_id <chr>,
#> #   experiment___ <chr>, explorer_url <chr>, feature_count <int>, filesize <dbl>, filetype <chr>,
#> #   is_primary_data <chr>, mean_genes_per_cell <dbl>, organism <chr>, organism_ontology_term_id <chr>, …
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
#> # A SingleCellExperiment-tibble abstraction: 100 × 57
#> # [90mFeatures=56239 | Cells=100 | Assays=counts[0m
#>    .cell                  dataset_id sample_id sample_ assay assay_ontology_term_id cell_count citation collection_id
#>    <chr>                  <chr>      <chr>     <chr>   <chr> <chr>                  <chr>      <chr>    <chr>        
#>  1 9c8fa5a8d2ae37179b579… 9f222629-… 9c8fa5a8… 9c8fa5… 10x … EFO:0009922            2282447    Publica… 6f6d381a-770…
#>  2 9c8fa5a8d2ae37179b579… 9f222629-… 9c8fa5a8… 9c8fa5… 10x … EFO:0009922            2282447    Publica… 6f6d381a-770…
#>  3 d0a8856647d20b1fa1e83… 9f222629-… d0a88566… d0a885… 10x … EFO:0011025            2282447    Publica… 6f6d381a-770…
#>  4 e4d7f8162faf68a85f61b… 9f222629-… e4d7f816… e4d7f8… 10x … EFO:0011025            2282447    Publica… 6f6d381a-770…
#>  5 9c8fa5a8d2ae37179b579… 9f222629-… 9c8fa5a8… 9c8fa5… 10x … EFO:0009922            2282447    Publica… 6f6d381a-770…
#>  6 a2459ad4272363e6eb775… 9f222629-… a2459ad4… a2459a… 10x … EFO:0011025            2282447    Publica… 6f6d381a-770…
#>  7 bfe624d44f7e5868cc22e… 9f222629-… bfe624d4… bfe624… 10x … EFO:0011025            2282447    Publica… 6f6d381a-770…
#>  8 270eb221dd0456cc06324… 9f222629-… 270eb221… 270eb2… 10x … EFO:0009899            2282447    Publica… 6f6d381a-770…
#>  9 4f067f7e5f960bc72b071… 9f222629-… 4f067f7e… 4f067f… 10x … EFO:0009922            2282447    Publica… 6f6d381a-770…
#> 10 e4d7f8162faf68a85f61b… 9f222629-… e4d7f816… e4d7f8… 10x … EFO:0011025            2282447    Publica… 6f6d381a-770…
#> # ℹ 90 more rows
#> # ℹ 48 more variables: dataset_version_id <chr>, default_embedding <chr>, development_stage <chr>,
#> #   development_stage_ontology_term_id <chr>, disease <chr>, disease_ontology_term_id <chr>, donor_id <chr>,
#> #   experiment___ <chr>, explorer_url <chr>, feature_count <int>, filesize <dbl>, filetype <chr>,
#> #   is_primary_data <chr>, mean_genes_per_cell <dbl>, organism <chr>, organism_ontology_term_id <chr>,
#> #   primary_cell_count <chr>, published_at <chr>, raw_data_location <chr>, revised_at <chr>, run_from_cell_id <chr>,
#> #   sample_heuristic <chr>, schema_version <chr>, self_reported_ethnicity <chr>, …
```

### Query metacell

The metadata includes a series of metacell aggregation levels, beginning
with 2, 4, 8, and so on. For example, the value of metacell_2 represents
a grouping of cells that can be split into two distinct metacells.

``` r
metacell_counts = 
   metadata |>
    dplyr::filter(
        self_reported_ethnicity == "African" &
        assay |> stringr::str_like("%10x%") &
        tissue == "lung parenchyma" &
        cell_type |> stringr::str_like("%CD4%")
    ) |>
    dplyr::filter(!is.na(metacell_2)) |>
    get_metacell(cell_aggregation = "metacell_2")
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Experiment.

metacell_counts
#> # A SingleCellExperiment-tibble abstraction: 713 × 40
#> # [90mFeatures=56239 | Cells=713 | Assays=counts[0m
#>    .cell        metacell_2 dataset_id sample_id assay assay_ontology_term_id development_stage development_stage_on…¹
#>    <chr>             <dbl> <chr>      <chr>     <chr> <chr>                  <chr>             <chr>                 
#>  1 d0a8856647d…         50 9f222629-… d0a88566… 10x … EFO:0011025            55-year-old huma… HsapDv:0000149        
#>  2 d0a8856647d…         34 9f222629-… d0a88566… 10x … EFO:0011025            55-year-old huma… HsapDv:0000149        
#>  3 d0a8856647d…         85 9f222629-… d0a88566… 10x … EFO:0011025            55-year-old huma… HsapDv:0000149        
#>  4 d0a8856647d…         13 9f222629-… d0a88566… 10x … EFO:0011025            55-year-old huma… HsapDv:0000149        
#>  5 d0a8856647d…         24 9f222629-… d0a88566… 10x … EFO:0011025            55-year-old huma… HsapDv:0000149        
#>  6 d0a8856647d…          7 9f222629-… d0a88566… 10x … EFO:0011025            55-year-old huma… HsapDv:0000149        
#>  7 d0a8856647d…         29 9f222629-… d0a88566… 10x … EFO:0011025            55-year-old huma… HsapDv:0000149        
#>  8 d0a8856647d…          3 9f222629-… d0a88566… 10x … EFO:0011025            55-year-old huma… HsapDv:0000149        
#>  9 d0a8856647d…          8 9f222629-… d0a88566… 10x … EFO:0011025            55-year-old huma… HsapDv:0000149        
#> 10 d0a8856647d…         52 9f222629-… d0a88566… 10x … EFO:0011025            55-year-old huma… HsapDv:0000149        
#> # ℹ 703 more rows
#> # ℹ abbreviated name: ¹​development_stage_ontology_term_id
#> # ℹ 32 more variables: disease <chr>, disease_ontology_term_id <chr>, donor_id <chr>, experiment___ <chr>,
#> #   explorer_url <chr>, feature_count <int>, is_primary_data <chr>, organism <chr>, organism_ontology_term_id <chr>,
#> #   published_at <chr>, raw_data_location <chr>, revised_at <chr>, sample_heuristic <chr>, schema_version <chr>,
#> #   self_reported_ethnicity <chr>, self_reported_ethnicity_ontology_term_id <chr>, sex <chr>,
#> #   sex_ontology_term_id <chr>, tissue <chr>, tissue_ontology_term_id <chr>, tissue_type <chr>, title <chr>, …
```

### Extract only a subset of genes

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
    get_single_cell_experiment(assays = "cpm", features = "PUM1")
#> ℹ Realising metadata.
#> ℹ Synchronising files
#> ℹ Reading files.
#> ℹ Compiling Experiment.

single_cell_counts
#> # A SingleCellExperiment-tibble abstraction: 1,800 × 83
#> # [90mFeatures=0 | Cells=1800 | Assays=cpm[0m
#>    .cell                       dataset_id observation_joinid sample_id cell_type cell_type_ontology_t…¹ sample_ assay
#>    <chr>                       <chr>      <chr>              <chr>     <chr>     <chr>                  <chr>   <chr>
#>  1 LAP87_ATCGTAGAGACTAGAT-1_d… 9f222629-… T-<6W?u>IS         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  2 LAP87_ATCTTCATCCTGTTAT-1_d… 9f222629-… (2|&?rOm9L         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  3 LAP87_GCTGGGTGTCATAGTC-1_d… 9f222629-… YK~|ZXL}LN         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  4 LAP87_GCTGGGTCAGCGTTTA-1_d… 9f222629-… 8(o)3D`Q%_         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  5 LAP87_TGCTTGCTCACCCTGT-1_d… 9f222629-… N=rI@P7%M#         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  6 LAP87_GACCCAGTCAAACGTC-1_d… 9f222629-… em#eUayTqZ         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  7 LAP87_TAGCACAGTTACGTAC-1_d… 9f222629-… N;1@x5XKt^         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  8 LAP87_AGCATCAAGTCGAAAT-1_d… 9f222629-… l1Cad2IZ8!         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#>  9 LAP87_CGTAATGGTTGGGTAG-1_d… 9f222629-… 26;xzkCsNW         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#> 10 LAP87_TACCGGGAGAAGCTCG-1_d… 9f222629-… 7XAbqq9MWA         9c8fa5a8… CD4-posi… CL:0000624             9c8fa5… 10x …
#> # ℹ 1,790 more rows
#> # ℹ abbreviated name: ¹​cell_type_ontology_term_id
#> # ℹ 75 more variables: assay_ontology_term_id <chr>, cell_count <chr>, citation <chr>, collection_id <chr>,
#> #   dataset_version_id <chr>, default_embedding <chr>, development_stage <chr>,
#> #   development_stage_ontology_term_id <chr>, disease <chr>, disease_ontology_term_id <chr>, donor_id <chr>,
#> #   experiment___ <chr>, explorer_url <chr>, feature_count <int>, filesize <dbl>, filetype <chr>,
#> #   is_primary_data <chr>, mean_genes_per_cell <dbl>, organism <chr>, organism_ontology_term_id <chr>, …
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
#> # A SingleCellExperiment-tibble abstraction: 2 × 83
#> # [90mFeatures=30867 | Cells=2 | Assays=counts[0m
#>   .cell dataset_id observation_joinid sample_id cell_type cell_type_ontology_t…¹ sample_ assay assay_ontology_term_id
#>   <chr> <chr>      <chr>              <chr>     <chr>     <chr>                  <chr>   <chr> <chr>                 
#> 1 GCTG… 218acb0f-… GA?EQwaM~g         33e161a6… CD8-posi… CL:0000625             33e161… 10x … EFO:0009899           
#> 2 ACAT… 218acb0f-… {QEUA;{iYj         118dde9a… CD8-posi… CL:0000625             118dde… 10x … EFO:0009899           
#> # ℹ abbreviated name: ¹​cell_type_ontology_term_id
#> # ℹ 74 more variables: cell_count <chr>, citation <chr>, collection_id <chr>, dataset_version_id <chr>,
#> #   default_embedding <chr>, development_stage <chr>, development_stage_ontology_term_id <chr>, disease <chr>,
#> #   disease_ontology_term_id <chr>, donor_id <chr>, experiment___ <chr>, explorer_url <chr>, feature_count <int>,
#> #   filesize <dbl>, filetype <chr>, is_primary_data <chr>, mean_genes_per_cell <dbl>, organism <chr>,
#> #   organism_ontology_term_id <chr>, primary_cell_count <chr>, published_at <chr>, raw_data_location <chr>,
#> #   revised_at <chr>, run_from_cell_id <chr>, sample_heuristic <chr>, schema_version <chr>, …
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

# Session Info

``` r
sessionInfo()
#> R version 4.4.1 (2024-06-14)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Red Hat Enterprise Linux 9.3 (Plow)
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
#> [1] ggplot2_3.5.1    cellNexus_1.6.11
#> 
#> loaded via a namespace (and not attached):
#>   [1] RcppAnnoy_0.0.22                splines_4.4.1                   later_1.3.2                    
#>   [4] filelock_1.0.3                  tibble_3.2.1                    polyclip_1.10-7                
#>   [7] preprocessCore_1.66.0           basilisk.utils_1.14.1           fastDummies_1.7.3              
#>  [10] lifecycle_1.0.4                 rprojroot_2.0.4                 globals_0.16.3                 
#>  [13] lattice_0.22-6                  MASS_7.3-61                     backports_1.5.0                
#>  [16] magrittr_2.0.3                  sass_0.4.9                      plotly_4.10.4                  
#>  [19] rmarkdown_2.26                  jquerylib_0.1.4                 yaml_2.3.10                    
#>  [22] httpuv_1.6.15                   Seurat_5.1.0                    sctransform_0.4.1              
#>  [25] spam_2.10-0                     sp_2.1-4                        spatstat.sparse_3.0-3          
#>  [28] reticulate_1.36.1               cowplot_1.1.3                   pbapply_1.7-2                  
#>  [31] DBI_1.2.2                       RColorBrewer_1.1-3              abind_1.4-5                    
#>  [34] zlibbioc_1.50.0                 Rtsne_0.17                      GenomicRanges_1.56.0           
#>  [37] purrr_1.0.2                     BiocGenerics_0.50.0             tidySingleCellExperiment_1.15.4
#>  [40] GenomeInfoDbData_1.2.12         IRanges_2.38.0                  S4Vectors_0.42.1               
#>  [43] ggrepel_0.9.5                   irlba_2.3.5.1                   listenv_0.9.1                  
#>  [46] spatstat.utils_3.0-5            goftest_1.2-3                   RSpectra_0.16-1                
#>  [49] spatstat.random_3.2-3           fitdistrplus_1.1-11             parallelly_1.37.1              
#>  [52] commonmark_1.9.1                leiden_0.4.3.1                  codetools_0.2-20               
#>  [55] DelayedArray_0.30.1             tidyselect_1.2.1                UCSC.utils_1.0.0               
#>  [58] matrixStats_1.3.0               stats4_4.4.1                    spatstat.explore_3.2-7         
#>  [61] duckdb_1.1.0                    jsonlite_1.8.8                  ellipsis_0.3.2                 
#>  [64] progressr_0.14.0                ggridges_0.5.6                  survival_3.7-0                 
#>  [67] tools_4.4.1                     ica_1.0-3                       Rcpp_1.0.13                    
#>  [70] glue_1.8.0                      gridExtra_2.3                   SparseArray_1.4.3              
#>  [73] xfun_0.48                       MatrixGenerics_1.16.0           GenomeInfoDb_1.40.0            
#>  [76] dplyr_1.1.4                     HDF5Array_1.32.1                withr_3.0.1                    
#>  [79] fastmap_1.1.1                   basilisk_1.14.3                 rhdf5filters_1.16.0            
#>  [82] fansi_1.0.6                     ttservice_0.4.0                 digest_0.6.35                  
#>  [85] R6_2.5.1                        mime_0.12                       colorspace_2.1-1               
#>  [88] scattermore_1.2                 tensor_1.5                      spatstat.data_3.0-4            
#>  [91] utf8_1.2.4                      tidyr_1.3.1                     generics_0.1.3                 
#>  [94] data.table_1.16.2               httr_1.4.7                      htmlwidgets_1.6.4              
#>  [97] S4Arrays_1.4.0                  uwot_0.2.2                      pkgconfig_2.0.3                
#> [100] gtable_0.3.5                    blob_1.2.4                      lmtest_0.9-40                  
#> [103] SingleCellExperiment_1.26.0     XVector_0.44.0                  brio_1.1.5                     
#> [106] htmltools_0.5.8.1               dotCall64_1.1-1                 SeuratObject_5.0.2             
#> [109] scales_1.3.0                    Biobase_2.64.0                  png_0.1-8                      
#> [112] knitr_1.48                      rstudioapi_0.16.0               tzdb_0.4.0                     
#> [115] reshape2_1.4.4                  checkmate_2.3.1                 nlme_3.1-166                   
#> [118] curl_5.2.1                      cachem_1.0.8                    zoo_1.8-12                     
#> [121] rhdf5_2.48.0                    stringr_1.5.1                   KernSmooth_2.23-24             
#> [124] parallel_4.4.1                  miniUI_0.1.1.1                  arrow_17.0.0                   
#> [127] zellkonverter_1.14.1            pillar_1.9.0                    grid_4.4.1                     
#> [130] vctrs_0.6.5                     RANN_2.6.2                      promises_1.3.0                 
#> [133] stringfish_0.16.0               dbplyr_2.5.0                    xtable_1.8-4                   
#> [136] cluster_2.1.6                   evaluate_1.0.1                  readr_2.1.5                    
#> [139] cli_3.6.3                       compiler_4.4.1                  rlang_1.1.4                    
#> [142] crayon_1.5.2                    future.apply_1.11.2             tidybulk_1.17.3                
#> [145] plyr_1.8.9                      stringi_1.8.4                   viridisLite_0.4.2              
#> [148] deldir_2.0-4                    assertthat_0.2.1                munsell_0.5.1                  
#> [151] lazyeval_0.2.2                  spatstat.geom_3.2-9             Matrix_1.7-0                   
#> [154] dir.expiry_1.12.0               RcppHNSW_0.6.0                  hms_1.1.3                      
#> [157] patchwork_1.2.0                 bit64_4.6.0-1                   future_1.33.2                  
#> [160] Rhdf5lib_1.26.0                 shiny_1.8.1.1                   highr_0.11                     
#> [163] SummarizedExperiment_1.34.0     ROCR_1.0-11                     igraph_2.1.1                   
#> [166] bslib_0.7.0                     RcppParallel_5.1.7              bit_4.0.5
```
