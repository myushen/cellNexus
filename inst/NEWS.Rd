\name{NEWS}
\title{News for Package \pkg{cellNexus}}

\section{News in version 0.99.31}{
\itemize{
    \item Added \code{donor_id} to harmonised metadata, enabling direct donor-level
    queries across datasets.
    \item Updated metadata parquet versions to \code{hca2024_v2.3.2}
    \item Expanded and corrected documentation for metadata columns on website.
    \item Removed internal columns (\code{sample_},
    \code{sample_heuristic}, \code{sample_chunk}, \code{cell_chunk}, 
    \code{sample_pseudobulk_chunk}) from documentation.
}}

\section{News in version 0.99.30}{
\itemize{
    \item \code{get_seurat()} now correctly handles multiple assays and non-default
    first assays.
    \itemize{
        \item When \code{assays = c("counts", "cpm")} (or any multi-assay request),
        all assays are transferred to the returned Seurat object: the first assay
        occupies the \code{"originalexp"} slot and every subsequent assay is added
        under its own name (e.g. \code{"cpm"}).
        \item When the first requested assay is not \code{"counts"} (e.g.
        \code{assays = "cpm"}), the conversion no longer errors with
        \emph{"No data in provided assay - counts"}. The first assay in the SCE is
        now passed explicitly to \code{SeuratObject::as.Seurat()}.
    }
    \item SCT assay: robust handling of samples for which SCT normalisation failed.
    \itemize{
        \item HTTP 404 responses for SCT h5ad files are now treated as
        \emph{"file not available"} rather than a fatal download error. The SCT
        sync path uses \code{ignore_not_found = TRUE} while all other assay syncs
        remain strict.
        \item A new internal helper \code{.skip_groups_missing_sct_file()} detects
        cached-file absence after sync, emits a \code{cellNexus says:} warning with
        the affected cell/sample counts and up to three example \code{sample_id}s,
        and instructs users to identify skipped samples post-retrieval via
        \code{dplyr::setdiff(unique(your_metadata\$sample_id),
        unique(SummarizedExperiment::colData(result)\$sample_id))}.
        \item When \code{assays = "sct"} only and \emph{every} sample is missing its
        SCT file, \code{get_single_cell_experiment()} aborts with a clear message.
        \item When \code{assays = c("counts", "sct")} (or any mixed-assay request)
        and \emph{every} sample is missing its SCT file, the function falls back
        gracefully to returning the non-SCT assays for all cells with a warning,
        rather than aborting.
    }
    \item Tests. Added coverage for \code{get_seurat()} multi-assay and
    non-default-assay paths.
}}

\section{News in version 0.99.27}{
\itemize{
    \item \code{get_pseudobulk()} now retains metadata columns that are constant within each
    \code{sample_id} \(\times\) \code{cell_type_unified_ensemble} combination, including
    user-added annotations, when building \code{colData}. Cell-level columns are dropped.
    \item \code{get_metacell()} uses the same once-per-query column selection for
    \code{sample_id} \(\times\) metacell-level annotations, so \code{cbind()} across
    files shares an identical \code{colData} schema.
    \item Added internal helpers \code{get_specific_annotation_columns()} and
    \code{keep_specific_annotation_columns()} to identify and keep columns functionally
    determined by selected key column(s). Works on local data frames and lazy DuckDB
    tables from \code{get_metadata()}.
    \item Tests. Added coverage for annotation-column detection and for preserving
    sample-grain user columns through \code{get_pseudobulk()}, using
    \code{SAMPLE_DATABASE_URL}.
}}

\section{News in version 0.99.10}{
\itemize{
    \item \code{get_pseudobulk()} gained \code{as_SummarizedExperiment} to optionally return a \code{SummarizedExperiment}. Coercion restores feature rownames, which are otherwise dropped by \code{as(x, "SummarizedExperiment")}.
    \item Tests. Added a test ensuring feature rownames are preserved when \code{as_SummarizedExperiment = TRUE}.
}}

\section{News in version 0.99.4}{
\itemize{
    \item Consistent features handling across retrieval functions. \code{get_pseudobulk()}, \code{get_single_cell_experiment()}, and \code{get_metacell()} now drop experiments/samples that do not contain all requested features and align rows to the requested feature set before merging. A warning indicates how many experiments were dropped.
    \item Stabilised feature subsetting during merge to avoid index/subscript errors and ensure identical row order across experiments prior to \code{cbind()}.
    \item Qualified functional calls (e.g., \code{purrr::map_lgl}) to avoid missing import issues.
    \item Documentation updates. Clarified the \code{features} parameter docs for all three functions to describe the drop-and-align behavior; regenerated help files.
    \item Tests. Added a test verifying pseudobulk subsetting to the specific feature \code{ENSG00000065485}.
}}

\section{News in version 0.1.0}{
\itemize{
    \item Added a NEWS file to track changes to the package.
}}
