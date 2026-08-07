# Functions that relate to the Seurat conversion

#' @importFrom methods as
#' @exportS3Method SeuratObject::as.sparse DelayedMatrix
as.sparse.DelayedMatrix <- function(x, ...) {
  # This is glue to ensure the SCE -> Seurat conversion works properly with
  # DelayedArray types
  as(x, "dgCMatrix")
}

#' Given a data frame of HCA metadata, returns a Seurat object corresponding to
#' the samples in that data frame
#'
#' @inheritDotParams get_single_cell_experiment
#' @param download_only Logical scalar. When `TRUE`, remote files are
#'   synchronised to `cache_directory` but not read or assembled. Returns
#'   `invisible(NULL)`. See [get_single_cell_experiment()] for details.
#' @export
#' @return A Seurat object containing the same data as a call to
#'   [get_single_cell_experiment()], or `invisible(NULL)` when
#'   `download_only = TRUE`. All requested assays are present in the
#'   returned object under their original names (e.g. `"counts"`, `"cpm"`).
#' @importFrom SummarizedExperiment assayNames assay
#' @examples
#' # Use the lightweight sample database URL (for fast checks during development only)
#' meta <- get_metadata(cloud_metadata = cellNexus::SAMPLE_DATABASE_URL) |> head(2)
#' seurat <- get_seurat(meta)
#' @references Shen, M., Y. Gao, N. Liu, D. Bhuva, M. Milton, J. Henao,
#'   J. Andrews, E. Yang, C. Zhan, N. Liu, S. Si, J. W. Hutchison,
#'   M. H. Shakeel, M. Morgan, A. T. Papenfuss, J. Iskander, J. M. Polo,
#'   and S. Mangiola. "cellNexus: Quality control, annotation, aggregation
#'   and analytical layers for the Human Cell Atlas data." bioRxiv (2026).
#'   doi:10.64898/2026.04.14.718336.
#' @source [Shen et al.,2026](https://www.biorxiv.org/content/10.64898/2026.04.14.718336v3)
get_seurat <- function(..., download_only = FALSE) {
  rlang::check_installed(c("Seurat", "SeuratObject"))

  sce <- get_single_cell_experiment(..., download_only = download_only)

  if (download_only) {
    return(invisible(NULL))
  }

  sce_assays <- assayNames(sce)
  first_assay <- sce_assays[1]
  seurat_obj <- SeuratObject::as.Seurat(sce, counts = first_assay, data = NULL)
  # as.Seurat() always names the first assay "originalexp"; replace it with a
  # copy under the real name to avoid the misleading label.
  if (first_assay != "originalexp") {
    assay_obj <- seurat_obj[["originalexp"]]
    SeuratObject::Key(assay_obj) <- paste0(first_assay, "_")
    seurat_obj[[first_assay]] <- assay_obj
    SeuratObject::DefaultAssay(seurat_obj) <- first_assay
    seurat_obj[["originalexp"]] <- NULL
  }
  # Attach any additional assays that were not part of the initial conversion.
  if (length(sce_assays) > 1) {
    for (assay_name in sce_assays[-1]) {
      seurat_obj[[assay_name]] <- SeuratObject::CreateAssayObject(
        counts = assay(sce, assay_name)
      )
    }
  }
  seurat_obj
}
