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
#' @export
#' @return A Seurat object containing the same data as a call to
#'   [get_single_cell_experiment()]
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
get_seurat <- function(...) {

  rlang::check_installed(c("Seurat", "SeuratObject"))

  get_single_cell_experiment(...) |>
    SeuratObject::as.Seurat(data = NULL)
}
