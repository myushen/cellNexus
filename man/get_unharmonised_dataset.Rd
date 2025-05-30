% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unharmonised.R
\name{get_unharmonised_dataset}
\alias{get_unharmonised_dataset}
\title{Returns unharmonised metadata for selected datasets.}
\source{
\href{https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3}{Mangiola et al.,2023}
}
\usage{
get_unharmonised_dataset(
  dataset_id,
  cells = NULL,
  conn = dbConnect(drv = duckdb(), read_only = TRUE),
  remote_url = UNHARMONISED_URL,
  cache_directory = get_default_cache_dir()
)
}
\arguments{
\item{dataset_id}{A character vector, where each entry is a dataset ID
obtained from the \verb{$file_id_cellNexus_single_cell} column of the table returned from
\code{\link[=get_metadata]{get_metadata()}}}

\item{cells}{An optional character vector of cell IDs. If provided, only
metadata for those cells will be returned.}

\item{conn}{An optional DuckDB connection object. If provided, it will re-use
the existing connection instead of opening a new one.}

\item{remote_url}{Optional character vector of length 1. An HTTP URL pointing
to the root URL under which all the unharmonised dataset files are located.}

\item{cache_directory}{Optional character vector of length 1. A file path on
your local system to a directory (not a file) that will be used to store
the unharmonised metadata files.}
}
\value{
A named list, where each name is a dataset file ID, and each value is
a "lazy data frame", ie a \code{tbl}.
}
\description{
Various metadata fields are \emph{not} common between datasets, so it does not
make sense for these to live in the main metadata table. This function is a
utility that allows easy fetching of this data if necessary.
}
\examples{
\dontrun{
dataset <- "838ea006-2369-4e2c-b426-b2a744a2b02b"
harmonised_meta <- get_metadata() |> 
    dplyr::filter(file_id_cellNexus_single_cell == dataset) |> dplyr::collect()
unharmonised_meta <- get_unharmonised_dataset(dataset)
unharmonised_tbl <- dplyr::collect(unharmonised_meta[[dataset]])
dplyr::left_join(harmonised_meta, unharmonised_tbl, by=c("file_id_cellNexus_single_cell", "cell_id"))
}
}
\references{
Mangiola, S., M. Milton, N. Ranathunga, C. S. N. Li-Wai-Suen,
A. Odainic, E. Yang, W. Hutchison et al. "A multi-organ map of the human
immune system across age, sex and ethnicity." bioRxiv (2023): 2023-06.
doi:10.1101/2023.06.08.542671.
}
