# Fast loop for README / figure-layout iteration.
#
# Knits temp_cellNexus.Rmd.orig -> temp_cellNexus.Rmd -> README.temp.md
# without running the full vignette or overwriting README.md.
#
# Usage (from the package root):
#   source("vignettes/temp_precompile.R")

proj_root <- rprojroot::find_package_root_file() |>
  normalizePath()
vig_dir <- file.path(proj_root, "vignettes")

# Chunk wd = vignettes/: plot files land next to the .Rmd (portable paths),
# static images use ../man/figures/ so R CMD build finds man/ from vignettes/.
knitr::opts_knit$set(
  root.dir = vig_dir,
  unnormalized.path = TRUE
)
knitr::knit(
  "vignettes/temp_cellNexus.Rmd.orig",
  "vignettes/temp_cellNexus.Rmd"
)
knitr::opts_knit$set(root.dir = NULL, unnormalized.path = NULL)

# README.temp.md lives at repo root: rewrite vignette-relative paths after render.
rewrite_readme_paths <- function(readme_path) {
  x <- readLines(readme_path, warn = FALSE, encoding = "UTF-8")
  x <- gsub("../man/figures/", "man/figures/", x, fixed = TRUE)
  x <- sub("(\\]\\(|src=\")plot-", "\\1vignettes/plot-", x, perl = TRUE)
  writeLines(x, readme_path, useBytes = FALSE)
}

rmarkdown::render(
  "vignettes/temp_cellNexus.Rmd",
  output_file   = "README.md",
  output_format = "github_document",
  output_dir    = proj_root,
  knit_root_dir = vig_dir
)
rewrite_readme_paths(file.path(proj_root, "README.temp.md"))
