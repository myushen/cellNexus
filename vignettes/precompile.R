# Pre-compute the cellNexus vignette and regenerate all derived outputs.
#
# Run this script whenever cellNexus.Rmd.orig is updated.
#  It executes all code once and produces three outputs:
#
#   1. vignettes/cellNexus.Rmd  — pre-computed vignette
#   2. vignettes/cellNexus.html — HTML vignette rendered from the above.
#   3. README.md                — GitHub-flavoured markdown rendered from the above.
#
# Usage (from the package root):
#   source("vignettes/precompile.R")

proj_root <- rprojroot::find_package_root_file() |>
  normalizePath()

# Pre-compute cellNexus.Rmd.
# knitr executes all code chunks and embeds the output as static markdown,
# so the resulting .Rmd contains no executable R code.
knitr::opts_knit$set(root.dir = proj_root)
knitr::knit("vignettes/cellNexus.Rmd.orig", "vignettes/cellNexus.Rmd")
knitr::opts_knit$set(root.dir = NULL)

# Render HTML vignette from the pre-computed file.
rmarkdown::render(
  "vignettes/cellNexus.Rmd",
  knit_root_dir = proj_root,
  output_format = BiocStyle::html_document(
    toc = TRUE,
    toc_float = TRUE
  )
)

# Render README.md from the pre-computed file.
rmarkdown::render(
  "vignettes/cellNexus.Rmd",
  output_file    = "README.md",
  output_format  = "github_document",
  output_dir     = proj_root,
  knit_root_dir  = proj_root
)
