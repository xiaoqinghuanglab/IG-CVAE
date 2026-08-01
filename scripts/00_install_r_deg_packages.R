#!/usr/bin/env Rscript

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  timeout = 1200
)

detected_cores <- parallel::detectCores(logical = FALSE)
ncpus <- if (is.na(detected_cores)) 1L else min(4L, max(1L, detected_cores))

all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)

if (length(file_arg) != 1L) {
  stop("Run this installer with Rscript scripts/00_install_r_deg_packages.R")
}

script_path <- normalizePath(
  sub("^--file=", "", file_arg),
  mustWork = TRUE
)
project_root <- normalizePath(
  file.path(dirname(script_path), ".."),
  mustWork = TRUE
)
local_library <- file.path(project_root, ".R", "library")

dir.create(local_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_library, .libPaths()))

cran_packages <- c(
  "data.table",
  "ggplot2",
  "ggrepel",
  "openxlsx",
  "statmod"
)

missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_cran)) {
  install.packages(
    missing_cran,
    lib = local_library,
    dependencies = NA,
    Ncpus = ncpus
  )
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    lib = local_library,
    Ncpus = ncpus
  )
}

if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install(
    "limma",
    lib = local_library,
    ask = FALSE,
    update = FALSE,
    Ncpus = ncpus
  )
}

required <- c("limma", cran_packages)

still_missing <- required[
  !vapply(
    required,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(still_missing)) {
  stop(
    "Package installation failed for: ",
    paste(still_missing, collapse = ", ")
  )
}

cat("\nDEG dependencies installed successfully.\n")
cat("R version: ", R.version.string, "\n", sep = "")
cat(
  "Bioconductor version: ",
  as.character(BiocManager::version()),
  "\n",
  sep = ""
)
cat("Project library: ", local_library, "\n", sep = "")

for (package in required) {
  cat(
    sprintf(
      "%-12s %s\n",
      package,
      as.character(utils::packageVersion(package))
    )
  )
}
