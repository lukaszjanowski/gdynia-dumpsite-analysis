# Resolve the repository root from this script's location. This makes
# source("D:/.../run_all.R") work even when getwd() points elsewhere.
.get_this_file <- function() {
  frames <- sys.frames()
  ofiles <- vapply(frames, function(x) {
    if (!is.null(x$ofile)) as.character(x$ofile) else NA_character_
  }, character(1))
  ofiles <- ofiles[!is.na(ofiles) & nzchar(ofiles)]
  if (length(ofiles)) return(tail(ofiles, 1L))
  NA_character_
}

.this_file <- .get_this_file()
if (!is.na(.this_file)) {
  .project_root <- dirname(normalizePath(.this_file, winslash = "/", mustWork = TRUE))
  setwd(.project_root)
}

if (!file.exists("00_config.R") || !dir.exists("R")) {
  stop(
    "Repository root could not be determined. Open the .Rproj file or run setwd() ",
    "to the repository folder before sourcing run_all.R. Current directory: ", getwd()
  )
}

stages <- c(
  "R/01_inventory.R",
  "R/02_harmonise_dems.R",
  "R/03_vertical_coregistration.R",
  "R/04_dod_and_volumes.R",
  "R/05_morphometry_acoustics.R",
  "R/06_depositional_patches.R",
  "R/07_boundary_accounting.R",
  "R/08_sensitivity_analysis.R",
  "R/10_publication_tables.R",
  "R/09_publication_figures.R",
  "R/11_reproducibility_report.R"
)

for (stage in stages) {
  message("\n=== ", stage, " ===")
  source(stage, encoding = "UTF-8", chdir = FALSE)
}
message("Pipeline v2 completed.")
