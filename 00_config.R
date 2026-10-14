required_packages <- c("terra", "sf", "data.table", "yaml", "ggplot2", "patchwork", "ggspatial", "scales", "cowplot")
missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(terra);library(sf);library(data.table)})
ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
CFG <- yaml::read_yaml(file.path(ROOT,"config","config.yml"))
OUT <- file.path(ROOT,"outputs")
OUTPUT_DIR <- OUT
DIRS <- c(rasters=file.path(OUT,"rasters"), vectors=file.path(OUT,"vectors"), tables=file.path(OUT,"tables"), figures=file.path(OUT,"figures"), logs=file.path(ROOT,"logs"))
invisible(lapply(DIRS, dir.create, recursive=TRUE, showWarnings=FALSE))
MAP_CRS <- CFG$project$analysis_crs
TARGET_RES <- CFG$analysis$target_resolution_m

# Vertical coordinate convention:
# "auto" detects whether input rasters use elevation or positive depth.
Z_MODE <- "auto"

LOD_VALUES <- unlist(CFG$analysis$lod_values_m)
PATCH_THRESHOLD_M <- CFG$analysis$primary_positive_threshold_m
MIN_PATCH_AREA_M2 <- CFG$analysis$primary_min_patch_area_m2
LARGE_PATCH_AREA_M2 <- CFG$analysis$large_patch_area_m2
CONNECT_DIRECTIONS <- CFG$analysis$connectivity_directions
COREG_APPLY_VERTICAL <- CFG$analysis$apply_vertical_coregistration
COREG_MAX_ABS_SHIFT_M <- CFG$analysis$max_vertical_shift_m
COREG_MIN_STABLE_CELLS <- CFG$analysis$minimum_stable_cells
# Compatibility aliases expected by the verified co-registration script.
MIN_STABLE_CELLS <- COREG_MIN_STABLE_CELLS
MAX_ABS_SHIFT_M <- COREG_MAX_ABS_SHIFT_M
options(stringsAsFactors=FALSE, scipen=999)

# Resolve required vector inputs from config, with GeoPackage/Shapefile fallback.
resolve_vector_input <- function(configured_path, stem, required = TRUE) {
  configured <- if (grepl("^([A-Za-z]:[/\\\\]|/)", configured_path)) configured_path else file.path(ROOT, configured_path)
  raw_dir <- file.path(ROOT, "data", "raw")
  candidates <- unique(c(
    configured,
    file.path(raw_dir, paste0(stem, ".gpkg")),
    file.path(raw_dir, paste0(stem, ".shp")),
    file.path(raw_dir, paste0(tolower(stem), ".gpkg")),
    file.path(raw_dir, paste0(tolower(stem), ".shp"))
  ))
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    if (required) stop("Vector input not found for " , stem, ". Checked:\n", paste(candidates, collapse = "\n"))
    return(NA_character_)
  }
  errors <- character()
  for (candidate in hit) {
    chosen <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
    test <- try(terra::vect(chosen), silent = TRUE)
    if (inherits(test, "try-error")) {
      errors <- c(errors, paste0(chosen, ": ", as.character(test)))
      next
    }
    if (nrow(test) < 1L) { errors <- c(errors, paste0(chosen, ": no features")); next }
    test_crs <- terra::crs(test)
    if (is.na(test_crs) || !nzchar(test_crs)) { errors <- c(errors, paste0(chosen, ": undefined CRS")); next }
    message("Using vector input: " , chosen)
    return(chosen)
  }
  stop("No readable vector input found for " , stem, ". Errors:\n", paste(errors, collapse = "\n"))
}

DUMPSITE_VECTOR <- resolve_vector_input(CFG$files$dumpsite, "DUMPSITE_VECTOR")
STABLE_VECTOR <- resolve_vector_input(CFG$files$stable, "STABLE_VECTOR")
# Optional analysis boundary.
# NULL means that the common valid coverage of the three DEMs defines
# the analysis domain.
ANALYSIS_VECTOR <- NULL

# Compatibility alias required by the co-registration module.
COREG_MAX_ABS_SHIFT <- COREG_MAX_ABS_SHIFT_M

LOD_PRIMARY <- CFG$analysis$primary_lod_m

# HOTFIX 04: reproducible 200 m dumpsite buffer
BOUNDARY_BUFFER_M <- CFG$analysis$external_zone_m
if (is.null(BOUNDARY_BUFFER_M) || !is.finite(BOUNDARY_BUFFER_M) || BOUNDARY_BUFFER_M <= 0) {
  stop("analysis.external_zone_m must be a positive number in config/config.yml")
}

# The buffer is a derived product, not a required raw input.
# Downstream scripts should create it from the official dumpsite boundary.
make_dumpsite_buffer <- function(dumpsite, width_m = BOUNDARY_BUFFER_M) {
  if (inherits(dumpsite, "SpatVector")) {
    if (terra::is.lonlat(dumpsite)) stop("Dumpsite must be projected before buffering.")
    return(terra::buffer(dumpsite, width = width_m))
  }
  if (inherits(dumpsite, "sf") || inherits(dumpsite, "sfc")) {
    if (sf::st_is_longlat(dumpsite)) stop("Dumpsite must be projected before buffering.")
    return(sf::st_buffer(dumpsite, dist = width_m))
  }
  stop("Unsupported dumpsite object class: ", paste(class(dumpsite), collapse = ", "))
}

# Compatibility path for scripts that still expect a BUFFER_VECTOR variable.
# The file is generated after stage 02 from outputs/vectors/dumpsite.gpkg.
BUFFER_VECTOR <- file.path(DIRS["vectors"], "dumpsite_buffer_200m.gpkg")
