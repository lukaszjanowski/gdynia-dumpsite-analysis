# ==============================================================================
# Gdynia dumpsite analysis
# Harmonisation and vertical co-registration of 2011, 2019 and 2023 DEMs
# ==============================================================================


# 1. Load configuration ---------------------------------------------------------

source(
  "00_config.R",
  encoding = "UTF-8"
)


# 2. Load paths identified during input QC -------------------------------------

paths_file <- file.path(
  OUT,
  "paths.rds"
)

if (!file.exists(paths_file)) {
  stop(
    "File paths.rds not found. Run 01_locate_and_qc.R first: ",
    paths_file
  )
}

paths <- readRDS(paths_file)


# 3. Validate required DEM paths ------------------------------------------------

required_dem_names <- c(
  "dem2011",
  "dem2019",
  "dem2023"
)

missing_dem_names <- setdiff(
  required_dem_names,
  names(paths)
)

if (length(missing_dem_names) > 0) {
  stop(
    "Missing DEM entries in paths.rds: ",
    paste(missing_dem_names, collapse = ", ")
  )
}

for (dem_name in required_dem_names) {
  if (!file.exists(paths[[dem_name]])) {
    stop(
      "DEM file not found for ",
      dem_name,
      ": ",
      paths[[dem_name]]
    )
  }
}


# 4. Function for reading and normalising DEMs ---------------------------------

read_dem <- function(file_path, survey_name) {

  message(
    "Reading ",
    survey_name,
    ": ",
    file_path
  )

  x <- terra::rast(file_path)

  if (terra::nlyr(x) > 1) {
    warning(
      survey_name,
      " contains more than one band. Only the first band will be used."
    )

    x <- x[[1]]
  }

  # Determine the central DEM value without loading the full raster into RAM.
  central_value <- terra::global(
    x,
    fun = "mean",
    na.rm = TRUE
  )[1, 1]

  if (!is.finite(central_value)) {
    stop(
      "Could not calculate a representative value for ",
      survey_name
    )
  }

  message(
    survey_name,
    " mean source value: ",
    round(central_value, 3)
  )

  # Convert positive water depth to negative seabed elevation if required.
  if (Z_MODE == "positive_depth") {

    x <- -x

    message(
      survey_name,
      ": sign changed from positive depth to seabed elevation."
    )

  } else if (
    Z_MODE == "auto" &&
    central_value > 0
  ) {

    x <- -x

    message(
      survey_name,
      ": positive values detected; sign changed automatically."
    )

  } else {

    message(
      survey_name,
      ": values retained without sign change."
    )
  }

  # Reproject only if required.
  if (!terra::same.crs(x, MAP_CRS)) {

    message(
      survey_name,
      ": projecting to ",
      MAP_CRS
    )

    x <- terra::project(
      x,
      MAP_CRS,
      method = "bilinear"
    )

  } else {

    message(
      survey_name,
      ": CRS already matches ",
      MAP_CRS
    )
  }

  names(x) <- survey_name

  return(x)
}


# 5. Read DEMs ------------------------------------------------------------------

d11_original <- read_dem(
  paths$dem2011,
  "DEM_2011"
)

d19_original <- read_dem(
  paths$dem2019,
  "DEM_2019"
)

d24_original <- read_dem(
  paths$dem2023,
  "DEM_2023"
)


# 6. Create the common 1 m reference grid --------------------------------------

# The 2019 DEM is used as the alignment reference.
reference_extent <- terra::ext(d19_original)

aligned_extent <- terra::ext(
  floor(reference_extent[1] / TARGET_RES) * TARGET_RES,
  ceiling(reference_extent[2] / TARGET_RES) * TARGET_RES,
  floor(reference_extent[3] / TARGET_RES) * TARGET_RES,
  ceiling(reference_extent[4] / TARGET_RES) * TARGET_RES
)

template <- terra::rast(
  aligned_extent,
  resolution = TARGET_RES,
  crs = MAP_CRS
)

names(template) <- "analysis_template"

message(
  "Reference grid created: ",
  terra::ncol(template),
  " columns x ",
  terra::nrow(template),
  " rows; resolution ",
  paste(terra::res(template), collapse = " x "),
  " m."
)


# 7. Align all DEMs to the common grid -----------------------------------------

align_dem <- function(x, survey_name) {

  message(
    "Aligning ",
    survey_name,
    " to the common grid."
  )

  y <- terra::resample(
    x,
    template,
    method = "bilinear"
  )

  names(y) <- survey_name

  return(y)
}

d11 <- align_dem(
  d11_original,
  "DEM_2011"
)

d19 <- align_dem(
  d19_original,
  "DEM_2019"
)

d24 <- align_dem(
  d24_original,
  "DEM_2023"
)


# 8. Read vector layers ---------------------------------------------------------

read_vector_checked <- function(
  file_path,
  layer_name,
  required = TRUE
) {

  if (
    is.null(file_path) ||
    length(file_path) == 0 ||
    is.na(file_path)
  ) {

    if (required) {
      stop(
        layer_name,
        " path is missing."
      )
    }

    return(NULL)
  }

  if (!file.exists(file_path)) {

    if (required) {
      stop(
        layer_name,
        " file not found: ",
        file_path
      )
    }

    return(NULL)
  }

  message(
    "Reading ",
    layer_name,
    ": ",
    file_path
  )

  v <- terra::vect(file_path)

  if (is.na(terra::crs(v))) {
    stop(
      layer_name,
      " does not have a defined coordinate reference system."
    )
  }

  if (!terra::same.crs(v, template)) {

    message(
      "Projecting ",
      layer_name,
      " to ",
      MAP_CRS
    )

    v <- terra::project(
      v,
      terra::crs(template)
    )
  }

  return(v)
}


dumpsite <- read_vector_checked(
  DUMPSITE_VECTOR,
  "dumpsite boundary",
  required = TRUE
)

stable <- read_vector_checked(
  STABLE_VECTOR,
  "stable seabed polygons",
  required = TRUE
)

analysis <- read_vector_checked(
  ANALYSIS_VECTOR,
  "analysis boundary",
  required = FALSE
)


# 9. Check vector geometry ------------------------------------------------------

if (nrow(dumpsite) < 1) {
  stop(
    "DUMPSITE_VECTOR does not contain any polygons."
  )
}

if (nrow(stable) < 1) {
  stop(
    "STABLE_VECTOR does not contain any polygons."
  )
}

message(
  "Dumpsite features: ",
  nrow(dumpsite)
)

message(
  "Stable-area features: ",
  nrow(stable)
)

message(
  "Total dumpsite area: ",
  round(
    sum(
      terra::expanse(
        dumpsite,
        unit = "km"
      )
    ),
    4
  ),
  " km2"
)

message(
  "Total stable-area extent: ",
  round(
    sum(
      terra::expanse(
        stable,
        unit = "km"
      )
    ),
    4
  ),
  " km2"
)


# 10. Create common valid-data mask --------------------------------------------

# All rasters already have identical geometry after resampling.
dem_stack <- c(
  d11,
  d19,
  d24
)

names(dem_stack) <- c(
  "DEM_2011",
  "DEM_2019",
  "DEM_2023_raw"
)

common_mask <- (
  !is.na(d11) &
  !is.na(d19) &
  !is.na(d24)
)

names(common_mask) <- "common_valid_coverage"

# Convert logical values to 1 and NA:
# 1 = all three DEMs contain valid data,
# NA = at least one DEM does not contain valid data.
common_mask <- terra::ifel(
  common_mask,
  1,
  NA
)

# If an explicit analysis polygon is provided, apply it additionally.
if (!is.null(analysis)) {

  message(
    "Applying the user-defined analysis boundary."
  )

  common_mask <- terra::mask(
    common_mask,
    analysis
  )
}

common <- terra::mask(
  dem_stack,
  common_mask
)

names(common) <- c(
  "DEM_2011",
  "DEM_2019",
  "DEM_2023_raw"
)

valid_cell_count <- terra::global(
  common_mask,
  fun = "sum",
  na.rm = TRUE
)[1, 1]

if (!is.finite(valid_cell_count)) {
  stop(
    "The three DEMs do not have a valid common coverage."
  )
}

valid_area_m2 <- valid_cell_count * prod(
  terra::res(common_mask)
)

message(
  "Common valid cells: ",
  format(
    valid_cell_count,
    big.mark = ",",
    scientific = FALSE
  )
)

message(
  "Common valid area: ",
  round(
    valid_area_m2 / 1e6,
    4
  ),
  " km2"
)


# 11. Save uncorrected common DEM stack ----------------------------------------

uncorrected_file <- file.path(
  DIRS["rasters"],
  "Gdynia_DEM_2011_2019_2023_aligned_uncorrected.tif"
)

terra::writeRaster(
  common,
  uncorrected_file,
  overwrite = TRUE,
  gdal = c(
    "COMPRESS=DEFLATE",
    "TILED=YES"
  )
)


# 12. Function for vertical-shift diagnostics ----------------------------------

calculate_shift <- function(
  early_dem,
  later_dem,
  stable_vector,
  pair_label
) {

  message(
    "Estimating vertical offset for ",
    pair_label
  )

  difference <- later_dem - early_dem

  difference_stable <- terra::mask(
    difference,
    stable_vector
  )

  stable_values <- terra::values(
    difference_stable,
    mat = FALSE,
    na.rm = TRUE
  )

  stable_values <- stable_values[
    is.finite(stable_values)
  ]

  number_of_cells <- length(stable_values)

  if (number_of_cells < MIN_STABLE_CELLS) {
    stop(
      "Too few valid stable-area cells for ",
      pair_label,
      ". Found ",
      number_of_cells,
      "; required at least ",
      MIN_STABLE_CELLS,
      "."
    )
  }

  median_shift <- stats::median(
    stable_values,
    na.rm = TRUE
  )

  mean_shift <- base::mean(
    stable_values,
    na.rm = TRUE
  )

  mad_value <- stats::mad(
    stable_values,
    center = median_shift,
    constant = 1,
    na.rm = TRUE
  )

  rmse_value <- sqrt(
    base::mean(
      stable_values^2,
      na.rm = TRUE
    )
  )

  standard_deviation <- stats::sd(
    stable_values,
    na.rm = TRUE
  )

  q025 <- as.numeric(
    stats::quantile(
      stable_values,
      probs = 0.025,
      na.rm = TRUE,
      names = FALSE
    )
  )

  q975 <- as.numeric(
    stats::quantile(
      stable_values,
      probs = 0.975,
      na.rm = TRUE,
      names = FALSE
    )
  )

  if (
    !is.finite(median_shift) ||
    !is.finite(rmse_value)
  ) {
    stop(
      "Non-finite co-registration statistics for ",
      pair_label
    )
  }

  if (
    abs(median_shift) >
    COREG_MAX_ABS_SHIFT
  ) {
    stop(
      "Vertical offset exceeds fail-safe threshold for ",
      pair_label,
      ". Estimated median offset: ",
      round(median_shift, 4),
      " m; maximum allowed: ",
      COREG_MAX_ABS_SHIFT,
      " m. Check vertical datums and DEM signs."
    )
  }

  result <- data.table::data.table(
    pair = pair_label,
    n = number_of_cells,
    median_shift = median_shift,
    mean_shift = mean_shift,
    mad = mad_value,
    sd = standard_deviation,
    rmse = rmse_value,
    q025 = q025,
    q975 = q975,
    minimum = min(
      stable_values,
      na.rm = TRUE
    ),
    maximum = max(
      stable_values,
      na.rm = TRUE
    ),
    status = "estimated"
  )

  return(result)
}


# 13. Estimate the 2011-2019 offset --------------------------------------------

shift_2011_2019 <- calculate_shift(
  common[[1]],
  common[[2]],
  stable,
  "2011-2019_before_correction"
)

message(
  "2011-2019 median vertical offset: ",
  round(
    shift_2011_2019$median_shift,
    4
  ),
  " m"
)


# 14. Correct the 2019 DEM ------------------------------------------------------

if (COREG_APPLY_VERTICAL) {

  dem_2019_corrected <- (
    common[[2]] -
    shift_2011_2019$median_shift
  )

  names(dem_2019_corrected) <- "DEM_2019"

} else {

  dem_2019_corrected <- common[[2]]

  names(dem_2019_corrected) <- "DEM_2019"

  warning(
    "The 2019 vertical correction was not applied because ",
    "COREG_APPLY_VERTICAL is FALSE."
  )
}


# 15. Recalculate the 2011-2019 residuals --------------------------------------

residual_2011_2019 <- calculate_shift(
  common[[1]],
  dem_2019_corrected,
  stable,
  "2011-2019_after_correction"
)


# 16. Estimate the 2019-2023 offset --------------------------------------------

# The corrected 2019 DEM must be used as the reference for the 2023 correction.

shift_2019_2023 <- calculate_shift(
  dem_2019_corrected,
  common[[3]],
  stable,
  "2019-2023_before_correction"
)

message(
  "2019-2023 median vertical offset: ",
  round(
    shift_2019_2023$median_shift,
    4
  ),
  " m"
)


# 17. Correct the 2023 DEM ------------------------------------------------------

if (COREG_APPLY_VERTICAL) {

  dem_2023_corrected <- (
    common[[3]] -
    shift_2019_2023$median_shift
  )

  names(dem_2023_corrected) <- "DEM_2023"

} else {

  dem_2023_corrected <- common[[3]]

  names(dem_2023_corrected) <- "DEM_2023"

  warning(
    "The 2023 vertical correction was not applied because ",
    "COREG_APPLY_VERTICAL is FALSE."
  )
}


# 18. Recalculate the 2019-2023 residuals --------------------------------------

residual_2019_2023 <- calculate_shift(
  dem_2019_corrected,
  dem_2023_corrected,
  stable,
  "2019-2023_after_correction"
)


# 19. Combine co-registration diagnostics --------------------------------------

coregistration_diagnostics <- data.table::rbindlist(
  list(
    shift_2011_2019,
    residual_2011_2019,
    shift_2019_2023,
    residual_2019_2023
  ),
  use.names = TRUE,
  fill = TRUE
)

coregistration_file <- file.path(
  DIRS["tables"],
  "coregistration_diagnostics.csv"
)

data.table::fwrite(
  coregistration_diagnostics,
  coregistration_file
)


# 20. Create final harmonised DEM stack ----------------------------------------

harmonised <- c(
  common[[1]],
  dem_2019_corrected,
  dem_2023_corrected
)

names(harmonised) <- c(
  "DEM_2011",
  "DEM_2019",
  "DEM_2023"
)

harmonised_file <- file.path(
  DIRS["rasters"],
  "Gdynia_DEM_2011_2019_2023_harmonized.tif"
)

terra::writeRaster(
  harmonised,
  harmonised_file,
  overwrite = TRUE,
  gdal = c(
    "COMPRESS=DEFLATE",
    "TILED=YES"
  )
)


# 21. Save the common-coverage mask --------------------------------------------

common_mask_file <- file.path(
  DIRS["rasters"],
  "Gdynia_common_coverage_2011_2019_2023.tif"
)

terra::writeRaster(
  common_mask,
  common_mask_file,
  overwrite = TRUE,
  datatype = "INT1U",
  NAflag = 255,
  gdal = c(
    "COMPRESS=DEFLATE",
    "TILED=YES"
  )
)


# 22. Save vector layers in GeoPackage format ----------------------------------

dumpsite_output <- file.path(
  DIRS["vectors"],
  "dumpsite.gpkg"
)

stable_output <- file.path(
  DIRS["vectors"],
  "stable_area.gpkg"
)

terra::writeVector(
  dumpsite,
  dumpsite_output,
  overwrite = TRUE
)

terra::writeVector(
  stable,
  stable_output,
  overwrite = TRUE
)

if (!is.null(analysis)) {

  analysis_output <- file.path(
    DIRS["vectors"],
    "analysis_area.gpkg"
  )

  terra::writeVector(
    analysis,
    analysis_output,
    overwrite = TRUE
  )
}


# 23. Save a compact processing summary ----------------------------------------

processing_summary <- data.table::data.table(
  setting = c(
    "target_resolution_m",
    "output_crs",
    "common_valid_cells",
    "common_valid_area_m2",
    "common_valid_area_km2",
    "coregistration_applied",
    "stable_area_features",
    "dumpsite_features"
  ),

  value = c(
    as.character(TARGET_RES),
    MAP_CRS,
    as.character(valid_cell_count),
    as.character(valid_area_m2),
    as.character(valid_area_m2 / 1e6),
    as.character(COREG_APPLY_VERTICAL),
    as.character(nrow(stable)),
    as.character(nrow(dumpsite))
  )
)

data.table::fwrite(
  processing_summary,
  file.path(
    DIRS["tables"],
    "harmonisation_summary.csv"
  )
)


# 24. Completion messages -------------------------------------------------------

message("")
message("Harmonisation and vertical co-registration completed.")
message("Final DEM stack: ", harmonised_file)
message("Common-coverage mask: ", common_mask_file)
message("Diagnostics: ", coregistration_file)
message("Dumpsite vector: ", dumpsite_output)
message("Stable-area vector: ", stable_output)
message("")

print(coregistration_diagnostics)