# ==============================================================================
# STEP 04b
# Object-based analysis of depositional patches, 2019-2023
#
# The script:
#   1. segments positive DoD patches using DoD > 0.35 m;
#   2. removes patches smaller than 500 m2;
#   3. classifies patches as:
#        - inside_official_dumpsite,
#        - outside_0_200m,
#        - outside_over_200m;
#   4. calculates patch area, volume, thickness, morphometry and relative BBS;
#   5. saves CSV, GeoPackage, labelled raster and diagnostic figures.
#
# Required existing products:
#   rasters/DoD_2019_2023_raw.tif
#   rasters/Gdynia_2023_morphometry.tif
#   rasters/BBS_2023_aligned_interpretive.tif
#   SHP/DUMPSITE_VECTOR.shp
#   outputs/vectors/dumpsite_buffer_200m.gpkg
# ============================================================================== 

source("00_config.R", encoding = "UTF-8")

# 1. Settings ------------------------------------------------------------------
PATCH_THRESHOLD_M <- 0.35
MIN_PATCH_AREA_M2 <- 500
LARGE_PATCH_AREA_M2 <- 1000
CONNECT_DIRECTIONS <- 8

buffer_file <- BUFFER_VECTOR

dod_file <- file.path(
  DIRS["rasters"],
  "DoD_2019_2023_raw.tif"
)

morphometry_file <- file.path(
  DIRS["rasters"],
  "Gdynia_2023_morphometry.tif"
)

bbs_file <- file.path(
  DIRS["rasters"],
  "BBS_2023_aligned_interpretive.tif"
)

required_files <- c(
  dod_file,
  morphometry_file,
  bbs_file,
  DUMPSITE_VECTOR,
  buffer_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Required files not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

message("Patch threshold: DoD > ", PATCH_THRESHOLD_M, " m")
message("Minimum patch area: ", MIN_PATCH_AREA_M2, " m2")

# 2. Read rasters ---------------------------------------------------------------
dod <- terra::rast(dod_file)[[1]]
names(dod) <- "DoD_2019_2023"

morphometry <- terra::rast(morphometry_file)

expected_morphometry <- c(
  "DEM_2023",
  "slope_deg",
  "TRI",
  "roughness",
  "BPI_11m",
  "BPI_51m"
)

if (!all(expected_morphometry %in% names(morphometry))) {
  stop(
    "Missing morphometry layers. Expected: ",
    paste(expected_morphometry, collapse = ", ")
  )
}

bbs <- terra::rast(bbs_file)[[1]]
names(bbs) <- "BBS_relative"

# Align auxiliary rasters only if necessary.
align_to_dod <- function(x, method = "bilinear") {
  if (!terra::same.crs(x, dod)) {
    x <- terra::project(x, terra::crs(dod), method = method)
  }
  if (!terra::compareGeom(x, dod, stopOnError = FALSE)) {
    x <- terra::resample(x, dod, method = method)
  }
  x
}

morphometry <- align_to_dod(morphometry, "bilinear")
bbs <- align_to_dod(bbs, "bilinear")

cell_area_m2 <- abs(prod(terra::res(dod)))
min_patch_cells <- ceiling(MIN_PATCH_AREA_M2 / cell_area_m2)

message("Raster cell area: ", cell_area_m2, " m2")
message("Minimum patch size: ", min_patch_cells, " cells")

# 3. Read and prepare vectors ---------------------------------------------------
read_polygon <- function(file_path, target, layer_name) {
  v <- terra::vect(file_path)
  if (nrow(v) < 1) stop(layer_name, " contains no features.")
  if (is.na(terra::crs(v))) stop(layer_name, " has no defined CRS.")
  if (!terra::same.crs(v, target)) {
    v <- terra::project(v, terra::crs(target))
  }
  v <- terra::makeValid(v)
  terra::aggregate(v)
}

official <- read_polygon(
  DUMPSITE_VECTOR,
  dod,
  "Official dumpsite boundary"
)

buffered <- read_polygon(
  buffer_file,
  dod,
  "Buffered dumpsite boundary"
)

# External 0-200 m ring.
ring_0_200m <- terra::erase(buffered, official)
ring_0_200m <- terra::makeValid(ring_0_200m)

if (nrow(ring_0_200m) < 1) {
  stop("The 0-200 m external ring could not be created.")
}

terra::writeVector(
  ring_0_200m,
  file.path(DIRS["vectors"], "outside_ring_0_200m.gpkg"),
  overwrite = TRUE
)

# Raster zone map: 1 inside, 2 external ring, 3 farther outside.
zone_raster <- terra::ifel(!is.na(dod), 3, NA)
zone_raster <- terra::rasterize(
  ring_0_200m,
  zone_raster,
  field = 2,
  background = zone_raster,
  touches = TRUE,
  update = TRUE
)
zone_raster <- terra::rasterize(
  official,
  zone_raster,
  field = 1,
  background = zone_raster,
  touches = TRUE,
  update = TRUE
)
names(zone_raster) <- "zone_code"

terra::writeRaster(
  zone_raster,
  file.path(DIRS["rasters"], "spatial_zones_official_0_200m_beyond.tif"),
  overwrite = TRUE,
  datatype = "INT1U",
  NAflag = 255,
  gdal = c("COMPRESS=DEFLATE", "TILED=YES")
)

# 4. Segment positive changes --------------------------------------------------
positive_mask <- terra::ifel(dod > PATCH_THRESHOLD_M, 1, NA)
names(positive_mask) <- "positive_change"

patch_id_raw <- terra::patches(
  positive_mask,
  directions = CONNECT_DIRECTIONS,
  zeroAsNA = TRUE
)
names(patch_id_raw) <- "patch_id_raw"

# Cell counts for each raw patch.
raw_frequency <- data.table::as.data.table(
  terra::freq(patch_id_raw, digits = 0)
)

if (nrow(raw_frequency) == 0) {
  stop("No positive patches found above the selected threshold.")
}

# terra::freq() normally returns: layer, value, count.
# Use column names explicitly; do not assume that the first two columns are
# patch ID and cell count.
if (!all(c("value", "count") %in% names(raw_frequency))) {
  stop(
    "Unexpected terra::freq() output columns: ",
    paste(names(raw_frequency), collapse = ", ")
  )
}

raw_frequency <- raw_frequency[, .(
  patch_id = as.integer(value),
  cell_count = as.numeric(count)
)]

raw_frequency[, area_m2 := cell_count * cell_area_m2]

keep_ids <- raw_frequency[
  area_m2 >= MIN_PATCH_AREA_M2,
  patch_id
]

if (length(keep_ids) == 0) {
  stop("No patches meet the minimum area criterion.")
}

message("Raw patches: ", nrow(raw_frequency))
message("Patches >= ", MIN_PATCH_AREA_M2, " m2: ", length(keep_ids))

# Keep original IDs, then renumber consecutively.
patch_id_filtered <- terra::ifel(
  patch_id_raw %in% keep_ids,
  patch_id_raw,
  NA
)

id_lookup <- data.table(
  old_id = sort(keep_ids),
  patch_id = seq_along(sort(keep_ids))
)

patch_id <- terra::classify(
  patch_id_filtered,
  as.matrix(id_lookup[, .(old_id, patch_id)]),
  others = NA
)
names(patch_id) <- "patch_id"

terra::writeRaster(
  patch_id,
  file.path(DIRS["rasters"], "depositional_patch_ids_2019_2023.tif"),
  overwrite = TRUE,
  datatype = "INT4S",
  gdal = c("COMPRESS=DEFLATE", "TILED=YES")
)

# 5. Convert patches to polygons -----------------------------------------------
patch_polygons <- terra::as.polygons(
  patch_id,
  aggregate = TRUE,
  values = TRUE,
  na.rm = TRUE
)

names(patch_polygons)[1] <- "patch_id"
patch_polygons$patch_id <- as.integer(patch_polygons$patch_id)
patch_polygons$area_m2 <- terra::expanse(patch_polygons, unit = "m")

# Centroids and coordinates.
centroids <- terra::centroids(patch_polygons, inside = TRUE)
centroid_xy <- terra::crds(centroids)
patch_polygons$centroid_x <- centroid_xy[, 1]
patch_polygons$centroid_y <- centroid_xy[, 2]

# 6. Patch-zone allocation ------------------------------------------------------
# Count cells of each zone in each patch. Majority-cell allocation is used.
zone_extract <- terra::extract(
  zone_raster,
  patch_polygons,
  cells = TRUE,
  touches = TRUE
)

zone_dt <- as.data.table(zone_extract)
setnames(zone_dt, old = names(zone_dt)[1], new = "polygon_row")
zone_dt <- zone_dt[!is.na(zone_code)]
zone_dt[, patch_id := patch_polygons$patch_id[polygon_row]]

zone_counts <- zone_dt[, .N, by = .(patch_id, zone_code)]
setorder(zone_counts, patch_id, -N)
zone_majority <- zone_counts[, .SD[1], by = patch_id]

zone_labels <- c(
  "1" = "inside_official_dumpsite",
  "2" = "outside_0_200m",
  "3" = "outside_over_200m"
)

zone_majority[, zone := unname(zone_labels[as.character(zone_code)])]

# Also calculate percentage of patch cells in each zone.
zone_wide <- data.table::dcast(
  zone_counts,
  patch_id ~ zone_code,
  value.var = "N",
  fill = 0
)

for (z in c("1", "2", "3")) {
  if (!z %in% names(zone_wide)) zone_wide[, (z) := 0]
}

setnames(
  zone_wide,
  c("1", "2", "3"),
  c("cells_inside", "cells_ring_0_200m", "cells_beyond_200m")
)

zone_wide[, cells_total := cells_inside + cells_ring_0_200m + cells_beyond_200m]
zone_wide[, pct_inside := 100 * cells_inside / cells_total]
zone_wide[, pct_ring_0_200m := 100 * cells_ring_0_200m / cells_total]
zone_wide[, pct_beyond_200m := 100 * cells_beyond_200m / cells_total]

# Boundary-crossing flag.
zone_wide[, crosses_official_boundary :=
  cells_inside > 0 & (cells_ring_0_200m + cells_beyond_200m) > 0]

# 7. Distance from each patch to the official boundary -------------------------
# Distance from centroid to the official polygon is used for a stable summary.
# It equals zero for centroids inside the polygon.
distance_matrix <- terra::distance(centroids, official)
centroid_distance_m <- apply(distance_matrix, 1, min, na.rm = TRUE)

# Signed distance: negative inside, positive outside.
inside_centroid <- lengths(terra::relate(centroids, official, relation = "within")) > 0
signed_distance_m <- centroid_distance_m
signed_distance_m[inside_centroid] <- -signed_distance_m[inside_centroid]

# 8. Zonal statistics helper ---------------------------------------------------
extract_stats <- function(raster_layer, polygons, prefix) {
  values_list <- terra::extract(
    raster_layer,
    polygons,
    fun = NULL,
    touches = TRUE
  )

  dt <- as.data.table(values_list)
  setnames(dt, old = names(dt)[1], new = "polygon_row")
  value_col <- setdiff(names(dt), "polygon_row")[1]
  setnames(dt, value_col, "value")
  dt[, patch_id := polygons$patch_id[polygon_row]]
  dt <- dt[is.finite(value)]

  out <- dt[, .(
    n_cells = .N,
    mean = mean(value),
    sd = stats::sd(value),
    median = stats::median(value),
    q25 = as.numeric(stats::quantile(value, 0.25)),
    q75 = as.numeric(stats::quantile(value, 0.75)),
    q98 = as.numeric(stats::quantile(value, 0.98)),
    minimum = min(value),
    maximum = max(value)
  ), by = patch_id]

  setnames(
    out,
    setdiff(names(out), "patch_id"),
    paste0(prefix, "_", setdiff(names(out), "patch_id"))
  )
  out
}

# 9. Calculate patch statistics -------------------------------------------------
patch_table <- data.table(
  patch_id = patch_polygons$patch_id,
  area_m2 = patch_polygons$area_m2,
  area_ha = patch_polygons$area_m2 / 10000,
  centroid_x = patch_polygons$centroid_x,
  centroid_y = patch_polygons$centroid_y,
  centroid_signed_distance_to_boundary_m = signed_distance_m
)

patch_table <- merge(patch_table, zone_majority[, .(patch_id, zone)], by = "patch_id", all.x = TRUE)
patch_table <- merge(patch_table, zone_wide, by = "patch_id", all.x = TRUE)

# DoD statistics.
dod_stats <- extract_stats(dod, patch_polygons, "dod_m")
patch_table <- merge(patch_table, dod_stats, by = "patch_id", all.x = TRUE)

# Positive volume calculated directly from segmented cells.
patch_table[, positive_volume_m3 := dod_m_mean * area_m2]

# Morphometry statistics.
for (layer_name in expected_morphometry) {
  stats_i <- extract_stats(
    morphometry[[layer_name]],
    patch_polygons,
    layer_name
  )
  patch_table <- merge(patch_table, stats_i, by = "patch_id", all.x = TRUE)
}

# Relative BBS statistics.
bbs_stats <- extract_stats(bbs, patch_polygons, "bbs_relative")
patch_table <- merge(patch_table, bbs_stats, by = "patch_id", all.x = TRUE)
patch_table[, bbs_relative_iqr := bbs_relative_q75 - bbs_relative_q25]

# Patch-size class and ranking.
patch_table[, size_class := data.table::fifelse(
  area_m2 >= LARGE_PATCH_AREA_M2,
  paste0("at_least_", LARGE_PATCH_AREA_M2, "_m2"),
  paste0(MIN_PATCH_AREA_M2, "_to_", LARGE_PATCH_AREA_M2 - 1, "_m2")
)]

setorder(patch_table, -positive_volume_m3)
patch_table[, volume_rank := seq_len(.N)]

# Dominant fields: two largest patches inside or crossing the official boundary.
dominant_ids <- patch_table[
  zone == "inside_official_dumpsite" | crosses_official_boundary == TRUE,
  head(patch_id, 2)
]
patch_table[, dominant_field := patch_id %in% dominant_ids]

# 10. Add attributes to polygons ------------------------------------------------
poly_attrs <- as.data.frame(
  patch_table[, .(
    patch_id,
    zone,
    area_m2,
    area_ha,
    positive_volume_m3,
    dod_m_mean,
    dod_m_maximum,
    centroid_signed_distance_to_boundary_m,
    crosses_official_boundary,
    dominant_field,
    volume_rank
  )]
)

patch_polygons <- merge(
  patch_polygons,
  poly_attrs,
  by = "patch_id",
  all.x = TRUE
)

# 11. Save outputs --------------------------------------------------------------
patch_csv <- file.path(
  DIRS["tables"],
  "depositional_patches_2019_2023.csv"
)

data.table::fwrite(patch_table, patch_csv)

patch_gpkg <- file.path(
  DIRS["vectors"],
  "depositional_patches_2019_2023.gpkg"
)

terra::writeVector(
  patch_polygons,
  patch_gpkg,
  overwrite = TRUE
)

# Zone summaries.
zone_summary <- patch_table[, .(
  patch_count = .N,
  total_area_m2 = sum(area_m2, na.rm = TRUE),
  total_positive_volume_m3 = sum(positive_volume_m3, na.rm = TRUE),
  median_patch_area_m2 = median(area_m2, na.rm = TRUE),
  median_patch_thickness_m = median(dod_m_mean, na.rm = TRUE),
  max_patch_thickness_m = max(dod_m_maximum, na.rm = TRUE)
), by = zone]

zone_summary[, pct_all_patch_area := 100 * total_area_m2 / sum(total_area_m2)]
zone_summary[, pct_all_patch_volume := 100 * total_positive_volume_m3 / sum(total_positive_volume_m3)]
setorder(zone_summary, -total_positive_volume_m3)

data.table::fwrite(
  zone_summary,
  file.path(DIRS["tables"], "depositional_patch_zone_summary_2019_2023.csv")
)

# Compact table for the two dominant fields.
dominant_fields <- patch_table[
  dominant_field == TRUE,
  .(
    patch_id,
    volume_rank,
    zone,
    area_m2,
    positive_volume_m3,
    dod_mean_m = dod_m_mean,
    dod_max_m = dod_m_maximum,
    slope_median_deg = slope_deg_median,
    slope_q98_deg = slope_deg_q98,
    TRI_median = TRI_median,
    roughness_median = roughness_median,
    BPI_11m_median = BPI_11m_median,
    BPI_51m_median = BPI_51m_median,
    bbs_relative_median,
    bbs_relative_iqr,
    crosses_official_boundary
  )
]

data.table::fwrite(
  dominant_fields,
  file.path(DIRS["tables"], "dominant_depositional_fields_2019_2023.csv")
)

# 12. Diagnostic maps -----------------------------------------------------------
map_file <- file.path(
  DIRS["figures"],
  "depositional_patches_2019_2023.png"
)

png(map_file, width = 3000, height = 2400, res = 300)
par(mar = c(4, 4, 4, 6))
plot(
  dod,
  range = c(-1, 3),
  col = hcl.colors(101, "Blue-Red 3", rev = TRUE),
  main = paste0(
    "Depositional patches 2019-2023 (DoD > ",
    PATCH_THRESHOLD_M,
    " m; area >= ",
    MIN_PATCH_AREA_M2,
    " m2)"
  )
)
lines(buffered, col = "grey35", lwd = 1.5, lty = 2)
lines(official, col = "black", lwd = 2.5)
lines(patch_polygons, col = "yellow", lwd = 1.3)
text(centroids, labels = patch_polygons$patch_id, cex = 0.55, col = "black")
legend(
  "topright",
  legend = c("official boundary", "200 m buffer", "retained patch"),
  col = c("black", "grey35", "yellow"),
  lwd = c(2.5, 1.5, 1.3),
  lty = c(1, 2, 1),
  bty = "n"
)
dev.off()

bbs_map_file <- file.path(
  DIRS["figures"],
  "BBS_2023_with_depositional_patches.png"
)

# terra::quantile() expects a multi-layer raster in some versions.
# Calculate plotting limits for the single BBS layer with terra::global().
bbs_limits <- as.numeric(
  terra::global(
    bbs,
    fun = function(x, ...) {
      stats::quantile(
        x,
        probs = c(0.02, 0.98),
        na.rm = TRUE,
        names = FALSE
      )
    },
    na.rm = TRUE
  )[1, ]
)

if (
  length(bbs_limits) != 2 ||
  any(!is.finite(bbs_limits)) ||
  bbs_limits[1] >= bbs_limits[2]
) {
  stop("Could not calculate valid BBS plotting limits.")
}

png(bbs_map_file, width = 3000, height = 2400, res = 300)
par(mar = c(4, 4, 4, 6))
plot(
  bbs,
  range = bbs_limits,
  col = hcl.colors(100, "Grays"),
  main = "Relative MBES backscatter and depositional patches, 2023"
)
lines(official, col = "red", lwd = 2.5)
lines(patch_polygons, col = "cyan", lwd = 1.4)
text(centroids, labels = patch_polygons$patch_id, cex = 0.55, col = "yellow")
dev.off()

# 13. Processing summary --------------------------------------------------------
summary_lines <- c(
  paste0("Threshold: DoD > ", PATCH_THRESHOLD_M, " m"),
  paste0("Minimum patch area: ", MIN_PATCH_AREA_M2, " m2"),
  paste0("Connectivity: ", CONNECT_DIRECTIONS, " directions"),
  paste0("Retained patches: ", nrow(patch_table)),
  paste0("Total retained area: ", round(sum(patch_table$area_m2) / 1e6, 4), " km2"),
  paste0("Total retained positive volume: ", round(sum(patch_table$positive_volume_m3), 0), " m3"),
  "",
  "Patch counts by zone:",
  capture.output(print(zone_summary)),
  "",
  "Two dominant depositional fields:",
  capture.output(print(dominant_fields))
)

writeLines(
  summary_lines,
  file.path(DIRS["tables"], "step_04b_summary.txt"),
  useBytes = TRUE
)

message("")
message("Step 04b completed successfully.")
message("Patch table: ", patch_csv)
message("Patch polygons: ", patch_gpkg)
message("Zone summary: depositional_patch_zone_summary_2019_2023.csv")
message("Dominant fields: dominant_depositional_fields_2019_2023.csv")
message("Map: ", map_file)
message("BBS map: ", bbs_map_file)
message("")
print(zone_summary)
print(dominant_fields)
