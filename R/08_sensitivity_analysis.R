# ==============================================================================
# 08_sensitivity_analysis.R
# Exact cell-based sensitivity analysis for the 2019-2023 depositional footprint
# ==============================================================================

source("00_config.R", encoding = "UTF-8")

THRESHOLDS_M <- unlist(CFG$analysis$positive_thresholds_m)
MIN_AREAS_M2 <- unlist(CFG$analysis$min_patch_areas_m2)
CONNECT_DIRECTIONS <- CFG$analysis$connectivity_directions
PRIMARY_THRESHOLD_M <- CFG$analysis$primary_positive_threshold_m
PRIMARY_MIN_AREA_M2 <- CFG$analysis$primary_min_patch_area_m2

f_dod <- file.path(DIRS["rasters"], "DoD_2019_2023_raw.tif")
f_off <- DUMPSITE_VECTOR
f_buf <- BUFFER_VECTOR

required_files <- c(f_dod, f_off, f_buf)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Required files not found:\n", paste(missing_files, collapse = "\n"))
}

message("Sensitivity thresholds: ", paste(THRESHOLDS_M, collapse = ", "), " m")
message("Minimum patch areas: ", paste(MIN_AREAS_M2, collapse = ", "), " m2")
message("Connectivity: ", CONNECT_DIRECTIONS, " directions")

dod <- terra::rast(f_dod)
if (terra::nlyr(dod) != 1L) stop("DoD input must contain exactly one layer.")
if (terra::is.lonlat(dod)) stop("DoD must use a projected CRS with metre units.")

official <- terra::vect(f_off)
buffered <- terra::vect(f_buf)
if (!terra::same.crs(official, dod)) official <- terra::project(official, terra::crs(dod))
if (!terra::same.crs(buffered, dod)) buffered <- terra::project(buffered, terra::crs(dod))

# Exact, mutually exclusive zones on the analytical grid.
official_r <- terra::rasterize(official, dod, field = 1, background = NA)
buffer_r <- terra::rasterize(buffered, dod, field = 1, background = NA)
zone_r <- terra::ifel(!is.na(official_r), 1,
             terra::ifel(!is.na(buffer_r), 2, 3))
zone_r <- terra::mask(zone_r, dod)
names(zone_r) <- "zone_id"
zone_labels <- c(
  `1` = "inside_official_dumpsite",
  `2` = "outside_0_200m",
  `3` = "outside_over_200m"
)

cell_area_m2 <- abs(prod(terra::res(dod)))
if (!is.finite(cell_area_m2) || cell_area_m2 <= 0) stop("Invalid raster cell area.")

safe_global_sum <- function(x) {
  value <- terra::global(x, "sum", na.rm = TRUE)[1, 1]
  if (is.na(value) || !is.finite(value)) 0 else as.numeric(value)
}

scenario_zone <- list()
scenario_counts <- list()
scenario_crossing <- list()
row_zone <- 0L
row_count <- 0L
row_cross <- 0L

for (threshold_m in THRESHOLDS_M) {
  message("Threshold ", threshold_m, " m")

  positive <- terra::ifel(dod > threshold_m, 1, NA)
  raw_patch <- terra::patches(
    positive,
    directions = CONNECT_DIRECTIONS,
    zeroAsNA = TRUE
  )

  patch_frequency <- terra::freq(raw_patch)
  if (is.null(patch_frequency) || !nrow(patch_frequency)) {
    warning("No positive patches at threshold ", threshold_m, " m.")
    next
  }
  # terra::freq() normally returns layer, value and count. Select by name
  # to remain robust across terra versions.
  value_col <- intersect(c("value", "ID", "id"), names(patch_frequency))[1]
  count_col <- intersect(c("count", "frequency", "freq"), names(patch_frequency))[1]
  if (is.na(value_col) || is.na(count_col)) {
    stop("Unexpected terra::freq() columns: ", paste(names(patch_frequency), collapse = ", "))
  }
  patch_frequency <- data.frame(
    patch_id = patch_frequency[[value_col]],
    cell_count = patch_frequency[[count_col]]
  )
  patch_frequency$area_m2 <- patch_frequency$cell_count * cell_area_m2
  all_patch_count <- nrow(patch_frequency)

  for (min_area_m2 in MIN_AREAS_M2) {
    keep_ids <- patch_frequency$patch_id[patch_frequency$area_m2 >= min_area_m2]
    keep_ids <- keep_ids[is.finite(keep_ids)]

    # Do not use `%in%` on a SpatRaster. In current terra versions the raster
    # method interprets the right-hand side as a raster matching request and
    # can fail with "[%in%] no matches supplied". `subst()` safely retains the
    # selected patch IDs and assigns NA to all remaining cells.
    if (!length(keep_ids)) {
      retained_patch <- terra::ifel(!is.na(raw_patch), NA, NA)
    } else if (length(keep_ids) == all_patch_count) {
      retained_patch <- raw_patch
    } else {
      retained_patch <- terra::subst(
        raw_patch,
        from = keep_ids,
        to = keep_ids,
        others = NA
      )
    }
    names(retained_patch) <- "patch_id"
    retained_mask <- !is.na(retained_patch)

    row_count <- row_count + 1L
    scenario_counts[[row_count]] <- data.table::data.table(
      threshold_m = threshold_m,
      min_patch_area_m2 = min_area_m2,
      all_patch_count = all_patch_count,
      retained_patch_count = length(keep_ids),
      removed_patch_count = all_patch_count - length(keep_ids)
    )

    for (zone_id in 1:3) {
      zone_mask <- retained_mask & zone_r == zone_id
      area_m2 <- safe_global_sum(terra::ifel(zone_mask, cell_area_m2, 0))
      positive_volume_m3 <- safe_global_sum(
        terra::ifel(zone_mask, dod * cell_area_m2, 0)
      )

      row_zone <- row_zone + 1L
      scenario_zone[[row_zone]] <- data.table::data.table(
        threshold_m = threshold_m,
        min_patch_area_m2 = min_area_m2,
        zone = unname(zone_labels[as.character(zone_id)]),
        retained_area_m2 = area_m2,
        positive_volume_m3 = positive_volume_m3
      )
    }

    # Identify patches occupying more than one management zone.
    if (length(keep_ids)) {
      cross <- try(terra::crosstab(c(retained_patch, zone_r), long = TRUE), silent = TRUE)
      if (!inherits(cross, "try-error") && !is.null(cross) && nrow(cross)) {
        cross_dt <- data.table::as.data.table(cross)
        data.table::setnames(cross_dt, names(cross_dt)[1:3], c("patch_id", "zone_id", "cell_count"))
        cross_dt <- cross_dt[!is.na(patch_id) & !is.na(zone_id) & cell_count > 0]
        crossing_ids <- cross_dt[, .(zone_count = data.table::uniqueN(zone_id)), by = patch_id][zone_count > 1, patch_id]
      } else {
        crossing_ids <- numeric()
      }
    } else {
      crossing_ids <- numeric()
    }

    row_cross <- row_cross + 1L
    scenario_crossing[[row_cross]] <- data.table::data.table(
      threshold_m = threshold_m,
      min_patch_area_m2 = min_area_m2,
      boundary_crossing_patch_count = length(crossing_ids),
      boundary_crossing_patch_ids = if (length(crossing_ids)) paste(crossing_ids, collapse = ";") else ""
    )
  }
}

zone_sensitivity <- data.table::rbindlist(scenario_zone, fill = TRUE)
patch_counts <- data.table::rbindlist(scenario_counts, fill = TRUE)
boundary_crossing <- data.table::rbindlist(scenario_crossing, fill = TRUE)

if (!nrow(zone_sensitivity)) stop("Sensitivity analysis produced no zone results.")

# Proportional allocation within every threshold/area scenario.
zone_sensitivity[, scenario_area_m2 := sum(retained_area_m2),
                 by = .(threshold_m, min_patch_area_m2)]
zone_sensitivity[, scenario_volume_m3 := sum(positive_volume_m3),
                 by = .(threshold_m, min_patch_area_m2)]
zone_sensitivity[, pct_scenario_area :=
                   data.table::fifelse(scenario_area_m2 > 0,
                                      100 * retained_area_m2 / scenario_area_m2, NA_real_)]
zone_sensitivity[, pct_scenario_volume :=
                   data.table::fifelse(scenario_volume_m3 > 0,
                                      100 * positive_volume_m3 / scenario_volume_m3, NA_real_)]

# Filter-retention table: compare each filtered scenario against min area = 0
# at the same DoD threshold. This is the missing `coverage` object in the
# truncated repository script.
totals <- zone_sensitivity[, .(
  retained_area_m2 = sum(retained_area_m2),
  retained_volume_m3 = sum(positive_volume_m3)
), by = .(threshold_m, min_patch_area_m2)]

baseline <- totals[min_patch_area_m2 == 0, .(
  threshold_m,
  all_area_m2 = retained_area_m2,
  all_volume_m3 = retained_volume_m3
)]

coverage <- merge(totals, baseline, by = "threshold_m", all.x = TRUE)
coverage[, retained_area_pct :=
           data.table::fifelse(all_area_m2 > 0,
                              100 * retained_area_m2 / all_area_m2, NA_real_)]
coverage[, retained_volume_pct :=
           data.table::fifelse(all_volume_m3 > 0,
                              100 * retained_volume_m3 / all_volume_m3, NA_real_)]
data.table::setorder(coverage, threshold_m, min_patch_area_m2)

primary <- zone_sensitivity[
  abs(threshold_m - PRIMARY_THRESHOLD_M) < 1e-9 &
    min_patch_area_m2 == PRIMARY_MIN_AREA_M2
]
if (!nrow(primary)) stop("Primary sensitivity scenario was not produced.")

# Write complete analytical outputs.
data.table::fwrite(
  zone_sensitivity,
  file.path(DIRS["tables"], "boundary_sensitivity_zone_summary.csv")
)
data.table::fwrite(
  patch_counts,
  file.path(DIRS["tables"], "boundary_sensitivity_patch_counts.csv")
)
data.table::fwrite(
  boundary_crossing,
  file.path(DIRS["tables"], "boundary_crossing_sensitivity.csv")
)
data.table::fwrite(
  primary,
  file.path(DIRS["tables"], "boundary_primary_sensitivity_summary.csv")
)
data.table::fwrite(
  coverage,
  file.path(DIRS["tables"], "boundary_patch_filter_coverage.csv")
)

# Diagnostic sensitivity plots.
p1 <- ggplot2::ggplot(
  zone_sensitivity[min_patch_area_m2 == PRIMARY_MIN_AREA_M2],
  ggplot2::aes(threshold_m, positive_volume_m3 / 1e6, colour = zone)
) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::labs(
    x = "Positive DoD threshold (m)",
    y = "Positive volume (million m3)",
    colour = "Zone"
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(
  file.path(DIRS["figures"], "boundary_threshold_sensitivity_volume.png"),
  p1, width = 8, height = 5, dpi = 300
)

p2 <- ggplot2::ggplot(
  zone_sensitivity[abs(threshold_m - PRIMARY_THRESHOLD_M) < 1e-9],
  ggplot2::aes(factor(min_patch_area_m2), positive_volume_m3 / 1e6, fill = zone)
) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::labs(
    x = "Minimum patch area (m2)",
    y = "Positive volume (million m3)",
    fill = "Zone"
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(
  file.path(DIRS["figures"], "boundary_min_area_sensitivity_volume.png"),
  p2, width = 8, height = 5, dpi = 300
)

message("")
message("Step 08 sensitivity analysis completed successfully.")
message("Scenarios: ", length(THRESHOLDS_M) * length(MIN_AREAS_M2))
message("Primary threshold: ", PRIMARY_THRESHOLD_M, " m")
message("Primary minimum patch area: ", PRIMARY_MIN_AREA_M2, " m2")
print(primary)
