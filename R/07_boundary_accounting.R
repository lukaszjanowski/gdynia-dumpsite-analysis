# ==============================================================================
# STEP 04c
# Exact cell-based allocation of depositional patches to spatial zones
# ==============================================================================

source(
  "00_config.R",
  encoding = "UTF-8"
)

dod <- terra::rast(
  file.path(
    DIRS["rasters"],
    "DoD_2019_2024_raw.tif"
  )
)[[1]]

patch_id <- terra::rast(
  file.path(
    DIRS["rasters"],
    "depositional_patch_ids_2019_2024.tif"
  )
)[[1]]

zone_raster <- terra::rast(
  file.path(
    DIRS["rasters"],
    "spatial_zones_official_0_200m_beyond.tif"
  )
)[[1]]

if (!terra::compareGeom(
  dod,
  patch_id,
  zone_raster,
  stopOnError = FALSE
)) {
  stop(
    "DoD, patch ID and zone rasters do not have matching geometry."
  )
}

cell_area_m2 <- abs(
  prod(
    terra::res(dod)
  )
)

# Extract only cells belonging to retained depositional patches
cell_stack <- c(
  patch_id,
  zone_raster,
  dod
)

names(cell_stack) <- c(
  "patch_id",
  "zone_code",
  "dod_m"
)

cell_data <- terra::as.data.frame(
  cell_stack,
  na.rm = TRUE
)

cell_data <- data.table::as.data.table(
  cell_data
)

cell_data <- cell_data[
  is.finite(patch_id) &
  is.finite(zone_code) &
  is.finite(dod_m) &
  dod_m > 0.35
]

cell_data[, patch_id := as.integer(patch_id)]
cell_data[, zone_code := as.integer(zone_code)]

zone_labels <- c(
  "1" = "inside_official_dumpsite",
  "2" = "outside_0_200m",
  "3" = "outside_over_200m"
)

cell_data[
  ,
  zone := unname(
    zone_labels[
      as.character(zone_code)
    ]
  )
]

# Exact patch-by-zone accounting
patch_zone_exact <- cell_data[
  ,
  .(
    cell_count = .N,
    area_m2 = .N * cell_area_m2,
    positive_volume_m3 = sum(
      dod_m * cell_area_m2,
      na.rm = TRUE
    ),
    mean_thickness_m = mean(
      dod_m,
      na.rm = TRUE
    ),
    median_thickness_m = median(
      dod_m,
      na.rm = TRUE
    ),
    maximum_thickness_m = max(
      dod_m,
      na.rm = TRUE
    )
  ),
  by = .(
    patch_id,
    zone
  )
]

data.table::setorder(
  patch_zone_exact,
  patch_id,
  zone
)

data.table::fwrite(
  patch_zone_exact,
  file.path(
    DIRS["tables"],
    "depositional_patch_zone_exact_2019_2024.csv"
  )
)

# Exact totals by zone
zone_exact <- patch_zone_exact[
  ,
  .(
    patch_count = data.table::uniqueN(
      patch_id
    ),
    total_cells = sum(
      cell_count
    ),
    total_area_m2 = sum(
      area_m2
    ),
    total_positive_volume_m3 = sum(
      positive_volume_m3
    ),
    mean_thickness_weighted_m =
      sum(positive_volume_m3) /
      sum(area_m2),
    maximum_thickness_m = max(
      maximum_thickness_m
    )
  ),
  by = zone
]

zone_exact[
  ,
  pct_total_area :=
    100 *
    total_area_m2 /
    sum(total_area_m2)
]

zone_exact[
  ,
  pct_total_positive_volume :=
    100 *
    total_positive_volume_m3 /
    sum(total_positive_volume_m3)
]

data.table::setorder(
  zone_exact,
  -total_positive_volume_m3
)

data.table::fwrite(
  zone_exact,
  file.path(
    DIRS["tables"],
    "depositional_zone_exact_summary_2019_2024.csv"
  )
)

# Exact split for boundary-crossing patches
crossing_patches <- patch_zone_exact[
  ,
  .(
    zone_count = data.table::uniqueN(zone),
    total_area_m2 = sum(area_m2),
    total_positive_volume_m3 =
      sum(positive_volume_m3)
  ),
  by = patch_id
][
  zone_count > 1
]

crossing_details <- merge(
  patch_zone_exact,
  crossing_patches[
    ,
    .(patch_id)
  ],
  by = "patch_id"
)

data.table::fwrite(
  crossing_details,
  file.path(
    DIRS["tables"],
    "boundary_crossing_patches_exact_2019_2024.csv"
  )
)

message("")
message("Exact cell-based zone accounting completed.")
message("")

print(zone_exact)

message("")
message("Boundary-crossing patch details:")
print(crossing_details)
