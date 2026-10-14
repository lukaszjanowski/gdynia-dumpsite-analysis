# ==============================================================================
# 06g_figure4_dominant_depositional_fields_v2_legend_below.R
# Figure 4: bathymetric change, 2019-2023, and two dominant depositional fields
# Main-text version: no labels or outlines for minor patches
# ==============================================================================

source("00_config.R", encoding = "UTF-8")

required <- c("terra", "sf", "ggplot2", "ggspatial", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))

library(terra)
library(sf)
library(ggplot2)
library(ggspatial)
library(scales)

# 1. Inputs --------------------------------------------------------------------
first_existing <- function(paths, label) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) stop(label, " not found. Checked:\n", paste(paths, collapse = "\n"))
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

dod_file <- first_existing(
  c(
    file.path(DIRS["rasters"], "DoD_2019_2023_raw.tif"),
    file.path(DIRS["rasters"], "DoD_2019_2023_raw.tif"),
    file.path(OUTPUT_DIR, "rasters", "DoD_2019_2023_raw.tif")
  ),
  "2019-2023 raw DoD"
)

patch_file <- first_existing(
  c(
    file.path(DIRS["vectors"], "depositional_patches_2019_2023.gpkg"),
    file.path(DIRS["vectors"], "depositional_patches_2019_2023.gpkg"),
    file.path(OUTPUT_DIR, "vectors", "depositional_patches_2019_2023.gpkg")
  ),
  "Depositional-patch vector"
)

boundary_file <- first_existing(
  c(
    file.path(DIRS["vectors"], "dumpsite.gpkg"),
    file.path(ROOT, CFG$files$dumpsite)
  ),
  "Official disposal-site boundary"
)

fig_dir <- if (exists("DIRS") && "figures" %in% names(DIRS)) {
  unname(DIRS["figures"])
} else {
  file.path(OUTPUT_DIR, "figures")
}
out_dir <- fig_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("DoD: ", dod_file)
message("Patches: ", patch_file)
message("Boundary: ", boundary_file)

# 2. Read data -----------------------------------------------------------------
dod <- rast(dod_file)
if (nlyr(dod) != 1L) stop("The DoD input must contain exactly one raster layer.")
names(dod) <- "2019-2023"

patches <- st_read(patch_file, quiet = TRUE)
boundary <- st_read(boundary_file, quiet = TRUE)
if (!nrow(patches)) stop("The depositional-patch vector is empty.")
if (is.na(st_crs(patches)) || is.na(st_crs(boundary))) {
  stop("Patch or boundary vector has no defined CRS.")
}

patches <- st_make_valid(st_transform(patches, crs(dod)))
boundary <- st_make_valid(st_transform(boundary, crs(dod)))

# 3. Select exactly the two dominant fields ------------------------------------
# Primary selection follows the existing dominant_field attribute. If an older
# file lacks that attribute, use the two largest polygons as a transparent and
# reproducible fallback.
logical_true <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

if ("dominant_field" %in% names(patches)) {
  dominant <- patches[logical_true(patches$dominant_field), ]
} else {
  dominant <- patches[0, ]
}

if (nrow(dominant) != 2L) {
  warning(
    "Expected two dominant_field polygons but found ", nrow(dominant),
    ". Falling back to the two largest patch geometries."
  )
  patch_area <- as.numeric(st_area(patches))
  dominant <- patches[order(patch_area, decreasing = TRUE)[seq_len(min(2L, nrow(patches)))], ]
}

if (nrow(dominant) != 2L) stop("Could not identify exactly two dominant depositional fields.")

# Dissolve multipart pieces within each dominant-field identity if possible.
id_field <- intersect(c("dominant_id", "field_id", "patch_id", "id"), names(dominant))[1]
if (!is.na(id_field) && length(unique(dominant[[id_field]])) == 2L) {
  dominant <- aggregate(dominant, by = list(dominant[[id_field]]), FUN = function(x) x[1])
}

dominant$map_class <- "Two dominant depositional fields"
boundary$map_class <- "Official disposal-site boundary"

message("Dominant fields selected: ", nrow(dominant))
message("Dominant total area: ", round(sum(as.numeric(st_area(dominant))) / 1e6, 3), " km2")

# 4. Common sign and colour convention -----------------------------------------
# DoD = later DEM minus earlier DEM. Positive values are accumulation/raising and
# are mapped in red. Negative values represent lowering and are mapped in blue.
abs_q <- as.numeric(global(dod, function(x, ...) {
  quantile(abs(x), 0.995, na.rm = TRUE, names = FALSE)
}, na.rm = TRUE)[1, 1])
DISPLAY_LIMIT_M <- max(3.0, ceiling(abs_q * 2) / 2)
message("Symmetric display range: -", DISPLAY_LIMIT_M, " to +", DISPLAY_LIMIT_M, " m")

break_values <- c(
  -DISPLAY_LIMIT_M, -2, -1, -0.5, -0.2, 0,
   0.2, 0.5, 1, 2, DISPLAY_LIMIT_M
)
break_values <- sort(unique(pmax(-DISPLAY_LIMIT_M, pmin(DISPLAY_LIMIT_M, break_values))))
base_colours <- c(
  "#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
  "#F7F7F7",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"
)
map_colours <- colorRampPalette(base_colours)(length(break_values))
map_positions <- rescale(break_values, from = c(-DISPLAY_LIMIT_M, DISPLAY_LIMIT_M))

# 5. Display data ---------------------------------------------------------------
TARGET_DISPLAY_CELLS <- 1400000
agg_factor <- max(1L, ceiling(sqrt(ncell(dod) / TARGET_DISPLAY_CELLS)))
message("Display aggregation factor: ", agg_factor)

dod_display <- if (agg_factor > 1L) {
  aggregate(dod, fact = agg_factor, fun = mean, na.rm = TRUE)
} else {
  dod
}

# Transform display copies to longitude/latitude to match Figures 1-3.
dod_ll <- project(dod_display, "EPSG:4326", method = "bilinear")
boundary_ll <- st_transform(boundary, 4326)
dominant_ll <- st_transform(dominant, 4326)

map_extent <- ext(dod_ll)
x_limits <- c(map_extent[1], map_extent[2])
y_limits <- c(map_extent[3], map_extent[4])

dod_df <- as.data.frame(dod_ll, xy = TRUE, na.rm = TRUE)
names(dod_df)[3] <- "change_m"

# 6. Plot ----------------------------------------------------------------------
lon_lab <- function(v) paste0(formatC(v, format = "f", digits = 2), "°E")
lat_lab <- function(v) paste0(formatC(v, format = "f", digits = 3), "°N")

legend_ticks <- pretty(c(-DISPLAY_LIMIT_M, DISPLAY_LIMIT_M), n = 7)
legend_ticks <- legend_ticks[legend_ticks >= -DISPLAY_LIMIT_M &
                               legend_ticks <= DISPLAY_LIMIT_M]
if (!0 %in% legend_ticks) legend_ticks <- sort(c(legend_ticks, 0))

FIELD_COLOUR <- "#00796B"
BOUNDARY_COLOUR <- "#111111"

p <- ggplot() +
  geom_raster(
    data = dod_df,
    aes(x = x, y = y, fill = change_m),
    interpolate = FALSE
  ) +
  # Only the two dominant fields are shown; minor patch outlines and labels are
  # deliberately omitted from the main-text figure.
  geom_sf(
    data = dominant_ll,
    aes(colour = map_class),
    fill = scales::alpha(FIELD_COLOUR, 0.09),
    linewidth = 0.85,
    inherit.aes = FALSE
  ) +
  geom_sf(
    data = boundary_ll,
    aes(colour = map_class),
    fill = NA,
    linewidth = 0.80,
    inherit.aes = FALSE
  ) +
  scale_fill_gradientn(
    colours = map_colours,
    values = map_positions,
    limits = c(-DISPLAY_LIMIT_M, DISPLAY_LIMIT_M),
    oob = squish,
    na.value = NA,
    name = "Elevation change (m)",
    breaks = legend_ticks,
    labels = function(x) sprintf("%+.1f", x),
    guide = guide_colourbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(62, "mm"),
      barheight = grid::unit(3.8, "mm"),
      ticks.colour = "#333333",
      frame.colour = "#555555",
      order = 1
    )
  ) +
  scale_colour_manual(
    name = NULL,
    breaks = c("Official disposal-site boundary", "Two dominant depositional fields"),
    values = c(
      "Official disposal-site boundary" = BOUNDARY_COLOUR,
      "Two dominant depositional fields" = FIELD_COLOUR
    ),
    guide = guide_legend(
      direction = "horizontal",
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        fill = c(NA, scales::alpha(FIELD_COLOUR, 0.09)),
        linewidth = c(0.9, 0.9),
        linetype = c("solid", "solid")
      ),
      order = 2
    )
  ) +
  scale_x_continuous(
    breaks = pretty(x_limits, n = 5),
    labels = lon_lab,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = pretty(y_limits, n = 5),
    labels = lat_lab,
    expand = c(0, 0)
  ) +
  coord_sf(
    xlim = x_limits,
    ylim = y_limits,
    crs = st_crs(4326),
    datum = st_crs(4326),
    expand = FALSE
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.18,
    text_cex = 0.55,
    line_width = 0.42,
    unit_category = "metric",
    pad_x = grid::unit(3, "mm"),
    pad_y = grid::unit(3, "mm")
  ) +
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = north_arrow_orienteering(text_size = 5.0),
    height = grid::unit(8, "mm"),
    width = grid::unit(8, "mm"),
    pad_x = grid::unit(3, "mm"),
    pad_y = grid::unit(3, "mm")
  ) +
  labs(
    title = "Bathymetric change, 2019-2023",
    subtitle = "Official boundary and two dominant depositional fields",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 9, base_family = "sans") +
  theme(
    panel.background = element_rect(fill = "#F4F9FC", colour = NA),
    panel.grid.major = element_line(colour = "#D8E3E8", linewidth = 0.22),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 8.7),
    axis.text = element_text(size = 7.3, colour = "#333333"),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 9.0, colour = "#444444", hjust = 0),
    plot.margin = margin(4, 5, 3, 5),
    plot.background = element_rect(fill = "white", colour = NA),
    # Both legends are outside the analytical map area.
    legend.position = "bottom",
    legend.justification = "center",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.title = element_text(face = "bold", size = 8.0, hjust = 0.5),
    legend.text = element_text(size = 6.8),
    legend.key.width = grid::unit(7, "mm"),
    legend.key.height = grid::unit(3.5, "mm"),
    legend.margin = margin(1, 2, 0, 2),
    legend.spacing.y = grid::unit(0.8, "mm"),
    legend.box.spacing = grid::unit(0.5, "mm")
  )

# 7. Export --------------------------------------------------------------------
output_stem <- file.path(
  out_dir,
  "Figure_4_bathymetric_change_2019_2023_dominant_fields_v2"
)

save_common <- list(
  plot = p,
  width = 180,
  height = 132,
  units = "mm",
  bg = "white",
  limitsize = FALSE
)

do.call(ggsave, c(list(filename = paste0(output_stem, ".png"), dpi = 300), save_common))
do.call(ggsave, c(list(
  filename = paste0(output_stem, ".tif"), dpi = 600, compression = "lzw"
), save_common))
do.call(ggsave, c(list(
  filename = paste0(output_stem, ".pdf"), device = grDevices::cairo_pdf
), save_common))

# 8. Reproducibility ------------------------------------------------------------
dominant_summary <- data.frame(
  field = seq_len(nrow(dominant)),
  area_m2 = as.numeric(st_area(dominant))
)
write.csv(
  dominant_summary,
  paste0(output_stem, "_dominant_fields.csv"),
  row.names = FALSE
)

record <- c(
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("DoD: ", dod_file),
  paste0("Patches: ", patch_file),
  paste0("Boundary: ", boundary_file),
  "Displayed survey interval: 2019-2023",
  "DoD convention: later DEM minus earlier DEM",
  "Positive values: accumulation/seabed raising; mapped in red",
  "Negative values: seabed lowering; mapped in blue",
  "Only two dominant depositional fields are shown; minor patches are omitted.",
  paste0("Common symmetric display limits: +/-", DISPLAY_LIMIT_M, " m"),
  paste0("Display aggregation factor: ", agg_factor),
  "Raster aggregation was used only for figure rendering."
)
writeLines(record, paste0(output_stem, "_reproducibility.txt"), useBytes = TRUE)

outputs <- c(
  paste0(output_stem, ".png"),
  paste0(output_stem, ".tif"),
  paste0(output_stem, ".pdf"),
  paste0(output_stem, "_dominant_fields.csv"),
  paste0(output_stem, "_reproducibility.txt")
)
check <- data.frame(
  file = outputs,
  exists = file.exists(outputs),
  size_mb = round(file.info(outputs)$size / 1024^2, 3)
)
print(check)
if (!all(check$exists)) stop("At least one Figure 4 output was not created.")

message("")
message("Figure 4 v2 completed successfully.")
message("Minor patch outlines and labels were omitted.")
message("PNG: ", outputs[1])
message("TIFF: ", outputs[2])
message("PDF: ", outputs[3])
message("")
