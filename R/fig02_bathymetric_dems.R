# ==============================================================================
# 06e_figure2_bathymetric_DEMs.R
# Figure 2: bathymetric DEMs with illuminated relief
# Survey years: 2011, 2019 and 2023
# ==============================================================================

source("00_config.R", encoding = "UTF-8")

required <- c("terra", "sf", "ggplot2", "patchwork", "ggspatial", "scales", "cowplot")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))

library(terra)
library(sf)
library(ggplot2)
library(patchwork)
library(ggspatial)
library(scales)

# 1. Files ---------------------------------------------------------------------
first_existing <- function(paths, label) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) stop(label, " not found. Checked:\n", paste(paths, collapse = "\n"))
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

dem_file <- first_existing(
  c(
    file.path(DIRS["rasters"], "Gdynia_DEM_2011_2019_2023_harmonized.tif"),
    file.path(DIRS["rasters"], "Gdynia_DEM_2011_2019_2023_harmonized.tif")
  ),
  "Harmonised three-layer DEM"
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

message("DEM stack: ", dem_file)
message("Official boundary: ", boundary_file)

# 2. Read and harmonise display domain -----------------------------------------
dems <- rast(dem_file)
if (nlyr(dems) != 3L) stop("The DEM file must contain exactly three layers.")
names(dems) <- c("2011", "2019", "2023")

# Use only the common valid footprint so all panels have identical data coverage.
common_mask <- !is.na(dems[[1]]) & !is.na(dems[[2]]) & !is.na(dems[[3]])
dems <- mask(dems, common_mask)
dems <- trim(dems)

boundary <- st_read(boundary_file, quiet = TRUE)
if (is.na(st_crs(boundary))) stop("The official boundary has no defined CRS.")
boundary <- st_make_valid(boundary)
boundary <- st_transform(boundary, crs(dems))

# 3. Common elevation scale -----------------------------------------------------
# Calculate robust limits jointly from all three surfaces.
quantile_fun <- function(x, ...) {
  quantile(x, probs = c(0.005, 0.995), na.rm = TRUE, names = FALSE)
}
q_layer <- global(dems, quantile_fun, na.rm = TRUE)
z_limits <- range(as.matrix(q_layer), na.rm = TRUE)

message("Shared colour range: ", round(z_limits[1], 2), " to ",
        round(z_limits[2], 2), " m")

# 4. Global Mapper-like lighting parameters -----------------------------------
# Based on the supplied Global Mapper settings:
#   light altitude 45 degrees
#   light azimuth 45 degrees
#   one light source
# User-requested vertical exaggeration: 5x
VERTICAL_EXAGGERATION <- 5.0
LIGHT_ALTITUDE <- 45
LIGHT_AZIMUTH <- 45

# RGB modulation parameters. Illumination is centred around 0.5. Shadows and
# highlights are deliberately restrained so bathymetric colours remain legible.
AMBIENT_LIGHT <- 0.68
DIRECT_LIGHT <- 0.52
MIN_BRIGHTNESS <- 0.48
MAX_BRIGHTNESS <- 1.18

# 5. Rendering resolution -------------------------------------------------------
# Aggregation applies only to the figure copy, never to analytical rasters.
TARGET_DISPLAY_CELLS <- 1100000
agg_factor <- max(1L, ceiling(sqrt(ncell(dems[[1]]) / TARGET_DISPLAY_CELLS)))
message("Display aggregation factor: ", agg_factor)

if (agg_factor > 1L) {
  dems_display <- aggregate(dems, fact = agg_factor, fun = mean, na.rm = TRUE)
} else {
  dems_display <- dems
}

# 6. Hillshade -----------------------------------------------------------------
# normalize=FALSE is essential. normalize=TRUE returns approximately 0-255 and
# caused the earlier grey artefact when interpreted as 0-1 alpha.
make_hillshade <- function(dem) {
  dem_exag <- dem * VERTICAL_EXAGGERATION
  slope <- terrain(dem_exag, v = "slope", unit = "radians", neighbors = 8)
  aspect <- terrain(dem_exag, v = "aspect", unit = "radians", neighbors = 8)
  hs <- shade(
    slope,
    aspect,
    angle = LIGHT_ALTITUDE,
    direction = LIGHT_AZIMUTH,
    normalize = FALSE
  )
  hs <- clamp(hs, lower = 0, upper = 1, values = TRUE)
  hs <- mask(hs, dem)
  names(hs) <- "illumination"
  hs
}

message("Calculating 5x hillshade...")
hill_display <- rast(lapply(1:3, function(i) make_hillshade(dems_display[[i]])))
names(hill_display) <- names(dems_display)

# 7. Reproject display copies to WGS84 -----------------------------------------
template_ll <- project(dems_display[[1]], "EPSG:4326", method = "bilinear")

dems_ll <- rast(lapply(1:3, function(i) {
  project(dems_display[[i]], template_ll, method = "bilinear")
}))
hill_ll <- rast(lapply(1:3, function(i) {
  project(hill_display[[i]], template_ll, method = "bilinear")
}))
names(dems_ll) <- names(dems)
names(hill_ll) <- names(dems)

boundary_ll <- st_transform(boundary, 4326)

map_extent <- ext(dems_ll[[1]])
x_limits <- c(map_extent[1], map_extent[2])
y_limits <- c(map_extent[3], map_extent[4])

# 8. Palette and direct RGB composition ----------------------------------------
# Hillshade is baked into RGB colours. No alpha raster is drawn, so NA cells can
# never create a grey rectangular overlay.
bathy_colours <- c(
  "#062B56", "#07577A", "#087F78", "#2A9862",
  "#86A93A", "#D2B33B", "#E89A63", "#F4B6A6", "#F2C6E6"
)
colour_function <- grDevices::colorRamp(bathy_colours, space = "Lab")

compose_rgb_df <- function(dem, illumination, year_label) {
  d <- as.data.frame(c(dem, illumination), xy = TRUE, na.rm = FALSE)
  names(d)[3:4] <- c("elevation", "illumination")

  # Keep only cells that are valid in the DEM and hillshade.
  d <- d[is.finite(d$elevation) & is.finite(d$illumination), , drop = FALSE]
  if (!nrow(d)) stop("No valid DEM-hillshade cells for year ", year_label)

  message(
    "  ", year_label, ": ", format(nrow(d), big.mark = ","),
    " cells; illumination ",
    round(min(d$illumination), 3), " to ",
    round(max(d$illumination), 3)
  )

  scaled_z <- (d$elevation - z_limits[1]) / (z_limits[2] - z_limits[1])
  scaled_z <- pmin(1, pmax(0, scaled_z))

  # colorRamp() may simplify its output in edge cases. Force an n x 3 matrix.
  base_rgb_raw <- colour_function(scaled_z)
  base_rgb <- matrix(
    as.numeric(base_rgb_raw),
    nrow = length(scaled_z),
    ncol = 3,
    byrow = FALSE
  )
  colnames(base_rgb) <- c("red", "green", "blue")

  # Global Mapper-like brightness modulation.
  brightness <- AMBIENT_LIGHT + DIRECT_LIGHT * d$illumination
  brightness <- pmin(MAX_BRIGHTNESS, pmax(MIN_BRIGHTNESS, brightness))

  # Explicit matrix multiplication by row, without dimension dropping.
  shaded_rgb <- base_rgb * matrix(
    brightness,
    nrow = length(brightness),
    ncol = 3
  )

  # Mild gamma adjustment and clipping to valid RGB values.
  shaded_rgb <- 255 * (shaded_rgb / 255) ^ 0.96
  shaded_rgb <- pmin(255, pmax(0, shaded_rgb))
  shaded_rgb <- matrix(
    shaded_rgb,
    nrow = nrow(d),
    ncol = 3,
    dimnames = list(NULL, c("red", "green", "blue"))
  )

  d$rgb <- grDevices::rgb(
    red = shaded_rgb[, "red"],
    green = shaded_rgb[, "green"],
    blue = shaded_rgb[, "blue"],
    maxColorValue = 255
  )

  d
}
message("Compositing illuminated RGB bathymetry...")
rgb_df <- lapply(1:3, function(i) {
  compose_rgb_df(dems_ll[[i]], hill_ll[[i]], names(dems)[i])
})

# 9. Shared style ---------------------------------------------------------------
lon_lab <- function(v) paste0(formatC(v, format = "f", digits = 2), "°E")
lat_lab <- function(v) paste0(formatC(v, format = "f", digits = 2), "°N")

panel_theme <- theme_minimal(base_size = 9, base_family = "sans") +
  theme(
    panel.background = element_rect(fill = "#F4F9FC", colour = NA),
    panel.grid.major = element_line(colour = "#D8E3E8", linewidth = 0.22),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7.0, colour = "#333333"),
    plot.title = element_text(face = "bold", size = 10.5, hjust = 0),
    plot.margin = margin(2, 2, 2, 2)
  )

make_panel <- function(i, letter, show_y = TRUE) {
  p <- ggplot() +
    geom_raster(
      data = rgb_df[[i]],
      aes(x = x, y = y, fill = rgb),
      interpolate = TRUE
    ) +
    scale_fill_identity() +
    geom_sf(
      data = boundary_ll,
      fill = NA,
      colour = "#B2182B",
      linewidth = 0.80,
      inherit.aes = FALSE
    ) +
    scale_x_continuous(
      breaks = pretty(x_limits, n = 3),
      labels = lon_lab,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = pretty(y_limits, n = 4),
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
      width_hint = 0.22,
      text_cex = 0.52,
      line_width = 0.38,
      unit_category = "metric",
      pad_x = grid::unit(2.5, "mm"),
      pad_y = grid::unit(2.5, "mm")
    ) +
    annotation_north_arrow(
      location = "tl",
      which_north = "true",
      style = north_arrow_orienteering(text_size = 4.7),
      height = grid::unit(7.5, "mm"),
      width = grid::unit(7.5, "mm"),
      pad_x = grid::unit(2.5, "mm"),
      pad_y = grid::unit(2.5, "mm")
    ) +
    labs(
      title = paste0(letter, "  ", names(dems)[i]),
      x = NULL,
      y = NULL
    ) +
    panel_theme

  if (!show_y) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  p
}

p1 <- make_panel(1, "a", TRUE)
p2 <- make_panel(2, "b", FALSE)
p3 <- make_panel(3, "c", FALSE)

# 10. Independent shared colour bar --------------------------------------------
# Because map colours are pre-composed RGB, use a small dummy plot solely to
# produce a correct unshaded elevation legend.
legend_data <- data.frame(
  x = seq(z_limits[1], z_limits[2], length.out = 300),
  y = 1
)

legend_plot <- ggplot(legend_data, aes(x = x, y = y, fill = x)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = bathy_colours,
    limits = z_limits,
    name = "Seabed elevation (m)",
    guide = guide_colourbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(58, "mm"),
      barheight = grid::unit(3.4, "mm"),
      ticks.colour = "#333333",
      frame.colour = "#555555"
    )
  ) +
  guides(fill = guide_colourbar()) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 8.5),
    legend.text = element_text(size = 7.2),
    plot.margin = margin(0, 0, 0, 0)
  )

# Retain only the legend area.
legend_only <- patchwork::wrap_elements(full = cowplot::get_legend(legend_plot))

# 11. Composition ---------------------------------------------------------------
map_row <- p1 + p2 + p3 + plot_layout(nrow = 1, widths = c(1, 1, 1))

final_plot <- map_row / legend_only +
  plot_layout(heights = c(1, 0.12)) +
  plot_annotation(
    title = "Harmonised bathymetric DEMs",
    subtitle = "Bathymetric surfaces after harmonisation and vertical co-registration",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = element_text(size = 9.2, colour = "#444444", hjust = 0),
      plot.margin = margin(4, 5, 2, 5),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

# 12. Export -------------------------------------------------------------------
# cowplot is needed only to reliably extract the shared legend.
if (!requireNamespace("cowplot", quietly = TRUE)) {
  stop("Install package 'cowplot' for the shared colour bar: install.packages('cowplot')")
}

output_stem <- file.path(
  out_dir,
  "Figure_2_bathymetric_DEMs_2011_2019_2023"
)

save_common <- list(
  plot = final_plot,
  width = 270,
  height = 104,
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

# 13. Reproducibility record ----------------------------------------------------
record <- c(
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("DEM stack: ", dem_file),
  paste0("Boundary: ", boundary_file),
  "Displayed survey years: 2011, 2019, 2023",
  paste0("Shared colour limits: ", paste(round(z_limits, 3), collapse = " to "), " m"),
  paste0("Vertical exaggeration: ", VERTICAL_EXAGGERATION),
  paste0("Light altitude: ", LIGHT_ALTITUDE, " degrees"),
  paste0("Light azimuth: ", LIGHT_AZIMUTH, " degrees"),
  paste0("Ambient light factor: ", AMBIENT_LIGHT),
  paste0("Direct light factor: ", DIRECT_LIGHT),
  paste0("Brightness range: ", MIN_BRIGHTNESS, " to ", MAX_BRIGHTNESS),
  paste0("Display aggregation factor: ", agg_factor),
  "Hillshade was composited directly into RGB colours; no alpha raster was used.",
  "The analytical DEM values and source resolution were not modified."
)
writeLines(record, paste0(output_stem, "_reproducibility.txt"), useBytes = TRUE)

outputs <- c(
  paste0(output_stem, ".png"),
  paste0(output_stem, ".tif"),
  paste0(output_stem, ".pdf"),
  paste0(output_stem, "_reproducibility.txt")
)
check <- data.frame(
  file = outputs,
  exists = file.exists(outputs),
  size_mb = round(file.info(outputs)$size / 1024^2, 3)
)
print(check)
if (!all(check$exists)) stop("At least one Figure 2 output was not created.")

message("")
message("Figure 2 Global Mapper-style v4 completed successfully.")
message("PNG: ", outputs[1])
message("TIFF: ", outputs[2])
message("PDF: ", outputs[3])
message("")
