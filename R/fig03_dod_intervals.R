# ==============================================================================
# 06f_figure3_DoD_intervals_final.R
# Figure 3: co-registered DEMs of Difference
# Positive values = seabed raising (red)
# Negative values = seabed lowering (blue)
# Intervals: 2011-2019 and 2019-2023
# ==============================================================================

source("00_config.R", encoding = "UTF-8")

required <- c("terra", "sf", "ggplot2", "patchwork", "ggspatial", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))

library(terra)
library(sf)
library(ggplot2)
library(patchwork)
library(ggspatial)
library(scales)

# 1. Source files ---------------------------------------------------------------
first_existing <- function(paths, label) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) stop(label, " not found. Checked:\n", paste(paths, collapse = "\n"))
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

dod_2011_2019_file <- first_existing(
  c(
    file.path(DIRS["rasters"], "DoD_2011_2019_raw.tif"),
    file.path(OUTPUT_DIR, "rasters", "DoD_2011_2019_raw.tif")
  ),
  "2011-2019 raw DoD"
)

dod_2019_2023_file <- first_existing(
  c(
    file.path(DIRS["rasters"], "DoD_2019_2023_raw.tif"),
    file.path(DIRS["rasters"], "DoD_2019_2023_raw.tif"),
    file.path(OUTPUT_DIR, "rasters", "DoD_2019_2023_raw.tif")
  ),
  "2019-2023 raw DoD"
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

message("DoD 2011-2019: ", dod_2011_2019_file)
message("DoD 2019-2023: ", dod_2019_2023_file)
message("Official boundary: ", boundary_file)

# 2. Read and align rasters -----------------------------------------------------
dod_11_19 <- rast(dod_2011_2019_file)
dod_19_23 <- rast(dod_2019_2023_file)

if (nlyr(dod_11_19) != 1L || nlyr(dod_19_23) != 1L) {
  stop("Each DoD input must contain exactly one raster layer.")
}

# Match the second interval to the first interval if grid metadata differ.
if (!compareGeom(dod_11_19, dod_19_23, stopOnError = FALSE)) {
  message("Aligning 2019-2023 DoD to the 2011-2019 grid...")
  dod_19_23 <- project(dod_19_23, dod_11_19, method = "bilinear")
}

names(dod_11_19) <- "2011-2019"
names(dod_19_23) <- "2019-2023"

# Use the common analytical footprint in both panels.
common_mask <- !is.na(dod_11_19) & !is.na(dod_19_23)
dod_11_19 <- mask(dod_11_19, common_mask)
dod_19_23 <- mask(dod_19_23, common_mask)

dod_stack <- c(dod_11_19, dod_19_23)
dod_stack <- trim(dod_stack)

boundary <- st_read(boundary_file, quiet = TRUE)
if (is.na(st_crs(boundary))) stop("The official boundary has no defined CRS.")
boundary <- st_make_valid(boundary)
boundary <- st_transform(boundary, crs(dod_stack))

# 3. Confirm sign convention ----------------------------------------------------
# DoDs are assumed to follow later DEM minus earlier DEM:
# positive = seabed raising; negative = seabed lowering.
interval_stats <- data.frame(
  interval = names(dod_stack),
  minimum_m = as.numeric(global(dod_stack, "min", na.rm = TRUE)[, 1]),
  maximum_m = as.numeric(global(dod_stack, "max", na.rm = TRUE)[, 1]),
  mean_m = as.numeric(global(dod_stack, "mean", na.rm = TRUE)[, 1])
)
print(interval_stats)

if (!any(interval_stats$minimum_m < 0) || !any(interval_stats$maximum_m > 0)) {
  warning("The DoD rasters do not contain both negative and positive values. Check sign convention.")
}

# 4. Common symmetric display scale --------------------------------------------
# Use a joint robust absolute limit and round upward to the next 0.5 m. The same
# symmetric scale is applied to both panels. A minimum of 3 m preserves the scale
# used in earlier Figure_2_DoD_intervals products.
abs_q_fun <- function(x, ...) {
  quantile(abs(x), probs = 0.995, na.rm = TRUE, names = FALSE)
}
robust_abs <- max(as.numeric(global(dod_stack, abs_q_fun, na.rm = TRUE)[, 1]), na.rm = TRUE)
DISPLAY_LIMIT_M <- max(3.0, ceiling(robust_abs * 2) / 2)

message("Common symmetric display range: -", DISPLAY_LIMIT_M,
        " to +", DISPLAY_LIMIT_M, " m")

# More colour resolution near zero while retaining the common +/- range.
colour_break_values <- c(
  -DISPLAY_LIMIT_M, -2, -1, -0.5, -0.2, 0,
   0.2, 0.5, 1, 2, DISPLAY_LIMIT_M
)
colour_break_values <- sort(unique(pmax(-DISPLAY_LIMIT_M,
                                         pmin(DISPLAY_LIMIT_M, colour_break_values))))

colour_break_colours <- c(
  "#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
  "#F7F7F7",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"
)
# Interpolate the full colour specification to the number of retained break values.
colour_break_colours <- grDevices::colorRampPalette(colour_break_colours)(
  length(colour_break_values)
)
colour_positions <- scales::rescale(
  colour_break_values,
  from = c(-DISPLAY_LIMIT_M, DISPLAY_LIMIT_M)
)

# 5. Display copies -------------------------------------------------------------
TARGET_DISPLAY_CELLS <- 1200000
agg_factor <- max(1L, ceiling(sqrt(ncell(dod_stack[[1]]) / TARGET_DISPLAY_CELLS)))
message("Display aggregation factor: ", agg_factor)

if (agg_factor > 1L) {
  dod_display <- aggregate(dod_stack, fact = agg_factor, fun = mean, na.rm = TRUE)
} else {
  dod_display <- dod_stack
}

# WGS84 display to match Figures 1 and 2.
template_ll <- project(dod_display[[1]], "EPSG:4326", method = "bilinear")
dod_ll <- rast(lapply(1:2, function(i) {
  project(dod_display[[i]], template_ll, method = "bilinear")
}))
names(dod_ll) <- names(dod_stack)
boundary_ll <- st_transform(boundary, 4326)

map_extent <- ext(dod_ll[[1]])
x_limits <- c(map_extent[1], map_extent[2])
y_limits <- c(map_extent[3], map_extent[4])

raster_to_df <- function(r) {
  d <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(d)[3] <- "change_m"
  d
}
dod_df <- lapply(1:2, function(i) raster_to_df(dod_ll[[i]]))

# 6. Shared map style -----------------------------------------------------------
lon_lab <- function(v) paste0(formatC(v, format = "f", digits = 2), "°E")
lat_lab <- function(v) paste0(formatC(v, format = "f", digits = 2), "°N")

panel_theme <- theme_minimal(base_size = 9, base_family = "sans") +
  theme(
    panel.background = element_rect(fill = "#F4F9FC", colour = NA),
    panel.grid.major = element_line(colour = "#D8E3E8", linewidth = 0.22),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7.1, colour = "#333333"),
    plot.title = element_text(face = "bold", size = 10.8, hjust = 0),
    plot.margin = margin(2, 3, 2, 3),
    legend.title = element_text(face = "bold", size = 8.5),
    legend.text = element_text(size = 7.1)
  )

legend_ticks <- pretty(c(-DISPLAY_LIMIT_M, DISPLAY_LIMIT_M), n = 7)
legend_ticks <- legend_ticks[legend_ticks >= -DISPLAY_LIMIT_M &
                               legend_ticks <= DISPLAY_LIMIT_M]
if (!0 %in% legend_ticks) legend_ticks <- sort(c(legend_ticks, 0))

make_panel <- function(i, letter, show_y = TRUE) {
  p <- ggplot() +
    geom_raster(
      data = dod_df[[i]],
      aes(x = x, y = y, fill = change_m),
      interpolate = FALSE
    ) +
    geom_sf(
      data = boundary_ll,
      fill = NA,
      colour = "black",
      linewidth = 0.78,
      inherit.aes = FALSE
    ) +
    scale_fill_gradientn(
      colours = colour_break_colours,
      values = colour_positions,
      limits = c(-DISPLAY_LIMIT_M, DISPLAY_LIMIT_M),
      oob = scales::squish,
      na.value = NA,
      name = "Elevation change (m)",
      breaks = legend_ticks,
      labels = function(x) sprintf("%+.1f", x),
      guide = guide_colourbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(64, "mm"),
        barheight = grid::unit(3.8, "mm"),
        ticks.colour = "#333333",
        frame.colour = "#555555"
      )
    ) +
    scale_x_continuous(
      breaks = pretty(x_limits, n = 4),
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
      width_hint = 0.19,
      text_cex = 0.53,
      line_width = 0.40,
      unit_category = "metric",
      pad_x = grid::unit(2.5, "mm"),
      pad_y = grid::unit(2.5, "mm")
    ) +
    annotation_north_arrow(
      location = "tl",
      which_north = "true",
      style = north_arrow_orienteering(text_size = 4.8),
      height = grid::unit(7.5, "mm"),
      width = grid::unit(7.5, "mm"),
      pad_x = grid::unit(2.5, "mm"),
      pad_y = grid::unit(2.5, "mm")
    ) +
    labs(
      title = paste0(letter, "  ", names(dod_ll)[i]),
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

# 7. Final layout ---------------------------------------------------------------
final_plot <- p1 + p2 +
  plot_layout(nrow = 1, widths = c(1, 1), guides = "collect") +
  plot_annotation(
    title = "Co-registered DEMs of Difference",
    subtitle = paste0(
      "Positive values indicate seabed raising; negative values indicate seabed lowering"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = element_text(size = 9.2, colour = "#444444", hjust = 0),
      plot.margin = margin(4, 5, 2, 5),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(0, 0, 0, 0)
  )

# 8. Export --------------------------------------------------------------------
output_stem <- file.path(
  out_dir,
  "Figure_3_DoD_intervals_2011_2019_2023_final"
)

save_common <- list(
  plot = final_plot,
  width = 220,
  height = 112,
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

# 9. Reproducibility ------------------------------------------------------------
write.csv(
  interval_stats,
  paste0(output_stem, "_interval_statistics.csv"),
  row.names = FALSE
)

record <- c(
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("2011-2019 DoD: ", dod_2011_2019_file),
  paste0("2019-2023 DoD: ", dod_2019_2023_file),
  paste0("Official boundary: ", boundary_file),
  "DoD convention: later DEM minus earlier DEM",
  "Positive values: seabed raising; mapped in red",
  "Negative values: seabed lowering; mapped in blue",
  paste0("Common symmetric display limits: +/-", DISPLAY_LIMIT_M, " m"),
  paste0("Robust joint absolute 99.5th percentile: ", round(robust_abs, 3), " m"),
  paste0("Display aggregation factor: ", agg_factor),
  "The analytical DoD values and source resolution were not modified."
)
writeLines(record, paste0(output_stem, "_reproducibility.txt"), useBytes = TRUE)

outputs <- c(
  paste0(output_stem, ".png"),
  paste0(output_stem, ".tif"),
  paste0(output_stem, ".pdf"),
  paste0(output_stem, "_interval_statistics.csv"),
  paste0(output_stem, "_reproducibility.txt")
)
check <- data.frame(
  file = outputs,
  exists = file.exists(outputs),
  size_mb = round(file.info(outputs)$size / 1024^2, 3)
)
print(check)
if (!all(check$exists)) stop("At least one Figure 3 output was not created.")

message("")
message("Figure 3 completed successfully.")
message("Colour direction confirmed: negative blue, positive red.")
message("PNG: ", outputs[1])
message("TIFF: ", outputs[2])
message("PDF: ", outputs[3])
message("")
