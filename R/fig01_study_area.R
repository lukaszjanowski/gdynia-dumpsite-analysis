# ==============================================================================
# 06d_fix_figure1_study_area_v8_asymmetric_context.R
# Final two-panel Figure 1
#
# Panel a: full PLB220005 regional setting + European locator inset
# Panel b: official dumpsite, 200 m buffer and stable-seabed polygons
#
# Explicit paths, no recursive searching, no global s2 geometry repair.
# ==============================================================================

source("00_config.R", encoding = "UTF-8")

# 1. Packages ------------------------------------------------------------------
required <- c(
  "sf", "ggplot2", "patchwork", "ggspatial",
  "rnaturalearth", "rnaturalearthdata"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Install missing packages first:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "), "))"
  )
}

library(sf)
library(ggplot2)
library(patchwork)
library(ggspatial)
library(rnaturalearth)

# 2. Repository paths ---------------------------------------------------------
files <- list(
  plb220005 = file.path(ROOT, CFG$files$plb220005),
  plh220032 = file.path(ROOT, "data", "raw", "PLH220032.gpkg"),
  coast = file.path(ROOT, CFG$files$coastline),
  dumpsite = file.path(ROOT, CFG$files$dumpsite),
  buffer = file.path(ROOT, CFG$files$dumpsite),
  stable = file.path(ROOT, CFG$files$stable)
)
missing_files <- unlist(files[c("plb220005", "plh220032", "coast", "dumpsite", "stable")])[!file.exists(unlist(files[c("plb220005", "plh220032", "coast", "dumpsite", "stable")]))]
if (length(missing_files)) stop("Missing Figure 1 source files:\n", paste(missing_files, collapse = "\n"))
fig_dir <- unname(DIRS["figures"])
out_dir <- fig_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 3. Read local data ------------------------------------------------------------
read_local <- function(path) {
  z <- suppressWarnings(st_read(path, quiet = TRUE, options = "ENCODING=UTF-8"))
  if (!nrow(z)) stop("Empty vector layer: ", path)
  if (is.na(st_crs(z))) stop("Undefined CRS: ", path)
  z <- st_zm(z, drop = TRUE, what = "ZM")

  old_s2 <- sf_use_s2()
  sf_use_s2(FALSE)
  z <- suppressWarnings(st_make_valid(z))
  sf_use_s2(old_s2)
  z
}

z <- lapply(files, read_local)
metric_crs <- st_crs(z$dumpsite)
if (is.na(metric_crs)) metric_crs <- st_crs(2180)

# Harmonise local data in the projected CRS before dissolving.
z_metric <- lapply(z, st_transform, crs = metric_crs)
z_metric$buffer <- st_buffer(st_union(z_metric$dumpsite), dist = CFG$analysis$external_zone_m)

union_to_ll <- function(object) {
  old_s2 <- sf_use_s2()
  sf_use_s2(FALSE)
  geom <- suppressWarnings(st_union(object))
  sf_use_s2(old_s2)
  st_transform(st_sf(id = 1L, geometry = geom), 4326)
}

plb_ll <- union_to_ll(z_metric$plb220005)
plh_ll <- union_to_ll(z_metric$plh220032)
dumpsite_ll <- union_to_ll(z_metric$dumpsite)
buffer_ll <- union_to_ll(z_metric$buffer)
stable_ll <- st_transform(z_metric$stable, 4326)
coast_ll <- st_transform(z_metric$coast, 4326)

# 4. Rectangular panel extents --------------------------------------------------
# Panel a: complete PLB220005 with 10 km westward and 25 km eastward context.
# west and east, respectively. The north-south extent follows the complete PLB220005 height.
# A very small vertical padding is retained only to prevent the boundary line
# from touching the map frame.
WEST_CONTEXT_M <- 10000
EAST_CONTEXT_M <- 25000
VERTICAL_PADDING_M <- 750

plb_bbox_metric <- st_bbox(
  st_union(z_metric$plb220005)
)

panel_a_bbox_metric <- st_bbox(
  c(
    xmin = as.numeric(plb_bbox_metric["xmin"]) - WEST_CONTEXT_M,
    ymin = as.numeric(plb_bbox_metric["ymin"]) - VERTICAL_PADDING_M,
    xmax = as.numeric(plb_bbox_metric["xmax"]) + EAST_CONTEXT_M,
    ymax = as.numeric(plb_bbox_metric["ymax"]) + VERTICAL_PADDING_M
  ),
  crs = metric_crs
)

panel_a_rect_metric <- st_as_sfc(panel_a_bbox_metric)
panel_a_rect_ll <- st_transform(panel_a_rect_metric, 4326)
panel_a_bbox_ll <- st_bbox(panel_a_rect_ll)

bbox_a <- c(
  xmin = as.numeric(panel_a_bbox_ll["xmin"]),
  ymin = as.numeric(panel_a_bbox_ll["ymin"]),
  xmax = as.numeric(panel_a_bbox_ll["xmax"]),
  ymax = as.numeric(panel_a_bbox_ll["ymax"])
)

regional_width_m <- as.numeric(
  panel_a_bbox_metric["xmax"] - panel_a_bbox_metric["xmin"]
)
regional_height_m <- as.numeric(
  panel_a_bbox_metric["ymax"] - panel_a_bbox_metric["ymin"]
)
regional_aspect <- regional_width_m / regional_height_m

# Panel b: metric 1.6 km margin around the 200 m buffer.
detail_metric <- st_buffer(st_union(z_metric$buffer), 1600)
detail_bbox_ll <- st_bbox(st_transform(st_as_sfc(st_bbox(detail_metric)), 4326))
bbox_b <- c(
  xmin = as.numeric(detail_bbox_ll["xmin"]),
  ymin = as.numeric(detail_bbox_ll["ymin"]),
  xmax = as.numeric(detail_bbox_ll["xmax"]),
  ymax = as.numeric(detail_bbox_ll["ymax"])
)

bbox_sf <- function(b) {
  st_sf(
    id = 1L,
    geometry = st_as_sfc(st_bbox(b, crs = st_crs(4326)))
  )
}

# Crop only for rendering efficiency; extents remain rectangular.
old_s2 <- sf_use_s2()
sf_use_s2(FALSE)
coast_a <- suppressWarnings(st_crop(coast_ll, st_bbox(bbox_sf(bbox_a))))
coast_b <- suppressWarnings(st_crop(coast_ll, st_bbox(bbox_sf(bbox_b))))
sf_use_s2(old_s2)

# 5. Safe Europe locator --------------------------------------------------------
# Select Europe before plotting and omit Russia. This avoids global multipart
# geometry artefacts and s2 errors while retaining a clear locator map.
world <- rnaturalearth::ne_countries(
  scale = "small", type = "countries", returnclass = "sf"
)
continent_field <- intersect(c("continent", "CONTINENT"), names(world))[1]
name_field <- intersect(c("admin", "name", "name_long", "sovereignt"), names(world))[1]
iso_field <- intersect(c("adm0_a3", "iso_a3", "ADM0_A3", "ISO_A3"), names(world))[1]
if (is.na(continent_field) || is.na(name_field)) {
  stop("Required Natural Earth country attributes were not found.")
}

europe <- world[
  (world[[continent_field]] == "Europe" |
     world[[name_field]] %in% c("Turkey", "Türkiye", "Cyprus")) &
    !world[[name_field]] %in% c("Russia", "Russian Federation"),
]
europe <- st_transform(europe, 4326)

if (!is.na(iso_field)) {
  poland <- europe[europe[[iso_field]] == "POL", ]
} else {
  poland <- europe[europe[[name_field]] %in% c("Poland", "Republic of Poland"), ]
}
if (!nrow(poland)) {
  poland <- europe[europe[[name_field]] %in% c("Poland", "Republic of Poland"), ]
}
if (!nrow(poland)) stop("Poland not found in Natural Earth data.")

site_point <- st_centroid(dumpsite_ll)
site_xy <- st_coordinates(site_point)[1, ]

p_inset <- ggplot() +
  geom_sf(data = europe, fill = "#E5E5E5", colour = "#858585", linewidth = 0.16) +
  geom_sf(data = poland, fill = "#7FA6C9", colour = "#3B668C", linewidth = 0.38) +
  geom_point(
    aes(x = site_xy[1], y = site_xy[2]),
    shape = 21, size = 2.2, stroke = 0.55,
    fill = "#D7191C", colour = "#7F0000"
  ) +
  coord_sf(
    xlim = c(-12, 35), ylim = c(34, 72),
    crs = st_crs(4326), datum = NA, expand = FALSE
  ) +
  labs(title = "European setting") +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 7.8, hjust = 0.5),
    panel.background = element_rect(fill = "white", colour = "#555555", linewidth = 0.4),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(2, 2, 2, 2)
  )

# 6. Helpers and shared styling -------------------------------------------------
lon_lab <- function(v) paste0(formatC(v, format = "f", digits = 2), "°E")
lat_lab <- function(v) paste0(formatC(v, format = "f", digits = 2), "°N")

map_theme <- theme_minimal(base_size = 9, base_family = "sans") +
  theme(
    panel.background = element_rect(fill = "#F4F9FC", colour = NA),
    panel.grid.major = element_line(colour = "#DBE5EA", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 8.3),
    axis.text = element_text(size = 7.0, colour = "#333333"),
    plot.title = element_text(face = "bold", size = 10.5),
    plot.subtitle = element_text(size = 7.8, colour = "#444444"),
    plot.margin = margin(2, 3, 2, 3),
    plot.background = element_rect(fill = "white", colour = NA)
  )

COL_SITE <- "#B2182B"
COL_BUFFER <- "#E6553A"
COL_STABLE <- "#1B7837"
COL_PLB <- "#6A51A3"
COL_PLH <- "#238B45"

# 7. Panel a -------------------------------------------------------------------
# Use direct cartographic labels rather than a separate legend. This avoids the
# hidden PLB220005 legend and reduces white space below the panel.
panel_a_base <- ggplot() +
  geom_sf(
    data = coast_a,
    fill = "#DED8CC",
    colour = "#6B6259",
    linewidth = 0.25
  ) +
  geom_sf(
    data = plb_ll,
    fill = "#E8E3F3",
    colour = COL_PLB,
    linewidth = 0.70,
    alpha = 0.80
  ) +
  geom_sf(
    data = plh_ll,
    fill = "#D9F0D3",
    colour = COL_PLH,
    linewidth = 0.70,
    alpha = 0.50
  ) +
  # Important: the disposal site is outline-only in panel a.
  # No red fill and no label are placed over the polygon.
  geom_sf(
    data = dumpsite_ll,
    aes(colour = "Gdynia disposal site"),
    fill = NA,
    linewidth = 1.10,
    show.legend = TRUE
  ) +
  annotate(
    "label",
    x = mean(st_bbox(plb_ll)[c("xmin", "xmax")]),
    y = st_bbox(plb_ll)["ymax"] - 0.035 *
      (st_bbox(plb_ll)["ymax"] - st_bbox(plb_ll)["ymin"]),
    label = "PLB220005 Zatoka Pucka",
    size = 3.0,
    fontface = "bold",
    colour = "#4B3B73",
    fill = scales::alpha("white", 0.78),
    linewidth = 0.12
  ) +
  annotate(
    "label",
    x = mean(st_bbox(plh_ll)[c("xmin", "xmax")]) - 0.2 *
      (bbox_a["xmax"] - bbox_a["xmin"]),
      (bbox_a["xmax"] - bbox_a["xmin"]),
    y = st_bbox(plh_ll)["ymin"] + 0.10 *
      (st_bbox(plh_ll)["ymax"] - st_bbox(plh_ll)["ymin"]) + 0.10 *
      (bbox_a["ymax"] - bbox_a["ymin"]),
    label = "PLH220032 Zatoka Pucka\ni Półwysep Helski",
    size = 2.8,
    fontface = "bold",
    colour = COL_PLH,
    fill = scales::alpha("white", 0.78),
    linewidth = 0.12
  ) +
  scale_colour_manual(
    name = NULL,
    breaks = "Gdynia disposal site",
    values = c("Gdynia disposal site" = COL_SITE)
  ) +
  scale_x_continuous(
    breaks = pretty(bbox_a[c("xmin", "xmax")], n = 5),
    labels = lon_lab
  ) +
  scale_y_continuous(
    breaks = pretty(bbox_a[c("ymin", "ymax")], n = 5),
    labels = lat_lab
  ) +
  coord_sf(
    xlim = bbox_a[c("xmin", "xmax")],
    ylim = bbox_a[c("ymin", "ymax")],
    crs = st_crs(4326),
    datum = st_crs(4326),
    expand = FALSE
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.20,
    text_cex = 0.57,
    line_width = 0.42,
    unit_category = "metric"
  ) +
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    pad_x = grid::unit(3, "mm"),
    pad_y = grid::unit(3, "mm"),
    style = north_arrow_orienteering(text_size = 5.3),
    height = grid::unit(9, "mm"),
    width = grid::unit(9, "mm")
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        fill = NA,
        linewidth = 1.10,
        linetype = "solid"
      )
    )
  ) +
  labs(
    title = "a  Regional setting",
    subtitle = "PLB220005 and PLH220032 Natura 2000 areas with regional context",
    x = NULL,
    y = NULL
  ) +
  map_theme +
  theme(
    legend.position = c(0.76, 0.08),
    legend.justification = c(0.5, 0.5),
    legend.direction = "horizontal",
    legend.background = element_rect(
      fill = scales::alpha("white", 0.88),
      colour = "#BDBDBD",
      linewidth = 0.25
    ),
    legend.key = element_blank(),
    legend.text = element_text(size = 6.8),
    legend.key.width = grid::unit(7, "mm"),
    legend.key.height = grid::unit(3.2, "mm"),
    legend.margin = margin(1.5, 3, 1.5, 3)
  )

panel_a <- panel_a_base +
  inset_element(
    p_inset,
    left = 0.65, bottom = 0.66, right = 0.985, top = 0.985,
    align_to = "panel", on_top = TRUE
  )

# 8. Panel b -------------------------------------------------------------------
# Restore the successful detailed panel from v3, with shorter legend labels and
# enough internal width to avoid clipping.
panel_b <- ggplot() +
  geom_sf(data = coast_b, fill = "#DED8CC", colour = "#6B6259", linewidth = 0.25) +
  geom_sf(
    data = buffer_ll,
    aes(colour = "200 m buffer"),
    fill = NA, linewidth = 0.75, linetype = "22"
  ) +
  geom_sf(
    data = dumpsite_ll,
    aes(colour = "Official boundary"),
    fill = scales::alpha(COL_SITE, 0.10), linewidth = 1.05
  ) +
  geom_sf(
    data = stable_ll,
    aes(colour = "Stable-seabed polygons"),
    fill = scales::alpha(COL_STABLE, 0.22), linewidth = 0.80
  ) +
  scale_colour_manual(
    name = NULL,
    breaks = c("200 m buffer", "Official boundary", "Stable-seabed polygons"),
    values = c(
      "200 m buffer" = COL_BUFFER,
      "Official boundary" = COL_SITE,
      "Stable-seabed polygons" = COL_STABLE
    )
  ) +
  scale_x_continuous(
    breaks = pretty(bbox_b[c("xmin", "xmax")], n = 4), labels = lon_lab
  ) +
  scale_y_continuous(
    breaks = pretty(bbox_b[c("ymin", "ymax")], n = 5), labels = lat_lab
  ) +
  coord_sf(
    xlim = bbox_b[c("xmin", "xmax")],
    ylim = bbox_b[c("ymin", "ymax")],
    crs = st_crs(4326), datum = st_crs(4326), expand = FALSE
  ) +
  annotation_scale(
    location = "bl", width_hint = 0.22,
    text_cex = 0.57, line_width = 0.42, unit_category = "metric"
  ) +
  annotation_north_arrow(
    location = "tl", which_north = "true",
    pad_x = grid::unit(3, "mm"), pad_y = grid::unit(3, "mm"),
    style = north_arrow_orienteering(text_size = 5.3),
    height = grid::unit(9, "mm"), width = grid::unit(9, "mm")
  ) +
  labs(
    title = "b  Analysis area",
    subtitle = "Official boundary, 200 m buffer and stable-seabed controls",
    x = NULL, y = NULL
  ) +
  guides(
    colour = guide_legend(
      nrow = 3, byrow = TRUE,
      override.aes = list(
        fill = NA,
        linewidth = c(0.9, 1.1, 0.9),
        linetype = c("22", "solid", "solid")
      )
    )
  ) +
  map_theme +
  theme(
    # Keep the analytical-layer legend inside panel b. This removes the large
    # white strip previously reserved below the map.
    legend.position = c(0.72, 0.105),
    legend.justification = c(0.5, 0.5),
    legend.direction = "vertical",
    legend.background = element_rect(
      fill = scales::alpha("white", 0.88),
      colour = "#BDBDBD",
      linewidth = 0.25
    ),
    legend.key = element_blank(),
    legend.text = element_text(size = 6.1),
    legend.key.width = grid::unit(6.0, "mm"),
    legend.key.height = grid::unit(3.2, "mm"),
    legend.spacing.x = grid::unit(1.0, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm"),
    legend.margin = margin(1.5, 2.5, 1.5, 2.5)
  )

# 9. Combine panels -------------------------------------------------------------
final_plot <- panel_a + panel_b +
  plot_layout(
    widths = c(1.12, 1.00),
    guides = "keep"
  ) +
  plot_annotation(
    title = "Location of the Gdynia dredged-material disposal site",
    subtitle = "Regional environmental context and detailed analytical boundaries",
    theme = theme(
      plot.title = element_text(face = "bold", size = 12.5, hjust = 0),
      plot.subtitle = element_text(size = 8.8, colour = "#444444", hjust = 0),
      plot.margin = margin(4, 5, 2, 5),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

# 10. Export -------------------------------------------------------------------
output_stem <- file.path(
  out_dir,
  "Figure_1_study_area_Gdynia_2011_2019_2023_two_panel_compact_v8"
)

save_common <- list(
  plot = final_plot,
  width = 235,
  height = 120,
  units = "mm",
  bg = "white",
  limitsize = FALSE
)

do.call(ggsave, c(list(filename = paste0(output_stem, ".png"), dpi = 300), save_common))
do.call(ggsave, c(list(
  filename = paste0(output_stem, ".tif"),
  dpi = 600, compression = "lzw"
), save_common))
do.call(ggsave, c(list(
  filename = paste0(output_stem, ".pdf"),
  device = grDevices::cairo_pdf
), save_common))

# 11. Check --------------------------------------------------------------------
outputs <- c(
  paste0(output_stem, ".png"),
  paste0(output_stem, ".tif"),
  paste0(output_stem, ".pdf")
)
check <- data.frame(
  file = outputs,
  exists = file.exists(outputs),
  size_mb = round(file.info(outputs)$size / 1024^2, 3)
)
print(check)
if (!all(check$exists)) stop("At least one output was not created.")

message("")
message("Compact two-panel Figure 1 v8 completed successfully.")
message("Panel-a western context: ", WEST_CONTEXT_M / 1000, " km")
message("Panel-a eastern context: ", EAST_CONTEXT_M / 1000, " km")
message("Panel-a vertical padding: ", VERTICAL_PADDING_M / 1000, " km north and south")
message("Panel-a physical aspect ratio: ", round(regional_aspect, 3))
message("PNG: ", outputs[1])
message("TIFF: ", outputs[2])
message("PDF: ", outputs[3])
message("")
