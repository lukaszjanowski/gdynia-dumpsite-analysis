source("00_config.R", encoding="UTF-8")
figure_scripts <- c("R/fig01_study_area.R","R/fig02_bathymetric_dems.R","R/fig03_dod_intervals.R","R/fig04_dominant_fields.R","R/figures_quantitative_and_supplementary.R")
for (f in figure_scripts) { message("Running ",f); source(f,encoding="UTF-8") }
