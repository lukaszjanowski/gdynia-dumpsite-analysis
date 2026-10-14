# ==============================================================================
# STEP 06 - Publication tables and figures
# Requires completed steps 03, 04, 04b/04c and 05.
# ============================================================================== 
source("00_config.R", encoding="UTF-8")

req <- c(
  file.path(DIRS["tables"], "dod_volume_summary_all_lod.csv"),
  file.path(DIRS["tables"], "vertical_offset_diagnostics.csv"),
  file.path(DIRS["tables"], "dominant_depositional_fields_2019_2023.csv"),
  file.path(DIRS["tables"], "boundary_primary_exact_summary.csv"),
  file.path(DIRS["rasters"], "Gdynia_DEM_2011_2019_2023_harmonized.tif"),
  file.path(DIRS["rasters"], "DoD_2011_2019_raw.tif"),
  file.path(DIRS["rasters"], "DoD_2019_2023_raw.tif")
)
miss <- req[!file.exists(req)]; if(length(miss)) stop("Missing required files:\n", paste(miss,collapse="\n"))

dodsum <- fread(req[1]); coreg <- fread(req[2]); fields <- fread(req[3]); boundary <- fread(req[4])

# Reconstruct the vertical offset from diagnostics if the column is absent.
# Standardise the actual co-registration diagnostic schema.
required_coreg_columns <- c(
  "pair",
  "n",
  "median_shift",
  "mean_shift",
  "mad",
  "rmse",
  "q025",
  "q975",
  "status"
)

missing_coreg_columns <- setdiff(
  required_coreg_columns,
  names(coreg)
)

if (length(missing_coreg_columns)) {
  stop(
    "Missing columns in vertical diagnostics: ",
    paste(missing_coreg_columns, collapse = ", ")
  )
}

# The pair field contains both the interval and correction stage.
coreg[, interval := sub(
  "_(before|after)_correction$",
  "",
  pair
)]

coreg[, stage := data.table::fifelse(
  grepl("_before_correction$", pair),
  "before_correction",
  data.table::fifelse(
    grepl("_after_correction$", pair),
    "after_correction",
    NA_character_
  )
)]

if (anyNA(coreg$stage)) {
  stop(
    "Could not determine correction stage from pair values: ",
    paste(coreg[is.na(stage), pair], collapse = ", ")
  )
}

# The applied vertical correction is the pre-correction median shift.
offset_lookup <- coreg[
  stage == "before_correction",
  .(
    applied_offset_m = median_shift[1]
  ),
  by = interval
]

coreg <- merge(
  coreg,
  offset_lookup,
  by = "interval",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(coreg$applied_offset_m)) {
  stop(
    "Applied vertical offset could not be assigned to every row."
  )
}

# Publication-facing column names.
coreg[, median_difference_m := median_shift]
coreg[, mean_difference_m := mean_shift]
coreg[, mad_m := mad]
coreg[, rmse_m := rmse]
coreg[, q025_m := q025]
coreg[, q975_m := q975]

t1 <- coreg[, .(
  pair = interval,
  stage,
  applied_offset_m,
  n,
  median_difference_m,
  mean_difference_m,
  mad_m,
  rmse_m,
  q025_m,
  q975_m,
  status
)]

t1[
  ,
  stage := factor(
    stage,
    levels = c(
      "before_correction",
      "after_correction"
    )
  )
]

t1[
  ,
  stage := factor(
    stage,
    levels = c(
      "before_correction",
      "after_correction"
    )
  )
]

data.table::setorder(
  t1,
  pair,
  stage
)

t1[
  ,
  stage := as.character(stage)
]

t1[
  ,
  stage := as.character(stage)
]
data.table::fwrite(
  t1,
  file.path(
    DIRS["tables"],
    "Table_1_coregistration.csv"
  )
)
fwrite(t1, file.path(DIRS["tables"], "Table_1_coregistration.csv"))

# Table 2: principal DoD metrics. 0.20 m for both, plus conservative 0.35 m for 2019-2023.
t2 <- dodsum[(pair=="2011_2019" & lod_m==0.20) | (pair=="2019_2023" & lod_m %in% c(0.20,0.35))]
t2 <- t2[zone %in% c("inside_dumpsite","outside_dumpsite","full_common_coverage")]
fwrite(t2, file.path(DIRS["tables"], "Table_2_DoD_principal_results.csv"))

# Table 3: exact boundary split at primary conservative criteria.
fwrite(boundary, file.path(DIRS["tables"], "Table_3_boundary_crossing.csv"))

# Table 4: dominant fields.
fwrite(fields, file.path(DIRS["tables"], "Table_4_dominant_fields.csv"))

# Figure 1: DEM time series.
dems <- rast(req[5]); names(dems) <- c("2011","2019","2023")
q_dem <- as.numeric(global(dems[[1]], function(x,...) quantile(x,c(.02,.98),na.rm=TRUE,names=FALSE), na.rm=TRUE)[1,])
png(file.path(DIRS["figures"], "Figure_1_DEM_time_series.png"), width=4200, height=1600, res=300)
par(mfrow=c(1,3), mar=c(3,3,3,5))
for(i in 1:3) { plot(dems[[i]], range=q_dem, col=hcl.colors(100,"Batlow"), main=paste("DEM",names(dems)[i]), axes=TRUE) }
dev.off()

# Figure 2: DoD panels with official boundary.
off <- vect(DUMPSITE_VECTOR); if(!same.crs(off,dems)) off <- project(off,crs(dems)); off <- aggregate(makeValid(off))
d11 <- rast(req[6]); d24 <- rast(req[7])
png(file.path(DIRS["figures"], "Figure_2_DoD_intervals.png"), width=3400, height=1700, res=300)
par(mfrow=c(1,2),mar=c(3,3,3,5))
plot(d11, range=c(-1,1.5), col=hcl.colors(101,"Blue-Red 3",rev=TRUE), main="2011-2019", axes=TRUE); lines(off,lwd=2)
plot(d24, range=c(-1,3), col=hcl.colors(101,"Blue-Red 3",rev=TRUE), main="2019-2023", axes=TRUE); lines(off,lwd=2)
dev.off()

# Figure 3: clean dominant-field and boundary map without all 88 labels.
pfile <- file.path(DIRS["vectors"], "depositional_patches_2019_2023.gpkg")
if(file.exists(pfile)) {
  pp <- vect(pfile)
  dom <- pp[pp$dominant_field %in% TRUE,]
  buf <- vect(BUFFER_VECTOR); if (!same.crs(buf, d24)) buf <- project(buf, crs(d24)); buf <- aggregate(makeValid(buf))
  png(file.path(DIRS["figures"], "Figure_3_dominant_fields_boundary.png"), width=2600,height=2200,res=300)
  par(mar=c(4,4,3,6)); plot(d24,range=c(-1,3),col=hcl.colors(101,"Blue-Red 3",rev=TRUE),main="Dominant depositional fields and boundary",axes=TRUE)
  lines(buf,col="grey35",lty=2,lwd=1.5); lines(off,lwd=2.5); lines(pp,col="gold",lwd=.5); lines(dom,col="cyan",lwd=3)
  legend("topright",c("official boundary","200 m buffer","retained patches","two dominant fields"),col=c("black","grey35","gold","cyan"),lty=c(1,2,1,1),lwd=c(2.5,1.5,.5,3),bty="n")
  dev.off()
}

# Figure 4: sensitivity plot already generated by step 05; create combined concise version.
sens <- fread(file.path(DIRS["tables"], "boundary_sensitivity_zone_summary.csv"))
p <- ggplot(sens[min_patch_area_m2==500], aes(threshold_m, positive_volume_m3/1e6, color=zone)) +
  geom_line(linewidth=1) + geom_point(size=2) + theme_bw() +
  labs(x="DoD threshold (m)",y="Positive volume (million m3)",color="Zone")
ggsave(file.path(DIRS["figures"], "Figure_4_sensitivity.png"),p,width=7.5,height=5,dpi=300)

# Captions text file.
captions <- c(
"Figure 1. Harmonised bathymetric DEMs for 2011, 2019, and 2023 on a common 1 m grid.",
"Figure 2. Co-registered DEMs of Difference for 2011-2019 and 2019-2023. Positive values denote seabed raising.",
"Figure 3. Conservative depositional patches (DoD > 0.35 m; area >= 500 m2), the two dominant fields, official boundary, and 200 m buffer.",
"Figure 4. Sensitivity of positive-volume estimates to the DoD threshold for patches >= 500 m2."
)
writeLines(captions,file.path(DIRS["tables"],"publication_figure_captions.txt"))
message("Step 10 publication tables completed successfully.")
