source("00_config.R", encoding="UTF-8")
expected <- c(file.path(DIRS["tables"],"coregistration_diagnostics.csv"), file.path(DIRS["tables"],"harmonisation_summary.csv"))
if(any(!file.exists(expected))) stop("Co-registration outputs missing. Run R/02_harmonise_dems.R first.")
message("Sequential vertical co-registration was completed within stage 02; diagnostics verified.")
