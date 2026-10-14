testthat::test_that("publication year is 2023", {cfg <- yaml::read_yaml("config/config.yml"); testthat::expect_true("dem2023" %in% names(cfg$files))})
testthat::test_that("final Figure 4 uses legend below", {x <- readLines("R/fig04_dominant_fields.R",warn=FALSE); testthat::expect_true(any(grepl('legend.position = "bottom"',x,fixed=TRUE)))})
testthat::test_that("core production stages exist", {testthat::expect_true(all(file.exists(c("R/02_03_harmonise_coregister.R","R/04_dod_and_volumes.R","R/08_sensitivity_analysis.R"))))})
