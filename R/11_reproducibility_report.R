source("00_config.R", encoding="UTF-8")
writeLines(capture.output(sessionInfo()),file.path(DIRS["logs"],"session_info.txt"))
files <- list.files(OUT,recursive=TRUE,full.names=TRUE)
fwrite(data.table(file=files,size_bytes=file.info(files)$size),file.path(DIRS["tables"],"output_inventory.csv"))
