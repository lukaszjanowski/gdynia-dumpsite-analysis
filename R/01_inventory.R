source("00_config.R", encoding = "UTF-8")

F <- CFG$files
required_names <- c("dem2011", "dem2019", "dem2023", "dumpsite", "stable")

paths <- lapply(F, function(x) {
  if (grepl("^([A-Za-z]:[/\\\\]|/)", x)) x else file.path(ROOT, x)
})
names(paths) <- names(F)

# Do not use `key = ...` inside data.table(); `key` is a reserved argument
# that requests setkeyv() instead of creating a column.
inventory <- data.table::data.table(
  input_name = names(paths),
  path = unname(unlist(paths, use.names = FALSE))
)

inventory[, required := input_name %in% required_names]
inventory[, exists := file.exists(path)]
inventory[, size_bytes := ifelse(exists, file.info(path)$size, NA_real_)]

inventory_file <- file.path(DIRS["tables"], "input_inventory.csv")
data.table::fwrite(inventory, inventory_file)

missing_required <- inventory[required == TRUE & exists == FALSE]
if (nrow(missing_required) > 0L) {
  stop(
    "Missing required input files:\n",
    paste0("- ", missing_required$input_name, ": ", missing_required$path, collapse = "\n"),
    "\n\nEdit config/config.yml or copy the files to data/raw/.\n",
    "Inventory written to: ", inventory_file
  )
}

saveRDS(paths, file.path(OUT, "paths.rds"))
message("Input inventory completed: ", inventory_file)
message("Path registry created: ", file.path(OUT, "paths.rds"))
print(inventory)
