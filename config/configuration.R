ref <- read.csv("dev/package_versions_reference.csv")
ref$version <- sapply(ref$pkg, \(p)
                      tryCatch(utils::packageDescription(p)$Version, error = \(e) NA))
write.csv(ref, "dev/package_versions_reference.csv", row.names = FALSE)
ref <- read.csv("package_versions_reference.csv")
ref <- ref[!is.na(ref$version), ]
renv::install(paste0(ref$pkg, "@", ref$version))


renv::snapshot()


lock <- jsonlite::read_json("renv.lock")$Packages
lock_df <- data.frame(
  pkg          = names(lock),
  lock_version = sapply(lock, `[[`, "Version")
)

check <- merge(ref, lock_df, by = "pkg", all.x = TRUE)
check$ok <- check$version == check$lock_version
check[!check$ok | is.na(check$ok), ]