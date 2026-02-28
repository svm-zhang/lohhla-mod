parse_file_path <- function(file) {
  if (!file.exists(file)) {
    print(paste("[ERROR] Cannot find the file provided: ", file, sep = ""))
    quit(status = 1)
  }
  file <- normalizePath(file, mustWork = TRUE)
}

parse_dir_path <- function(dir, create) {
  if (!dir.exists(dir)) {
    if (create != TRUE) {
      print(paste("[ERROR] Cannot find the file provided: ", dir, sep = ""))
    } else {
      dir.create(dir, recursive = TRUE)
      dir <- normalizePath(dir, mustWork = TRUE)
    }
  }
}
