if (tolower(Sys.info()[["user"]]) == "jardang") {
  .libPaths(c("E:/Documents/jardang/r_packages/win-library/4.4.1", .libPaths()))
  DATA_DIR <- "X:/Documents/JARDANG/data"
  REPO_DIR <- "C:/Users/jardang/Documents/facts-emissions-across-network"
  OUTPUT_DIR <- "X:/Documents/JARDANG/facts_emissions_output"
} else if (tolower(Sys.info()[["user"]]) == "jota_"){
  DATA_DIR <- "C:/Users/jota_/Documents/NBB_data/"
  REPO_DIR <- local({
    for (i in sys.nframe():1) {
      fn <- sys.frame(i)$ofile
      if (!is.null(fn)) return(normalizePath(dirname(fn), winslash = "/"))
    }
    normalizePath(getwd(), winslash = "/")
  })
  OUTPUT_DIR <- "C:/Users/jota_/Documents/facts_emissions_output"
} else {
  stop("Define directories for this user.")
}

PROC_DATA <- file.path(DATA_DIR, "processed")
RAW_DATA <- file.path(DATA_DIR, "raw")
INT_DATA <- file.path(DATA_DIR, "intermediate")
