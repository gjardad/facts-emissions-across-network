###############################################################################
# preprocess/split_proxy_rdata.R
#
# PURPOSE
#   One-off script to split the combined deployment_proxy_avg.RData (which
#   contains both proxy_avg and proxy_list) into two separate files.
#   Run once on RMD after run_subsampled_en.R was executed before the
#   two-file save convention was in place.
#
# RUNS ON: RMD
###############################################################################

if (tolower(Sys.info()[["user"]]) == "jardang") {
  REPO_DIR <- "C:/Users/jardang/Documents/facts-emissions-across-network"
} else if (tolower(Sys.info()[["user"]]) == "jota_") {
  REPO_DIR <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
                        error = function(e) normalizePath(getwd(), winslash = "/"))
  while (!file.exists(file.path(REPO_DIR, "paths.R"))) REPO_DIR <- dirname(REPO_DIR)
} else {
  stop("Define REPO_DIR for this user.")
}
source(file.path(REPO_DIR, "paths.R"))

cat("Loading combined deployment_proxy_avg.RData...\n")
load(file.path(PROC_DATA, "deployment_proxy_avg.RData"))

save(proxy_avg,  file = file.path(PROC_DATA, "deployment_proxy_avg.RData"))
cat("Saved deployment_proxy_avg.RData  —", nrow(proxy_avg), "rows\n")

save(proxy_list, file = file.path(PROC_DATA, "deployment_proxy_list.RData"))
cat("Saved deployment_proxy_list.RData —", length(proxy_list), "draws\n")
