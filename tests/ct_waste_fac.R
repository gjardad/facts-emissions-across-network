d <- read.csv("C:/Users/jota_/Documents/NBB_data/raw/Climate TRACE/BEL/DATA/power/electricity-generation_emissions_sources_v5_4_1.csv", stringsAsFactors=FALSE)
fac <- unique(d[, c("source_id","source_name","source_type")])
w2 <- fac[grepl("waste", fac$source_type, ignore.case=TRUE),]
cat("Waste source_type facilities in CT electricity:\n")
for(i in seq_len(nrow(w2))) {
  cat(sprintf("  %s | type=%s | %s\n", w2$source_id[i], w2$source_type[i], w2$source_name[i]))
}
