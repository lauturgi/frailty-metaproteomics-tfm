# Load libraries
library(KEGGREST)
library(dplyr)

#core_pep_kegg_db_old <- read.csv("/mnt/d/proteomica/fragilidad/Pepfunk/database/core_pep_kegg_db.csv")

# ==================================
# KO from KEGG IDs
# ==================================
# Working directory
work_path <- getwd()

print("Reading core_pep_kegg_db.csv ...")
# Read core_pep_kegg_db.csv
core_pep_kegg_db <- read.csv(paste0(work_path, "/core_pep_kegg_db.csv")) %>%
  select(-X)
print("Done.")

# Get all unique KEGG entries
kegg <- as.data.frame(unique(core_pep_kegg_db$KEGG))
colnames(kegg) <- c("KEGG")
print(paste("Number unique kegg:", nrow(kegg)))
write.csv(kegg, file = paste0(work_path, "/unique_kegg.csv"))
print("unique_kegg.csv saved.")