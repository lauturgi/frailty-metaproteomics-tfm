# Load libraries
library(dplyr)
library(data.table)

# ==================================
# Get unique KEGG IDs
# ==================================
# Working directory
work_path <- getwd()

# Read core_pep_kegg_db.csv
core_pep_kegg_db <- fread(paste0(work_path, "/psva/core_pep_kegg_db.csv"),
                                 sep = ",", header = TRUE) %>%
  select(Peptide, Protein, UniRef90, KEGG)
nrow(core_pep_kegg_db)
# 73998314

# Get all unique KEGG entries
kegg <- as.data.frame(unique(core_pep_kegg_db$KEGG))
colnames(kegg) <- c("KEGG")
print(paste("Number unique kegg:", nrow(kegg)))
# Number unique kegg: 549153

# Save unique_kegg.csv
write.csv(kegg, file = paste0(work_path, "/psva/unique_kegg.csv"))
print("unique_kegg.csv saved.")