# Load libraries
library(dplyr)

# ==================================
# Number KO associated to each peptide
# ==================================

# Working directory
work_path <- getwd()

# Load kegg_ko
load(paste0(work_path, "/kegg_ko.RData"))
     
print("Reading core_pep_kegg_db.csv ...")
# Read core_pep_kegg_db.csv
core_pep_kegg_db <- read.csv(paste0(work_path, "/core_pep_kegg_db.csv")) %>%
  select(-X)
print("Done.")
     
# Join core_pep_kegg_db and kegg_ko by KEGG
core_pep_ko_db <- inner_join(core_pep_kegg_db, kegg_ko, by = "KEGG")

# Count how often each peptide is associated with each KEGG term by summing
# the occurrences per term
print(paste("Counting how often each peptide is associated with each KO",
            "(peptides shared by multiple proteins) ... "))
count_core_pep_ko_db <- core_pep_ko_db %>% count(Peptide, KO)
print("Done.")

# Remove ko: suffix
count_core_pep_ko_db$KO <- mapply(function(x) substr(x, 4, nchar(x)),
                                  count_core_pep_ko_db$KO)

# Save peptide-KO-count dataset as csv file
write.csv(count_core_pep_ko_db, file = paste0(work_path,
                                              "/count_core_pep_ko_db.csv"))
print("count_core_pep_ko_db.csv saved.")