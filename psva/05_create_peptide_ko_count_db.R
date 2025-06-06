# Load libraries
library(dplyr)

# ==================================
# Number KO associated to each peptide
# ==================================

# Working directory
work_path <- getwd()

# Load kegg_ko
load(paste0(work_path, "/psva/kegg_ko.RData"))
     
# Read core_pep_kegg_db.csv
core_pep_kegg_db <- read.csv(paste0(work_path, "/psva/core_pep_kegg_db.csv")) %>%
  select(Peptide, Protein, ID, KEGG)
nrow(core_pep_kegg_db)
# 
     
# Join core_pep_kegg_db and kegg_ko by KEGG
core_pep_ko_db <- inner_join(core_pep_kegg_db, kegg_ko, by = "KEGG")
nrow(core_pep_ko_db)
#

# Count how often each peptide is associated with each KEGG term 
count_core_pep_ko_db <- core_pep_ko_db %>% count(Peptide, KO)
nrow(count_core_pep_ko_db)
# 

# Remove ko: suffix
count_core_pep_ko_db$KO <- mapply(function(x) substr(x, 4, nchar(x)),
                                  count_core_pep_ko_db$KO)

# Save peptide-KO-count dataset as csv file
write.csv(count_core_pep_ko_db, file = paste0(work_path,
                                              "/psva/count_core_pep_ko_db.csv"))
print("count_core_pep_ko_db.csv saved.")