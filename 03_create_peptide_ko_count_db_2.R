# Load libraries
library(dplyr)

# ==================================
# Number KO associated to each peptide
# ==================================

# Working directory
work_path <- getwd()

# Read core_pep_kegg_db.csv
core_pep_kegg_db <- read.csv(paste0(work_path, "/core_pep_kegg_db.csv")) %>%
  select(-X)

# Count how often each peptide is associated with each KEGG term by summing
# the occurrences per term
count_core_pep_ko_db <- core_pep_ko_db %>% count(Peptide, KO)

# Remove ko: suffix
count_core_pep_ko_db$KO <- mapply(function(x) substr(x, 4, nchar(x)),
                                  count_core_pep_ko_db$KO)

# Save peptide-KO-count dataset as csv file
write.csv(count_core_pep_ko_db, file = paste0(work_path,
                                              "/count_core_pep_ko_db_2.csv"))
print("count_core_pep_ko_db_2.csv saved.")