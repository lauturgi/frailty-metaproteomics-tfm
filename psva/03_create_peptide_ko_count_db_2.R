# Load libraries
library(dplyr)
library(data.table)

# ==================================
# Number KO associated to each peptide
# ==================================

# Working directory
work_path <- getwd()

# Read core_pep_kegg_db.csv
core_pep_kegg_db <- fread(paste0(work_path, "/psva/core_pep_kegg_db_2.csv"),
                          sep = ",", header = TRUE) %>%
  select(Peptide, Protein, ID, KO)

nrow(core_pep_kegg_db)
# 477498868

# Count how often each peptide is associated with each KO
count_core_pep_ko_db <- core_pep_kegg_db %>% count(Peptide, KO)

nrow(count_core_pep_ko_db)
# 330446378

# Remove ko: suffix
count_core_pep_ko_db$KO <- mapply(function(x) substr(x, 4, nchar(x)),
                                  count_core_pep_ko_db$KO)

# Rename columns
count_core_pep_ko_db <- count_core_pep_ko_db %>%
  setNames(c("pep", "ko", "count"))

# Get information on number of proteins for a peptide
count_core_pep_ko_db <- count_core_pep_ko_db %>% group_by(pep) %>% 
  summarize(total = sum(count))  %>%
  merge(., count_core_pep_ko_db, by='pep', all.y=T) %>%
  mutate(prop = count/total) %>% select(pep, ko, prop)

nrow(count_core_pep_ko_db)
# 330446378

# Unique names for duplicate peptides
count_core_pep_ko_db$newpep_name <- make.names(count_core_pep_ko_db$pep,
                                               unique=T)

# Save count_core_pep_ko_db
save(count_core_pep_ko_db,
     file = paste0(work_path, "/psva/data/count_core_pep_ko_db_2.RData"))

# Save peptide-KO-count dataset as csv file
write.csv(count_core_pep_ko_db, file = paste0(work_path,
                                              "/psva/count_core_pep_ko_db_2.csv"))
print("count_core_pep_ko_db_2.csv saved.")