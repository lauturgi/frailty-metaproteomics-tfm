# Load libraries
library(dplyr)
library(tibble)

# Working directory
work_path <- getwd()

# Data path
data_path <- paste0(work_path, "/psva/data/")

# Load pep_exp
load(paste0(data_path, "pep_exp.RData"))
norm_pep_exp <- pep_exp

nrow(norm_pep_exp)
# 699

# Load count_core_pep_ko_db_2
load(paste0(data_path, "count_core_pep_ko_db_2.RData"))
nrow(count_core_pep_ko_db)
# 330446378

# Join peptides and core_pep_kegg by "pep" to get ko, newpep_name and prop
matched_core_pep_ko <- norm_pep_exp %>%
  rownames_to_column(var="pep") %>%
  merge(., count_core_pep_ko_db, by="pep")
nrow(matched_core_pep_ko)
# 1139

# Number unique peptides
length(unique(matched_core_pep_ko[,"pep"]))
# 556

# All peptides from norm_pep_exp
all_pep <- rownames_to_column(norm_pep_exp, var = "pep")

# Get those not in matched_core_pep_ko$pep
unmatched_pep <- all_pep[!all_pep$pep %in% matched_core_pep_ko$pep, ]

# Save unmatched_pep
save(unmatched_pep, file = paste0(data_path, "unmatched_pep.RData"))

# Normalize intensities by proportion of functional annotation
matched_core_pep_ko[, grep("FT", names(matched_core_pep_ko))] <- 
  lapply(matched_core_pep_ko[, grep("FT", names(matched_core_pep_ko))], 
         function(x) log2(x * matched_core_pep_ko$prop))

# Save matched_core_pep_ko
save(matched_core_pep_ko, file = paste0(data_path, "matched_core_pep_ko.RData"))

# Set rownames to "newpep_name" and select keep "Intensity" columns
pep_kegg_exp <- matched_core_pep_ko %>%
  column_to_rownames(var="newpep_name") %>%
  select(contains("FT"))
nrow(pep_kegg_exp)
# 1139

# Save pep_kegg_exp
save(pep_kegg_exp, file = paste0(data_path, "pep_kegg_exp.RData"))
