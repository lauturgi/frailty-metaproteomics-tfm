# Load libraries
library(dplyr)
library(tibble)

# Working directory
work_path <- getwd()

# Load norm_pep_exp
#load(paste0(work_path, "/data/norm_pep_exp.RData"))

load(paste0(work_path, "/data/pep_exp.RData"))
norm_pep_exp <- pep_exp

# Read old core peptide db
#core_pep_kegg_old <- read.delim2(paste0("/mnt/d/proteomica/fragilidad/Pepfunk/",
#                                        "database/core_pep_kegg_db.csv"),
#                                 sep = ",", header=F) %>%
#  setNames(c("pep", "kegg", "count"))

# Read core peptide db
core_pep_kegg <- read.delim2(paste0(work_path,"/count_core_pep_ko_db.csv"), 
                             sep=",", header=T) %>%
  select(-X) %>%
  setNames(c("pep", "kegg", "count"))

# Get information on number of proteins for a peptide
core_pep_kegg <- core_pep_kegg %>% group_by(pep) %>% 
  summarize(total = sum(count))  %>%
  merge(., core_pep_kegg, by='pep', all.y=T) %>%
  mutate(prop = count/total) %>% select(pep, kegg, prop)

# Unique names for duplicate peptides
core_pep_kegg$newpep_name <- make.names(core_pep_kegg$pep, unique=T)

# Save core_pep_kegg
save(core_pep_kegg, file = paste0(work_path, "/data/core_pep_kegg.RData"))

# Join peptides and core_pep_kegg by "pep" to get ko, newpep_name and prop
matched_core_pep_kegg <- norm_pep_exp %>%
  rownames_to_column(var="pep") %>%
  merge(., core_pep_kegg, by="pep") 

# All peptides from norm_pep_exp
all_pep <- rownames_to_column(norm_pep_exp, var = "pep")

# Get those not in matched_core_pep_kegg$pep
unmatched_pep <- all_pep[!all_pep$pep %in% matched_core_pep_kegg$pep, ]

# Save unmatched_pep
save(unmatched_pep, file = paste0(work_path, "/data/unmatched_pep.RData"))

# Normalise intensities by proportion of functional annotation
matched_core_pep_kegg[, grep("Intensity", names(matched_core_pep_kegg))] <- 
  lapply(matched_core_pep_kegg[, grep("Intensity",
                                      names(matched_core_pep_kegg))], 
         function(x) log2(1 + (x * matched_core_pep_kegg$prop)))

# Save matched_core_pep_kegg
save(matched_core_pep_kegg, file = paste0(work_path,
                                          "/data/matched_core_pep_kegg.RData"))

# Set rownames to "newpep_name" and select keep "Intensity" columns
pep_kegg_exp <- matched_core_pep_kegg %>%
  column_to_rownames(var="newpep_name") %>%
  select(contains("Intensity"))

# Remove ".Intensity" from column names
colnames(pep_kegg_exp) <- sapply(strsplit(colnames(pep_kegg_exp), 
                                           "\\."), `[`, 1)

# Save pep_kegg_exp
save(pep_kegg_exp, file = paste0(work_path, "/data/pep_kegg_exp.RData"))
