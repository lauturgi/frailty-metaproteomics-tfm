# Load libraries
library(dplyr)
library(tibble)

# Working directory
work_path <- getwd()

# Choose option:
# A: MaxLFQ intensities
# B: Raw intensities

psva_option <- "A"
if (psva_option == "A"){
  data_path <- paste0(work_path, "/psva/data/maxlfq/")
  load(paste0(data_path, "pep_exp.RData"))
  norm_pep_exp <- pep_exp
} else if (de_option == "B"){
  data_path <- paste0(work_path, "/psva/data/lfq/")
  load(paste0(data_path, "norm_pep_exp.RData"))
} else {
  stop("No valid option.")
}

nrow(norm_pep_exp)
# A: 699
# B: 699

# Load count_core_pep_ko_db_2
load(paste0(work_path, "/psva/data/count_core_pep_ko_db_2.RData"))
nrow(count_core_pep_ko_db)
# 330446378

# Read old core peptide db
#core_pep_kegg_old <- read.delim2(paste0("/mnt/d/proteomica/fragilidad/Pepfunk/",
#                                        "database/core_pep_kegg_db.csv"),
#                                 sep = ",", header=F) %>%
#  setNames(c("pep", "kegg", "count"))


# Join peptides and core_pep_kegg by "pep" to get ko, newpep_name and prop
matched_core_pep_ko <- norm_pep_exp %>%
  rownames_to_column(var="pep") %>%
  merge(., count_core_pep_ko_db, by="pep")
nrow(matched_core_pep_ko)
# A: 1139
# B: 1139

# All peptides from norm_pep_exp
all_pep <- rownames_to_column(norm_pep_exp, var = "pep")

# Get those not in matched_core_pep_ko$pep
unmatched_pep <- all_pep[!all_pep$pep %in% matched_core_pep_ko$pep, ]

# Save unmatched_pep
save(unmatched_pep, file = paste0(data_path, "unmatched_pep.RData"))

# Normalise intensities by proportion of functional annotation
matched_core_pep_ko[, grep("Intensity", names(matched_core_pep_ko))] <- 
  lapply(matched_core_pep_ko[, grep("Intensity",
                                      names(matched_core_pep_ko))], 
         function(x) log2(1 + (x * matched_core_pep_ko$prop)))

# Save matched_core_pep_ko
save(matched_core_pep_ko, file = paste0(data_path, "matched_core_pep_ko.RData"))

# Set rownames to "newpep_name" and select keep "Intensity" columns
pep_kegg_exp <- matched_core_pep_ko %>%
  column_to_rownames(var="newpep_name") %>%
  select(contains("Intensity"))
nrow(pep_kegg_exp)
# A: 1139
# B: 1139

# Remove ".Intensity" from column names
colnames(pep_kegg_exp) <- sapply(strsplit(colnames(pep_kegg_exp), 
                                           "\\."), `[`, 1)

# Save pep_kegg_exp
save(pep_kegg_exp, file = paste0(data_path, "pep_kegg_exp.RData"))
print("pep_kegg_exp.RData saved.")
