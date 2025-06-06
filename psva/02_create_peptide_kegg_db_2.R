# Load libraries
library(data.table)
library(dplyr)
library(tidyr)

# Working directory
work_path <- getwd()

# UHGP annotation file path
#uhgp_path <- "/mnt/d/proteomica/fragilidad/Pepfunk/database/uhgp-90"

# Read UHGP annotation file
uhgp_ann <- fread(file.path(paste0(work_path, "/psva/uhgp-90_eggNOG.tsv")),
                  sep = "\t", header = TRUE) %>%
  select(1, 12) %>%
  filter(KEGG_ko != "-")  # Filter UHGP proteins without KEGG KO
colnames(uhgp_ann) <- c("Protein", "KO")
nrow(uhgp_ann)
# 4851159

# Separate rows with more than one KO in several rows
uhgp_ann <- uhgp_ann %>%
  separate_rows(KO, sep = ",")
nrow(uhgp_ann)
# 5379202

# Read uhgp_pep
Pep2Prot <- fread(file.path(work_path, "/psva/uhgp_pep.tsv"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)
nrow(Pep2Prot)
# 1027319075

# Read diamond_uhgp_out
Prot2UR <- fread(file.path(work_path, "/psva/diamond_uhgp_best_hits_out.tsv"),
                 header = FALSE) %>%
  select(1, 3) %>%
  rename(Protein = V1, ID = V3)
nrow(Prot2UR)
# 12317685

# Join by "Protein" to get UniRef90
Pep2Prot2UR <- left_join(Pep2Prot, Prot2UR, by = "Protein",
                         relationship = "many-to-many")
nrow(Pep2Prot2UR)
# 1027319075

# Join by "Protein" to get KO
Pep2Prot2UR2KEGG <- left_join(Pep2Prot2UR, uhgp_ann, by = "Protein",
                              relationship = "many-to-many")
nrow(Pep2Prot2UR2KEGG)
# 1079321301

# Numer missing KO
sum(is.na(Pep2Prot2UR2KEGG$KO))
# 601822433
# Number proteins with KO
nrow(unique(Pep2Prot2UR2KEGG[!is.na(Pep2Prot2UR2KEGG$KO), "Protein"]))
# 4850365
# Number proteins without KO
nrow(unique(Pep2Prot2UR2KEGG[is.na(Pep2Prot2UR2KEGG$KO), "Protein"]))
# 8946788

# Number missing UniRef90
sum(is.na(Pep2Prot2UR2KEGG$ID))
# 54138006
# Number proteins with UniRef90
nrow(unique(Pep2Prot2UR2KEGG[!is.na(Pep2Prot2UR2KEGG$ID), "Protein"]))
# 12310128
# Number proteins without UniRef90
nrow(unique(Pep2Prot2UR2KEGG[is.na(Pep2Prot2UR2KEGG$ID), "Protein"]))
# 1487025

# Number proteins with KO but not UniRef90
nrow(unique(Pep2Prot2UR2KEGG[!is.na(Pep2Prot2UR2KEGG$KO) &
                               is.na(Pep2Prot2UR2KEGG$ID), "Protein"]))
# 344
# Number of proteins with KO and UniRef90
nrow(unique(Pep2Prot2UR2KEGG[!is.na(Pep2Prot2UR2KEGG$KO) &
                               !is.na(Pep2Prot2UR2KEGG$ID), "Protein"]))
# 4850021

# Filter peptides without KO annotation
Pep2Prot2UR2KEGG <- Pep2Prot2UR2KEGG %>% filter(!is.na(KO))
nrow(Pep2Prot2UR2KEGG)
# 477498868
nrow(unique(Pep2Prot2UR2KEGG[, "Protein"]))
# 4850365

# Save peptide-protein-UniRef90-KEGG dataset as csv file
write.csv(Pep2Prot2UR2KEGG, file = paste0(work_path,
                                          "/psva/core_pep_kegg_db_2.csv"))
print("core_pep_kegg_db_2.csv saved.")