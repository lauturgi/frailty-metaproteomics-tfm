# Load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


# ==================================
# KEGG annotation
# ==================================
# Select proteins with KEGG annotation
# ==================================
# Working directory
work_path <- getwd()

# Load KEGG_UniRef90.idmapping.dat
UR2KEGG <- fread(file.path(work_path, "/KEGG_UniRef90.idmapping.dat"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

sum(lengths(UR2KEGG$UniRef90) > 1)
# 0

sum(lengths(UR2KEGG$KEGG) > 1)
# 357028

# Unnest KEGG
UR2KEGG <- unnest(UR2KEGG, c(UniRef90, KEGG))
nrow(UR2KEGG)
# 24617890

# Read uhgp_pep
#Pep2Prot <- fread(file.path(work_path, "/uhgp_pep.tsv"), header = FALSE) %>%
#  rename(Peptide = V1, Protein = V2)
#nrow(Pep2Prot)
# 1027319075

# Read uhgp_pep fractions
Pep2Protaa <- fread(file.path(work_path, "uhgp_pepaa"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)
Pep2Protab <- fread(file.path(work_path, "uhgp_pepab"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)
Pep2Protac <- fread(file.path(work_path, "uhgp_pepac"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)
Pep2Protad <- fread(file.path(work_path, "uhgp_pepad"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)
Pep2Protae <- fread(file.path(work_path, "uhgp_pepae"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)
Pep2Protaf <- fread(file.path(work_path, "uhgp_pepaf"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)
Pep2Protag <- fread(file.path(work_path, "uhgp_pepag"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)
Pep2Protah <- fread(file.path(work_path, "uhgp_pepah"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)
Pep2Protai <- fread(file.path(work_path, "uhgp_pepai"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)
Pep2Protaj <- fread(file.path(work_path, "uhgp_pepaj"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)
Pep2Protak <- fread(file.path(work_path, "uhgp_pepak"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)

# Read diamond_uhgp_out
Prot2UR <- fread(file.path(work_path, "/diamond_uhgp_best_hits_out.tsv"),
                      header = FALSE) %>%
  select(1, 3) %>%
  rename(Protein = V1, ID = V3)
nrow(Prot2UR)
# 12317685

# Join by values of the "Protein"
#Pep2Prot2UR <- left_join(Pep2Prot, Prot2UR, by = "Protein",
#                           relationship = "many-to-many")
#nrow(Pep2Prot2UR)
# 1027319075

Pep2Prot2URaa <- left_join(Pep2Protaa, Prot2UR, by = "Protein",
                           relationship = "many-to-many")
Pep2Prot2URab <- left_join(Pep2Protab, Prot2UR, by = "Protein",
                           relationship = "many-to-many")
Pep2Prot2URac <- left_join(Pep2Protac, Prot2UR, by = "Protein",
                           relationship = "many-to-many")
Pep2Prot2URad <- left_join(Pep2Protad, Prot2UR, by = "Protein",
                           relationship = "many-to-many")
Pep2Prot2URae <- left_join(Pep2Protae, Prot2UR, by = "Protein",
                           relationship = "many-to-many")
Pep2Prot2URaf <- left_join(Pep2Protaf, Prot2UR, by = "Protein",
                           relationship = "many-to-many")
Pep2Prot2URag <- left_join(Pep2Protag, Prot2UR, by = "Protein",
                           relationship = "many-to-many")
Pep2Prot2URah <- left_join(Pep2Protah, Prot2UR, by = "Protein",
                           relationship = "many-to-many")
Pep2Prot2URai <- left_join(Pep2Protai, Prot2UR, by = "Protein",
                           relationship = "many-to-many")
Pep2Prot2URaj <- left_join(Pep2Protaj, Prot2UR, by = "Protein",
                           relationship = "many-to-many")
Pep2Prot2URak <- left_join(Pep2Protak, Prot2UR, by = "Protein",
                           relationship = "many-to-many")
# Join by "UniRef90"
#Pep2Prot2UR2KEGG <- left_join(Pep2Prot2UR, UR2KEGG,  by = c("ID" = "UniRef90"),
#                              relationship = "many-to-many")
#nrow(Pep2Prot2UR2KEGG)
# 1032315537

Pep2Prot2UR2KEGGaa <- left_join(Pep2Prot2URaa, UR2KEGG,
                                by = c("ID" = "UniRef90"),
                                relationship = "many-to-many")
Pep2Prot2UR2KEGGab <- left_join(Pep2Prot2URab, UR2KEGG,
                                by = c("ID" = "UniRef90"),
                                relationship = "many-to-many")
Pep2Prot2UR2KEGGac <- left_join(Pep2Prot2URac, UR2KEGG,
                                by = c("ID" = "UniRef90"),
                                relationship = "many-to-many")
Pep2Prot2UR2KEGGad <- left_join(Pep2Prot2URad, UR2KEGG,
                                by = c("ID" = "UniRef90"),
                                relationship = "many-to-many")
Pep2Prot2UR2KEGGae <- left_join(Pep2Prot2URae, UR2KEGG,
                                by = c("ID" = "UniRef90"),
                                relationship = "many-to-many")
Pep2Prot2UR2KEGGaf <- left_join(Pep2Prot2URaf, UR2KEGG,
                                by = c("ID" = "UniRef90"),
                                relationship = "many-to-many")
Pep2Prot2UR2KEGGag <- left_join(Pep2Prot2URag, UR2KEGG,
                                by = c("ID" = "UniRef90"),
                                relationship = "many-to-many")
Pep2Prot2UR2KEGGah <- left_join(Pep2Prot2URah, UR2KEGG,
                                by = c("ID" = "UniRef90"),
                                relationship = "many-to-many")
Pep2Prot2UR2KEGGai <- left_join(Pep2Prot2URai, UR2KEGG,
                                by = c("ID" = "UniRef90"),
                                relationship = "many-to-many")
Pep2Prot2UR2KEGGaj <- left_join(Pep2Prot2URaj, UR2KEGG,
                                by = c("ID" = "UniRef90"),
                                relationship = "many-to-many")
Pep2Prot2UR2KEGGak <- left_join(Pep2Prot2URak, UR2KEGG,
                                by = c("ID" = "UniRef90"),
                                relationship = "many-to-many")

# Filter peptides without KEGG annotation
#Pep2Prot2UR2KEGG <- Pep2Prot2UR2KEGG %>% filter(!is.na(KEGG))
#nrow(Pep2Prot2UR2KEGG)
# 73998314
Pep2Prot2UR2KEGGaa <- Pep2Prot2UR2KEGGaa %>% filter(!is.na(KEGG))
Pep2Prot2UR2KEGGab <- Pep2Prot2UR2KEGGab %>% filter(!is.na(KEGG))
Pep2Prot2UR2KEGGac <- Pep2Prot2UR2KEGGac %>% filter(!is.na(KEGG))
Pep2Prot2UR2KEGGad <- Pep2Prot2UR2KEGGad %>% filter(!is.na(KEGG))
Pep2Prot2UR2KEGGae <- Pep2Prot2UR2KEGGae %>% filter(!is.na(KEGG))
Pep2Prot2UR2KEGGaf <- Pep2Prot2UR2KEGGaf %>% filter(!is.na(KEGG))
Pep2Prot2UR2KEGGag <- Pep2Prot2UR2KEGGag %>% filter(!is.na(KEGG))
Pep2Prot2UR2KEGGah <- Pep2Prot2UR2KEGGah %>% filter(!is.na(KEGG))
Pep2Prot2UR2KEGGai <- Pep2Prot2UR2KEGGai %>% filter(!is.na(KEGG))
Pep2Prot2UR2KEGGaj <- Pep2Prot2UR2KEGGaj %>% filter(!is.na(KEGG))
Pep2Prot2UR2KEGGak <- Pep2Prot2UR2KEGGak %>% filter(!is.na(KEGG))

# Bind fractions
Pep2Prot2UR2KEGG <- bind_rows(Pep2Prot2UR2KEGGaa, Pep2Prot2UR2KEGGab,
                              Pep2Prot2UR2KEGGac, Pep2Prot2UR2KEGGad, 
                              Pep2Prot2UR2KEGGae, Pep2Prot2UR2KEGGaf, 
                              Pep2Prot2UR2KEGGag, Pep2Prot2UR2KEGGah, 
                              Pep2Prot2UR2KEGGai, Pep2Prot2UR2KEGGaj, 
                              Pep2Prot2UR2KEGGak)

Pep2Prot2UR2KEGG <- select(Pep2Prot2UR2KEGG, -ID) #get rid of ID column

# Save peptide-protein-UniRef90-KEGG dataset as csv file
write.csv(Pep2Prot2UR2KEGG, file = paste0(work_path, "/core_pep_kegg_db.csv"))
