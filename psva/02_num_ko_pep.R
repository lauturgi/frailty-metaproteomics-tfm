library(data.table)
library(dplyr)
library(tidyr)

work_path <- getwd()

Pep2Prot2UR2KEGG <- fread(paste0(work_path,"/core_pep_kegg_db_2.csv"),
                          sep = ",", header = TRUE) %>%
  select(Peptide, Protein, ID, KO)

count_ko_pep <- Pep2Prot2UR2KEGG %>%
  group_by(Peptide, KO) %>%
  summarise(KO_count = n())

save(count_ko_pep, file = paste0(work_path, "/data/count_ko_pep.RData"))

total_ko_pep <- count_ko_pep %>%
  group_by(Peptide) %>%
  summarise(Total_KO_count = sum(KO_count))

save(total_ko_pep, file = paste0(work_path, "/data/total_ko_pep.RData"))

total_ko_perc <- total_ko_pep %>%
  count(Total_KO_count) %>%                     
  mutate(percent = 100 * n / sum(n))

save(total_ko_perc, file = paste0(work_path, "/psva/data/total_ko_perc.RData"))

dist_ko_pep <- count_ko_pep %>%
  group_by(Peptide) %>%
  summarise(KO_dist = n_distinct(KO))

save(dist_ko_pep, file = paste0(work_path, "/data/dist_ko_pep.RData"))

