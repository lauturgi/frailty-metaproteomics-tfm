# Load libraries
library(plyr)
library(dplyr)

# Working directory
work_path <- getwd()

# Load kegg database
load(paste0(work_path, "/psva/data/kegg_ko_l3_4.RData"))

# Load count_core_pep_ko_db_2
load(paste0(work_path, "/psva/data/count_core_pep_ko_db_2.RData"))

# Group ko by pathways
l4_3 <- dlply(ko_df %>% select(L3_DESC, L4),
                      .(L3_DESC))
# Number ko sets
length(l4_3)
# 560

# Get peptide sets joining by ko
kegg_pepsets <- lapply(l4_3, function(df) {
  subset <- count_core_pep_ko_db[count_core_pep_ko_db$ko %in% df$L4, ]
  subset$newpep_name
})
length(kegg_pepsets)
# 560

# Number of non-empty ko sets
sum(sapply(kegg_pepsets, function(x) length(x) > 0))
# 538

# Save kegg_pepsets
save(kegg_pepsets, file = paste0(work_path, "/psva/data/kegg_pepsets.RData"))
