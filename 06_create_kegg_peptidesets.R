# Load libraries
library(plyr)
library(dplyr)

# Working directory
work_path <- getwd()

# Load kegg database
load(paste0(work_path, "/data/kegg_ko_l3_4.RData"))

# Load core_pep_kegg
load(paste0(work_path, "/data/core_pep_kegg.RData"))

# Unique pathways vector
# l3 <- ko_df$L3_DESC %>% unique()

# Group ko by pathways
l4_3 <- dlply(ko_df %>% select(L3_DESC, L4),
                      .(L3_DESC))

# Get peptide sets joining kegg and pathway_kegg by kegg
kegg_pepsets <- lapply(l4_3, function(df) {
  subset <- core_pep_kegg[core_pep_kegg$kegg %in% df$L4, ]
  subset$newpep_name
})

# Save kegg_pepsets
save(kegg_pepsets, file = paste0(work_path, "/data/kegg_pepsets.RData"))
