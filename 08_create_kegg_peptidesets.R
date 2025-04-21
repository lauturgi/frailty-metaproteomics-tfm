# Load libraries
library(plyr)
library(dplyr)

# Working directory
work_path <- getwd()

# Load kegg database
load(paste0(work_path, "/data/kegg_ko_path.RData"))

# Load core_pep_kegg
load(paste0(work_path, "/data/core_pep_kegg.RData"))

# Unique pathways vector
pathways <- kegg_ko_path$PATHWAY_DESCR %>% unique()

# Group ko by pathways
pathway_kegg <- dlply(kegg_ko_path %>% select(PATHWAY_DESCR, KEGG),
                      .(PATHWAY_DESCR))

# Get peptide sets joining kegg and pathway_kegg by kegg
kegg_pepsets <- lapply(pathway_kegg, function(df) {
  subset <- core_pep_kegg[core_pep_kegg$kegg %in% df$KEGG, ]
  subset$newpep_name
})

# Save kegg_pepsets
save(kegg_pepsets, file = paste0(work_path, "/data/kegg_pepsets.RData"))

