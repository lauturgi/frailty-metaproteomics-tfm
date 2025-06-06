# Load libraries
library(KEGGREST)
library(dplyr)

# ==================================
# KO from KEGG IDs
# ==================================
# Working directory
work_path <- getwd()

# Read unique KEGG
kegg <- read.csv(paste0(work_path, "/psva/unique_kegg.csv")) %>%
  select(-X)

# Get unique organisms
kegg$org <- sapply(strsplit(kegg$KEGG, ":"), `[`, 1)
orgs <- unique(kegg$org)

# List to add data.frames
kegg_ko_list <- list()

# For each organism, get ko associated to each gene
for (org in orgs) {
  Sys.sleep(1)
  res <- keggLink("ko", org)
  kegg_ids <- names(res)  # Get KEGG
  ko_ids <- unname(res)  # Get KO
  df <- data.frame(KEGG = kegg_ids, KO = ko_ids)
  kegg_ko_list[[length(kegg_ko_list) + 1]] <- df
}

# Combine all data.frames
kegg_ko <- do.call(rbind, kegg_ko_list)

save(kegg_ko, file = paste0(work_path, "/psva/data/kegg_ko.RData"))