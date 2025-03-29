# ==================================
# Peptide enrichment analysis
# ==================================
#renv::install("DESeq2")
#renv::install("GSVA")
# Load libraries
library(DESeq2)
library(GSVA)
library(limma)
library(dplyr)
# Load functions
source("filter_valids.R")


# Working directory
work_path <- getwd()

# Load kegg database
load(paste0(work_path, "/data/kegg_ko_path.RData"))

# Load core peptide db
core_pep_kegg <- read.delim2(paste0(work_path,"/data/core_pep_kegg_db.csv"), 
                             sep=",", header=F, col.names = c("pep", "kegg",
                                                              "count", "eval"))

# Get information on number of proteins for a peptide
core_pep_kegg <- core_pep_kegg %>% group_by(pep) %>% 
  summarize(total = sum(count))  %>%
  merge(., core_pep_kegg, by='pep', all.y=T) %>%
  mutate(prop = count/total) %>% select(pep, kegg, prop)

# Update pep if it is a duplicate
core_pep_kegg$newpep_name <- make.names(core_pep_kegg$pep, unique=T)