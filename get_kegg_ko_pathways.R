# Load libraries
library(KEGGREST)
library(tidyr)

# ==================================
# Build KEGG database
# ==================================
# 1. KEGG pathways
# ==================================
# Working directory
work_path <- getwd()

# Get ko functional orthologs
kegg_ko <- as.data.frame(keggList("ko"))
kegg_ko$KEGG <- row.names(kegg_ko)
row.names(kegg_ko) <- NULL
kegg_ko <- kegg_ko[, c(2,1)]
colnames(kegg_ko)[2] <- "KEGG_DESCR"

# Get pathways
kegg_path <- as.data.frame(keggList("pathway"))
kegg_path$PATHWAY <- row.names(kegg_path)
row.names(kegg_path) <- NULL
kegg_path <- kegg_path[, c(2,1)]
colnames(kegg_path)[2] <- "PATHWAY_DESCR"

# Get link between ko and pathways
kegg_ko_path <- keggLink("pathway", "ko")
kegg_ko_path <- data.frame(
  KEGG = names(kegg_ko_path),
  PATHWAY = unname(kegg_ko_path)
)

# Get kegg term (ko:K00001 --> K00001)
kegg_ko_path$KEGG <- sapply(strsplit(kegg_ko_path$KEGG, ":"), `[`, 2)

# Each kegg term has two entries for the same pathway: "path:map" and "path:ko" 
# Keep only one and the corresponding ID
kegg_ko_path <- kegg_ko_path[grepl("^path:ko", kegg_ko_path$PATHWAY), ]
kegg_ko_path$PATHWAY <- gsub("^path:ko", "", kegg_ko_path$PATHWAY)

# Get kegg description
kegg_ko_path <- merge(kegg_ko_path, kegg_ko, by = "KEGG")

# Get pathway description
kegg_path$PATHWAY <- gsub("^map", "", kegg_path$PATHWAY)
kegg_ko_path <- merge(kegg_ko_path, kegg_path, by = "PATHWAY")
kegg_ko_path <- kegg_ko_path[, c("KEGG", "KEGG_DESCR", "PATHWAY",
                                 "PATHWAY_DESCR")]  # Reorder columns
kegg_ko_path$PATHWAY_DESCR <- paste0(kegg_ko_path$PATHWAY_DESCR, " [PATH:ko",
                                     kegg_ko_path$PATHWAY , "]")  # from path

# ==================================
# 2. KEGG bride
# ==================================
# Get bride functional hierarchies
kegg_brite <- as.data.frame(keggList("brite", "ko"))
kegg_brite$PATHWAY <- row.names(kegg_brite)
kegg_brite$PATHWAY <- gsub("^ko", "", kegg_brite$PATHWAY)
row.names(kegg_brite) <- NULL
kegg_brite <- kegg_brite[, c(2,1)]
colnames(kegg_brite)[2] <- "PATHWAY_DESCR"

# Get link between ko and bride
kegg_ko_brite <- keggLink("brite", "ko")
kegg_ko_brite <- data.frame(
  KEGG = names(kegg_ko_brite),
  PATHWAY = unname(kegg_ko_brite)
)

# Keep only ko brite hierarchies (br:ko)
kegg_ko_brite <- kegg_ko_brite[grepl("^br:ko", kegg_ko_brite$PATHWAY), ]
kegg_ko_brite$KEGG <- gsub("^ko:", "", kegg_ko_brite$KEGG)
kegg_ko_brite$PATHWAY <- gsub("^br:ko", "", kegg_ko_brite$PATHWAY)

# Remove orthologs and modules
kegg_brite <- kegg_brite[!(kegg_brite$PATHWAY %in%
                             c("00001", "00002", "00003")), ]

# Get bride description
kegg_ko_brite <- merge(kegg_ko_brite, kegg_brite, by = "PATHWAY")

# Get kegg description
kegg_ko_brite <- merge(kegg_ko_brite, kegg_ko, by = "KEGG")
kegg_ko_brite$PATHWAY_DESCR <- paste0(kegg_ko_brite$PATHWAY_DESCR, " [BR:ko",
                                      kegg_ko_brite$PATHWAY , "]")  # from bride
kegg_ko_brite <- kegg_ko_brite[, c("KEGG", "KEGG_DESCR", "PATHWAY",
                                   "PATHWAY_DESCR")] # reorder columns

# ==================================
# 3. Combine br and path
# ==================================
kegg_ko_path <- rbind(kegg_ko_path, kegg_ko_brite)

save_path <- paste0(work_path,"/data/kegg_ko_path.RData")
save(kegg_ko_path, file = save_path)

# ==================================
# KEGG pathways from FISABIO
# ==================================
# Load KEGG pathways
kegg_path <- "/mnt/d/proteomica/fragilidad/Pepfunk/kegg/"
kegg_ko_path_2 <- paste0(kegg_path, "kegg.pathways.tsv")
kegg_ko_path_2 <- read.csv2(kegg_ko_path_2, sep = "\t")

## There are several pathways for each kegg term separated by "|"

# Separate in multiple rows by "|"
kegg_ko_path_2 <- kegg_ko_path_2 %>%
  separate_rows(PATHWAYS, PATHWAYS_DESCRIPTION, sep = "\\|")

# Check 
head(kegg_ko_path_2)

# Compare with kegg_L3.txt (old) from pepfunk
pepfunk_path <- "/mnt/d/proteomica/fragilidad/Pepfunk/database/"
kegg_old <- read.csv2(paste0(pepfunk_path,"kegg_L3.txt"), sep = "\t", 
                      header = FALSE, colClasses = "character")
colnames(kegg_old) <- c("KEGG", "KEGG_DESCRIPTION", "PATHWAYS",
                        "PATHWAYS_DESCRIPTION")

kegg_check <- merge(kegg_ko_path_2, kegg_old, by.x = c("KEGG", "PATHWAYS"), 
                    by.y = c("KEGG", "PATHWAYS"), all = TRUE)

kegg_check <- kegg_check[is.na(kegg_check$PATHWAYS_DESCRIPTION.x),]
