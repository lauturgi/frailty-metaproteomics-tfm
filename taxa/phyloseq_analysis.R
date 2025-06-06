# Install phyloseq in R environment from Bioconductor
#renv::install("bioc::phyloseq")

# Load necessary libraries
library("phyloseq")
library("dplyr")
library("tibble")
library("DESeq2")
library("ggplot2")

# Working directory
work_path <- getwd()

# ==================================
# Create taxonomy table
# ==================================

# Load Unipept LCA output
unipept_path <- paste0(work_path,"/data/unipept_lca/")
list.files(unipept_path)
for (file in list.files(unipept_path)) {
  if (substr(file, 1, 19) == "unipept_lca_no_dup_") {
    load(paste0(unipept_path, file))
  }
}

# List of LCA dataframes
unipept_dfs <- ls(pattern = "unipept_lca_")

# Create an empty matrix with 0 rows and 7 columns for taxonomy
taxmat <- matrix(ncol = 7, nrow = 0)

# Loop over the dataframes
for (i in 1:length(unipept_dfs)) {
  replicate <- unipept_dfs[i]
  df <- get(unipept_dfs[i])
  if (i == 1) {
    new_rows <- df %>%
      group_by(peptide, superkingdom_name, phylum_name, class_name, order_name,
               family_name, genus_name, species_name) %>%
      summarize(count = n())
    new_rows <- as.data.frame(new_rows)
    taxmat <- rbind(taxmat, new_rows[, c("superkingdom_name",
                                         "phylum_name", "class_name",
                                         "order_name", "family_name",
                                         "genus_name", "species_name")])
    rownames(taxmat) <- new_rows$peptide
  }
  else {
    for (j in 1:nrow(df)) {
      if (df[j, "peptide"] %in% rownames(taxmat)) {
        next
      }
      else {
        new_row <- df[j, c("superkingdom_name", "phylum_name", "class_name",
                           "order_name", "family_name", "genus_name",
                           "species_name")]
        taxmat <- rbind(taxmat, new_row)
        rownames(taxmat)[nrow(taxmat)] <- df[j, "peptide"]
      }
    }
  }
}

# Replace empty strings with NA
taxmat <- apply(taxmat, c(1, 2), function(x) ifelse(x == "", NA, x))

# Assign column names
colnames(taxmat) <- c("superkingdom", "phylum", "class", "order", "family",
                      "genus", "species")

# Filter peptides without rank
#taxmat <- taxmat[!is.na(taxmat[, "superkingdom"]), ]

# ==================================
# Create OTU table
# ==================================

# Loop over the dataframes
#for (i in 1:length(unipept_dfs)) {
#  # Get replicate and dataframe from sample i
#  replicate <- substr(unipept_dfs[i], 13, 15)
#  df <- get(unipept_dfs[i])
#  
#  # Count the occurrences of each peptide in the sample
#  peptide_counts <- table(df$peptide)
#  
#  if (i == 1) {
#    otumat <- as.matrix(peptide_counts)  # Create otumat from peptide_counts
#    colnames(otumat)[i] <- replicate  # Rename column with replicate number
#  }
#  else {
#    new_col <- rep(NA, nrow(otumat))  # Create a new column with NA values
#    otumat <- cbind(otumat, new_col)  # Add the new column to the matrix
#    colnames(otumat)[ncol(otumat)] <- replicate  # Set the column name
#    
#    # Loop through each peptide and update the matrix
#    for (peptide in names(peptide_counts)) {
#      if (peptide %in% rownames(otumat)) {
#        # Assign count to the peptide row
#        otumat[peptide, replicate] <- peptide_counts[peptide]
#      } else {
#        # Add a new row for the new peptide
#        new_row <- rep(NA, ncol(otumat))
#        otumat <- rbind(otumat, new_row)
#        rownames(otumat)[nrow(otumat)] <- peptide
#        otumat[peptide, replicate] <- peptide_counts[peptide]
#      }
#    }
#  }
#}

# Load metadata
load(paste0(work_path, "/data/metadata.RData"))

# Path to combined_peptide.tsv
fragpipe_path <- paste0("/mnt/d/proteomica/fragilidad/datos/ProteinIdent/",
                        "MSfraggerSPbacteria/DatosProtIdentCombinados/")
# Load combined_peptide.tsv
peptide_file <- paste0(fragpipe_path,"combined_peptide.tsv")
peptides  <- read.delim(peptide_file, row.names = 1) %>% as.data.frame() %>% 
  select(contains("Spectral")) %>%  # Spectral columns
  select(contains(metadata_ft$replicate))

# Change colnames to include the group each replicate belongs to
replicates_ft <- metadata_ft[metadata_ft$Fragilidad == "FT", "replicate"]
replicates_nft <- metadata_ft[metadata_ft$Fragilidad == "NFT", "replicate"]
for (i in seq_along(colnames(peptides))) {
  colname <- colnames(peptides)[i]  # Get column name
  
  if (substr(colname, 5, 7) %in% replicates_ft) {
    colnames(peptides)[i] <- paste0("FT", substr(colname, 4, 7))
  } else {
    colnames(peptides)[i] <- paste0("NFT", substr(colname, 4, 7))
  }
}

# Create otumat from peptide intensities
otumat <- as.matrix(peptides)  # Create otumat from peptide_counts

# "I" > "L" in peptides
rownames(otumat) <- gsub("I", "L", rownames(otumat))

# Match peptides from taxmat and otumat
otus <- list()

# Compare rownames from taxmat and otumat
for (otu in rownames(taxmat)) {
  # Search all rownames in otumat that contains otu
  match_peptides <- grep(otu, rownames(otumat), value = TRUE)
  
  # Keep in dictionary
  if (length(match_peptides) > 0) {
    otus[[otu]] <- match_peptides
  }
}

length(otus)
sum(sapply(otus, length) > 1)

###############
# DOBLE-CHECK #  
###############
# Sum peptide intensities that belong to the same OTU
# Convert 'otus' into a data frame 
otus <- stack(otus) %>% rename(peptide = values, otu = ind)

# Convert 'otumat' to a data frame
otumat <- as.data.frame(otumat) 
otumat$peptide <- rownames(otumat)  # Add peptide names as a new column

# Merge otumat with otus by peptide
otumat <- inner_join(otumat, otus, by = "peptide")

# Group by otu and sum values
otumat <- otumat %>%
  select(-peptide) %>%  # Remove peptide column
  group_by(otu) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)))  # Sum

# Convert  back to matrix with otu as rownames
otumat <- otumat %>%
  column_to_rownames(var = "otu") %>%
  as.matrix()


# ==================================
# Build Phyloseq object
# ==================================

OTU <- otu_table(otumat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)


# ==================================
# Merge sample data
# ==================================
sampledata <- data.frame(metadata_ft$Fragilidad)
rownames(sampledata) <- metadata_ft$replicate
for (i in seq_along(rownames(sampledata))) {
  rowname <- rownames(sampledata)[i]  # Get row name
  
  if (rowname %in% replicates_ft) {
    rownames(sampledata)[i] <- paste0("FT_", rowname)
  } else {
    rownames(sampledata)[i] <- paste0("NFT_", rowname)
  }
}
colnames(sampledata)[1] <- "frailty"
sampledata <- sample_data(sampledata)


physeq <- phyloseq(OTU, TAX, sampledata)


sample_data(physeq)
tax_table(physeq)
sample_variables(physeq)

# !! Quedarme con el taxon que al menos tiene 3 peptidos representados
# ==================================
# Normalise by DESeq2
# ==================================
diagdds = phyloseq_to_deseq2(kostic, ~ DIAGNOSIS)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
# ==================================
# 
# ==================================

physeq.bact <- subset_taxa(physeq, superkingdom == "Bacteria")
physeq.bact.glom = tax_glom(physeq.bact, "phylum")
plot_bar(physeq.bact.glom, fill = "phylum")
plot_bar(physeq, 
         x = "frailty", 
         fill = "superkingdom") + 
  theme(axis.text.x = element_text(size = 4))
plot_heatmap(physeq, taxa.label="superkingdom")
