# Install phyloseq in R environment from Bioconductor
#renv::install("bioc::phyloseq")

# Load necessary libraries
library("phyloseq")
library("dplyr")

# Load Unipept LCA output
data_path <- paste0(getwd(),"/data")
for (file in list.files(data_path) {
  if (substring(file, 1, 12) == "unipept_lca_") {
    load(paste0(work_path, "/data/", file))
  }
}

# ==================================
# Create OTU table
# ==================================

# List of LCA dataframes
unipept_dfs <- ls(pattern = "unipept_lca_")

# Create an empty list to store counts of peptides in each sample
otu_list <- list()

# Loop over the dataframes
for (i in 1:length(unipept_dfs)) {
  # Get replicate and dataframe from sample i
  replicate <- substr(unipept_dfs[i], 13, 15)
  df <- get(unipept_dfs[i])
  
  # Count the occurrences of each peptide in the sample
  peptide_counts <- table(df$peptide)
  
  if (i == 1) {
    otumat <- as.matrix(peptide_counts)  # Create otumat from peptide_counts
    colnames(otumat)[i] <- replicate  # Rename column with replicate number
  }
  else {
    new_col <- rep(NA, nrow(otumat))  # Create a new column with NA values
    otumat <- cbind(otumat, new_col)  # Add the new column to the matrix
    colnames(otumat)[ncol(otumat)] <- replicate  # Set the column name
    
    # Loop through each peptide and update the matrix
    for (peptide in names(peptide_counts)) {
      if (peptide %in% rownames(otumat)) {
        # Assign count to the peptide row
        otumat[peptide, replicate] <- peptide_counts[peptide]
      } else {
        # Add a new row for the new peptide
        new_row <- rep(NA, ncol(otumat))
        otumat <- rbind(otumat, new_row)
        rownames(otumat)[nrow(otumat)] <- peptide
        otumat[peptide, replicate] <- peptide_counts[peptide]
      }
    }
  }
}

# ==================================
# Create taxonomy table
# ==================================
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


# ==================================
# Build Phyloseq object
# ==================================

OTU <- otu_table(otumat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
physeq <- phyloseq(OTU, TAX)

physeq


# ==================================
# Merge sample data
# ==================================

sampledata <- data.frame(metadata_ft$Fragilidad)
rownames(sampledata) <- metadata_ft$replicate
colnames(sampledata)[1] <- "frailty"
sampledata <- sample_data(sampledata)

physeq <- merge_phyloseq(physeq, sampledata)


# ==================================
# 
# ==================================

plot_bar(physeq, x = "frailty", fill = "superkingdom")
plot_heatmap(physeq, taxa.label="superkingdom")
