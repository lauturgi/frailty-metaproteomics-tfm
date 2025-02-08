# Install phyloseq in R environment from Bioconductor
#renv::install("bioc::phyloseq")

# Load necessary libraries
library("phyloseq")
library("dplyr")

# ==================================
# Build Phyloseq object
# ==================================

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
  replicate <- unipept_dfs[i]
  df <- get(unipept_dfs[i])
  
  # Count the occurrences of each peptide in the sample
  peptide_counts <- table(df$peptide)
  
  if (i == 1) {
    # Create OTU matrix from peptide counts
    otu_table <- as.matrix(peptide_counts)
    # Rename column with replicate number
    colnames(otu_table)[i] <- replicate
  }
  else {
    new_col <- rep(NA, nrow(otu_table))  # Create a new column with NA values
    otu_table <- cbind(otu_table, new_col)  # Add the new column to the matrix
    colnames(otu_table)[ncol(otu_table)] <- replicate  # Set the column name
    
    # Loop through each peptide and update the matrix
    for (peptide in names(peptide_counts)) {
      if (peptide %in% rownames(otu_table)) {
        otu_table[peptide, replicate] <- peptide_counts[peptide]
      } else {
        # Add a new row for the new peptide
        new_row <- rep(NA, ncol(otu_table))
        otu_table <- rbind(otu_table, new_row)
        rownames(otu_table)[nrow(otu_table)] <- peptide
        otu_table[peptide, replicate] <- peptide_counts[peptide]
      }
    }
  }
}

# ==================================
# Create taxonomy table
# ==================================
# Create an empty matrix with 0 rows and 7 columns for taxonomy
tax_table <- matrix(ncol = 7, nrow = 0)

# Assign column names
colnames(tax_table) <- c("superkingdom", "phylum", "class", "order", "family",
                         "genus", "species")

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
    tax_table <- rbind(tax_table, new_rows[, c("superkingdom_name",
                                               "phylum_name", "class_name",
                                               "order_name", "family_name",
                                               "genus_name", "species_name")])
    rownames(tax_table) <- new_rows$peptide
  }
  else {
    for (j in 1:nrow(df)) {
      if (df[j, "peptide"] %in% rownames(tax_table)) {
        next
      }
      else {
        new_row <- df[j, c("superkingdom_name", "phylum_name", "class_name",
                           "order_name", "family_name", "genus_name",
                           "species_name")]
        tax_table <- rbind(tax_table, new_row)
        rownames(tax_table)[nrow(tax_table)] <- df[j, "peptide"]
      }
    }
  }
}
  
  