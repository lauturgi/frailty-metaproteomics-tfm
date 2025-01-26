# renv::install("cli")
library(FragPipeAnalystR)
library(ggplot2)
library(limma)
library(vsn)
#renv::install("nortest")
library(nortest)
library(stringr)

# ==================================
# 1. Making SummarizedExperiment
# ==================================

# Read combined_protein
combined_protein_df <- read.table(combined_protein,
           header = TRUE,
           fill = TRUE, # to fill any missing data
           sep = "\t",
           quote = "",
           comment.char = "",
           blank.lines.skip = F,
           check.names = F
)

# See column names
colnames(combined_protein_df)

# Replace 0 by NA
combined_protein_df[ , !colnames(combined_protein_df) %in% "Protein ID"] <- lapply(combined_protein_df[ , !colnames(combined_protein_df) %in% "Protein ID"], function(x) replace(x, x == 0, NA))


# Select column "Protein ID" and those that refer to MaxLFQ intensities
protein_id_col <- which(colnames(combined_protein_df) == "Protein ID")
lfq_col <- grep("MaxLFQ", colnames(combined_protein_df))

# Filter "MaxLFQ" columns to include only those matching replicates in metadata_ft$replicate
lfq_col <- lfq_col[sapply(lfq_col, function(idx) {
  any(grepl(paste(metadata_ft$replicate, collapse = "|"), colnames(combined_protein_df)[idx]))
})]

# Combine the indices of "Protein ID" and "MaxLFQ" columns
to_keep <- c(protein_id_col, lfq_col)

# Subset the data frame to keep only the selected columns
combined_protein_filtered_df <- combined_protein_df[, to_keep]

# Rename columns containing "MaxLFQ"
colnames(combined_protein_filtered_df) <- sapply(colnames(combined_protein_filtered_df), function(col) {
  if (grepl("MaxLFQ", col)) {
    strsplit(col, " ")[[1]][1]  # Split by space and keep the first substring
  } else {
    col  # Keep column name unchanged if it doesn't contain "MaxLFQ"
  }
})

# Check selected column names
colnames(combined_protein_filtered_df)

# Get MaxLFQ column names from combined_protein_filtered_df
sample_cols <- grep("_", colnames(combined_protein_filtered_df), value = TRUE)

# Get replicate from each column name
replicate_map <- data.frame(
  replicate = str_extract(sample_cols, "\\d+"),  # Extract numeric part of column names
  sample_name = sample_cols                     # Corresponding full column names
)

# Match replicate order from metadata and get sample_name
metadata_ft$sample_name <- replicate_map$sample_name[match(metadata_ft$replicate, replicate_map$replicate)]

# Save protein MaxLFQ intensities as csv
output_file <- "combined_protein_naguider.csv"
write.csv(combined_protein_filtered_df, file = output_file, row.names = TRUE)

# Save sample names and group as csv
output_file <- "metadata_ft_naguider.csv"
write.csv(metadata_ft[, c("sample_name", "Fragilidad")], file = output_file, row.names = TRUE)

