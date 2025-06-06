# Load libraries
library(FragPipeAnalystR)
library(ggplot2)
library(limma)
library(vsn)
#library(nortest)
library(dplyr)
library(SummarizedExperiment)
#library(reshape2)
library(ComplexHeatmap)
library(UniprotR)
library(impute)

# Working directory
work_path <- getwd()

# Load functions
source(paste0(work_path, "/functions/create_bar_plot.R"))
source(paste0(work_path, "/functions/create_density_plot.R"))
source(paste0(work_path, "/functions/create_pca_plot.R"))
source(paste0(work_path, "/functions/create_lfq_boxplot.R"))

# ==================================
# 1. Making SummarizedExperiment
# ==================================

# Path to combined_protein.tsv and experiment_annotation.tsv
fragpipe_path <- paste0("/mnt/d/proteomica/fragilidad/datos/ProteinIdent/",
                        "MSfraggerSPbacteria/DatosProtIdentCombinados/")

combined_protein <- paste0(fragpipe_path, "combined_protein.tsv")
experiment_ann <- paste0(fragpipe_path, "experiment_annotation.tsv")

# Read quant table
prot_quant <- read.table(combined_protein,
                         header = T,
                         sep = "\t",
                         quote = "",
                         comment.char = "",
                         blank.lines.skip = F,
                         check.names = F
)

# Delete contam rows
prot_quant <- prot_quant[!grepl("contam", prot_quant$Protein),]

# Take first identifier per row and make unique names. If there is no name, then
# ID will be taken
prot_uniq <- prot_quant %>%
  mutate(
    name = get("Gene"),
    ID = get("Protein ID"),
    name = make.unique(ifelse(name == "" | is.na(name), ID, name))
  )

# Set rownames
rownames(prot_uniq) <- prot_uniq$ID

# Select Intensity columns
lfq_col <- grep("Intensity", colnames(prot_uniq))
lfq_col <- lfq_col[!grepl("MaxLFQ", lfq_col)]
prot_lfq <- prot_uniq[, lfq_col]

# Replace 0 by NA
prot_lfq[prot_lfq == 0] <- NA

# Read annotation table
exp_anno <- read.table(experiment_ann,
                       header = T,
                       sep = "\t",
                       stringsAsFactors = F)

exp_anno$label <- paste(exp_anno$sample, "Intensity", sep = " ")

# Set rownames
rownames(exp_anno) <- exp_anno$label

# Match column names quant with label from annotation
matched <- match(
  make.names(exp_anno$label),
  make.names(colnames(prot_lfq))
)

# Check if labels in annotation match with column names in quant
if (any(is.na(matched))) {
  print(make.names(exp_anno$label))
  print(make.names(colnames(prot_lfq)))
}

# Set rownames from annotation to sample name
rownames(exp_anno) <- exp_anno$sample_name
# Set colnames matched from quant to sample name
colnames(prot_lfq)[matched] <- exp_anno$sample_name
# Keep column name not NA and reorder
prot_lfq <- prot_lfq[, !is.na(colnames(prot_lfq))][rownames(exp_anno)]

# Create rowData
row_data <- prot_uniq[, -lfq_col]
rownames(row_data) <- prot_uniq$ID

# Make SummarizedExperiment
se <- SummarizedExperiment(
  assays = as.matrix(prot_lfq),
  colData = exp_anno,
  rowData = row_data,
  metadata = list("log2transform"=F, "lfq_type"="Intensity",
                  "level"="protein")
)

# Check class of se object
class(se)

# Check log2, exp, lfq_type and level
metadata(se)
# $log2transform
# FALSE

#$lfq_type
# "Intensity"

#$level
# "protein"

# Check number of rows and columns
dim(se)
# 6854  274

# Phenotypic variables
colData(se)
names(colData(se))
# "file" "sample"      "sample_name" "condition"   "replicate"   "label"

# Check if more than one sample per case
table(colData(se)$sample_name)

# Intensity matrix
head(assay(se))

# Row names matrix
head(colnames(se))

# ==================================
# 2. Filtering frailty cohort
# ==================================

# Load metadatos
load(paste0(work_path, "/data/metadata.RData"))

# Keep only data from FT cohort
se <- se[, colData(se)$replicate %in% metadata_ft$replicate]
unique(colData(se)$replicate)

# Check number of rows and columns
dim(se)
# 6854  203

# ==================================
# 3. Adding phenotypic data
# ==================================

# Create an aligned dataframe with the order of `replicate` in colData(se)
replicate_se <- colData(se)$replicate
replicate_metadata <- metadata_ft$replicate
metadata_ft <- metadata_ft[match(replicate_se, replicate_metadata), ]

# Add phenotypic variables to colData
colData(se)$frailty <- metadata_ft$frailty
colData(se)$sex <- metadata_ft$sex
colData(se)$education <- metadata_ft$education
colData(se)$alcohol <- metadata_ft$alcohol
colData(se)$tobacco <- metadata_ft$tobacco
colData(se)$diabetes <- metadata_ft$diabetes
colData(se)$chf <- metadata_ft$chf
colData(se)$af <- metadata_ft$af
colData(se)$hipfracture <- metadata_ft$hipfracture
colData(se)$depression <- metadata_ft$depression
colData(se)$osteoarthritis <- metadata_ft$osteoarthritis
colData(se)$sarcopenia <- metadata_ft$sarcopenia
colData(se)$ilef <- metadata_ft$ilef
colData(se)$age <- metadata_ft$age
colData(se)$medas <- metadata_ft$medas
colData(se)$energy <- metadata_ft$energy
colData(se)$bmi <- metadata_ft$bmi

# Check
head(colData(se))

# ==================================
# 4. Filtering
# ==================================
# 4.1. By Missing
# ==================================

# Get intensities matrix from SummarizedExperiment
se_assay <- assay(se)

# Number of proteins
nrow(se_assay)
# 6854

# Number of NA values in assay
table(is.na(se_assay))
# FALSE    TRUE 
# 148253 1243109 
1243109/(1243109+148253)
# 0.8934476

# Number of proteins with all NA
proteins_to_keep <- apply(se_assay, 1, function(x) !all(is.na(x)))
sum(!proteins_to_keep)
# 454 

# Filter proteins with all NA
se_assay <- se_assay[proteins_to_keep, ]
se <- se[proteins_to_keep,]

# Save SummarizedExperiment
save_path <- paste0(work_path,"/data/lfq/se.RData")
save(se, file = save_path)

# Get colData from SummarizedExperiment
col_data <- as.data.frame(colData(se))

# Number of proteins after filtering those with all NA
nrow(se_assay)
# 6400

# Number of NA values in assay after filtering those with all NA
table(is.na(se_assay))
# FALSE    TRUE 
# 148253 1150947 
1150947/(1150947+148253)
# 0.885889
# Number of proteins detected per sample
num_prot_sample <- apply(se_assay, 2, function(x) sum(!is.na(x)))
num_prot_sample <- as.data.frame(num_prot_sample)
colnames(num_prot_sample) <- "Proteins"
num_prot_sample$Sample <- substr(rownames(num_prot_sample), 1, 7)
rownames(num_prot_sample) <- NULL
# Add frailty group
num_prot_sample <- merge(num_prot_sample, col_data[, c("sample", "frailty")], 
                         by.x = "Sample", by.y = "sample")

# Summary total proteins detected per sample
summary(num_prot_sample)
# Sample             Proteins      frailty  
# Length:203         Min.   :  40.0   NFT:138  
# Class :character   1st Qu.: 520.5   FT : 65  
# Mode  :character   Median : 720.0            
# Mean   : 730.3            
# 3rd Qu.: 913.0            
# Max.   :2054.0            

# Check proteins detected in samples with a low detection
rownames(se_assay[!is.na(se_assay[,"YDD_248"]), ])

# Density plot number of proteins per sample in FT vs NFT
p <- ggplot(num_prot_sample, aes(x = Proteins, fill = frailty)) +
  geom_density(alpha = 0.4, linewidth = 0.2) +
  labs(title = "Number of proteins per sample", x = "Number of proteins", 
       y = "Density") +
  scale_fill_manual(values = c("FT" = "blue", "NFT" = "red")) +
  theme_minimal()

save_path <- paste0(work_path, "/plots/preprocessing/lfq/",
                    "density_plot_num_prot_sample.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

# Number of missing per sample
num_na_sample <- apply(se_assay, 2, function(x) sum(is.na(x)))
num_na_sample <- as.data.frame(num_na_sample)
colnames(num_na_sample) <- "NAs"
num_na_sample$Sample <- substr(rownames(num_na_sample), 1, 7)
rownames(num_na_sample) <- NULL
# Add frailty group
num_na_sample <- merge(num_na_sample, col_data[, c("sample", "frailty")], 
                       by.x = "Sample", by.y = "sample")

# Density plot number of missing per sample in FT vs NFT
p <- ggplot(num_na_sample, aes(x = NAs, fill = frailty)) +
  geom_density(alpha = 0.4, linewidth = 0.2) +
  labs(title = "Number of NAs per sample", x = "Number of NAs", y = "Density") +
  scale_fill_manual(values = c("FT" = "blue", "NFT" = "red")) +
  theme_minimal()

save_path <- paste0(work_path, "/plots/preprocessing/lfq/", 
                    "density_plot_num_na_sample.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

# Make binary assay (1 = valid value, 0 = missing value)
missval <- ifelse(is.na(se_assay), 0, 1)

# Convert matrix to long format dataframe
df <- melt(missval)

# Rename columns
colnames(df) <- c("protein", "sample", "intensity")

# Add frailty group
df <- merge(df, col_data[, c("sample", "frailty")], by = "sample")

# Plot heatmap
p <- ggplot(df, aes(sample, protein, fill = intensity)) + 
  geom_tile() +
  facet_grid(. ~ frailty, scales = "free_x", space = "free_x") +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 2, angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_gradient(low = "white", high = "black")

save_path <- paste0(work_path,"/plots/preprocessing/lfq/heatmap_na.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

# Complex Heatmap
# Clustering rows and columns by binary distance
# Column split by frailty groups
p <- Heatmap(missval, name = "NA", col = c("white", "black"),
             show_row_dend = FALSE, 
             column_names_gp = gpar(fontsize = 4),
             show_row_names = FALSE, 
             # cluster_columns = FALSE,
             # cluster_rows = FALSE,
             column_split = colData(se)$frailty, 
             clustering_distance_rows = "binary", 
             clustering_distance_columns = "binary", 
             use_raster = FALSE)

save_path <- paste0(work_path,"/plots/preprocessing/lfq/heatmap_clust_na.png")
png(save_path, width = 8, height = 6, units = "in", res = 300)
draw(p)
dev.off()

# Annotation of columns by frailty group
column_ann <- HeatmapAnnotation(frailty = colData(se)$frailty,
                                col = list(frailty = c("FT" = "red", 
                                                       "NFT" = "blue")))
p <- Heatmap(missval, name = "NA", col = c("white", "black"), 
             show_row_dend = FALSE, column_names_gp = gpar(fontsize = 4),
             show_row_names = FALSE, 
             #cluster_columns = FALSE, 
             #cluster_rows = FALSE,
             clustering_distance_rows = "binary", 
             clustering_distance_columns = "binary",
             top_annotation = column_ann,
             use_raster = FALSE)

save_path <- paste0(work_path,"/plots/preprocessing/lfq/heatmap_clust_na_ann.png")
png(save_path, width = 8, height = 6, units = "in", res = 300)
draw(p)
dev.off()

# Sample proportion for each protein
presence_proportion <- rowSums(!is.na(se_assay)) / ncol(se_assay)
sample_percent <- presence_proportion * 100
sample_percent <- data.frame(SamplePercentage = sample_percent)
p <- create_density_plot(sample_percent, "SamplePercentage")

save_path <- paste0(work_path, "/plots/preprocessing/lfq/",
                    "density_plot_SamplePercentage.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

# Number of proteins present in more than 25%
sum(presence_proportion > 0.25)
# 897

# Number of proteins present in more than 50%
sum(presence_proportion > 0.5)
# 403

# Number of proteins present in more than 75%
sum(presence_proportion > 0.75)
# 171

# Number of proteins present in 100%
sum(presence_proportion == 1)
# 9

# Proteins present in 100% of samples
rownames(se_assay[presence_proportion == 1,])
# P00761 --> trypsin
# P04264 --> Keratin, type II cytoskeletal 1
# P06702 --> Protein S100-A9: regulation of inflammatory processes and
# immune response
# A6KYK2
# P0DTE7
# P35527
# P94316
# Q8RQP4
# Q9UGM3

# Proteins in more than 75% of samples
rownames(se_assay[presence_proportion >= 0.75,])

# Proteins in less than 50% of samples
rownames(se_assay[presence_proportion < 0.5,])

# Keep proteins with minimum fraction of valid values in at least one condition
min_fraction <- 0.5
fraction_ft <- apply(se_assay[,colData(se)$frailty == "FT"], 1,
                     function(x) sum(!is.na(x)) / length(x))
fraction_nft <- apply(se_assay[, colData(se)$frailty == "NFT"], 1,
                      function(x) sum(!is.na(x)) / length(x))

proteins_to_keep <- (fraction_ft >= min_fraction) | (fraction_nft >= 
                                                       min_fraction)

# Number of proteins present in >= 50% of samples in at least one condition
sum(proteins_to_keep) 
# 450

# Filter SummarizedExperiment
se_filt_miss <- se[proteins_to_keep, ]

# Get protein intensities matrix from filtered SummarizedExperiment
se_filt_miss_assay <- assay(se_filt_miss)

# Number of proteins after filtering those present in <50% of samples in both
# conditions
nrow(se_filt_miss_assay)
# 450

# Number of proteins detected per sample
num_prot_sample <- apply(se_filt_miss_assay, 2, function(x) sum(!is.na(x)))
num_prot_sample <- as.data.frame(num_prot_sample)
colnames(num_prot_sample) <- "Proteins"
num_prot_sample$Sample <- substr(rownames(num_prot_sample), 1, 7)
rownames(num_prot_sample) <- NULL
# Add frailty group
num_prot_sample <- merge(num_prot_sample, col_data[, c("sample", "frailty")], 
                         by.x = "Sample", by.y = "sample")

# Summary total proteins detected per sample
summary(num_prot_sample)
# Sample             Proteins     frailty  
# Length:203         Min.   : 35.0   NFT:138  
# Class :character   1st Qu.:263.0   FT : 65  
# Mode  :character   Median :338.0            
# Mean   :315.2            
# 3rd Qu.:385.0            
# Max.   :446.0            

# Density plot number of proteins per sample in FT vs NFT
p <- ggplot(num_prot_sample, aes(x = Proteins, fill = frailty)) +
  geom_density(alpha = 0.4, linewidth = 0.2) +
  labs(title = "Number of proteins per sample", x = "Number of proteins", 
       y = "Density") +
  scale_fill_manual(values = c("FT" = "blue", "NFT" = "red")) +
  theme_minimal()

save_path <- paste0(work_path, "/plots/preprocessing/lfq/density_plot_",
                    "num_prot_sample_filt_min_fract.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

# Number of NA values in assay
table(is.na(se_filt_miss_assay))
# FALSE  TRUE 
# 63977 27373 

# Number of missing per sample
num_na_sample <- apply(se_filt_miss_assay, 2, function(x) sum(is.na(x)))
num_na_sample <- as.data.frame(num_na_sample)
colnames(num_na_sample) <- "NAs"
num_na_sample$Sample <- substr(rownames(num_na_sample), 1, 7)
rownames(num_na_sample) <- NULL
# Add frailty group
num_na_sample <- merge(num_na_sample, col_data[, c("sample", "frailty")], 
                       by.x = "Sample", by.y = "sample")

# Density plot number of missing per sample in FT vs NFT
p <- ggplot(num_na_sample, aes(x = NAs, fill = frailty)) +
  geom_density(alpha = 0.4, linewidth = 0.2) +
  labs(title = "Number of NAs per sample", x = "Number of NAs", y = "Density") +
  scale_fill_manual(values = c("FT" = "blue", "NFT" = "red")) +
  theme_minimal()

save_path <- paste0(work_path, "/plots/preprocessing/lfq/density_plot_",
                    "num_na_sample_filt_min_fract.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

# Make binary assay (1 = valid value, 0 = missing value)
missval <- ifelse(is.na(se_filt_miss_assay), 0, 1)

# Convert matrix to long format dataframe
df <- melt(missval)

# Rename columns
colnames(df) <- c("protein", "sample", "intensity")

# Add frailty group
df <- merge(df, col_data[, c("sample", "frailty")], by = "sample")

# Plot heatmap
p <- ggplot(df, aes(sample, protein, fill = intensity)) + 
  geom_tile() +
  facet_grid(. ~ frailty, scales = "free_x", space = "free_x") +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 2, angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_gradient(low = "white", high = "black")

save_path <- paste0(work_path, "/plots/preprocessing/lfq/",
                    "heatmap_na_filt_min_fract.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

# Complex Heatmap
# Clustering rows and columns by binary distance
# Column split by frailty groups
p <- Heatmap(missval, name = "NA", col = c("white", "black"),
             show_row_dend = FALSE, 
             column_names_gp = gpar(fontsize = 4), show_row_names = FALSE,
             # cluster_columns = FALSE,
             # cluster_rows = FALSE,
             column_split = col_data$frailty,
             clustering_distance_rows = "binary",
             clustering_distance_columns = "binary")

save_path <- paste0(work_path, "/plots/preprocessing/lfq/heatmap_clust_",
                    "na_filt_min_fract.png")
png(save_path, width = 8, height = 6, units = "in", res = 300)
draw(p)
dev.off()

# Annotation of columns by frailty group
column_ann <- HeatmapAnnotation(frailty = col_data$frailty,
                                col = list(frailty = c("FT" = "red", 
                                                       "NFT" = "blue")))
p <- Heatmap(missval, name = "NA", col = c("white", "black"),
             show_row_dend = FALSE,
             column_names_gp = gpar(fontsize = 4), show_row_names = FALSE,
             # cluster_columns = FALSE,
             # cluster_rows = FALSE,
             clustering_distance_rows = "binary",
             clustering_distance_columns = "binary",
             top_annotation = column_ann)

save_path <- paste0(work_path,"/plots/preprocessing/lfq/heatmap_clust_",
                    "na_ann_filt_min_fract.png")
png(save_path, width = 8, height = 6, units = "in", res = 300)
draw(p)
dev.off()

# ==================================
# 4.2. By organism
# ==================================
# Get lineage from UniProt
lineage <- GetProteinAnnontate(rownames(se_filt_miss_assay),
                               columns = c("lineage"))
lineage_df <- as.data.frame(lineage)
lineage_df$protein_id <- rownames(se_filt_miss_assay)
colnames(lineage_df)[colnames(lineage_df) == "lineage_df"] <- "lineage"
# Get superkingdom from lineage
lineage_df$superkingdom <- sapply(strsplit(lineage_df$lineage, ","),
                                  function(x) x[2])
lineage_df$genus <- ifelse(
  grepl("Eukaryota", lineage_df$superkingdom),
  sapply(strsplit(lineage_df$lineage, ","), function(x) tail(x, 1)),
  NA
)


# Get protein_id from bacteria proteins
proteins_to_keep <- lineage_df[grepl("Bacteria", lineage_df$superkingdom) |
                                 grepl("Homo", lineage_df$genus), 
                               "protein_id"]
# Filter SummarizedExperiment
se_filt_org <- se_filt_miss[proteins_to_keep, ]

# Number of proteins after filtering no bacteria proteins
nrow(se_filt_org)
# 442

se_filt_org_assay <- assay(se_filt_org)

# ==================================
# 4.3. (NO) By abundance
# ==================================
# - Ribosome subunits
# - Elongation factor Tu
# - Keratin
# ==================================
ribosub <- rowData(se_filt_org)[grep("ribosomal",
                                      rowData(se_filt_org)$Description),
                                 c("Protein ID", "Organism", "Description")]
ribosub <- as.data.frame(ribosub)
eftu <- rowData(se_filt_org)[grep("Elongation factor Tu",
                                   rowData(se_filt_org)$Description),
                              c("Protein ID", "Organism", "Description")]
eftu <- as.data.frame(eftu)
keratin <- rowData(se_filt_org)[grep("Keratin",
                                      rowData(se_filt_org)$Description),
                                 c("Protein ID", "Organism", "Description")]
keratin <- as.data.frame(keratin)

#se_filtered <- se_filt_bact[-grep("ribosomal", 
#                                  rowData(se_filt_bact)$Description), ]
#se_filtered <- se_filtered[-grep("Elongation factor Tu", 
#                                 rowData(se_filtered)$Description), ]
#se_filtered <- se_filtered[-grep("Keratin", 
#                                 rowData(se_filtered)$Description), ]

se_filtered <- se_filt_org
se_filtered_assay <- se_filt_org_assay

nrow(se_filtered_assay)
# 442

# Save SummarizedExperiment filtered
save_path <- paste0(work_path,"/data/lfq/se_filtered.RData")
save(se_filtered, file = save_path)

# Boxplot lfq after filtering
#res <- create_lfq_boxplot(se_filtered, col_data, lfq_type = "Intensity")
#print(res$plot_samples)
#print(res$plot_groups)
#print(res$lfq_summary)

# Transform se_filtered assay to long
se_filtered_df <- as.data.frame(se_filtered_assay)
se_filtered_df$protein_id <- rownames(se_filtered_df)
rownames(se_filtered_df) <- NULL
se_filtered_long <- melt(se_filtered_df)
colnames(se_filtered_long) <- c("protein_id", "Samples", "Intensity")

# Add metadata
se_filtered_long <- as.data.frame(merge(se_filtered_long, 
                                        col_data[,c("sample", "frailty")], 
                                        by.x = "Samples", by.y = "sample"))

summary(se_filtered_long$Intensity)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
# 192     42321     90238    284668    211653 643911360     26842 

# Plot intensity vs samples
p <- ggplot(se_filtered_long, aes(x = Samples, y = Intensity, fill = frailty)) + 
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    legend.position = c(0.95, 0.95),
    legend.background = element_rect(fill = alpha("white", 0.5))) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple")) +
  ylim(0, 250000) +
  facet_wrap(~frailty, scales = "free_x", ncol = 1)

save_path <- paste0(work_path,"/plots/preprocessing/lfq/boxplot_lfq_samples.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

# Plot intensity vs frailty group
p <- ggplot(se_filtered_long, aes(x = frailty, y = Intensity, fill = frailty)) + 
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    legend.position = c(0.95, 0.95),
    legend.background = element_rect(fill = alpha("white", 0.5))) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple")) +
  ylim(0, 250000)

save_path <- paste0(work_path,"/plots/preprocessing/lfq/boxplot_lfq_groups.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

# Intensities density plot
p <- ggplot(se_filtered_long, aes(x = Intensity)) +
  geom_density(fill = "blue", alpha = 0.4) +
  theme_minimal() + 
  xlim(0, 700000)

save_path <- paste0(work_path,"/plots/preprocessing/lfq/density_plot_lfq.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

# Get proteins with intensities higher than percentile 95%
#Q95 <- quantile(se_filtered_long$LFQ, 0.95, na.rm = TRUE)
#se_filtered_long_q95 <- se_filtered_long[se_filtered_long$LFQ > Q95, ]
#se_filtered_long_q95 <- se_filtered_long_q95[!is.na(se_filtered_long_q95), ]
#prot_q95 <-  se_filtered_long_q95 %>% arrange(desc(LFQ))
#prot_q95 <- unique(prot_q95$protein_id)
#prot_q95_ann <- GetProteinAnnontate(prot_q95,columns = c("gene_primary",
#                                                         "organism_name",
#                                                          "protein_name",
#                                                          "cc_function"))
#write.csv2(prot_q95_ann, "prot_q95_ann.csv")


# Intensities density plot per sample
save_path <- paste0(work_path, "/plots/preprocessing/lfq/densities_plot_lfq.png")
png(save_path, width = 8, height = 6, units = "in", res = 300)
plotDensities(se_filtered_assay, legend = FALSE)
dev.off()

#for (i in 1:nrow(se_filtered_assay)){
#  plotDensities(se_filtered_assay[, i], legend = TRUE)
#}

# Standard deviation (sd) and mean are calculated row-wise from the expression
# matrix. The scatterplot of these versus each other to verify whether there is
# a dependence of the sd on the mean. The red line running median estimator
# (window-width 10%). If there is no variance-mean dependence, the line should
# be aprox. horizontal.
save_path <- paste0(work_path, "/plots/preprocessing/lfq/scatterplot.png")
png(save_path, width = 8, height = 6, units = "in", res = 300)
meanSdPlot(se_filtered_assay)
dev.off()

#plotMA(se_filtered_assay)

# Q-Q sample quantiles vs theoretical quantiles for a normal distribution
for (i in 1:ncol(se_filtered_assay)) {
  qqnorm(se_filtered_assay[, i], main = paste("QQ Plot -", 
                                              colnames(se_filtered_assay)[i]))
  qqline(se_filtered_assay[, i], col = "red")
}

# Normal contrast
p_values <- apply(se_filtered_assay, 2,
                  function(column) shapiro.test(column)$p.value)
norm_results <- p_values >= 0.05 # TRUE if normal, FALSE if not
summary(norm_results)
#    Mode   FALSE
#logical     203

# ==================================
# 5. Normalization
# ==================================
# 5.1. vsn
# ==================================

fit <- vsnMatrix(se_filtered_assay)
se_norm <- se_filtered
assay(se_norm) <- predict(fit, se_filtered_assay)
se_norm_assay <- assay(se_norm)

# Save SummarizedExperiment normalized by vsn
save_path <- paste0(work_path,"/data/lfq/se_norm.RData")
save(se_norm, file = save_path)

# Boxplot lfq after vsn
res <- create_lfq_boxplot(se_norm, col_data, lfq_type = "Intensity")
p <- res$plot_samples
save_path <- paste0(work_path,"/plots/preprocessing/lfq/boxplot_lfq_",
                    "samples_vsn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

p <- res$plot_groups
save_path <- paste0(work_path,"/plots/preprocessing/lfq/boxplot_lfq_",
                    "groups_vsn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

print(res$lfq_summary)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 7.579  15.350  16.430  16.580  17.626  27.863   26842  

# Intensities density plot
se_norm_long <- res$se_long
p <- ggplot(se_norm_long, aes(x = Intensity)) +
  geom_density(fill = "blue", alpha = 0.4) +
  theme_minimal()

save_path <- paste0(work_path,"/plots/preprocessing/lfq/density_plot_",
                    "lfq_vsn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

# Intensities density plot per sample
save_path <- paste0(work_path, "/plots/preprocessing/lfq/densities_plot_",
                    "lfq_vsn.png")
png(save_path, width = 8, height = 6, units = "in", res = 300)
plotDensities(se_norm_assay, legend = FALSE)
dev.off()

#for (i in 1:nrow(se_norm_assay)){
#  plotDensities(se_norm_assay[, i], legend = TRUE)
#}

# Scatterplot
save_path <- paste0(work_path, "/plots/preprocessing/lfq/scatterplot_vsn.png")
png(save_path, width = 8, height = 6, units = "in", res = 300)
meanSdPlot(se_norm_assay)
dev.off()

# Normal contrast
p_values <- apply(se_norm_assay, 2,
                  function(column) shapiro.test(column)$p.value)
norm_results <- p_values >= 0.05 # TRUE if normal, FALSE if not
summary(norm_results)
#    Mode   FALSE    TRUE 
#logical     189      14 

# ==================================
# 6. Imputation
# ==================================
# 6.1. No imputation
# ==================================
se_no_imp <- se_norm

# Save SummarizedExperiment no imputated
save_path <- paste0(work_path,"/data/lfq/se_no_imp.RData")
save(se_no_imp, file = save_path)

# ==================================
# 6.2. Impute missing values
# ==================================
# 6.2.1. Perseus
# ==================================
# Missing values are replaced with random values generated from a shifted and
# scaled normal distribution based on the existing data
se_perseus <- manual_impute(se_no_imp)
se_perseus_assay <- assay(se_perseus)

# Plot PCA
p <- create_pca_plot(se_perseus, n_top_loadings = 5)
save_path <- paste0(work_path,"/plots/preprocessing/lfq/pca_plot_perseus.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

# Boxplot lfq after perseus
res <- create_lfq_boxplot(se_perseus, col_data, lfq_type = "Intensity")
p <- res$plot_samples
save_path <- paste0(work_path,"/plots/preprocessing/lfq/boxplot_lfq_",
                    "samples_perseus.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

p <- res$plot_groups
save_path <- paste0(work_path,"/plots/preprocessing/lfq/boxplot_lfq_",
                    "groups_perseus.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

print(res$lfq_summary)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7.579  13.829  15.516  15.598  17.075  27.863 

# Intensities density plot
se_perseus_long <- res$se_long
p <- ggplot(se_perseus_long, aes(x = Intensity)) +
  geom_density(fill = "blue", alpha = 0.4) +
  theme_minimal()

save_path <- paste0(work_path,"/plots/preprocessing/lfq/density_plot_lfq_",
                    "perseus.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

# Intensities density plot per sample
save_path <- paste0(work_path, "/plots/preprocessing/lfq/densities_plot_lfq_",
                    "perseus.png")
png(save_path, width = 8, height = 6, units = "in", res = 300)
plotDensities(se_perseus_assay, legend = FALSE)
dev.off()

# Save SummarizedExperiment imputed by perseus
save_path <- paste0(work_path,"/data/lfq/se_perseus.RData")
save(se_perseus, file = save_path)


# ==================================
# 6.2.2. KNN
# ==================================
# Delete 10 samples with more than 80% missing values
se_knn <- se_no_imp
se_knn <- se_knn[, colMeans(is.na(assay(se_knn))) <= 0.8]

# Get colData from SummarizedExperiment
col_data <- as.data.frame(colData(se_knn))

# Impute using k-nearest neighbors
assay(se_knn) <- impute.knn(assay(se_knn), rowmax = 0.6)$data
se_knn_assay <- assay(se_knn)

# Plot PCA
p <- create_pca_plot(se_knn, n_top_loadings = 5)
save_path <- paste0(work_path,"/plots/preprocessing/lfq/pca_plot_knn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

# Boxplot lfq after knn
res <- create_lfq_boxplot(se_knn, col_data, lfq_type = "Intensity")
p <- res$plot_samples
save_path <- paste0(work_path,"/plots/preprocessing/lfq/boxplot_lfq_",
                    "samples_knn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

p <- res$plot_groups
save_path <- paste0(work_path,"/plots/preprocessing/lfq/boxplot_lfq_",
                    "groups_knn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

print(res$lfq_summary)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7.579  15.253  16.059  16.324  17.179  27.863

# Intensities density plot
se_knn_long <- res$se_long
p <- ggplot(se_knn_long, aes(x = Intensity)) +
  geom_density(fill = "blue", alpha = 0.4) +
  theme_minimal()

save_path <- paste0(work_path,"/plots/preprocessing/lfq/density_plot_",
                    "lfq_knn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)

# Intensities density plot per sample
save_path <- paste0(work_path, "/plots/preprocessing/lfq/densities_plot_",
                    "lfq_knn.png")
png(save_path, width = 8, height = 6, units = "in", res = 300)
plotDensities(se_knn_assay, legend = FALSE)
dev.off()

# Save SummarizedExperiment imputed by knn
save_path <- paste0(work_path,"/data/lfq/se_knn.RData")
save(se_knn, file = save_path)
