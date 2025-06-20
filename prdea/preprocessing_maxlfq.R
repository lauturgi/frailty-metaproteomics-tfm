# Load libraries
library(FragPipeAnalystR)
library(ggplot2)
library(limma)
library(vsn)
library(dplyr)
library(SummarizedExperiment)
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
combined_protein <- paste0(work_path, "/prdea/data/combined_protein.tsv")
experiment_ann <- paste0(work_path, "/prdea/data/experiment_annotation.tsv")

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

# Select MaxLFQ columns
lfq_col <- grep("MaxLFQ", colnames(prot_uniq))
prot_maxlfq <- prot_uniq[, lfq_col]

# Replace 0 by NA
prot_maxlfq[prot_maxlfq == 0] <- NA

# Read annotation table
exp_anno <- read.table(experiment_ann,
                       header = T,
                       sep = "\t",
                       stringsAsFactors = F)

exp_anno$label <- paste(exp_anno$sample, "MaxLFQ.Intensity", sep = " ")

# Set rownames
rownames(exp_anno) <- exp_anno$label

# Match column names quant with label from annotation
matched <- match(
  make.names(exp_anno$label),
  make.names(colnames(prot_maxlfq))
  )

# Check if labels in annotation match with column names in quant
if (any(is.na(matched))) {
  print(make.names(exp_anno$label))
  print(make.names(colnames(prot_maxlfq)))
}

# Set rownames from annotation to sample name
rownames(exp_anno) <- exp_anno$sample_name
# Set colnames matched from quant to sample name
colnames(prot_maxlfq)[matched] <- exp_anno$sample_name
# Keep column name not NA and reorder
prot_maxlfq <- prot_maxlfq[, !is.na(colnames(prot_maxlfq))][rownames(exp_anno)]

# Create rowData
row_data <- prot_uniq[, -lfq_col]
rownames(row_data) <- prot_uniq$ID

# Make SummarizedExperiment (SE)
se <- SummarizedExperiment(
  assays = as.matrix(prot_maxlfq),
  colData = exp_anno,
  rowData = row_data,
  metadata = list("log2transform"=F, "lfq_type"="MaxLFQ",
                  "level"="protein")
)

# Check class of se object
class(se)

# Check log2, exp, lfq_type and level
metadata(se)
# $log2transform
# FALSE

#$lfq_type
# "MaxLFQ"

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

# Get intensities matrix from SE
se_assay <- assay(se)

# Number of proteins
nrow(se_assay)
# 6854

# Number of NA values in assay
na_table <- table(is.na(se_assay))
na_table
# FALSE    TRUE 
# 61378 1329984

prop_na <- na_table["TRUE"] / sum(na_table)
prop_na
# 0.9558864

# Proteins with at least one intensity value
proteins_to_keep <- apply(se_assay, 1, function(x) !all(is.na(x)))

# Number of proteins to keep
sum(proteins_to_keep)
# 2902

# Number of proteins with all NA
sum(!proteins_to_keep)
# 3952 

# Filter proteins with all NA
se <- se[proteins_to_keep,]
se_assay <- se_assay[proteins_to_keep, ]

# Save SE
save_path <- paste0(work_path,"/prdea/data/maxlfq/se.RData")
save(se, file = save_path)

# Get colData from SE
col_data <- as.data.frame(colData(se))

# Check number of proteins after filtering those with all NA
nrow(se_assay)
# 2902

# Number of NA values in assay after filtering those with all NA
na_table <- table(is.na(se_assay))
na_table
# FALSE   TRUE 
# 61378 527728

# NA proportion
prop_na <- na_table["TRUE"] / sum(na_table)
prop_na
# 0.8958116

# Number of proteins detected per sample
num_prot_sample <- apply(se_assay, 2, function(x) sum(!is.na(x)))
num_prot_sample <- as.data.frame(num_prot_sample)
colnames(num_prot_sample) <- "proteins"
num_prot_sample$sample <- rownames(num_prot_sample)
rownames(num_prot_sample) <- NULL
# Add frailty group
num_prot_sample <- merge(num_prot_sample, col_data[, c("sample", "frailty")], 
                         by = "sample")

# Summary total proteins detected per sample
summary(num_prot_sample)
#Sample             Proteins      frailty  
#Length:203         Min.   :  14.0   NFT:138  
#Class :character   1st Qu.: 183.5   FT : 65  
#Mode  :character   Median : 288.0            
#                   Mean   : 302.4            
#                   3rd Qu.: 378.5            
#                   Max.   :1101.0  

# Check proteins detected in samples with a low detection
rownames(se_assay[!is.na(se_assay[,"YDD_248"]), ])

# Density plot number of proteins per sample in FT vs NFT
p <- ggplot(num_prot_sample, aes(x = proteins, fill = frailty)) +
  geom_density(linewidth = 0.4, alpha = 0.3) +
  labs(title = "Distribution of proteins per sample",
       x = "Number of proteins", 
       y = "Density") +
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail")) +
  theme_minimal(base_size = 12) +
  coord_cartesian(ylim = c(0, 0.003)) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/",
                    "density_plot_num_prot_sample.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Mean number of proteins in frail
mean(num_prot_sample[num_prot_sample$frailty == "FT", "proteins"])
# 275.9692
# Standard deviation number of proteins in frail
sd(num_prot_sample[num_prot_sample$frailty == "FT", "proteins"])
# 176.3123

# Mean number of proteins in non-frail
mean(num_prot_sample[num_prot_sample$frailty == "NFT", "proteins"])
# 314.7826
# Standard deviation number of proteins in non-frail
sd(num_prot_sample[num_prot_sample$frailty == "NFT", "proteins"])
# 188.8102

# Wilcox test number of proteins FT vs NFT
wilcox.test(proteins ~ frailty, data = num_prot_sample)
# Wilcoxon rank sum test with continuity correction
# data:  Proteins by frailty
# W = 5127, p-value = 0.1004
# alternative hypothesis: true location shift is not equal to 0

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

save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/",
                    "heatmap_na.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Complex Heatmap
# Clustering rows and columns by binary distance
# Column split by frailty groups
p <- Heatmap(missval, name = "NA", col = c("white", "black"),
             show_row_dend = FALSE, 
             column_names_gp = gpar(fontsize = 4),
             show_row_names = FALSE, 
             # cluster_columns = FALSE,
             # cluster_rows = FALSE,
             column_split = col_data$frailty, 
             clustering_distance_rows = "binary", 
             clustering_distance_columns = "binary", 
             use_raster = FALSE)

save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/",
                    "heatmap_clust_na.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
draw(p)
dev.off()

# Annotation of columns by frailty group
column_ann <- HeatmapAnnotation(
  frailty = col_data$frailty,
  col = list(frailty = c("FT" = "red", "NFT" = "blue")),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

p <- Heatmap(missval, col = c("white", "black"),
             show_row_dend = FALSE, column_names_gp = gpar(fontsize = 4),
             show_row_names = FALSE, 
             show_column_names = FALSE,
             #cluster_columns = FALSE, 
             #cluster_rows = FALSE,
             clustering_distance_rows = "binary", 
             clustering_distance_columns = "binary",
             top_annotation = column_ann,
             use_raster = FALSE,
             show_heatmap_legend = FALSE
)
# Heatmap legend
legend_1 <- Legend(at = c(0, 1),
                   labels = c("Missing value", "Valid value"),
                   legend_gp = gpar(fill = c("white", "black")),
                   title = NULL)
# Annotation legend
legend_2 <- Legend(at = c("FT", "NFT"),
                   labels = c("Frail", "Non-Frail"),
                   legend_gp = gpar(fill = c("red", "blue")),
                   title = NULL)

# Combine legends
legends <- packLegend(legend_2, legend_1, direction = "horizontal")

# Plot heatmap
save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/",
                    "heatmap_clust_na_ann.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
draw(p, heatmap_legend_list = legends, heatmap_legend_side = "top")
dev.off()

# Sample proportion for each protein
sample_proportion <- rowSums(!is.na(se_assay)) / ncol(se_assay)
sample_perce <- sample_proportion * 100
sample_perce <- data.frame(percentage = sample_perce)
p <- create_density_plot(sample_perce, "percentage",
                         "Sample proportion (%) per protein", 
                         "Distribution of protein detection across samples")

save_path <- paste0(work_path, "/prdea/plots/preprocessing/maxlfq/",
                    "density_plot_sample_proportion.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Number of proteins present in more than 25%
sum(sample_proportion > 0.25)
# 320

# Number of proteins present in more than 50%
sum(sample_proportion > 0.5)
# 159

# Number of proteins present in more than 75%
sum(sample_proportion > 0.75)
# 72

# Number of proteins present in 100%
sum(sample_proportion == 1)
# 3

# Proteins present in 100% of samples
rownames(se_assay[sample_proportion == 1,])
# P00761; trypsin
# P04264; Keratin, type II cytoskeletal 1
# P06702; Protein S100-A9: regulation of inflammatory processes and immune response

# Proteins in more than 75% of samples
rownames(se_assay[sample_proportion >= 0.75,])

# Proteins in less than 50% of samples
rownames(se_assay[sample_proportion < 0.5,])


# Keep proteins with minimum fraction of valid values in at least one condition
min_fraction <- 0.5
fraction_ft <- apply(se_assay[, col_data$frailty == "FT"], 1,
                     function(x) sum(!is.na(x)) / length(x))
fraction_nft <- apply(se_assay[, col_data$frailty == "NFT"], 1,
                      function(x) sum(!is.na(x)) / length(x))

proteins_to_keep <- (fraction_ft >= min_fraction) | (fraction_nft >= 
                                                       min_fraction)

# Number of proteins present in >= 50% of samples in at least one condition
sum(proteins_to_keep) 
# 185

# Filter SE
se_filt_miss <- se[proteins_to_keep, ]

# Get protein intensity matrix
se_filt_miss_assay <- assay(se_filt_miss)

# Check number of proteins after filtering by sample proportion
nrow(se_filt_miss_assay)
# 185

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

# Get protein_id from bacteria and homo
proteins_to_keep <- lineage_df[grepl("Bacteria", lineage_df$superkingdom) |
                                 grepl("Homo", lineage_df$genus), "protein_id"]

# Filter SE
se_filtered <- se_filt_miss[proteins_to_keep, ]

# Save SE filtered
save_path <- paste0(work_path,"/prdea/data/maxlfq/se_filtered.RData")
save(se_filtered, file = save_path)

# Get protein intensity matrix
se_filtered_assay <- assay(se_filtered)

# Check number of proteins after filtering no bacteria or homo
nrow(se_filtered_assay)
# 182

# Number of NA values in assay after filtering by sample proportion and organism
na_table <- table(is.na(se_filtered_assay))
na_table
# FALSE  TRUE 
# 25454 11492

# NA proportion
prop_na <- na_table["TRUE"] / sum(na_table)
prop_na
# 0.3110486

# Number of proteins detected per sample
num_prot_sample <- apply(se_filtered_assay, 2, function(x) sum(!is.na(x)))
num_prot_sample <- as.data.frame(num_prot_sample)
colnames(num_prot_sample) <- "proteins"
num_prot_sample$sample <- rownames(num_prot_sample)
rownames(num_prot_sample) <- NULL
# Add frailty group
num_prot_sample <- merge(num_prot_sample, col_data[, c("sample", "frailty")], 
                         by = "sample")

# Summary total proteins detected per sample
summary(num_prot_sample)
# sample             proteins     frailty  
# Length:203         Min.   : 13.0   NFT:138  
# Class :character   1st Qu.:101.5   FT : 65  
# Mode  :character   Median :137.0            
#                    Mean   :125.4            
#                    3rd Qu.:159.0            
#                    Max.   :180.0            

# Density plot number of proteins per sample in FT vs NFT
p <- ggplot(num_prot_sample, aes(x = proteins, fill = frailty)) +
  geom_density(linewidth = 0.4, alpha = 0.3) +
  labs(title = "Distribution of proteins per sample",
       x = "Number of proteins", 
       y = "Density") +
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail")) +
  theme_minimal(base_size = 12) +
  coord_cartesian(ylim = c(0, 0.015)) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

save_path <- paste0(work_path, "/prdea/plots/preprocessing/maxlfq/",
                    "density_plot_num_prot_sample_filt.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Mean number of proteins in frail
mean(num_prot_sample[num_prot_sample$frailty == "FT", "proteins"])
# 118.3385
# Standard deviation number of proteins in frail
sd(num_prot_sample[num_prot_sample$frailty == "FT", "proteins"])
# 38.23826

# Mean number of proteins in non-frail
mean(num_prot_sample[num_prot_sample$frailty == "NFT", "proteins"])
# 128.7101
# Standard deviation number of proteins in non-frail
sd(num_prot_sample[num_prot_sample$frailty == "NFT", "proteins"])
# 43.28402

# Wilcox test number of proteins FT vs NFT
wilcox.test(proteins ~ frailty, data = num_prot_sample)
# Wilcoxon rank sum test with continuity correction
# data:  Proteins by frailty
# W = 5435, p-value = 0.01503
# alternative hypothesis: true location shift is not equal to 0

# Make binary assay
missval <- ifelse(is.na(se_filtered_assay), 0, 1)

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

save_path <- paste0(work_path, "/prdea/plots/preprocessing/maxlfq/",
                    "heatmap_na_filt.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Complex Heatmap
# Clustering rows and columns by binary distance
# Column split by frailty groups
p <- Heatmap(missval, name = "NA", col = c("white", "black"),
             show_row_dend = FALSE, 
             column_names_gp = gpar(fontsize = 4),
             show_row_names = FALSE,
             # cluster_columns = FALSE,
             # cluster_rows = FALSE,
             column_split = col_data$frailty,
             clustering_distance_rows = "binary",
             clustering_distance_columns = "binary",
             use_raster = FALSE)

save_path <- paste0(work_path, "/prdea/plots/preprocessing/maxlfq/",
                    "heatmap_clust_na_filt_min.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
draw(p)
dev.off()

# Annotation of columns by frailty group
column_ann <- HeatmapAnnotation(
  frailty = col_data$frailty,
  col = list(frailty = c("FT" = "red", "NFT" = "blue")),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

p <- Heatmap(missval, col = c("white", "black"),
             show_row_dend = FALSE,
             column_names_gp = gpar(fontsize = 4),
             show_row_names = FALSE,
             show_column_names = FALSE,
             # cluster_columns = FALSE,
             # cluster_rows = FALSE,
             clustering_distance_rows = "binary",
             clustering_distance_columns = "binary",
             top_annotation = column_ann,
             use_raster = FALSE,
             show_heatmap_legend = FALSE
             )

# Heatmap legend
legend_1 <- Legend(at = c(0, 1),
                   labels = c("Missing value", "Valid value"),
                   legend_gp = gpar(fill = c("white", "black")),
                   title = NULL)
# Annotation legend
legend_2 <- Legend(at = c("FT", "NFT"),
                   labels = c("Frail", "Non-Frail"),
                   legend_gp = gpar(fill = c("red", "blue")),
                   title = NULL)

# Combine legends
legends <- packLegend(legend_2, legend_1, direction = "horizontal")

# Plot heatmap
save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/heatmap_clust_",
                    "na_ann_filt.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
draw(p, heatmap_legend_list = legends, heatmap_legend_side = "top")
dev.off()

# Transform se_filtered assay to long
se_filtered_df <- as.data.frame(se_filtered_assay)
se_filtered_df$protein_id <- rownames(se_filtered_df)
rownames(se_filtered_df) <- NULL
se_filtered_long <- melt(se_filtered_df)
colnames(se_filtered_long) <- c("protein_id", "sample", "maxlfq")

# Add metadata
se_filtered_long <- as.data.frame(merge(se_filtered_long, 
                                        col_data[,c("sample", "frailty")], 
                                        by = "sample"))

summary(se_filtered_long$maxlfq)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# 6773    55858    87660   164843   152299 79194328    11492

sd(se_filtered_long$maxlfq, na.rm = TRUE)
# 615188.3

# Plot maxlfq vs sample
p <- ggplot(se_filtered_long, aes(x = sample, y = maxlfq, fill = frailty)) + 
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Samples", 
       y = "Protein MaxLFQ intensities") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(), #element_text(angle = 90, size = 4.5)
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail")) +
  facet_wrap("frailty", scales = "free_x", ncol = 1) +
  ylim(0, 1000000)

save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/",
                    "boxplot_maxlfq_samples.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Plot maxlfq vs frailty group
p <- ggplot(se_filtered_long, aes(x = frailty, y = maxlfq, fill = frailty)) + 
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Frailty group", 
       y = "Protein MaxLFQ intensities") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 8),
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail"))

save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/",
                    "boxplot_maxlfq_groups.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Intensities density plot
p <- ggplot(se_filtered_long, aes(x = maxlfq)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  labs(title = "Distribution of protein MaxLFQ intensities",
       x = "Protein MaxLFQ intensities", 
       y = "Density") +
  theme_minimal(base_size = 12)

save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/",
                    "density_plot_maxlfq.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Get proteins with intensities higher than percentile 99
q99 <- quantile(se_filtered_long$maxlfq, 0.99, na.rm = TRUE)
se_filtered_long_q99 <- se_filtered_long[se_filtered_long$maxlfq > q99, ]
se_filtered_long_q99 <- se_filtered_long_q99[!is.na(se_filtered_long_q99), ]
prot_q99 <-  se_filtered_long_q99 %>% arrange(desc(maxlfq))
prot_q99 <- unique(prot_q99$protein_id)
prot_q99_ann <- GetProteinAnnontate(prot_q99,columns = c("gene_primary",
                                                         "organism_name",
                                                          "protein_name",
                                                          "cc_function"))
write.csv2(prot_q99_ann, paste0(work_path, "/prdea/output/prot_q99_ann.csv"))

# Intensities density plot per sample
plotDensities(se_filtered_assay, legend = FALSE)
for (i in 1:nrow(se_filtered_assay)){
  plotDensities(se_filtered_assay[, i], legend = TRUE)
}

# Compute row-wise mean and sd. Scatterplot to check mean-variance dependence.
# Red line: median estimator (window-width 10%). If horizontal, no dependence.
save_path <- paste0(work_path, "/prdea/plots/preprocessing/maxlfq/",
                    "scatterplot.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
meanSdPlot(se_filtered_assay)
dev.off()

# Q-Q sample vs theoretical quantiles normal distribution
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
#    Mode   FALSE    TRUE 
#logical     202      1

# ==================================
# 5. Normalization
# ==================================
# 5.1. vsn
# ==================================
fit <- vsnMatrix(se_filtered_assay)
se_norm <- se_filtered
assay(se_norm) <- predict(fit, se_filtered_assay)
se_norm_assay <- assay(se_norm)

# Save SE normalized by vsn
save_path <- paste0(work_path,"/prdea/data/maxlfq/se_norm.RData")
save(se_norm, file = save_path)

# Boxplot maxlfq after vsn
res <- create_lfq_boxplot(se_norm_assay, col_data)
p <- res$plot_samples
save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/",
                    "boxplot_maxlfq_samples_vsn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- res$plot_groups
save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/",
                    "boxplot_maxlfq_groups_vsn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

print(res$lfq_summary)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 11.98   15.88   16.50   16.66   17.24   24.92    11492 

# Intensities density plot
se_norm_long <- res$se_long

p <- ggplot(se_norm_long, aes(x = maxlfq)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  labs(title = "Distribution of protein MaxLFQ intensities",
       x = "Protein MaxLFQ intensities", 
       y = "Density") +
  theme_minimal(base_size = 12)

save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/density_plot_",
                    "maxlfq_vsn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Intensities density plot per sample
plotDensities(se_norm_assay, legend = FALSE)
for (i in 1:nrow(se_norm_assay)){
  plotDensities(se_norm_assay[, i], legend = TRUE)
}

# Scatterplot
save_path <- paste0(work_path, "/prdea/plots/preprocessing/maxlfq/",
                    "scatterplot_vsn.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
meanSdPlot(se_norm_assay)
dev.off()

# Normal contrast
p_values <- apply(se_norm_assay, 2,
                  function(column) shapiro.test(column)$p.value)
norm_results <- p_values >= 0.05 # TRUE if normal, FALSE if not
summary(norm_results)
#    Mode   FALSE    TRUE 
#logical     169      34 

# ==================================
# 5.2. Log 2 transform
# ==================================
# Just log2
se_log <- se_filtered
assay(se_log) <- log2(assay(se_log))
se_log_assay <- assay(se_log)

# Save SE log2
save_path <- paste0(work_path,"/prdea/data/maxlfq/se_log.RData")
save(se_log, file = save_path)

# Boxplot maxlfq after log2
res <- create_lfq_boxplot(se_log_assay, col_data)
p <- res$plot_samples
save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/boxplot_maxlfq_",
                    "samples_log2.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- res$plot_groups
save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/boxplot_maxlfq_",
                    "groups_log2.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

print(res$lfq_summary)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 12.73   15.77   16.42   16.60   17.22   26.24   11492 

print(res$lfq_sd)
# 1.20735

# Intensities density plot
se_log_long <- res$se_long

p <- ggplot(se_log_long, aes(x = maxlfq)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  labs(title = "Distribution of protein MaxLFQ intensities",
       x = "Protein MaxLFQ intensities", 
       y = "Density") +
  theme_minimal(base_size = 12)

save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/density_plot_",
                    "maxlfq_log2.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Intensities density plot per sample
plotDensities(se_log_assay, legend = FALSE)

# Scatterplot
save_path <- paste0(work_path, "/prdea/plots/preprocessing/maxlfq/",
                    "scatterplot_log2.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
meanSdPlot(se_log_assay)
dev.off()

# ==================================
# 6. Imputation
# ==================================
# 6.1. No imputation
# ==================================
se_no_imp <- se_log

# Save SE no imputated
save_path <- paste0(work_path,"/prdea/data/maxlfq/se_no_imp.RData")
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
save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/",
                    "pca_plot_perseus.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Boxplot maxlfq after perseus
res <- create_lfq_boxplot(se_perseus_assay, col_data)
p <- res$plot_samples
save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/",
                    "boxplot_maxlfq_samples_perseus.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- res$plot_groups
save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/",
                    "boxplot_maxlfq_groups_perseus.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

print(res$lfq_summary)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 11.54   14.81   15.87   15.94   16.82   26.24 

# Intensities density plot
se_perseus_long <- res$se_long
p <- ggplot(se_perseus_long, aes(x = maxlfq)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  labs(title = "Distribution of protein MaxLFQ intensities",
       x = "Protein MaxLFQ intensities", 
       y = "Density") +
  theme_minimal(base_size = 12)

save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/",
                    "density_plot_maxlfq_perseus.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Intensities density plot per sample
plotDensities(se_perseus_assay, legend = FALSE)

# Save SE imputed by perseus
save_path <- paste0(work_path,"/prdea/data/maxlfq/se_perseus.RData")
save(se_perseus, file = save_path)


# ==================================
# 6.2.2. KNN
# ==================================
# Delete samples with more than 80% missing values
se_knn <- se_no_imp
se_knn <- se_knn[, colMeans(is.na(assay(se_knn))) <= 0.8]

# Get colData from KNN SE
col_data <- as.data.frame(colData(se_knn))

# Impute using k-nearest neighbors
assay(se_knn) <- impute.knn(assay(se_knn), rowmax = 0.6)$data
se_knn_assay <- assay(se_knn)

# Plot PCA
p <- create_pca_plot(se_knn, n_top_loadings = 5)
save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/pca_plot_knn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Boxplot maxlfq after knn
res <- create_lfq_boxplot(se_knn_assay, col_data)
p <- res$plot_samples
save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/boxplot_maxlfq_",
                    "samples_knn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- res$plot_groups
save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/boxplot_maxlfq_",
                    "groups_knn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

print(res$lfq_summary)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 12.73   15.74   16.32   16.50   17.05   26.24 

# Intensities density plot
se_knn_long <- res$se_long
p <- ggplot(se_knn_long, aes(x = maxlfq)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  labs(title = "Distribution of protein MaxLFQ intensities",
       x = "Protein MaxLFQ intensities", 
       y = "Density") +
  theme_minimal(base_size = 12)

save_path <- paste0(work_path,"/prdea/plots/preprocessing/maxlfq/density_plot_",
                    "maxlfq_knn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Intensities density plot per sample
plotDensities(se_knn_assay, legend = FALSE)

# Save SE imputed by knn
save_path <- paste0(work_path,"/prdea/data/maxlfq/se_knn.RData")
save(se_knn, file = save_path)
