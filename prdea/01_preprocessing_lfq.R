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
source(paste0(work_path, "/functions/build_summarized_experiment.R"))
source(paste0(work_path, "/functions/create_bar_plot.R"))
source(paste0(work_path, "/functions/create_density_plot.R"))
source(paste0(work_path, "/functions/create_pca_plot.R"))
source(paste0(work_path, "/functions/create_lfq_boxplot.R"))

# Data, plot and output paths
data_path <- paste0(work_path, "/prdea/data/")
plot_path <- paste0(work_path, "/prdea/plots/")
output_path <- paste0(work_path, "/prdea/output/")

# Path to combined_protein.tsv and experiment_annotation.tsv
combined_protein <- paste0(data_path, "combined_protein.tsv")
experiment_ann <- paste0(data_path, "experiment_annotation.tsv")

# ==================================
# 1. Preprocessing MaxLFQ intensities
# ==================================

# Make SummarizedExperiment from maxlfq
se_maxlfq <- build_summarized_experiment(combined_file = combined_protein,
                                  experiment_ann = experiment_ann)

# Check class of se_maxlfq object
class(se_maxlfq)

# Check log2, exp, lfq_type and level
metadata(se_maxlfq)
#$log2transform
# FALSE
#$lfq_type
# "MaxLFQ"
#$level
# "protein"

# Check number of rows and columns
dim(se_maxlfq)
# 6854  274

# Phenotypic variables
colData(se_maxlfq)
names(colData(se_maxlfq))
# "file" "sample"      "sample_name" "condition"   "replicate"   "label"

# Check if more than one sample per case
table(colData(se_maxlfq)$sample_name)

# Intensity matrix
head(assay(se_maxlfq))

# Row names matrix
head(colnames(se_maxlfq))

# ==================================
# 1.1. Filtering frailty cohort
# ==================================

# Load metadatos
load(paste0(work_path, "/data/metadata.RData"))

# Keep only data from FT cohort
se_maxlfq <- se_maxlfq[, colData(se_maxlfq)$replicate %in%
                         metadata_ft$replicate]
unique(colData(se_maxlfq)$replicate)

# Check number of rows and columns
dim(se_maxlfq)
# 6854  203

# ==================================
# 1.2. Adding phenotypic data
# ==================================

# Create an aligned dataframe with the order of `replicate` in colData(se)
replicate_se <- colData(se_maxlfq)$replicate
replicate_metadata <- metadata_ft$replicate
metadata_ft <- metadata_ft[match(replicate_se, replicate_metadata), ]

# Add phenotypic variables to colData
colData(se_maxlfq)$frailty <- metadata_ft$frailty
colData(se_maxlfq)$sex <- metadata_ft$sex
colData(se_maxlfq)$education <- metadata_ft$education
colData(se_maxlfq)$alcohol <- metadata_ft$alcohol
colData(se_maxlfq)$tobacco <- metadata_ft$tobacco
colData(se_maxlfq)$diabetes <- metadata_ft$diabetes
colData(se_maxlfq)$chf <- metadata_ft$chf
colData(se_maxlfq)$af <- metadata_ft$af
colData(se_maxlfq)$hipfracture <- metadata_ft$hipfracture
colData(se_maxlfq)$depression <- metadata_ft$depression
colData(se_maxlfq)$osteoarthritis <- metadata_ft$osteoarthritis
colData(se_maxlfq)$sarcopenia <- metadata_ft$sarcopenia
colData(se_maxlfq)$ilef <- metadata_ft$ilef
colData(se_maxlfq)$age <- metadata_ft$age
colData(se_maxlfq)$medas <- metadata_ft$medas
colData(se_maxlfq)$energy <- metadata_ft$energy
colData(se_maxlfq)$bmi <- metadata_ft$bmi

# Check
head(colData(se_maxlfq))

# ==================================
# 1.3. Filtering
# ==================================
# 1.3.1. By Missing
# ==================================

# Get intensities matrix from SE
se_maxlfq_assay <- assay(se_maxlfq)

# Number of proteins
nrow(se_maxlfq_assay)
# 6854

# Number of NA values in assay
na_table <- table(is.na(se_maxlfq_assay))
na_table
# FALSE    TRUE 
# 61378 1329984

prop_na <- na_table["TRUE"] / sum(na_table)
prop_na
# 0.9558864

# Proteins with at least one intensity value
proteins_to_keep <- apply(se_maxlfq_assay, 1, function(x) !all(is.na(x)))

# Number of proteins to keep
sum(proteins_to_keep)
# 2902

# Number of proteins with all NA
sum(!proteins_to_keep)
# 3952 

# Filter proteins with all NA
se_maxlfq <- se_maxlfq[proteins_to_keep,]
se_maxlfq_assay <- se_maxlfq_assay[proteins_to_keep, ]

# Save SE
save_path <- paste0(data_path,"se_maxlfq.RData")
save(se_maxlfq, file = save_path)

# Get colData from SE
col_data <- as.data.frame(colData(se_maxlfq))

# Check number of proteins after filtering those with all NA
nrow(se_maxlfq_assay)
# 2902

# Number of NA values in assay after filtering those with all NA
na_table <- table(is.na(se_maxlfq_assay))
na_table
# FALSE   TRUE 
# 61378 527728

# NA proportion
prop_na <- na_table["TRUE"] / sum(na_table)
prop_na
# 0.8958116

# Number of proteins detected per sample
num_prot_sample <- apply(se_maxlfq_assay, 2, function(x) sum(!is.na(x)))
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
rownames(se_maxlfq_assay[!is.na(se_maxlfq_assay[,"YDD_248"]), ])

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

save_path <- paste0(plot_path, "preprocessing/density_plot_maxlfq_num_prot_",
                    "sample.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- ggplot(num_prot_sample, aes(x = frailty, y = proteins, fill = frailty)) + 
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, outlier.alpha = 0.5) +
  coord_cartesian(ylim = c(0, 2010)) +
  labs(x = "Frailty group", 
       y = "Number of proteins") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail"))

save_path <- paste0(plot_path, "preprocessing/boxplot_maxlfq_num_prot_",
                    "sample.png")
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
missval <- ifelse(is.na(se_maxlfq_assay), 0, 1)

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

save_path <- paste0(plot_path,"preprocessing/heatmap_maxlfq_na.png")
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

save_path <- paste0(plot_path,"preprocessing/heatmap_maxlfq_clust_na.png")
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
save_path <- paste0(plot_path,"preprocessing/heatmap_maxlfq_clust_na_ann.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
draw(p, heatmap_legend_list = legends, heatmap_legend_side = "top")
dev.off()

# Sample proportion for each protein
sample_proportion <- rowSums(!is.na(se_maxlfq_assay)) / ncol(se_maxlfq_assay)
sample_perce <- sample_proportion * 100
sample_perce <- data.frame(percentage = sample_perce)
p <- create_density_plot(sample_perce, "percentage",
                         "Sample proportion (%) per protein", 
                         "Distribution of protein detection across samples")

save_path <- paste0(plot_path, "preprocessing/density_plot_maxlfq_sample_",
                    "proportion.png")
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
rownames(se_maxlfq_assay[sample_proportion == 1,])
# P00761; trypsin
# P04264; Keratin, type II cytoskeletal 1
# P06702; Protein S100-A9: regulation of inflammatory processes and immune response

# Proteins in more than 75% of samples
rownames(se_maxlfq_assay[sample_proportion >= 0.75,])

# Proteins in less than 50% of samples
rownames(se_maxlfq_assay[sample_proportion < 0.5,])


# Keep proteins with minimum fraction of valid values in at least one condition
min_fraction <- 0.5
fraction_ft <- apply(se_maxlfq_assay[, col_data$frailty == "FT"], 1,
                     function(x) sum(!is.na(x)) / length(x))
fraction_nft <- apply(se_maxlfq_assay[, col_data$frailty == "NFT"], 1,
                      function(x) sum(!is.na(x)) / length(x))

proteins_to_keep <- (fraction_ft >= min_fraction) | (fraction_nft >= 
                                                       min_fraction)

# Number of proteins present in >= 50% of samples in at least one condition
sum(proteins_to_keep) 
# 185

# Filter SE
se_maxlfq_filt_miss <- se_maxlfq[proteins_to_keep, ]

# Get protein intensity matrix
se_maxlfq_filt_miss_assay <- assay(se_maxlfq_filt_miss)

# Check number of proteins after filtering by sample proportion
nrow(se_maxlfq_filt_miss_assay)
# 185

# ==================================
# 1.3.2. By organism
# ==================================
# Get lineage from UniProt
lineage <- GetProteinAnnontate(rownames(se_maxlfq_filt_miss_assay),
                               columns = c("lineage"))
lineage_df <- as.data.frame(lineage)
lineage_df$protein_id <- rownames(se_maxlfq_filt_miss_assay)
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
se_maxlfq_filt <- se_maxlfq_filt_miss[proteins_to_keep, ]

# Save SE filtered
save_path <- paste0(data_path,"se_maxlfq_filt.RData")
save(se_maxlfq_filt, file = save_path)

# Get protein intensity matrix
se_maxlfq_filt_assay <- assay(se_maxlfq_filt)

# Check number of proteins after filtering no bacteria or homo
nrow(se_maxlfq_filt_assay)
# 182

# Number of NA values in assay after filtering by sample proportion and organism
na_table <- table(is.na(se_maxlfq_filt_assay))
na_table
# FALSE  TRUE 
# 25454 11492

# NA proportion
prop_na <- na_table["TRUE"] / sum(na_table)
prop_na
# 0.3110486

# Number of proteins detected per sample
num_prot_sample <- apply(se_maxlfq_filt_assay, 2, function(x) sum(!is.na(x)))
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

save_path <- paste0(plot_path, "preprocessing/density_plot_maxlfq_num_prot_",
                    "sample_filt.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- ggplot(num_prot_sample, aes(x = frailty, y = proteins, fill = frailty)) + 
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, outlier.alpha = 0.5) +
  coord_cartesian(ylim = c(0, 450)) +
  labs(x = "Frailty group", 
       y = "Number of proteins") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail"))

save_path <- paste0(plot_path, "preprocessing/boxplot_maxlfq_num_prot_sample_",
                    "filt.png")
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
missval <- ifelse(is.na(se_maxlfq_filt_assay), 0, 1)

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

save_path <- paste0(plot_path, "preprocessing/heatmap_maxlfq_na_filt.png")
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

save_path <- paste0(plot_path, "preprocessing/heatmap_maxlfq_clust_na_filt_",
                    "min.png")
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
save_path <- paste0(plot_path,"preprocessing/heatmap_maxlfq_clust_na_ann_",
                    "filt.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
draw(p, heatmap_legend_list = legends, heatmap_legend_side = "top")
dev.off()

# Transform se_maxlfq_filt assay to long
se_maxlfq_filt_df <- as.data.frame(se_maxlfq_filt_assay)
se_maxlfq_filt_df$protein_id <- rownames(se_maxlfq_filt_df)
rownames(se_maxlfq_filt_df) <- NULL
se_maxlfq_filt_long <- melt(se_maxlfq_filt_df)
colnames(se_maxlfq_filt_long) <- c("protein_id", "sample", "maxlfq")

# Add metadata
se_maxlfq_filt_long <- as.data.frame(merge(se_maxlfq_filt_long,
                                           col_data[,c("sample", "frailty")],
                                           by = "sample"))

summary(se_maxlfq_filt_long$maxlfq)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# 6773    55858    87660   164843   152299 79194328    11492

sd(se_maxlfq_filt_long$maxlfq, na.rm = TRUE)
# 615188.3

# Plot maxlfq vs sample
p <- ggplot(se_maxlfq_filt_long, aes(x = sample, y = maxlfq, fill = frailty)) + 
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

save_path <- paste0(plot_path,"preprocessing/boxplot_maxlfq_samples.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Plot maxlfq vs frailty group
p <- ggplot(se_maxlfq_filt_long, aes(x = frailty, y = maxlfq, fill = frailty)) + 
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

save_path <- paste0(plot_path,"preprocessing/boxplot_maxlfq_groups.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Intensities density plot
p <- ggplot(se_maxlfq_filt_long, aes(x = maxlfq)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  coord_cartesian(ylim = c(0, 0.000025)) +
  labs(title = "Distribution of protein MaxLFQ intensities",
       x = "Protein MaxLFQ intensities", 
       y = "Density") +
  theme_minimal(base_size = 12)

save_path <- paste0(plot_path,"preprocessing/density_plot_maxlfq.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Get proteins with intensities higher than percentile 99
q99 <- quantile(se_maxlfq_filt_long$maxlfq, 0.99, na.rm = TRUE)
se_maxlfq_filt_long_q99 <- se_maxlfq_filt_long[se_maxlfq_filt_long$maxlfq > q99,]
se_maxlfq_filt_long_q99 <- se_maxlfq_filt_long_q99[!is.na(se_maxlfq_filt_long_q99),]
prot_q99 <-  se_maxlfq_filt_long_q99 %>% arrange(desc(maxlfq))
prot_q99 <- unique(prot_q99$protein_id)
prot_q99_ann <- GetProteinAnnontate(prot_q99,columns = c("gene_primary",
                                                         "organism_name",
                                                          "protein_name",
                                                          "cc_function"))
write.csv2(prot_q99_ann, paste0(output_path, "prot_q99_ann_maxlfq.csv"))

# Intensities density plot per sample
plotDensities(se_maxlfq_filt_assay, legend = FALSE)
for (i in 1:nrow(se_maxlfq_filt_assay)){
  plotDensities(se_maxlfq_filt_assay[, i], legend = TRUE)
}

# Compute row-wise mean and sd. Scatterplot to check mean-variance dependence.
# Red line: median estimator (window-width 10%). If horizontal, no dependence.
save_path <- paste0(plot_path, "preprocessing/scatterplot_maxlfq.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
meanSdPlot(se_maxlfq_filt_assay)
dev.off()

# Q-Q sample vs theoretical quantiles normal distribution
for (i in 1:ncol(se_maxlfq_filt_assay)) {
  qqnorm(se_maxlfq_filt_assay[, i],
         main = paste("QQ Plot -", colnames(se_maxlfq_filt_assay)[i]))
  qqline(se_maxlfq_filt_assay[, i], col = "red")
}

# Normal contrast
p_values <- apply(se_maxlfq_filt_assay, 2,
                  function(column) shapiro.test(column)$p.value)
norm_results <- p_values >= 0.05 # TRUE if normal, FALSE if not
summary(norm_results)
#    Mode   FALSE    TRUE 
#logical     202      1

# ==================================
# 1.4. Normalization
# ==================================
# 1.4.1. vsn
# ==================================
fit <- vsnMatrix(se_maxlfq_filt_assay)
se_maxlfq_norm <- se_maxlfq_filt
assay(se_maxlfq_norm) <- predict(fit, se_maxlfq_filt_assay)
se_maxlfq_norm_assay <- assay(se_maxlfq_norm)

# Save SE normalized by vsn
save_path <- paste0(data_path,"se_maxlfq_norm.RData")
save(se_maxlfq_norm, file = save_path)

# Boxplot maxlfq after vsn
res <- create_lfq_boxplot(se_maxlfq_norm_assay, col_data)
p <- res$plot_samples
save_path <- paste0(plot_path,"preprocessing/boxplot_maxlfq_samples_vsn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- res$plot_groups
save_path <- paste0(plot_path,"preprocessing/boxplot_maxlfq_groups_vsn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

print(res$lfq_summary)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 11.98   15.88   16.50   16.66   17.24   24.92    11492 

# Intensities density plot
se_maxlfq_norm_long <- res$se_long
p <- ggplot(se_maxlfq_norm_long, aes(x = maxlfq)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  labs(title = "Distribution of protein MaxLFQ intensities",
       x = "Protein MaxLFQ intensities", 
       y = "Density") +
  theme_minimal(base_size = 12)

save_path <- paste0(plot_path,"preprocessing/density_plot_maxlfq_vsn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Intensities density plot per sample
plotDensities(se_maxlfq_norm_assay, legend = FALSE)
for (i in 1:nrow(se_maxlfq_norm_assay)){
  plotDensities(se_maxlfq_norm_assay[, i], legend = TRUE)
}

# Scatterplot
save_path <- paste0(plot_path, "preprocessing/scatterplot_maxlfq_vsn.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
meanSdPlot(se_maxlfq_norm_assay)
dev.off()

# Normal contrast
p_values <- apply(se_maxlfq_norm_assay, 2,
                  function(column) shapiro.test(column)$p.value)
norm_results <- p_values >= 0.05 # TRUE if normal, FALSE if not
summary(norm_results)
#    Mode   FALSE    TRUE 
#logical     169      34 

# ==================================
# 1.4.2. Log 2 transform
# ==================================
# Just log2
se_maxlfq_log <- se_maxlfq_filt
assay(se_maxlfq_log) <- log2(assay(se_maxlfq_log))
se_maxlfq_log_assay <- assay(se_maxlfq_log)

# Save SE log2
save_path <- paste0(data_path,"se_maxlfq_log.RData")
save(se_maxlfq_log, file = save_path)

# Boxplot maxlfq after log2
res <- create_lfq_boxplot(se_maxlfq_log_assay, col_data, y_lim = c(7.5, 28),
                          y_lab = "Protein log2 MaxLFQ intensity")
p <- res$plot_samples
save_path <- paste0(plot_path,"preprocessing/boxplot_maxlfq_samples_log2.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- res$plot_groups
save_path <- paste0(plot_path,"preprocessing/boxplot_maxlfq_groups_log2.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

print(res$lfq_summary)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 12.73   15.77   16.42   16.60   17.22   26.24   11492 

print(res$lfq_sd)
# 1.20735

# Intensities density plot
se_maxlfq_log_long <- res$se_long

p <- ggplot(se_maxlfq_log_long, aes(x = maxlfq)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  coord_cartesian(xlim = c(7.5, 27.9), ylim = c(0, 0.4)) +
  labs(title = "Distribution of protein MaxLFQ intensities",
       x = "Protein log2 MaxLFQ intensity", 
       y = "Density") +
  theme_minimal(base_size = 12)

save_path <- paste0(plot_path, "preprocessing/density_plot_maxlfq_log2.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Intensities density plot per sample
plotDensities(se_maxlfq_log_assay, legend = FALSE)

# Scatterplot
save_path <- paste0(plot_path, "preprocessing/scatterplot_maxlfq_log2.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
meanSdPlot(se_maxlfq_log_assay)
dev.off()

# ==================================
# 1.5. Imputation
# ==================================
# 1.5.1. No imputation
# ==================================
se_maxlfq_no_imp <- se_maxlfq_log

# Save SE no imputated
save_path <- paste0(data_path,"se_maxlfq_no_imp.RData")
save(se_maxlfq_no_imp, file = save_path)

# ==================================
# 1.5.2. Perseus
# ==================================
# Missing values are replaced with random values generated from a shifted and
# scaled normal distribution based on the existing data
se_maxlfq_perse <- manual_impute(se_maxlfq_no_imp)
se_maxlfq_perse_assay <- assay(se_maxlfq_perse)

# Plot PCA
p <- create_pca_plot(se_maxlfq_perse, n_top_loadings = 5)
save_path <- paste0(plot_path,"preprocessing/pca_maxlfq_perseus.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Boxplot maxlfq after perseus
res <- create_lfq_boxplot(se_maxlfq_perse_assay, col_data)
p <- res$plot_samples
save_path <- paste0(plot_path,"preprocessing/boxplot_maxlfq_samples_",
                    "perseus.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- res$plot_groups
save_path <- paste0(plot_path,"preprocessing/boxplot_maxlfq_groups_perseus.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

print(res$lfq_summary)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 11.54   14.81   15.87   15.94   16.82   26.24 

# Intensities density plot
se_maxlfq_perse_long <- res$se_long
p <- ggplot(se_maxlfq_perse_long, aes(x = maxlfq)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  labs(title = "Distribution of protein MaxLFQ intensities",
       x = "Protein MaxLFQ intensities", 
       y = "Density") +
  theme_minimal(base_size = 12)

save_path <- paste0(plot_path,"preprocessing/density_plot_maxlfq_perseus.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Intensities density plot per sample
plotDensities(se_maxlfq_perse_assay, legend = FALSE)

# Save SE imputed by perseus
save_path <- paste0(data_path,"se_maxlfq_perse.RData")
save(se_maxlfq_perse, file = save_path)


# ==================================
# 1.5.3. KNN
# ==================================
# Delete samples with more than 80% missing values
se_maxlfq_knn <- se_maxlfq_no_imp
se_maxlfq_knn <- se_maxlfq_knn[, colMeans(is.na(assay(se_maxlfq_knn))) <= 0.8]

# Get colData from KNN SE
col_data <- as.data.frame(colData(se_maxlfq_knn))

# Impute using k-nearest neighbors
assay(se_maxlfq_knn) <- impute.knn(assay(se_maxlfq_knn), rowmax = 0.6)$data
se_maxlfq_knn_assay <- assay(se_maxlfq_knn)

# Plot PCA
p <- create_pca_plot(se_maxlfq_knn, n_top_loadings = 5)
save_path <- paste0(plot_path,"preprocessing/pca_maxlfq_knn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Boxplot maxlfq after knn
res <- create_lfq_boxplot(se_maxlfq_knn_assay, col_data)
p <- res$plot_samples
save_path <- paste0(plot_path,"preprocessing/boxplot_maxlfq_samples_knn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- res$plot_groups
save_path <- paste0(plot_path,"preprocessing/boxplot_maxlfq_groups_knn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

print(res$lfq_summary)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 12.73   15.74   16.32   16.50   17.05   26.24 

# Intensities density plot
se_maxlfq_knn_long <- res$se_long
p <- ggplot(se_maxlfq_knn_long, aes(x = maxlfq)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  labs(title = "Distribution of protein MaxLFQ intensities",
       x = "Protein MaxLFQ intensities", 
       y = "Density") +
  theme_minimal(base_size = 12)

save_path <- paste0(plot_path,"preprocessing/density_plot_maxlfq_knn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Intensities density plot per sample
plotDensities(se_maxlfq_knn_assay, legend = FALSE)

# Save SE imputed by knn
save_path <- paste0(data_path,"se_maxlfq_knn.RData")
save(se_maxlfq_knn, file = save_path)

# ==================================
# 2. Preprocessing Top-N intensities
# ==================================
# Make SummarizedExperiment from lfq
se_lfq <- build_summarized_experiment(combined_file = combined_protein,
                                      experiment_ann = experiment_ann,
                                      lfq_type = "LFQ")

# Check class of se object
class(se_lfq)

# Check log2, exp, lfq_type and level
metadata(se_lfq)
#$log2transform
# FALSE

#$lfq_type
# "Intensity"

#$level
# "protein"

# Check number of rows and columns
dim(se_lfq)
# 6854  274

# Phenotypic variables
colData(se_lfq)
names(colData(se_lfq))
# "file" "sample"      "sample_name" "condition"   "replicate"   "label"

# Check if more than one sample per case
table(colData(se_lfq)$sample_name)

# Intensity matrix
head(assay(se_lfq))

# Row names matrix
head(colnames(se_lfq))

# ==================================
# 2.1. Filtering frailty cohort
# ==================================

# Keep only data from FT cohort
se_lfq <- se_lfq[, colData(se_lfq)$replicate %in% metadata_ft$replicate]
unique(colData(se_lfq)$replicate)

# Check number of rows and columns
dim(se_lfq)
# 6854  203

# ==================================
# 2.2. Adding phenotypic data
# ==================================

# Add phenotypic variables to colData
colData(se_lfq)$frailty <- metadata_ft$frailty
colData(se_lfq)$sex <- metadata_ft$sex
colData(se_lfq)$education <- metadata_ft$education
colData(se_lfq)$alcohol <- metadata_ft$alcohol
colData(se_lfq)$tobacco <- metadata_ft$tobacco
colData(se_lfq)$diabetes <- metadata_ft$diabetes
colData(se_lfq)$chf <- metadata_ft$chf
colData(se_lfq)$af <- metadata_ft$af
colData(se_lfq)$hipfracture <- metadata_ft$hipfracture
colData(se_lfq)$depression <- metadata_ft$depression
colData(se_lfq)$osteoarthritis <- metadata_ft$osteoarthritis
colData(se_lfq)$sarcopenia <- metadata_ft$sarcopenia
colData(se_lfq)$ilef <- metadata_ft$ilef
colData(se_lfq)$age <- metadata_ft$age
colData(se_lfq)$medas <- metadata_ft$medas
colData(se_lfq)$energy <- metadata_ft$energy
colData(se_lfq)$bmi <- metadata_ft$bmi

# Check
head(colData(se_lfq))

# ==================================
# 2.3. Filtering
# ==================================
# 2.3.1. By Missing
# ==================================

# Get intensities matrix from SE
se_lfq_assay <- assay(se_lfq)

# Number of proteins
nrow(se_lfq_assay)
# 6854

# Number of NA values in assay
na_table <- table(is.na(se_lfq_assay))
na_table
# FALSE    TRUE 
# 148253 1243109

prop_na <- na_table["TRUE"] / sum(na_table)
prop_na
# 0.8934476

# Proteins with at least one intensity value
proteins_to_keep <- apply(se_lfq_assay, 1, function(x) !all(is.na(x)))

# Number of proteins to keep
sum(proteins_to_keep)
# 6400

# Number of proteins with all NA
sum(!proteins_to_keep)
# 454 

# Filter proteins with all NA
se_lfq <- se_lfq[proteins_to_keep,]
se_lfq_assay <- se_lfq_assay[proteins_to_keep, ]

# Save SE
save_path <- paste0(data_path,"se_lfq.RData")
save(se_lfq, file = save_path)

# Get colData from SE
col_data <- as.data.frame(colData(se_lfq))

# Check number of proteins after filtering those with all NA
nrow(se_lfq_assay)
# 6400

# Number of NA values in assay after filtering those with all NA
na_table <- table(is.na(se_lfq_assay))
na_table
# FALSE    TRUE 
# 148253 1150947

# NA proportion
prop_na <- na_table["TRUE"] / sum(na_table)
prop_na
# 0.885889

# Number of proteins detected per sample
num_prot_sample <- apply(se_lfq_assay, 2, function(x) sum(!is.na(x)))
num_prot_sample <- as.data.frame(num_prot_sample)
colnames(num_prot_sample) <- "proteins"
num_prot_sample$sample <- rownames(num_prot_sample)
rownames(num_prot_sample) <- NULL
# Add frailty group
num_prot_sample <- merge(num_prot_sample, col_data[, c("sample", "frailty")], 
                         by = "sample")

# Summary total proteins detected per sample
summary(num_prot_sample)
# Sample             Proteins      frailty  
# Length:203         Min.   :  40.0   NFT:138  
# Class :character   1st Qu.: 520.5   FT : 65  
# Mode  :character   Median : 720.0            
#                    Mean   : 730.3            
#                    3rd Qu.: 913.0            
#                    Max.   :2054.0            

# Check proteins detected in samples with a low detection
rownames(se_lfq_assay[!is.na(se_lfq_assay[,"YDD_248"]), ])

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

save_path <- paste0(plot_path, "preprocessing/density_plot_lfq_num_prot_",
                    "sample.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- ggplot(num_prot_sample, aes(x = frailty, y = proteins, fill = frailty)) + 
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, outlier.alpha = 0.5) +
  coord_cartesian(ylim = c(0, 2010)) +
  labs(x = "Frailty group", 
       y = "Number of proteins") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail"))

save_path <- paste0(plot_path, "preprocessing/boxplot_lfq_num_prot_sample.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)


# Mean number of proteins in frail
mean(num_prot_sample[num_prot_sample$frailty == "FT", "proteins"])
# 660.2308
# Standard deviation number of proteins in frail
sd(num_prot_sample[num_prot_sample$frailty == "FT", "proteins"])
# 335.1852

# Mean number of proteins in non-frail
mean(num_prot_sample[num_prot_sample$frailty == "NFT", "proteins"])
# 763.3188
# Standard deviation number of proteins in non-frail
sd(num_prot_sample[num_prot_sample$frailty == "NFT", "proteins"])
# 361.1816

# Wilcox test number of proteins FT vs NFT
wilcox.test(proteins ~ frailty, data = num_prot_sample)
# Wilcoxon rank sum test with continuity correction
# data:  proteins by frailty
# W = 5356.5, p-value = 0.02571
# alternative hypothesis: true location shift is not equal to 0

# Make binary assay (1 = valid value, 0 = missing value)
missval <- ifelse(is.na(se_lfq_assay), 0, 1)

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

save_path <- paste0(plot_path,"preprocessing/heatmap_lfq_na.png")
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

save_path <- paste0(plot_path,"preprocessing/heatmap_lfq_clust_na.png")
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
save_path <- paste0(plot_path,"preprocessing/heatmap_lfq_clust_na_ann.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
draw(p, heatmap_legend_list = legends, heatmap_legend_side = "top")
dev.off()

# Sample proportion for each protein
sample_proportion <- rowSums(!is.na(se_lfq_assay)) / ncol(se_lfq_assay)
sample_perce <- sample_proportion * 100
sample_perce <- data.frame(percentage = sample_perce)
p <- create_density_plot(sample_perce, "percentage",
                         "Sample proportion (%) per protein",
                         "Distribution of protein detection across samples")

save_path <- paste0(plot_path, "preprocessing/density_plot_lfq_sample_",
                    "proportion.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Number of proteins present in more than 25%
sum(sample_proportion > 0.25)
# 897

# Number of proteins present in more than 50%
sum(sample_proportion > 0.5)
# 403

# Number of proteins present in more than 75%
sum(sample_proportion > 0.75)
# 171

# Number of proteins present in 100%
sum(sample_proportion == 1)
# 9

# Proteins present in 100% of samples
rownames(se_lfq_assay[sample_proportion == 1,])

# Proteins in more than 75% of samples
rownames(se_lfq_assay[sample_proportion >= 0.75,])

# Proteins in less than 50% of samples
rownames(se_lfq_assay[sample_proportion < 0.5,])

# Keep proteins with minimum fraction of valid values in at least one condition
min_fraction <- 0.5
fraction_ft <- apply(se_lfq_assay[,col_data$frailty == "FT"], 1,
                     function(x) sum(!is.na(x)) / length(x))
fraction_nft <- apply(se_lfq_assay[, col_data$frailty == "NFT"], 1,
                      function(x) sum(!is.na(x)) / length(x))

proteins_to_keep <- (fraction_ft >= min_fraction) | (fraction_nft >= 
                                                       min_fraction)

# Number of proteins present in >= 50% of samples in at least one condition
sum(proteins_to_keep) 
# 450

# Filter SE
se_lfq_filt_miss <- se_lfq[proteins_to_keep, ]

# Get protein intensities matrix
se_lfq_filt_miss_assay <- assay(se_lfq_filt_miss)

# Check number of proteins after filtering by sample proportion
nrow(se_lfq_filt_miss_assay)
# 450

# ==================================
# 2.3.2. By organism
# ==================================
# Get lineage from UniProt
lineage <- GetProteinAnnontate(rownames(se_lfq_filt_miss_assay),
                               columns = c("lineage"))
lineage_df <- as.data.frame(lineage)
lineage_df$protein_id <- rownames(se_lfq_filt_miss_assay)
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
se_lfq_filt <- se_lfq_filt_miss[proteins_to_keep, ]

# Save SE filtered
save_path <- paste0(data_path,"se_lfq_filt.RData")
save(se_lfq_filt, file = save_path)

# Get  protein intensity matrix
se_lfq_filt_assay <- assay(se_lfq_filt)

# Check number of proteins after filtering no bacteria or homo
nrow(se_lfq_filt_assay)
# 442

# Number of NA values in assay after filtering by sample proportion and organism
na_table <- table(is.na(se_lfq_filt_assay))
na_table
# FALSE  TRUE 
# 62884 26842 

# NA proportion
prop_na <- na_table["TRUE"] / sum(na_table)
prop_na
# 0.2992552

# Number of proteins detected per sample
num_prot_sample <- apply(se_lfq_filt_assay, 2, function(x) sum(!is.na(x)))
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
# Length:203         Min.   : 34.0   NFT:138  
# Class :character   1st Qu.:257.5   FT : 65  
# Mode  :character   Median :331.0            
#                    Mean   :309.8            
#                    3rd Qu.:379.5            
#                    Max.   :438.0     

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

save_path <- paste0(plot_path, "preprocessing/density_plot_lfq_num_prot_",
                    "sample_filt.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- ggplot(num_prot_sample, aes(x = frailty, y = proteins, fill = frailty)) + 
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, outlier.alpha = 0.5) +
  coord_cartesian(ylim = c(0, 450)) +
  labs(x = "Frailty group", 
       y = "Number of proteins") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail"))

save_path <- paste0(plot_path, "preprocessing/boxplot_lfq_num_prot_sample_",
                    "filt.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Mean number of proteins in frail
mean(num_prot_sample[num_prot_sample$frailty == "FT", "proteins"])
# 292.8615
# Standard deviation number of proteins in frail
sd(num_prot_sample[num_prot_sample$frailty == "FT", "proteins"])
# 83.60059

# Mean number of proteins in non-frail
mean(num_prot_sample[num_prot_sample$frailty == "NFT", "proteins"])
# 317.7391
# Standard deviation number of proteins in non-frail
sd(num_prot_sample[num_prot_sample$frailty == "NFT", "proteins"])
# 92.70598

# Wilcox test number of proteins FT vs NFT
wilcox.test(proteins ~ frailty, data = num_prot_sample)
# Wilcoxon rank sum test with continuity correction
# data:  proteins by frailty
# W = 5500.5, p-value = 0.00934
# alternative hypothesis: true location shift is not equal to 0

# Make binary assay
missval <- ifelse(is.na(se_lfq_filt_assay), 0, 1)

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

save_path <- paste0(plot_path, "preprocessing/heatmap_lfq_na_filt.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Complex Heatmap
# Clustering rows and columns by binary distance
# Column split by frailty groups
p <- Heatmap(missval, col = c("white", "black"),
             show_row_dend = FALSE, 
             column_names_gp = gpar(fontsize = 4),
             show_row_names = FALSE,
             # cluster_columns = FALSE,
             # cluster_rows = FALSE,
             column_split = col_data$frailty,
             clustering_distance_rows = "binary",
             clustering_distance_columns = "binary",
             use_raster = FALSE)

save_path <- paste0(plot_path, "preprocessing/heatmap_lfq_clust_na_filt.png")
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
save_path <- paste0(plot_path,"preprocessing/heatmap_lfq_clust_na_ann_filt.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
draw(p, heatmap_legend_list = legends, heatmap_legend_side = "top")
dev.off()

# Transform se_lfq_filt assay to long
se_lfq_filt_df <- as.data.frame(se_lfq_filt_assay)
se_lfq_filt_df$protein_id <- rownames(se_lfq_filt_df)
rownames(se_lfq_filt_df) <- NULL
se_lfq_filt_long <- melt(se_lfq_filt_df)
colnames(se_lfq_filt_long) <- c("protein_id", "sample", "intensity")

# Add metadata
se_lfq_filt_long <- as.data.frame(merge(se_lfq_filt_long, 
                                        col_data[,c("sample", "frailty")], 
                                        by = "sample"))

summary(se_lfq_filt_long$intensity)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
# 192     42321     90238    284668    211653 643911360     26842 

sd(se_lfq_filt_long$intensity, na.rm = TRUE)
# 2956390

# Plot intensity vs sample
p <- ggplot(se_lfq_filt_long, aes(x = sample, y = intensity, fill = frailty)) + 
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Samples", 
       y = "Protein Top-N intensities") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(), #element_text(angle = 90, size = 4.5)
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail")) +
  facet_wrap("frailty", scales = "free_x", ncol = 1) +
  ylim(0, 1000000)

save_path <- paste0(plot_path,"preprocessing/boxplot_lfq_samples.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Plot intensity vs frailty group
p <- ggplot(se_lfq_filt_long, aes(x = frailty, y = intensity, fill = frailty)) + 
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Frailty group", 
       y = "Protein Top-N intensities") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 8),
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail"))

save_path <- paste0(plot_path,"preprocessing/boxplot_lfq_groups.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Intensities density plot
p <- ggplot(se_lfq_filt_long, aes(x = intensity)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  coord_cartesian(ylim = c(0, 0.000025)) +
  labs(title = "Distribution of protein Top-N intensities",
       x = "Protein Top-N intensities", 
       y = "Density") +
  theme_minimal(base_size = 12)

save_path <- paste0(plot_path,"preprocessing/density_plot_lfq.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Get proteins with intensities higher than percentile 99
q99 <- quantile(se_lfq_filt_long$intensity, 0.99, na.rm = TRUE)
se_lfq_filt_long_q99 <- se_lfq_filt_long[se_lfq_filt_long$intensity > q99,]
se_lfq_filt_long_q99 <- se_lfq_filt_long_q99[!is.na(se_lfq_filt_long_q99),]
prot_q99 <-  se_lfq_filt_long_q99 %>% arrange(desc(intensity))
prot_q99 <- unique(prot_q99$protein_id)
prot_q99_ann <- GetProteinAnnontate(prot_q99,columns = c("gene_primary",
                                                         "organism_name",
                                                         "protein_name",
                                                         "cc_function"))
write.csv2(prot_q99_ann, paste0(output_path, "prot_q99_ann_lfq.csv"))

# Intensities density plot per sample
plotDensities(se_lfq_filt_assay, legend = FALSE)
for (i in 1:nrow(se_lfq_filt_assay)){
  plotDensities(se_lfq_filt_assay[, i], legend = TRUE)
}

# Scatterplot
save_path <- paste0(plot_path, "preprocessing/scatterplot_lfq.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
meanSdPlot(se_lfq_filt_assay)
dev.off()

# Q-Q
for (i in 1:ncol(se_lfq_filt_assay)) {
  qqnorm(se_lfq_filt_assay[, i],
         main = paste("QQ Plot -", colnames(se_lfq_filt_assay)[i]))
  qqline(se_lfq_filt_assay[, i], col = "red")
}

# Normal contrast
p_values <- apply(se_lfq_filt_assay, 2,
                  function(column) shapiro.test(column)$p.value)
norm_results <- p_values >= 0.05 # TRUE if normal, FALSE if not
summary(norm_results)
#    Mode   FALSE
#logical     203

# ==================================
# 2.4. Normalization
# ==================================
# 2.4.1. vsn
# ==================================

fit <- vsnMatrix(se_lfq_filt_assay)
se_lfq_norm <- se_lfq_filt
assay(se_lfq_norm) <- predict(fit, se_lfq_filt_assay)
se_lfq_norm_assay <- assay(se_lfq_norm)

# Save SE normalized by vsn
save_path <- paste0(data_path,"se_lfq_norm.RData")
save(se_lfq_norm, file = save_path)


# Boxplot lfq after vsn
res <- create_lfq_boxplot(se_lfq_norm_assay, col_data, lfq_type = "intensity", 
                          y_lim = c(7.5, 28),
                          y_lab = "Protein VSN Top-N intensity")
p <- res$plot_samples
save_path <- paste0(plot_path,"preprocessing/boxplot_lfq_samples_vsn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- res$plot_groups
save_path <- paste0(plot_path,"preprocessing/boxplot_lfq_groups_vsn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

print(res$lfq_summary)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 7.579  15.350  16.430  16.580  17.626  27.863   26842

print(res$lfq_sd)
# 1.667724

# Intensities density plot
se_lfq_norm_long <- res$se_long
p <- ggplot(se_lfq_norm_long, aes(x = intensity)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  coord_cartesian(xlim = c(7.5, 27.9), ylim = c(0, 0.4)) +
  labs(title = "Distribution of protein Top-N intensities",
       x = "Protein VSN Top-N intensity", 
       y = "Density") +
  theme_minimal(base_size = 12)

save_path <- paste0(plot_path,"preprocessing/density_plot_lfq_vsn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Intensities density plot per sample
plotDensities(se_lfq_norm_assay, legend = FALSE)
for (i in 1:nrow(se_lfq_norm_assay)){
  plotDensities(se_lfq_norm_assay[, i], legend = TRUE)
}

# Scatterplot
save_path <- paste0(plot_path, "preprocessing/scatterplot_lfq_vsn.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
meanSdPlot(se_lfq_norm_assay)
dev.off()

# Normal contrast
p_values <- apply(se_lfq_norm_assay, 2,
                  function(column) shapiro.test(column)$p.value)
norm_results <- p_values >= 0.05 # TRUE if normal, FALSE if not
summary(norm_results)
#    Mode   FALSE    TRUE 
#logical     189      14 

# ==================================
# 2.5. Imputation
# ==================================
# 2.5.1. No imputation
# ==================================
se_lfq_no_imp <- se_lfq_norm

# Save SE no imputated
save_path <- paste0(data_path,"se_lfq_no_imp.RData")
save(se_lfq_no_imp, file = save_path)

# ==================================
# 2.5.2. Perseus
# ==================================
# Missing values are replaced with random values generated from a shifted and
# scaled normal distribution based on the existing data
se_lfq_perse <- manual_impute(se_lfq_no_imp)
se_lfq_perse_assay <- assay(se_lfq_perse)

# Plot PCA
p <- create_pca_plot(se_lfq_perse, n_top_loadings = 5)
save_path <- paste0(plot_path,"preprocessing/pca_lfq_perseus.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Boxplot lfq after perseus
res <- create_lfq_boxplot(se_lfq_perse_assay, col_data, lfq_type = "intensity",
                          y_lab = "Protein Top-N intensities")
p <- res$plot_samples
save_path <- paste0(plot_path,"preprocessing/boxplot_lfq_samples_perseus.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- res$plot_groups
save_path <- paste0(plot_path,"preprocessing/boxplot_lfq_groups_perseus.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

print(res$lfq_summary)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7.579  13.829  15.516  15.598  17.075  27.863 

# Intensities density plot
se_lfq_perse_long <- res$se_long
p <- ggplot(se_lfq_perse_long, aes(x = intensity)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  labs(title = "Distribution of protein Top-N intensities",
       x = "Protein Top-N intensities", 
       y = "Density") +
  theme_minimal(base_size = 12)

save_path <- paste0(plot_path,"preprocessing/density_plot_lfq_perseus.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Intensities density plot per sample
plotDensities(se_lfq_perse_assay, legend = FALSE)

# Save SE imputed by perseus
save_path <- paste0(data_path,"se_lfq_perse.RData")
save(se_lfq_perse, file = save_path)


# ==================================
# 2.5.3. KNN
# ==================================
# Delete samples with more than 80% missing values
se_lfq_knn <- se_lfq_no_imp
se_lfq_knn <- se_lfq_knn[, colMeans(is.na(assay(se_lfq_knn))) <= 0.8]

# Impute using k-nearest neighbors
assay(se_lfq_knn) <- impute.knn(assay(se_lfq_knn), rowmax = 0.6)$data
se_lfq_knn_assay <- assay(se_lfq_knn)

# Plot PCA
p <- create_pca_plot(se_lfq_knn, n_top_loadings = 5)
save_path <- paste0(plot_path,"preprocessing/pca_lfq_knn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Boxplot lfq after knn
res <- create_lfq_boxplot(se_lfq_knn_assay, col_data, lfq_type = "intensity",
                          y_lab = "Protein Top-N intensities")
p <- res$plot_samples
save_path <- paste0(plot_path,"preprocessing/boxplot_lfq_samples_knn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- res$plot_groups
save_path <- paste0(plot_path,"preprocessing/boxplot_lfq_groups_knn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

print(res$lfq_summary)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7.579  15.253  16.059  16.324  17.179  27.863

# Intensities density plot
se_lfq_knn_long <- res$se_long
p <- ggplot(se_lfq_knn_long, aes(x = intensity)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  labs(title = "Distribution of protein Top-N intensities",
       x = "Protein Top-N intensities", 
       y = "Density") +
  theme_minimal(base_size = 12)

save_path <- paste0(plot_path,"preprocessing/density_plot_lfq_knn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Intensities density plot per sample
plotDensities(se_lfq_knn_assay, legend = FALSE)

# Save SE imputed by knn
save_path <- paste0(data_path,"se_lfq_knn.RData")
save(se_lfq_knn, file = save_path)
