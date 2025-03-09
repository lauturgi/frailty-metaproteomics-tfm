#renv::install("nortest")
library(FragPipeAnalystR)
library(ggplot2)
library(limma)
library(vsn)
library(nortest)
library(dplyr)
library(SummarizedExperiment)
library(reshape2)
library(ComplexHeatmap)
library(UniprotR)

#renv::install("bioc::impute")
#renv::install()

# ==================================
# 1. Making SummarizedExperiment
# ==================================

# Working directory
work_path <- getwd()

# Path to combined_protein.tsv and experiment_annotation.tsv

fragpipe_path <- "/mnt/d/proteomica/fragilidad/datos/ProteinIdent/MSfraggerSPbacteria/DatosProtIdentCombinados/"

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

# Take the first identifier per row and make unique names.
# If there is no name, the ID will be taken.
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

# Select rowData
row_data <- prot_uniq[, -lfq_col]
rownames(row_data) <- prot_uniq$ID

# Make SummarizedExperiment
se <- SummarizedExperiment(
  assays = as.matrix(prot_maxlfq),
  colData = exp_anno,
  rowData = row_data,
  metadata = list("log2transform"=F, "lfq_type"="MaxLFQ",
                  "level"="protein")
)
# Make summarizedExperiment
# By default, make_se_from_files applies a log2 transformation
#se <- make_se_from_files(combined_protein, experiment_ann, type = "LFQ",
#                         level = "protein", lfq_type = "MaxLFQ")

# Check class of se object
class(se)

# Check log2 (TRUE), exp (LFQ), lfq_type (MaxLFQ), level (protein)
metadata(se)

# Check number of rows and columns
dim(se)

# Phenotypic variables
colData(se)
names(colData(se))

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
load(paste0(work_path, "/data/metadatos.RData"))

# Keep only data from FT cohort
metadata_ft$replicate <- substring(metadata_ft$ID.FISABIO, 3, 5)
head(metadata_ft$replicate)
se <- se[, colData(se)$replicate %in% metadata_ft$replicate]
unique(colData(se)$replicate)

# Check number of rows and columns
dim(se)

# ==================================
# 3. Adding phenotypic data
# ==================================

# Create an aligned dataframe with the order of `replicate` in colData(se)
replicate_se <- colData(se)$replicate
replicate_metadata <- metadata_ft$replicate

metadata_ft <- metadata_ft[match(replicate_se, replicate_metadata), ]

# Add phenotypic variables to colData
colData(se)$frailty <- factor(metadata_ft$Fragilidad, levels = c("NFT", "FT"))
colData(se)$sex <- metadata_ft$w23sex
# Group alcohol consumption in 3 levels in 0: Rare | 1: Monthly | 2: Weekly
# 1: A diario o casi a diario --> 2
# 2: 5-6 días por semana --> 2
# 3: 3-4 días por semana --> 2
# 4: 1-2 días por semana --> 2
# 5: 2-3 días en un mes --> 1
# 6: 1 vez al mes --> 1
# 7: Menos de una vez al mes --> 1
# 8: No en los últimos 12 meses, he dejado de tomar alcohol --> 0
# 9: Nunca o solamente unos sorbos para probar a lo largo de toda la vida --> 0
colData(se)$alcohol <- factor(ifelse(metadata_ft$w23alcohol %in% 
                                       c("1", "2", "3", "4"), "2",
                              ifelse(metadata_ft$w23alcohol %in% 
                                       c("5", "6", "7"), "1",
                                     ifelse(metadata_ft$w23alcohol %in% 
                                              c("8", "9"), "0", NA)
                                     )
                              ), levels = c("0", "1", "2"))

colData(se)$tobacco <- metadata_ft$w23tobacco
colData(se)$diabetes <- metadata_ft$w23diabetes
colData(se)$chf <- metadata_ft$w23chf
colData(se)$depression <- metadata_ft$w23depression
colData(se)$osteoarthritis <- metadata_ft$w23osteoarthritis
colData(se)$sarcopenia <- metadata_ft$w23sarcopenia

colData(se)$bmi <- metadata_ft$w23bmi
colData(se)$energy <- metadata_ft$w23energy

# ILEF calculated from SPPB score
colData(se)$ilef <- as.factor(ifelse(metadata_ft$w23sppbscore > 9, "0",
                                     ifelse(metadata_ft$w23sppbscore <= 9, "1", 
                                            NA)
                                     )
                              )


# Check
head(colData(se))

# Plot alcohol 3 levels
create_bar_plot(colData(se), "alcohol", title = paste("Bar plot for", 
                                                      "alcohol"), 
                paste0(work_path,"/plots/bar_plot_", "alcohol_3",".png"))
# Plot ILEF 
create_bar_plot(colData(se), "ilef", title = paste("Bar plot for", "ILEF"), 
                paste0(work_path,"/plots/bar_plot_", "ILEF",".png"))

# Save SummarizedExperiment
save_path <- paste0(work_path,"/data/se.RData")
save(se, file = save_path)

# ==================================
# 4. Filtering
# ==================================
# 4.1. By Missing
# ==================================

# Get protein intensities matrix from SummarizedExperiment
se_assay <- assay(se)

# Number of proteins
nrow(se_assay)
# 6854

# Number of NA values in assay
table(is.na(se_assay))
# FALSE    TRUE 
# 61378 1329984

# Number of proteins with all NA
proteins_to_keep <- apply(se_assay, 1, function(x) !all(is.na(x)))
sum(!proteins_to_keep)
# 3952 

# Filter proteins with all NA
se_assay <- se_assay[proteins_to_keep, ]
se <- se[proteins_to_keep,]

# Number of proteins after filtering those with all NA
nrow(se_assay)
# 2902

# Number of NA values in assay after filtering those with all NA
table(is.na(se_assay))
# FALSE   TRUE 
# 61378 527728 

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
ggplot(num_prot_sample, aes(x = Proteins, fill = frailty)) +
  geom_density(alpha = 0.4, linewidth = 0.2) +
  labs(title = "Number of proteins per sample", x = "Number of proteins", 
       y = "Density") +
  scale_fill_manual(values = c("FT" = "blue", "NFT" = "red")) +
  theme_minimal()


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
ggplot(num_na_sample, aes(x = NAs, fill = frailty)) +
  geom_density(alpha = 0.4, linewidth = 0.2) +
  labs(title = "Number of NAs per sample", x = "Number of NAs", y = "Density") +
  scale_fill_manual(values = c("FT" = "blue", "NFT" = "red")) +
  theme_minimal()


# Make binary assay (1 = valid value, 0 = missing value)
missval <- ifelse(is.na(se_assay), 0, 1)
colnames(missval) <- substr(colnames(missval), 1, 7)

# Convert matrix to long format dataframe
df <- melt(missval)

# Rename columns
colnames(df) <- c("protein", "sample", "intensity")

# Rename samples
df$sample <- substr(df$sample, 1, 7)

# Add frailty group
col_data <- as.data.frame(colData(se))
df <- merge(df, col_data[, c("sample", "frailty")], by = "sample")

# Plot heatmap
ggplot(df, aes(sample, protein, fill = intensity)) + 
  geom_tile() +
  facet_grid(. ~ frailty, scales = "free_x", space = "free_x") +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 2, angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_gradient(low = "white", high = "black")

# Complex Heatmap
# Clustering rows and columns by binary distance
# Column split by frailty groups
Heatmap(missval, name = "NA", col = c("white", "black"), show_row_dend = FALSE, 
        column_names_gp = gpar(fontsize = 4), show_row_names = FALSE, 
        # cluster_columns = FALSE,
        # cluster_rows = FALSE,
        column_split = colData(se)$frailty, 
        clustering_distance_rows = "binary", 
        clustering_distance_columns = "binary", 
        use_raster = FALSE)

# Annotation of columns by frailty group
column_ann <- HeatmapAnnotation(frailty = colData(se)$frailty,
                                col = list(frailty = c("FT" = "red", 
                                                       "NFT" = "blue")))
Heatmap(missval, name = "NA", col = c("white", "black"), show_row_dend = FALSE, 
        column_names_gp = gpar(fontsize = 4), show_row_names = FALSE, 
        # cluster_columns = FALSE,
        # cluster_rows = FALSE,
        clustering_distance_rows = "binary", 
        clustering_distance_columns = "binary",
        top_annotation = column_ann, 
        use_raster = FALSE)

# Filter missing by sample proportion

# Sample proportion for each protein
presence_proportion <- rowSums(!is.na(se_assay)) / ncol(se_assay)
sample_percent <- presence_proportion * 100
sample_percent <- data.frame(SamplePercentage = sample_percent)
create_density_plot(sample_percent, "SamplePercentage")

# Number of proteins present in more than 25%
sum(presence_proportion > 0.25)
# 320

# Number of proteins present in more than 50%
sum(presence_proportion > 0.5)
# 159

# Number of proteins present in more than 75%
sum(presence_proportion > 0.75)
# 72

# Number of proteins present in 100%
sum(presence_proportion == 1)
# 3

# Proteins present in 100% of samples
rownames(se_assay[presence_proportion == 1,])
# P00761 --> trypsin
# P04264 --> Keratin, type II cytoskeletal 1
# P06702 --> Protein S100-A9: regulation of inflammatory processes and
# immune response

# Proteins in more than 75% of samples
rownames(se_assay[presence_proportion >= 0.75,])

# Proteins in less than 50% of samples
rownames(se_assay[presence_proportion < 0.5,])

#########################################################
# ANALIZAR PROTEINAS PRESENTES EN < 50% DE LAS MUESTRAS #
#########################################################

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
# 185

# Filter SummarizedExperiment
se_filtered <- se[proteins_to_keep, ]

# Get protein intensities matrix from filtered SummarizedExperiment
se_filtered_assay <- assay(se_filtered)

# Number of proteins after filtering those present in <50% of samples in both
# conditions
nrow(se_filtered_assay)
# 185

# Number of proteins detected per sample
num_prot_sample <- apply(se_filtered_assay, 2, function(x) sum(!is.na(x)))
num_prot_sample <- as.data.frame(num_prot_sample)
colnames(num_prot_sample) <- "Proteins"
num_prot_sample$Sample <- substr(rownames(num_prot_sample), 1, 7)
rownames(num_prot_sample) <- NULL
# Add frailty group
num_prot_sample <- merge(num_prot_sample, col_data[, c("sample", "frailty")], 
                         by.x = "Sample", by.y = "sample")

# Summary total proteins detected per sample
summary(num_prot_sample)
#Sample             Proteins     frailty  
#Length:203         Min.   : 14.0   NFT:138  
#Class :character   1st Qu.:103.0   FT : 65  
#Mode  :character   Median :140.0            
#                   Mean   :127.6            
#                   3rd Qu.:161.5            
#                   Max.   :182.0  

# Density plot number of proteins per sample in FT vs NFT
ggplot(num_prot_sample, aes(x = Proteins, fill = frailty)) +
  geom_density(alpha = 0.4, linewidth = 0.2) +
  labs(title = "Number of proteins per sample", x = "Number of proteins", 
       y = "Density") +
  scale_fill_manual(values = c("FT" = "blue", "NFT" = "red")) +
  theme_minimal()


# Number of NA values in assay
table(is.na(se_filtered_assay))
# FALSE  TRUE 
# 25894 11661 

# Number of missing per sample
num_na_sample <- apply(se_filtered_assay, 2, function(x) sum(is.na(x)))
num_na_sample <- as.data.frame(num_na_sample)
colnames(num_na_sample) <- "NAs"
num_na_sample$Sample <- substr(rownames(num_na_sample), 1, 7)
rownames(num_na_sample) <- NULL
# Add frailty group
num_na_sample <- merge(num_na_sample, col_data[, c("sample", "frailty")], 
                       by.x = "Sample", by.y = "sample")

# Density plot number of missing per sample in FT vs NFT
ggplot(num_na_sample, aes(x = NAs, fill = frailty)) +
  geom_density(alpha = 0.4, linewidth = 0.2) +
  labs(title = "Number of NAs per sample", x = "Number of NAs", y = "Density") +
  scale_fill_manual(values = c("FT" = "blue", "NFT" = "red")) +
  theme_minimal()

# Make binary assay (1 = valid value, 0 = missing value)
missval <- ifelse(is.na(se_filtered_assay), 0, 1)
colnames(missval) <- substr(colnames(missval), 1, 7)

# Convert matrix to long format dataframe
df <- melt(missval)

# Rename columns
colnames(df) <- c("protein", "sample", "intensity")

# Rename samples
df$sample <- substr(df$sample, 1, 7)

# Add frailty group
col_data <- as.data.frame(colData(se))
df <- merge(df, col_data[, c("sample", "frailty")], by = "sample")

# Plot heatmap
ggplot(df, aes(sample, protein, fill = intensity)) + 
  geom_tile() +
  facet_grid(. ~ frailty, scales = "free_x", space = "free_x") +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 2, angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_gradient(low = "white", high = "black")

# Complex Heatmap
# Clustering rows and columns by binary distance
# Column split by frailty groups
Heatmap(missval, name = "NA", col = c("white", "black"), show_row_dend = FALSE, 
        column_names_gp = gpar(fontsize = 4), show_row_names = FALSE, 
        # cluster_columns = FALSE,
        # cluster_rows = FALSE,
        column_split = colData(se_filtered)$frailty, 
        clustering_distance_rows = "binary", 
        clustering_distance_columns = "binary")

# Annotation of columns by frailty group
column_ann <- HeatmapAnnotation(frailty = colData(se_filtered)$frailty,
                                col = list(frailty = c("FT" = "red", 
                                                       "NFT" = "blue")))
Heatmap(missval, name = "NA", col = c("white", "black"), show_row_dend = FALSE, 
        column_names_gp = gpar(fontsize = 4), show_row_names = FALSE, 
        # cluster_columns = FALSE,
        # cluster_rows = FALSE,
        clustering_distance_rows = "binary", 
        clustering_distance_columns = "binary",
        top_annotation = column_ann)

# Save SummarizedExperiment filtered
save_path <- paste0(work_path,"/data/se_filtered.RData")
save(se_filtered, file = save_path)

# ==================================
# 4.2. By organism
# ==================================
# Get lineage from UniProt
lineage_df <- GetProteinAnnontate(rownames(se_filtered_assay),
                                  columns = c("lineage"))
lineage_df <- as.data.frame(lineage_df)
lineage_df$protein_id <- rownames(se_filtered_assay)
colnames(lineage_df)[colnames(lineage_df) == "lineage_df"] <- "lineage"
# Get superkingdom from lineage
lineage_df$superkingdom <- sapply(strsplit(lineage_df$lineage, ","), function(x) x[2])

# Get protein_id from bacteria proteins
proteins_to_keep <- lineage_df[grep("Bacteria", lineage_df$superkingdom) , 
                               "protein_id"]
# Filter SummarizedExperiment
se_filt_bact <- se_filtered[proteins_to_keep, ]

# Number of proteins after filtering no bacteria proteins
nrow(se_filt_bact)
# 163

se_filt_bact_assay <- assay(se_filt_bact)

# ==================================
# 5. Normalization
# ==================================
# Boxplot
se_filt_bact_df <- as.data.frame(se_filt_bact_assay)
se_filt_bact_df$protein_id <- rownames(se_filt_bact_df)
rownames(se_filt_bact_df) <- NULL
se_filt_bact_long <- melt(se_filt_bact_df)
colnames(se_filt_bact_long) <- c("protein_id", "Samples", "MaxLFQ")
head(se_filt_bact_long)
se_filt_bact_long$Samples <- substr(se_filt_bact_long$Samples, 1, 7)
se_filt_bact_long <- as.data.frame(merge(se_filt_bact_long, 
                                        colData(se_filt_bact)[,c("sample", 
                                                                "frailty")], 
                                        by.x = "Samples", by.y = "sample"))
# Summary MaxLFQ to set limit to y axis
summary(se_filt_bact_long$MaxLFQ)

# Plot MaxLFQ vs samples
ggplot(se_filt_bact_long, aes(x = Samples, y = MaxLFQ, fill = frailty)) + 
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    legend.position = c(0.95, 0.95),
    legend.background = element_rect(fill = alpha("white", 0.5))) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple")) +
  ylim(0, 250000) +
  facet_wrap(~frailty, scales = "free_x", ncol = 1)

# Plot MaxLFQ vs frailty group
ggplot(se_filt_bact_long, aes(x = frailty, y = MaxLFQ, fill = frailty)) + 
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    legend.position = c(0.95, 0.95),
    legend.background = element_rect(fill = alpha("white", 0.5))) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple")) +
  ylim(0, 250000)

# Intensities density plot
ggplot(se_filt_bact_long, aes(x = MaxLFQ)) +
  geom_density(fill = "blue", alpha = 0.4) +
  theme_minimal()# + 
  #xlim(0, 500000)

# Quantile 95%
q99 <- quantile(se_filt_bact_long$MaxLFQ, 0.99, na.rm = TRUE)

# Get proteins with intensities higher than q95
#se_filt_long_q99 <- se_filtered_long[se_filtered_long$MaxLFQ > q99, ]
#se_filt_long_q99 <- se_filt_long_q99[!is.na(se_filt_long_q99), ]
#unique(se_filt_long_q99$protein_id)

# Intensities density plot per sample
plotDensities(se_filtered_assay, legend = FALSE)

for (i in 1:nrow(se_filtered_assay)){
  plotDensities(se_filtered_assay[, i], legend = TRUE)
}

# Standard deviation (sd) and mean are calculated row-wise from the expression
# matrix. The scatterplot of these versus each other to verify whether there is
# a dependence # of the sd on the mean. The red line running median estimator
# (window-width 10%). If there is no variance-mean dependence, the line should
# be aprox. horizontal.
meanSdPlot(se_filtered_assay)

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
#    Mode   FALSE    TRUE 
#logical     172      31 

#norm_results <- data.frame(
#  sample = substring(colnames(se_filtered_assay), 1, 7),
#  p_value = p_values,
#  norm = ifelse(norm_results, "1", "0")
#)

#norm_results$frailty <- colData(se_filtered)$frailty

# Percentage normal and not normal
#sum(norm_results$norm == "0")/nrow(norm_results)
#sum(norm_results$norm == "1")/nrow(norm_results)

# Percentage of not normal and normal in frailty patients
#sum(norm_results$norm == "0" & norm_results$frailty == "FT")/
#  sum(norm_results$frailty == "FT")
#sum(norm_results$norm == "1" & norm_results$frailty == "FT")/
#  sum(norm_results$frailty == "FT")

# Percentage of not normal and normal in non-frailty patients
#sum(norm_results$norm == "0" & norm_results$frailty == "NFT")/
#  sum(norm_results$frailty == "NFT")
#sum(norm_results$norm == "1" & norm_results$frailty == "NFT")/
#  sum(norm_results$frailty == "NFT")

# Normalization
# se_norm <- VSN_normalization(se_filtered)
fit = vsn2(se_filtered)
se_norm = predict(fit, se_filtered)
se_norm_assay <- assay(se_norm)

# Intensities density plot per sample
plotDensities(se_norm_assay, legend = FALSE)

for (i in 1:nrow(se_norm_assay)){
  plotDensities(se_norm_assay[, i], legend = TRUE)
}

meanSdPlot(se_norm_assay)

# Normal contrast
p_values <- apply(se_norm_assay, 2,
                  function(column) shapiro.test(column)$p.value)
norm_results <- p_values >= 0.05 # TRUE if normal, FALSE if not
summary(norm_results)
#    Mode   FALSE    TRUE 
#logical     172      31 

# ==================================
# 6. Imputation
# ==================================
# 6.1. No imputation
# ==================================
se_no_imp <- se_filtered

# Save SummarizedExperiment no imputated
save_path <- paste0(work_path,"/data/se_no_imp.RData")
save(se_no_imp, file = save_path)

# ==================================
# 6.2. Impute missing values
# ==================================
# 6.2.1. Perseus
# ==================================
# Missing values are replaced with random values generated from a shifted and
# scaled normal distribution based on the existing data
se_perseus <- manual_impute(se_no_imp)

# Plot PCA
create_pca_plot(se_perseus, n_top_loadings = 5)

# ==================================
# 6.2. KNN
# ==================================
se_no_imp_prepro <- se_no_imp[, colMeans(is.na(assay(se_no_imp))) <= 0.8]
se_knn <- impute(se_no_imp_prepro, fun = "knn")
plot_pca(se_knn, x = 1, y = 2, indicate = c("frailty"), point_size = 8, 
         interactive = TRUE)

# Save SummarizedExperiment imputed
save_path <- paste0(work_path,"/data/se_perseus.RData")
save(se_perseus, file = save_path)
