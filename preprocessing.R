#renv::install("nortest")
library(FragPipeAnalystR)
library(ggplot2)
library(limma)
library(vsn)
library(nortest)
library(dplyr)
library(SummarizedExperiment)
library(reshape2)
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

# Make summarizedExperiment
# By default, make_se_from_files applies a log2 transformation
se <- make_se_from_files(combined_protein, experiment_ann, type = "LFQ",
                         level = "protein", lfq_type = "MaxLFQ")

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
colData(se)$alcohol <- factor(ifelse(metadata_ft$w23alcohol %in% c("1", "2", "3", "4"), "2",
                              ifelse(metadata_ft$w23alcohol %in% c("5", "6", "7"), "1",
                                     ifelse(metadata_ft$w23alcohol %in% c("8", "9"), "0", NA)
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
# 4. Missing
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

# Filter proteins with all NA
se_assay <- se_assay[apply(se_assay, 1, function(x) !all(is.na(x))), ]

# Number of proteins after filtering those with all NA
nrow(se_assay)
# 2902

# Number of NA values in assay after filtering those with all NA
table(is.na(se_assay))
# FALSE   TRUE 
# 61378 527728 

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

# ==================================
# 4.1. Filtering missing
# ==================================
#############
# REVISAR!! #
#############
# Percentage of proteins detected per sample
presence_proportion <- rowSums(!is.na(se_assay)) / ncol(se_assay)
sample_percent <- presence_proportion * 100
sample_percent <- data.frame(SamplePercentage = sample_percent)
create_density_plot(sample_percent, "SamplePercentage")

# Number of proteins present in more than 25%
sum(presence_proportion > 0.25)

# Number of proteins present in more than 50%
sum(presence_proportion > 0.5)

# Number of proteins present in more than 75%
sum(presence_proportion > 0.75)

# Keep proteins with a minimum fraction of valid values in at least one condition
min_fraction <- 0.5
fraction_ft <- apply(data_matrix[,colData(se)$frailty == "FT"], 1,
                     function(x) sum(!is.na(x)) / length(x))
fraction_nft <- apply(data_matrix[, colData(se)$frailty == "NFT"], 1,
                      function(x) sum(!is.na(x)) / length(x))

proteins_to_keep <- (fraction_ft >= min_fraction) | (fraction_nft >= 
                                                       min_fraction)
sum(proteins_to_keep) # 185

# Filter SummarizedExperiment
se_filtered <- se[proteins_to_keep, ]
data_matrix <- assay(se_filtered)

dim(se_filtered)

sum(is.na(assay(se_filtered)))

# Save SummarizedExperiment filtered
save_path <- paste0(work_path,"/data/se_filtered.RData")
save(se_filtered, file = save_path)

table(is.na(assay(se_filtered)))

# ==================================
# 5. Normalization
# ==================================

# Intensities density plot
summary(data_matrix)
plotDensities(data_matrix, legend = FALSE)
meanSdPlot(data_matrix)
#plotMA(data_matrix)

#se_norm <- VSN_normalization(se_filtered)
#data_matrix_norm <- assay(se_norm)
#plotDensities(data_matrix_norm, legend = FALSE)
#meanSdPlot(data_matrix_norm)

# Q-Q sample quantiles vs theoretical quantiles for a normal distribution
for (i in 1:ncol(data_matrix)) {
  qqnorm(data_matrix[, i], main = paste("QQ Plot -", colnames(data_matrix)[i]))
  qqline(data_matrix[, i], col = "red")
}

# Normal contrast
p_values <- apply(data_matrix, 2, function(column) shapiro.test(column)$p.value)
norm_results <- p_values >= 0.05 # TRUE if normal, FALSE if not
summary(norm_results)
norm_results <- data.frame(
  sample = substring(colnames(data_matrix), 1, 7),
  p_value = p_values,
  norm = ifelse(norm_results, "1", "0")
)

norm_results$frailty <- colData(se_filtered)$frailty

# Percentage normal and not normal
sum(norm_results$norm == "0")/nrow(norm_results)
sum(norm_results$norm == "1")/nrow(norm_results)

# Percentage of not normal and normal in frailty patients
sum(norm_results$norm == "0" & norm_results$frailty == "FT")/
  sum(norm_results$frailty == "FT")
sum(norm_results$norm == "1" & norm_results$frailty == "FT")/
  sum(norm_results$frailty == "FT")

# Percentage of not normal and normal in non-frailty patients
sum(norm_results$norm == "0" & norm_results$frailty == "NFT")/
  sum(norm_results$frailty == "NFT")
sum(norm_results$norm == "1" & norm_results$frailty == "NFT")/
  sum(norm_results$frailty == "NFT")

se_no_imp <- se_filtered

# Save SummarizedExperiment no imputated
save_path <- paste0(work_path,"/data/se_no_imp.RData")
save(se_no_imp, file = save_path)

# ==================================
# 6. Imputation
# ==================================
# Impute missing values
# ==================================
# 6.1. Perseus
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
