#renv::install("nortest")
library(FragPipeAnalystR)
library(ggplot2)
library(limma)
library(vsn)
library(nortest)

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

sum(is.na(assay(se))) # 1793735


# Save SummarizedExperiment
save_path <- paste0(work_path,"/data/se.RData")
save(se, file = save_path)

# Check class of se object
class(se)

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

# Keep only data from FT cohort --> ¡¡Existe función en FragpipeAnalystR para eliminar muestras!!
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
                                     ifelse(metadata_ft$w23sppbscore <= 9, "1", NA)
                                     )
                              )


# Check
head(colData(se))

# Plot alcohol 3 levels
create_bar_plot(colData(se), "alcohol", title = paste("Bar plot for", "alcohol"), paste0(work_path,"/plots/bar_plot_", "alcohol_3",".png"))
# Plot ILEF 
create_bar_plot(colData(se), "ilef", title = paste("Bar plot for", "ILEF"), paste0(work_path,"/plots/bar_plot_", "ILEF",".png"))


# ==================================
# 4. Filtering missing
# ==================================

# Get protein intensities matrix from SummarizedExperiment
data_matrix <- assay(se)

# Percentage of proteins detected per sample
presence_proportion <- rowSums(!is.na(data_matrix)) / ncol(data_matrix)
sample_percent <- presence_proportion * 100
sample_percent <- data.frame(SamplePercentage = sample_percent)
create_density_plot(sample_percent, "SamplePercentage")

# Number of proteins present in more than 25%
sum(presence_proportion > 0.25)

# Number of proteins present in more than 50%
sum(presence_proportion > 0.5)

# Number of proteins present in more than 75%
sum(presence_proportion > 0.75)

# Retain proteins with a minimum fraction of valid values per condition
min_fraction <- 0.5
fraction_ft <- apply(data_matrix[,colData(se)$frailty == "FT"], 1,
                     function(x) sum(!is.na(x)) / length(x))
fraction_nft <- apply(data_matrix[, colData(se)$frailty == "NFT"], 1,
                      function(x) sum(!is.na(x)) / length(x))

# Keep proteins with at least 50% of valid values in at least one condition
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

# ==================================
# 5. Normalization
# ==================================

# Intensities density plot
summary(data_matrix)
plotDensities(data_matrix, legend = FALSE)
meanSdPlot(data_matrix)
#plotMA(data_matrix)

# Q-Q sample quantiles vs theoretical quantiles for a normal distribution

for (i in 1:ncol(data_matrix)) {
  qqnorm(data_matrix[, i], main = paste("QQ Plot -", colnames(data_matrix)[i])) # Usar nombre de la columna
  qqline(data_matrix[, i], col = "red")
}

# Normal contrast
p_values <- apply(data_matrix, 2, function(column) pearson.test(column)$p.value)
normality_results <- p_values > 0.05 # TRUE if normal, FALSE if not
summary_table <- data.frame(
  sample = colnames(data_matrix),
  p_value = p_values,
  norm = ifelse(normality_results, "Yes", "No")
)
summary(normality_results)

#se_norm <- VSN_normalization(se_filtered)
#data_matrix_norm <- assay(se_norm)
#plotDensities(data_matrix_norm, legend = FALSE)
#meanSdPlot(data_matrix_norm)


se_no_imp <- se_filtered

# Save SummarizedExperiment no imputated
save_path <- paste0(work_path,"/data/se_no_imp.RData")
save(se_no_imp, file = save_path)

# ==================================
# 6. Imputation
# ==================================

# Impute missing values. Missing values are replaced with random values
# generated from a shifted and scaled normal distribution based on the existing
# data
se_perseus <- manual_impute(se_no_imp)

# Save SummarizedExperiment imputed
save_path <- paste0(work_path,"/data/se_perseus.RData")
save(se_perseus, file = save_path)
