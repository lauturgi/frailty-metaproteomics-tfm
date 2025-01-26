#renv::install("bioc::ComplexHeatmap")
#renv::install("bioc::limma")
## on bash: sudo apt-get install libnetcdf-dev
#renv::install("bioc::MSnbase")
#renv::install("bioc::SummarizedExperiment")
#renv::install("bioc::cmapR")
#renv::install("bioc::ConsensusClusterPlus")
#renv::install("Nesvilab/FragPipeAnalystR")
#renv::install("bioc::DEP")

library(FragPipeAnalystR)
library(ggplot2)
library(limma)
library(DEP)

# Load SummarizedExperiment
work_path <- getwd()
load(paste0(work_path, "/data/se.RData"))

# Number of proteins per sample
plot_feature_numbers(se)

# Missing value heatmap
plot_missval_heatmap(se)

# Plot correlation heatmap
png("heatmap_plot_ft.png", width = 2000, height = 1600, res = 300)
plot_correlation_heatmap(se, font_size = 2, indicate = c("frailty"))
dev.off()

# Plot PCA
## x: indicate PC to represent in x axis
## y: indicate PC to represent in y axis
## indicate: select metadata used to differentiate samples. Only condition, 
## since there are many samples to be represented with different shapes
## (limit of 6)
plot_pca(se, x = 1, y = 2, indicate = c("frailty"), point_size = 8, 
         interactive = TRUE)

# ¡¡It is adviced to first remove proteins with too many missing values!!
## https://rdrr.io/bioc/DEP/src/R/functions.R. 2 options:
## - filter_missval: Filters to retain proteins with a limited number of missing
## values in at least one condition. Parameters:
##    - se: SummarizedExperiment object
##    - thr: Integer threshold defining the maximum missing values per condition
## - filter_proteins: Multiple filtering options based on missing values
##    - se: SummarizedExperiment object
##    - type: The type of filtering to apply:
##       - "complete": Retains proteins with no missing values across all samples
##       - "condition": Retains proteins with a limited number of missing values
##                      in at least one condition (controlled by thr)
##       - "fraction": Retains proteins with a minimum fraction of valid values
##                     across all samples (controlled by min)
##    - thr: Threshold for the number of missing values per condition
##    - min: Minimum fraction of valid values

#se_filtered <- filter_proteins(se, type = "fraction", min = 0.25)

#table(se$condition)

plot_missval_heatmap(se_filtered)

## --> plot_feature(se, "A1A0H0")

# Impute missing values. Missing values are replaced with random values
# generated from a shifted and scaled normal distribution based on the existing
# data
se_imputed <- manual_impute(se_filtered)


plot_pca(se_imputed, x = 1, y = 2, indicate = c("frailty"), point_size = 8, 
         interactive = TRUE)

plot_correlation_heatmap(se_imputed, font_size = 2, indicate = c("frailty"))


# Limma test FragpipeAnalystR
## test_limma function performs differential expression analysis. Parameters:
## - se: SummarizedExperiment object
## - type: Type of contrasts to test: "control", "all" or "manual". Default: "control"
## - control: Control condition for comparison (if type == "control")
## - test: Specific contrasts to test if type == "manual". Vector with contrasts
##         in the form condition1_vs_condition2
## - design_formula: Formula used to create the design matrix for the linear
##                   model. Default: ~ 0 + condition, meaning no intercept and a
##                   condition-based model
## - paired: Logical indicating if the samples are paired. If TRUE, it includes
##           a replicate variable in the model formula for paired analysis
de <- test_limma(se_filtered, type = "manual", test = "FT_vs_NFT")
head(rowData(de)[, c("FT_vs_NFT_diff", "FT_vs_NFT_p.val", "FT_vs_NFT_p.adj")])


# Get significant proteins
## add_rejections function adds columns to SummarizedExperiment to indicate
## which proteins are significantly differentially expressed based on adjusted
## p-values and log2 fold changes. Parameters:
## - alpha (default 0.05): The adjusted p-value threshold for significance
## - lfc (default 1): The absolute log2 fold change threshold for significance

de_sig <- add_rejections(de)
tail(colnames(rowData(de_sig)))
head(rowData(de_sig)[, c("FT_vs_NFT_diff", "FT_vs_NFT_p.val", "FT_vs_NFT_p.adj", "significant")])
sum(rowData(de_sig)$FT_vs_NFT_significant == TRUE)
plot_volcano(de_sig, "FT_vs_NFT")
