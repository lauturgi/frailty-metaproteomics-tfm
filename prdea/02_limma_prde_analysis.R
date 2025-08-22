# Load libraries
library(SummarizedExperiment)
library(ggplot2)
library(limma)
library(dplyr)
library(UniprotR)

# ==================================
# Limma DE analysis
# ==================================
# MaxLFQ preprocessing:
# ==================================
# - MaxLFQ,
# - log2transformed,
# - NA in at least 1 condition < 0.5,
# - only bacteria and human proteins,
# - no imputation
# ==================================
# Top-N preprocessing:
# ==================================
# - Norm intensities by VSN,
# - NA in at least 1 condition < 0.5,
# - only bacteria and human proteins,
# - no imputation
# ==================================

# Working directory
work_path <- getwd()

# Load functions
source(paste0(work_path, "/functions/create_volcano_plot.R"))

# Data, plot and output paths
data_path <- paste0(work_path, "/prdea/data/")
plot_path <- paste0(work_path, "/prdea/plots/limma_analysis/")
output_path <- paste0(work_path, "/prdea/output/sig_prot_no_imp_")

# Load MaxLFQ SE
load(paste0(data_path, "se_maxlfq_no_imp.RData"))

# Load LFQ SE
load(paste0(data_path, "se_lfq_no_imp.RData"))

# ==================================
# Frailty
# ==================================
# A) MaxLFQ intensities
# ==================================
# Build model matrix considering only frailty
design = model.matrix(~ colData(se_maxlfq_no_imp)$frailty)
colnames(design) = c("constant", "frailty")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_maxlfq_no_imp), design)

# Estimated coefficients for each protein
head(coef(fit))

# Moderate t-test of the differential expression by means of Bayes' empirical
# moderation of the standard errors to a global value
fit1 = eBayes(fit)

# p-values adjusted using the Benjamini-Hochberg method
diff_1_max <- topTable(fit1, coef=2, adjust="BH", sort.by="P",
                       number=nrow(assay(se_maxlfq_no_imp)))

# Volcano plot
create_volcano_plot(diff_1_max, save_path = paste0(plot_path, "volcano_plot_", 
                                                   "maxlfq_frail.png"),
                    title = "Differentially expressed proteins - Frailty")

diff_1_max <- diff_1_max[diff_1_max$adj.P.Val < 0.05, ]
nrow(diff_1_max)
# 10
rownames(diff_1_max)
# "P06702" "A6KYJ0" "A6L792" "A9KRZ4" "P05109" "Q8A9M2" "P94360" "A6KYJ4" "C4KZP0" "A9KJL3"

# Get annotation from Uniprot: gen, organism, protein, function, keywords
# and go terms
prot_ann <- GetProteinAnnontate(rownames(diff_1_max),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_1_max_ann <- merge(diff_1_max, prot_ann, by = "row.names")
row.names(diff_1_max_ann) <- diff_1_max_ann$Row.names
diff_1_max_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_1_max_ann <- diff_1_max_ann[order(diff_1_max_ann$adj.P.Val,
                                       -abs(diff_1_max_ann$logFC)), ]

# Save sig_prot_no_imp_maxlfq_ft.csv
write.csv(diff_1_max_ann, file = paste0(output_path, "maxlfq_ft.csv"),
          row.names = TRUE)

sel0 = which("A9KJL3"==rownames(assay(se_maxlfq_no_imp)))
df0 = data.frame(frailty=colData(se_maxlfq_no_imp)[,c("frailty")],
                 expression=assay(se_maxlfq_no_imp)[sel0,])
ggplot(df0,aes(x=frailty,y=expression)) + geom_boxplot()

# ==================================
# B) Top-N intensities
# ==================================
# Build model matrix considering only frailty
design = model.matrix(~ colData(se_lfq_no_imp)$frailty)
colnames(design) = c("constant", "frailty")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_lfq_no_imp), design)

# Estimated coefficients for each protein
head(coef(fit))

# Moderate t-test of the differential expression by means of Bayes' empirical
# moderation of the standard errors to a global value
fit1 = eBayes(fit)

# p-values adjusted using the Benjamini-Hochberg method
diff_1_lfq <- topTable(fit1, coef=2, adjust="BH", sort.by="P",
                       number=nrow(assay(se_lfq_no_imp)))

# Volcano plot
create_volcano_plot(diff_1_lfq, save_path = paste0(plot_path, "volcano_plot_", 
                                                   "lfq_frail.png"),
                    title = "Differentially expressed proteins - Frailty")

diff_1_lfq <- diff_1_lfq[diff_1_lfq$adj.P.Val < 0.05, ]
nrow(diff_1_lfq)
# 24
rownames(diff_1_lfq)
# "P06702" "A9KJL3" "P05109" "P33656" "Q8A015" "A9KMF6" "B2UYT8" "Q9Z9L6" "A6L4M1" "A6L0U5" "A6LPS3" "Q05650" "A9KNK6" "A9KNC4" "A6L7J5" "Q5LHW2" "Q5L8B5" "Q8A477" "Q5WLM8" "A9KRZ1" "A6L792" "C4Z2R3" "A9KJI4" "A6L2R5"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_1_lfq),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_1_lfq_ann <- merge(diff_1_lfq, prot_ann, by = "row.names")
row.names(diff_1_lfq_ann) <- diff_1_lfq_ann$Row.names
diff_1_lfq_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_1_lfq_ann <- diff_1_lfq_ann[order(diff_1_lfq_ann$adj.P.Val,
                                       -abs(diff_1_lfq_ann$logFC)), ]

# Save sig_prot_no_imp_lfq_ft.csv
write.csv(diff_1_lfq_ann, file = paste0(output_path, "lfq_ft.csv"),
          row.names = TRUE)

sel0 = which("A9KJL3"==rownames(assay(se_lfq_no_imp)))
df0 = data.frame(frailty=colData(se_lfq_no_imp)[,c("frailty")],
                 expression=assay(se_lfq_no_imp)[sel0,])
ggplot(df0,aes(x=frailty,y=expression)) + geom_boxplot()

# ==================================
# Two predictor variables:
# ==================================
#   Fraily and sex
# ==================================
#   A) MaxLFQ intensities
# ==================================
# Build model matrix considering frailty, sex and interaction
design = model.matrix(~ colData(se_maxlfq_no_imp)$frailty*
                        colData(se_maxlfq_no_imp)$sex)
colnames(design) = c("constant", "frailty", "sex", "frailty:sex")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_maxlfq_no_imp), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_2_max <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_maxlfq_no_imp)))
diff_3_max <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_maxlfq_no_imp)))
diff_4_max <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_maxlfq_no_imp)))

create_volcano_plot(diff_2_max, save_path = paste0(plot_path, "volcano_plot_",
                                               "maxlfq_frail_sex__frail.png"),
                     title = paste("Differentially expressed proteins",
                                   "- Frailty"))
create_volcano_plot(diff_3_max, save_path = paste0(plot_path, "volcano_plot_",
                                               "maxlfq_frail_sex__sex.png"),
                    title = "Differentially expressed proteins - Sex")
create_volcano_plot(diff_4_max, save_path = paste0(plot_path, "volcano_plot_",
                                               "maxlfq_frail_sex__frail_sex.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Sex"))

diff_2_max <- diff_2_max[diff_2_max$adj.P.Val < 0.05, ]
nrow(diff_2_max)
# 30
rownames(diff_2_max)
# "A6KYJ0" "C4KZP0" "A6L792" "A9KJJ5" "Q5L923" "A6KYJ7" "P05109" "A6KYJ4" "P95544" "A9KRZ4" "Q5L9B6" "C4ZBD3" "Q5L8C5" "Q8A9M2" "C4ZD46" "A8YXK9" "B9E9L7" "Q8A1A2" "A6KYH8" "A9KKU0" "P55990" "C4ZI85" "Q59199" "C4ZBG1" "P06702" "A6KYJ6" "A6KYH1" "C4ZF71" "Q18CF4" "C4Z0Q6"

diff_3_max <- diff_3_max[diff_3_max$adj.P.Val < 0.05, ]
nrow(diff_3_max)
# 0
rownames(diff_3_max)

diff_4_max <- diff_4_max[diff_4_max$adj.P.Val < 0.05, ]
nrow(diff_4_max)
# 1
rownames(diff_4_max)
# "C4Z0Q6"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_4_max),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_4_max_ann <- merge(diff_4_max, prot_ann, by = "row.names")
row.names(diff_4_max_ann) <- diff_4_max_ann$Row.names
diff_4_max_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_4_max_ann <- diff_4_max_ann[order(diff_4_max_ann$adj.P.Val,
                                 -abs(diff_4_max_ann$logFC)), ]

# Save sig_prot_no_imp_maxlfq_ft_sex.csv
write.csv(diff_4_max_ann, file = paste0(output_path, "maxlfq_ft_sex.csv"),
          row.names = TRUE)

sel0 = which("C4Z0Q6"==rownames(assay(se_maxlfq_no_imp)))
df0 = data.frame(frailty=colData(se_maxlfq_no_imp)[,c("frailty")], 
                 sex=colData(se_maxlfq_no_imp)[,c("sex")], 
                 expression=assay(se_maxlfq_no_imp)[sel0,])
df0$frailty_sex <- paste(df0$frailty, df0$sex, sep = "_")
ggplot(df0,aes(x=frailty_sex,y=expression)) + geom_boxplot()

# ==================================
#   B) Top-N intensities
# ==================================
# Build model matrix considering frailty, sex and interaction
design = model.matrix(~ colData(se_lfq_no_imp)$frailty*
                        colData(se_lfq_no_imp)$sex)
colnames(design) = c("constant", "frailty", "sex", "frailty:sex")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_lfq_no_imp), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_2_lfq <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                       number=nrow(assay(se_lfq_no_imp)))
diff_3_lfq <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                       number=nrow(assay(se_lfq_no_imp)))
diff_4_lfq <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                       number=nrow(assay(se_lfq_no_imp)))

create_volcano_plot(diff_2_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                   "lfq_frail_sex__frail.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty"))
create_volcano_plot(diff_3_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                   "lfq_frail_sex__sex.png"),
                    title = "Differentially expressed proteins - Sex")
create_volcano_plot(diff_4_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                   "lfq_frail_sex__frail_sex.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Sex"))

diff_2_lfq <- diff_2_lfq[diff_2_lfq$adj.P.Val < 0.05, ]
nrow(diff_2_lfq)
# 3
rownames(diff_2_lfq)
# "A9KNC4" "Q05650" "A6LPS3"

diff_3_lfq <- diff_3_lfq[diff_3_lfq$adj.P.Val < 0.05, ]
nrow(diff_3_lfq)
# 0
rownames(diff_3_lfq)

diff_4_lfq <- diff_4_lfq[diff_4_lfq$adj.P.Val < 0.05, ]
nrow(diff_4_lfq)
# 0
rownames(diff_4_lfq)

# ==================================
#   Frailty and education
# ==================================
#   A) MaxLFQ intensities
# ==================================
# Build model matrix considering frailty, education and interaction
se_maxlfq_no_imp_educa <- se_maxlfq_no_imp[, !is.na(colData(se_maxlfq_no_imp)$education)]
dim(se_maxlfq_no_imp_educa)
# 182 201
design = model.matrix(~ colData(se_maxlfq_no_imp_educa)$frailty*
                        colData(se_maxlfq_no_imp_educa)$education)
colnames(design) = c("constant", "frailty", "secondary", "university",
                     "frailty:secondary", "frailty:university")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_maxlfq_no_imp_educa), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 5 coefficients per each column of model matrix. Adjusted using BH
diff_5_max <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_maxlfq_no_imp_educa)))
diff_6_max <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_maxlfq_no_imp_educa)))
diff_7_max <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_maxlfq_no_imp_educa)))
diff_8_max <- topTable(fit1,coef=5, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_maxlfq_no_imp_educa)))
diff_9_max <- topTable(fit1,coef=6, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_maxlfq_no_imp_educa)))

create_volcano_plot(diff_5_max, save_path = paste0(plot_path, "volcano_plot_",
                                               "maxlfq_frail_educa__frail.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty"))
create_volcano_plot(diff_6_max, save_path = paste0(plot_path, "volcano_plot_",
                                               "maxlfq_frail_educa__seco.png"),
                    title = paste("Differentially expressed proteins -",
                                  "Secondary"))
create_volcano_plot(diff_7_max, save_path = paste0(plot_path, "volcano_plot_",
                                               "maxlfq_frail_educa__uni.png"),
                    title = paste("Differentially expressed proteins",
                                  "- University"))
create_volcano_plot(diff_8_max, save_path = paste0(plot_path, "volcano_plot_",
                                               "maxlfq_frail_educa__frail_",
                                               "seco.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Secondary"))
create_volcano_plot(diff_9_max, save_path = paste0(plot_path, "volcano_plot_",
                                               "maxlfq_frail_educa__frail_",
                                               "uni.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:University"))

diff_5_max <- diff_5_max[diff_5_max$adj.P.Val < 0.05, ]
nrow(diff_5_max)
# 5
rownames(diff_5_max)
# "P06702" "A6KYJ4" "P05109" "A6L048" "A6L792"

diff_6_max <- diff_6_max[diff_6_max$adj.P.Val < 0.05, ]
nrow(diff_6_max)
# 0
rownames(diff_6_max)

diff_7_max <- diff_7_max[diff_7_max$adj.P.Val < 0.05, ]
nrow(diff_7_max)
# 0
rownames(diff_7_max)

diff_8_max <- diff_8_max[diff_8_max$adj.P.Val < 0.05, ]
nrow(diff_8_max)
# 0
rownames(diff_8_max)

diff_9_max <- diff_9_max[diff_9_max$adj.P.Val < 0.05, ]
nrow(diff_9_max)
# 0
rownames(diff_9_max)

# ==================================
#   B) Top-N intensities
# ==================================
# Build model matrix considering frailty, education and interaction
se_lfq_no_imp_educa <- se_lfq_no_imp[, !is.na(colData(se_lfq_no_imp)$education)]
dim(se_lfq_no_imp_educa)
# 442 201
design = model.matrix(~ colData(se_lfq_no_imp_educa)$frailty*
                        colData(se_lfq_no_imp_educa)$education)
colnames(design) = c("constant", "frailty", "secondary", "university",
                     "frailty:secondary", "frailty:university")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_lfq_no_imp_educa), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 5 coefficients per each column of model matrix. Adjusted using BH
diff_5_lfq <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                       number=nrow(assay(se_lfq_no_imp_educa)))
diff_6_lfq <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                       number=nrow(assay(se_lfq_no_imp_educa)))
diff_7_lfq <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                       number=nrow(assay(se_lfq_no_imp_educa)))
diff_8_lfq <- topTable(fit1,coef=5, adjust ="BH", sort.by="P",
                       number=nrow(assay(se_lfq_no_imp_educa)))
diff_9_lfq <- topTable(fit1,coef=6, adjust ="BH", sort.by="P",
                       number=nrow(assay(se_lfq_no_imp_educa)))

create_volcano_plot(diff_5_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                   "lfq_frail_educa__frail.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty"))
create_volcano_plot(diff_6_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                   "lfq_frail_educa__seco.png"),
                    title = paste("Differentially expressed proteins -",
                                  "Secondary"))
create_volcano_plot(diff_7_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                   "lfq_frail_educa__uni.png"),
                    title = paste("Differentially expressed proteins",
                                  "- University"))
create_volcano_plot(diff_8_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                   "lfq_frail_educa__frail_",
                                                   "seco.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Secondary"))
create_volcano_plot(diff_9_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                   "lfq_frail_educa__frail_",
                                                   "uni.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:University"))

diff_5_lfq <- diff_5_lfq[diff_5_lfq$adj.P.Val < 0.05, ]
nrow(diff_5_lfq)
# 4
rownames(diff_5_lfq)
# "P06702" "P05109" "Q8A015" "A2RC28"

diff_6_lfq <- diff_6_lfq[diff_6_lfq$adj.P.Val < 0.05, ]
nrow(diff_6_lfq)
# 0
rownames(diff_6_lfq)

diff_7_lfq <- diff_7_lfq[diff_7_lfq$adj.P.Val < 0.05, ]
nrow(diff_7_lfq)
# 0
rownames(diff_7_lfq)

diff_8_lfq <- diff_8_lfq[diff_8_lfq$adj.P.Val < 0.05, ]
nrow(diff_8_lfq)
# 0
rownames(diff_8_lfq)

diff_9_lfq <- diff_9_lfq[diff_9_lfq$adj.P.Val < 0.05, ]
nrow(diff_9_lfq)
#  0
rownames(diff_9_lfq)

# ==================================
#   Frailty and tobacco
# ==================================
#   A) MaxLFQ intensities
# ==================================
# Build model matrix considering frailty, tobacco and interaction
se_maxlfq_no_imp_tobacco <- se_maxlfq_no_imp[, !is.na(colData(se_maxlfq_no_imp)$tobacco)]
dim(se_maxlfq_no_imp_tobacco)
# 182 201
design = model.matrix(~ colData(se_maxlfq_no_imp_tobacco)$frailty*
                        colData(se_maxlfq_no_imp_tobacco)$tobacco)
colnames(design) = c("constant", "frailty", "former", "current",
                     "frailty:former", "frailty:current")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_maxlfq_no_imp_tobacco), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 6 coefficients per each column of the model matrix. Adjusted using BH
diff_10_max <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_tobacco)))
diff_11_max <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_tobacco)))
diff_12_max <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_tobacco)))
diff_13_max <- topTable(fit1,coef=5, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_tobacco)))
diff_14_max <- topTable(fit1,coef=6, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_tobacco)))

create_volcano_plot(diff_10_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_toba__frail.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_11_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_toba__former.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Tobacco former"))
create_volcano_plot(diff_12_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_toba__current.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Tobacco current"))
create_volcano_plot(diff_13_max, save_path = paste0(plot_path, "/volcano_plot_",
                                                "maxlfq_frail_toba__frail_",
                                                "former.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Tobacco former"))
create_volcano_plot(diff_14_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_toba__frail_",
                                                "current.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Tobacco current"))

diff_10_max <- diff_10_max[diff_10_max$adj.P.Val < 0.05, ]
nrow(diff_10_max)
# 0
rownames(diff_10_max)

diff_11_max <- diff_11_max[diff_11_max$adj.P.Val < 0.05, ]
nrow(diff_11_max)
# 0
rownames(diff_11_max)

diff_12_max <- diff_12_max[diff_12_max$adj.P.Val < 0.05, ]
nrow(diff_12_max)
# 0
rownames(diff_12_max)

diff_13_max <- diff_13_max[diff_13_max$adj.P.Val < 0.05, ]
nrow(diff_13_max)
# 0
rownames(diff_13_max)

diff_14_max <- diff_14_max[diff_14_max$adj.P.Val < 0.05, ]
nrow(diff_14_max)
# 79
rownames(diff_14_max)
# "A6KYH6" "P0C2E7" "C4ZBT7" "C4ZBL1" "A6KYK2" "C4ZB90" "A6L1X1" "A6KYJ7" "A9VT65" "C4Z2R7" "A9KJL3" "A9KRZ2" "A6LFQ4" "P22983" "A6KYH8" "C4Z2R9" "C4ZBU3" "C4ZD46" "A8YXK9" "C4ZBS5" "C4ZBD3" "C4Z2T1" "C4Z3R4" "A9KK92" "Q59309" "A6KYK6" "A9KK94" "A6TH53" "A6KYI1" "C4ZBG1" "Q8A9M2" "P0C0G6" "A6KXA0" "C4ZBT1" "A6KYH0" "P94316" "A6KYJ8" "Q8RQP4" "A9KJI8" "C4Z2T8" "C4ZBS1" "A9KKU0" "Q1Q899" "C4Z2V8" "A6KYJ2" "Q88XX2" "A6L048" "Q8A1A2" "B9E9L7" "A6L1L8" "C0R090" "C4ZBD5" "A9KRZ3" "A6KXL2" "A6KYK3" "B8I7Y6" "C4Z1J4" "A6L792" "P95544" "C4KZP0" "A6TWI7" "Q5L8C5" "A9KJI0" "A6KYI8" "A6KYK9" "P19543" "A6KYJ6" "Q5L9B6" "A6L0V1" "A6L903" "P02768" "Q5L923" "C4ZBS3" "Q18CF4" "A4XI37" "C4ZF71" "A9KJJ5" "A5I7J4" "P24295"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_14_max),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_14_max_ann <- merge(diff_14_max, prot_ann, by = "row.names")
row.names(diff_14_max_ann) <- diff_14_max_ann$Row.names
diff_14_max_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_14_max_ann <- diff_14_max_ann[order(diff_14_max_ann$adj.P.Val,
                                         -abs(diff_14_max_ann$logFC)), ]

# Save sig_prot_no_imp_maxlfq_ft_tobacco_current.csv
write.csv(diff_14_max_ann, file = paste0(output_path, "maxlfq_ft_tobacco_",
                                         "current.csv"),
          row.names = TRUE)

sel0 = which("P0C2E7"==rownames(assay(se_maxlfq_no_imp_tobacco)))
df0 = data.frame(frailty=colData(se_maxlfq_no_imp_tobacco)[,c("frailty")], 
                 tobacco=colData(se_maxlfq_no_imp_tobacco)[,c("tobacco")], 
                 expression=assay(se_maxlfq_no_imp_tobacco)[sel0,])
df0$frailty_tobacco <- paste(df0$frailty, df0$tobacco, sep = "_")
ggplot(df0,aes(x=frailty_tobacco,y=expression)) + geom_boxplot()

# ==================================
#   B) Top-N intensities
# ==================================
# Build model matrix considering frailty, tobacco and interaction
se_lfq_no_imp_tobacco <- se_lfq_no_imp[, !is.na(colData(se_lfq_no_imp)$tobacco)]
dim(se_lfq_no_imp_tobacco)
# 442 201
design = model.matrix(~ colData(se_lfq_no_imp_tobacco)$frailty*
                        colData(se_lfq_no_imp_tobacco)$tobacco)
colnames(design) = c("constant", "frailty", "former", "current",
                     "frailty:former", "frailty:current")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_lfq_no_imp_tobacco), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 6 coefficients per each column of the model matrix. Adjusted using BH
diff_10_lfq <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_tobacco)))
diff_11_lfq <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_tobacco)))
diff_12_lfq <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_tobacco)))
diff_13_lfq <- topTable(fit1,coef=5, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_tobacco)))
diff_14_lfq <- topTable(fit1,coef=6, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_tobacco)))

create_volcano_plot(diff_10_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_toba__frail.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_11_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_toba__former.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Tobacco former"))
create_volcano_plot(diff_12_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_toba__current.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Tobacco current"))
create_volcano_plot(diff_13_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_toba__frail_",
                                                    "former.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Tobacco former"))
create_volcano_plot(diff_14_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_toba__frail_",
                                                    "current.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Tobacco current"))

diff_10_lfq <- diff_10_lfq[diff_10_lfq$adj.P.Val < 0.05, ]
nrow(diff_10_lfq)
# 3
rownames(diff_10_lfq)
# "A2RC28" "Q8A015" "A9KJL3"

diff_11_lfq <- diff_11_lfq[diff_11_lfq$adj.P.Val < 0.05, ]
nrow(diff_11_lfq)
# 0
rownames(diff_11_lfq)

diff_12_lfq <- diff_12_lfq[diff_12_lfq$adj.P.Val < 0.05, ]
nrow(diff_12_lfq)
# 0
rownames(diff_12_lfq)

diff_13_lfq <- diff_13_lfq[diff_13_lfq$adj.P.Val < 0.05, ]
nrow(diff_13_lfq)
# 0
rownames(diff_13_lfq)

diff_14_lfq <- diff_14_lfq[diff_14_lfq$adj.P.Val < 0.05, ]
nrow(diff_14_lfq)
# 0
rownames(diff_14_lfq)

# ==================================
#   Frailty and alcohol
# ==================================
#   A) MaxLFQ intensities
# ==================================
# Build model matrix considering frailty, alcohol and interaction
se_maxlfq_no_imp_alcohol <- se_maxlfq_no_imp[, !is.na(colData(se_maxlfq_no_imp)$alcohol)]
dim(se_maxlfq_no_imp_alcohol)
# 182 195
design = model.matrix(~ colData(se_maxlfq_no_imp_alcohol)$frailty*
                        colData(se_maxlfq_no_imp_alcohol)$alcohol)
colnames(design) = c("constant", "frailty", "alcohol_monthly", "alcohol_weekly",
                     "frailty:alcohol_monthly", "frailty:alcohol_weekly")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_maxlfq_no_imp_alcohol), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 6 coefficients per each column of model matrix. Adjusted using BH
diff_15_max <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_alcohol)))
diff_16_max <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_alcohol)))
diff_17_max <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_alcohol)))
diff_18_max <- topTable(fit1,coef=5, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_alcohol)))
diff_19_max <- topTable(fit1,coef=6, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_alcohol)))

create_volcano_plot(diff_15_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_alco__frail.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_16_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_alco__monthly.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Alcohol monthly"))
create_volcano_plot(diff_17_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_alco__weekly.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Alcohol weekly"))
create_volcano_plot(diff_18_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_alco__frail_",
                                                "monthly.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Alcohol monthly"))
create_volcano_plot(diff_19_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_alco__frail_",
                                                "weekly.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Alcohol weekly"))

diff_15_max <- diff_15_max[diff_15_max$adj.P.Val < 0.05, ]
nrow(diff_15_max) 
# 0
rownames(diff_15_max)

diff_16_max <- diff_16_max[diff_16_max$adj.P.Val < 0.05, ]
nrow(diff_16_max)
# 0
rownames(diff_16_max)

diff_17_max <- diff_17_max[diff_17_max$adj.P.Val < 0.05, ]
nrow(diff_17_max)
# 1
rownames(diff_17_max)
# P55259

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_17_max),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_17_max_ann <- merge(diff_17_max, prot_ann, by = "row.names")
row.names(diff_17_max_ann) <- diff_17_max_ann$Row.names
diff_17_max_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_17_max_ann <- diff_17_max_ann[order(diff_17_max_ann$adj.P.Val,
                                         -abs(diff_17_max_ann$logFC)), ]

# Save sig_prot_no_imp_maxlfq_alcohol_weekly.csv
write.csv(diff_17_max_ann, file = paste0(output_path,
                                         "maxlfq_alcohol_weekly.csv"),
          row.names = TRUE)

diff_18_max <- diff_18_max[diff_18_max$adj.P.Val < 0.05, ]
nrow(diff_18_max)
# 0
rownames(diff_18_max)

diff_19_max <- diff_19_max[diff_19_max$adj.P.Val < 0.05, ]
nrow(diff_19_max) 
# 1
rownames(diff_19_max)
# "A6TH53"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_19_max),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_19_max_ann <- merge(diff_19_max, prot_ann, by = "row.names")
row.names(diff_19_max_ann) <- diff_19_max_ann$Row.names
diff_19_max_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_19_max_ann <- diff_19_max_ann[order(diff_19_max_ann$adj.P.Val,
                                         -abs(diff_19_max_ann$logFC)), ]

# Save sig_prot_no_imp_maxlfq_ft_alcohol_weekly.csv
write.csv(diff_19_max_ann, file = paste0(output_path,
                                         "maxlfq_ft_alcohol_weekly.csv"),
          row.names = TRUE)


sel0 = which("A6TH53"==rownames(assay(se_maxlfq_no_imp_alcohol)))
df0 = data.frame(frailty=colData(se_maxlfq_no_imp_alcohol)[,c("frailty")], 
                 alcohol=colData(se_maxlfq_no_imp_alcohol)[,c("alcohol")], 
                 expression=assay(se_maxlfq_no_imp_alcohol)[sel0,])
df0$frailty_alcohol <- paste(df0$frailty, df0$alcohol, sep = "_")
ggplot(df0,aes(x=frailty_alcohol,y=expression)) + geom_boxplot()

# ==================================
#   B) Top-N intensities
# ==================================
# Build model matrix considering frailty, alcohol and interaction
se_lfq_no_imp_alcohol <- se_lfq_no_imp[, !is.na(colData(se_lfq_no_imp)$alcohol)]
dim(se_lfq_no_imp_alcohol)
# 442 195
design = model.matrix(~ colData(se_lfq_no_imp_alcohol)$frailty*
                        colData(se_lfq_no_imp_alcohol)$alcohol)
colnames(design) = c("constant", "frailty", "alcohol_monthly", "alcohol_weekly",
                     "frailty:alcohol_monthly", "frailty:alcohol_weekly")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_lfq_no_imp_alcohol), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 6 coefficients per each column of model matrix. Adjusted using BH
diff_15_lfq <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_alcohol)))
diff_16_lfq <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_alcohol)))
diff_17_lfq <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_alcohol)))
diff_18_lfq <- topTable(fit1,coef=5, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_alcohol)))
diff_19_lfq <- topTable(fit1,coef=6, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_alcohol)))

create_volcano_plot(diff_15_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_alco__frail.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_16_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_alco__",
                                                    "monthly.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Alcohol monthly"))
create_volcano_plot(diff_17_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frailty_alco__",
                                                    "weekly.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Alcohol weekly"))
create_volcano_plot(diff_18_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_alco__frail_",
                                                    "monthly.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Alcohol monthly"))
create_volcano_plot(diff_19_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_alco__frail_",
                                                    "weekly.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Alcohol weekly"))

diff_15_lfq <- diff_15_lfq[diff_15_lfq$adj.P.Val < 0.05, ]
nrow(diff_15_lfq) 
# 4
rownames(diff_15_lfq)
# "B2UYT8" "Q5L8B5" "A9KJL3" "Q8A477"

diff_16_lfq <- diff_16_lfq[diff_16_lfq$adj.P.Val < 0.05, ]
nrow(diff_16_lfq)
# 0
rownames(diff_16_lfq)

diff_17_lfq <- diff_17_lfq[diff_17_lfq$adj.P.Val < 0.05, ]
nrow(diff_17_lfq)
# 0
rownames(diff_17_lfq)

diff_18_lfq <- diff_18_lfq[diff_18_lfq$adj.P.Val < 0.05, ]
nrow(diff_18_lfq)
# 0
rownames(diff_18_lfq)

diff_19_lfq <- diff_19_lfq[diff_19_lfq$adj.P.Val < 0.05, ]
nrow(diff_19_lfq) 
# 0
rownames(diff_19_lfq)

# ==================================
#   Frailty and diabetes
# ==================================
#   A) MaxLFQ intensities
# ==================================
# Build model matrix considering frailty, diabetes and interaction
se_maxlfq_no_imp_diabetes <- se_maxlfq_no_imp[, !is.na(colData(se_maxlfq_no_imp)$diabetes)]
dim(se_maxlfq_no_imp_diabetes)
# 182 201
design = model.matrix(~ colData(se_maxlfq_no_imp_diabetes)$frailty*
                        colData(se_maxlfq_no_imp_diabetes)$diabetes)
colnames(design) = c("constant", "frailty", "diabetes", "frailty:diabetes")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_maxlfq_no_imp_diabetes), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_20_max <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_maxlfq_no_imp_diabetes)))
diff_21_max <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_maxlfq_no_imp_diabetes)))
diff_22_max <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_maxlfq_no_imp_diabetes)))

create_volcano_plot(diff_20_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_diabe__frail.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_21_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_diabe__diabe.png"),
                    title = "Differentially expressed proteins - Diabetes")
create_volcano_plot(diff_22_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_diabe__frail_",
                                                "diabe.png"),
                    title = paste("Differentially expressed proteins -",
                                  "Frailty:Diabetes"))

diff_20_max <- diff_20_max[diff_20_max$adj.P.Val < 0.05, ]
nrow(diff_20_max)
# 1
rownames(diff_20_max)
# "A6KYJ0"

diff_21_max <- diff_21_max[diff_21_max$adj.P.Val < 0.05, ]
nrow(diff_21_max)
# 0
rownames(diff_21_max)

diff_22_max <- diff_22_max[diff_22_max$adj.P.Val < 0.05, ]
nrow(diff_22_max)
# 0
rownames(diff_22_max)

# ==================================
#   B) Top-N intensities
# ==================================
# Build model matrix considering frailty, diabetes and interaction
se_lfq_no_imp_diabetes <- se_lfq_no_imp[, !is.na(colData(se_lfq_no_imp)$diabetes)]
dim(se_lfq_no_imp_diabetes)
# 442 201
design = model.matrix(~ colData(se_lfq_no_imp_diabetes)$frailty*
                        colData(se_lfq_no_imp_diabetes)$diabetes)
colnames(design) = c("constant", "frailty", "diabetes", "frailty:diabetes")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_lfq_no_imp_diabetes), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_20_lfq <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_diabetes)))
diff_21_lfq <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_diabetes)))
diff_22_lfq <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_diabetes)))

create_volcano_plot(diff_20_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_diabe__frail.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_21_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_diabe__diabe.png"),
                    title = "Differentially expressed proteins - Diabetes")
create_volcano_plot(diff_22_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_diabe__frail_",
                                                    "diabe.png"),
                    title = paste("Differentially expressed proteins -",
                                  "Frailty:Diabetes"))

diff_20_lfq <- diff_20_lfq[diff_20_lfq$adj.P.Val < 0.05, ]
nrow(diff_20_lfq)
# 0
rownames(diff_20_lfq)

diff_21_lfq <- diff_21_lfq[diff_21_lfq$adj.P.Val < 0.05, ]
nrow(diff_21_lfq)
# 0
rownames(diff_21_lfq)

diff_22_lfq <- diff_22_lfq[diff_22_lfq$adj.P.Val < 0.05, ]
nrow(diff_22_lfq)
# 0
rownames(diff_22_lfq)

# ==================================
#   Frailty and chf
# ==================================
#   A) MaxLFQ intensities
# ==================================
# Build model matrix considering frailty, chf and interaction
se_maxlfq_no_imp_chf <- se_maxlfq_no_imp[, !is.na(colData(se_maxlfq_no_imp)$chf)]
dim(se_maxlfq_no_imp_chf)
# 182 199
design = model.matrix(~ colData(se_maxlfq_no_imp_chf)$frailty*
                        colData(se_maxlfq_no_imp_chf)$chf)
colnames(design) = c("constant", "frailty", "chf", "frailty:chf")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_maxlfq_no_imp_chf), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_23_max <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_chf)))
diff_24_max <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_chf)))
diff_25_max <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_chf)))

create_volcano_plot(diff_23_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_chf__frail.",
                                                    "png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_24_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_chf__chf.",
                                                    "png"),
                    title = "Differentially expressed proteins - Chf")
create_volcano_plot(diff_25_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_chf__frail_",
                                                    "chf.png"),
                    title = paste("Differentially expressed proteins -",
                                  "Frailty:Chf"))

diff_23_max <- diff_23_max[diff_23_max$adj.P.Val < 0.05, ]
nrow(diff_23_max)
# 5
rownames(diff_23_max)
# "P06702" "A6KYJ0" "P05109" "A9KJL3" "A6L792"

diff_24_max <- diff_24_max[diff_24_max$adj.P.Val < 0.05, ]
nrow(diff_24_max)
# 0
rownames(diff_24_max)

diff_25_max <- diff_25_max[diff_25_max$adj.P.Val < 0.05, ]
nrow(diff_25_max)
# 0
rownames(diff_25_max)

# ==================================
#   B) Top-N intensities
# ==================================
# Build model matrix considering frailty, chf and interaction
se_lfq_no_imp_chf <- se_lfq_no_imp[, !is.na(colData(se_lfq_no_imp)$chf)]
dim(se_lfq_no_imp_chf)
# 442 199
design = model.matrix(~ colData(se_lfq_no_imp_chf)$frailty*
                        colData(se_lfq_no_imp_chf)$chf)
colnames(design) = c("constant", "frailty", "chf", "frailty:chf")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_lfq_no_imp_chf), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_23_lfq <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_chf)))
diff_24_lfq <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_chf)))
diff_25_lfq <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_chf)))

create_volcano_plot(diff_23_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_chf__frail.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_24_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_chf__chf.png"),
                    title = "Differentially expressed proteins - Chf")
create_volcano_plot(diff_25_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_chf__frail_chf.",
                                                    "png"),
                    title = paste("Differentially expressed proteins -",
                                  "Frailty:Chf"))

diff_23_lfq <- diff_23_lfq[diff_23_lfq$adj.P.Val < 0.05, ]
nrow(diff_23_lfq)
# 22
rownames(diff_23_lfq)
# "P06702" "P05109" "A9KJL3" "A6LPS3" "A6L4M1" "Q9Z9L6" "A6L0U5" "A9KMF6" "Q05650" "B2UYT8" "A9KNK6" "P33656" "Q5LHW2" "Q8A015" "A6L7J5" "Q5L8B5" "A6L792" "A9KRZ1" "A2RC28" "Q5LH68" "Q5WLM8" "Q8A477"

diff_24_lfq <- diff_24_lfq[diff_24_lfq$adj.P.Val < 0.05, ]
nrow(diff_24_lfq)
# 1
rownames(diff_24_lfq)
# "A8F4S6"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_24_lfq),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_24_lfq_ann <- merge(diff_24_lfq, prot_ann, by = "row.names")
row.names(diff_24_lfq_ann) <- diff_24_lfq_ann$Row.names
diff_24_lfq_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_24_lfq_ann <- diff_24_lfq_ann[order(diff_24_lfq_ann$adj.P.Val,
                                         -abs(diff_24_lfq_ann$logFC)), ]

# Save sig_prot_no_imp_lfq_chf.csv
write.csv(diff_24_lfq_ann, file = paste0(output_path, "lfq_chf.csv"),
          row.names = TRUE)

diff_25_lfq <- diff_25_lfq[diff_25_lfq$adj.P.Val < 0.05, ]
nrow(diff_25_lfq)
# 0
rownames(diff_25_lfq)

# ==================================
#   Frailty and af
# ==================================
#   A) MaxLFQ intensities
# ==================================
# Build model matrix considering frailty, af and interaction
se_maxlfq_no_imp_af <- se_maxlfq_no_imp[, !is.na(colData(se_maxlfq_no_imp)$af)]
dim(se_maxlfq_no_imp_af)
# 182 201
design = model.matrix(~ colData(se_maxlfq_no_imp_af)$frailty*
                        colData(se_maxlfq_no_imp_af)$af)
colnames(design) = c("constant", "frailty", "af", "frailty:af")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_maxlfq_no_imp_af), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_26_max <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_af)))
diff_27_max <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_af)))
diff_28_max <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_af)))

create_volcano_plot(diff_26_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_af__frail.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_27_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_af__af.png"),
                    title = "Differentially expressed proteins - AF")
create_volcano_plot(diff_28_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_af__frail_af.png"),
                    title = paste("Differentially expressed proteins -",
                                  "Frailty:AF"))

diff_26_max <- diff_26_max[diff_26_max$adj.P.Val < 0.05, ]
nrow(diff_26_max)
# 6
rownames(diff_26_max)
# "P06702" "A6L792" "A6KYJ0" "Q8A9M2" "A9KRZ4" "P05109"

diff_27_max <- diff_27_max[diff_27_max$adj.P.Val < 0.05, ]
nrow(diff_27_max)
# 0
rownames(diff_27_max)

diff_28_max <- diff_28_max[diff_28_max$adj.P.Val < 0.05, ]
nrow(diff_28_max)
# 0
rownames(diff_28_max)

# ==================================
#   B) Top-N intensities
# ==================================
# Build model matrix considering frailty, af and interaction
se_lfq_no_imp_af <- se_lfq_no_imp[, !is.na(colData(se_lfq_no_imp)$af)]
dim(se_lfq_no_imp_af)
# 442 201
design = model.matrix(~ colData(se_lfq_no_imp_af)$frailty*
                        colData(se_lfq_no_imp_af)$af)
colnames(design) = c("constant", "frailty", "af", "frailty:af")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_lfq_no_imp_af), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_26_lfq <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_af)))
diff_27_lfq <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_af)))
diff_28_lfq <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_af)))

create_volcano_plot(diff_26_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_af__frail.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_27_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_af__af.png"),
                    title = "Differentially expressed proteins - AF")
create_volcano_plot(diff_28_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_af__frail_af.png"),
                    title = paste("Differentially expressed proteins -",
                                  "Frailty:AF"))

diff_26_lfq <- diff_26_lfq[diff_26_lfq$adj.P.Val < 0.05, ]
nrow(diff_26_lfq)
# 25
rownames(diff_26_lfq)
# "A9KJL3" "P06702" "P05109" "Q8A015" "Q05650" "A6L4M1" "A6L7J5" "A9KMF6" "A6L0U5" "Q9Z9L6" "A6KYG9" "B2UYT8" "A6LPS3" "P33656" "A6L792" "Q5LHW2" "Q5L8B5" "A9KRZ1" "A9KNK6" "Q8A477" "Q5WLM8" "A6L2R5" "Q5LH68" "A9KJI4" "C4Z2R3"

diff_27_lfq <- diff_27_lfq[diff_27_lfq$adj.P.Val < 0.05, ]
nrow(diff_27_lfq)
# 1
rownames(diff_27_lfq)
# "Q2T0I7"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_27_lfq),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_27_lfq_ann <- merge(diff_27_lfq, prot_ann, by = "row.names")
row.names(diff_27_lfq_ann) <- diff_27_lfq_ann$Row.names
diff_27_lfq_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_27_lfq_ann <- diff_27_lfq_ann[order(diff_27_lfq_ann$adj.P.Val,
                                         -abs(diff_27_lfq_ann$logFC)), ]

# Save sig_prot_no_imp_lfq_af.csv
write.csv(diff_27_lfq_ann, file = paste0(output_path, "lfq_af.csv"),
          row.names = TRUE)

diff_28_lfq <- diff_28_lfq[diff_28_lfq$adj.P.Val < 0.05, ]
nrow(diff_28_lfq)
# 2
rownames(diff_28_lfq)
# "Q2T0I7" "A6KYG9"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_28_lfq),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_28_lfq_ann <- merge(diff_28_lfq, prot_ann, by = "row.names")
row.names(diff_28_lfq_ann) <- diff_28_lfq_ann$Row.names
diff_28_lfq_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_28_lfq_ann <- diff_28_lfq_ann[order(diff_28_lfq_ann$adj.P.Val,
                                         -abs(diff_28_lfq_ann$logFC)), ]

# Save sig_prot_no_imp_lfq_ft_af.csv
write.csv(diff_28_lfq_ann, file = paste0(output_path, "lfq_ft_af.csv"),
          row.names = TRUE)


# ==================================
#   Frailty and osteoarthritis
# ==================================
#   A) MaxLFQ intensities
# ==================================
# Build model matrix considering frailty, osteoarthritis and interaction
se_maxlfq_no_imp_osteo <- se_maxlfq_no_imp[, !is.na(colData(se_maxlfq_no_imp)$osteoarthritis)]
dim(se_maxlfq_no_imp_osteo)
# 182 196
design = model.matrix(~ colData(se_maxlfq_no_imp_osteo)$frailty*
                        colData(se_maxlfq_no_imp_osteo)$osteoarthritis)
colnames(design) = c("constant", "frailty", "osteoarthritis",
                     "frailty:osteoarthritis")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_maxlfq_no_imp_osteo), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_29_max <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_osteo)))
diff_30_max <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_osteo)))
diff_31_max <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_osteo)))

create_volcano_plot(diff_29_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_osteo__",
                                                "frail.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_30_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_osteo__",
                                                "osteo.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Osteoarthritis"))
create_volcano_plot(diff_31_max, save_path = paste0(plot_path, "volcano_plot_",
                                                "maxlfq_frail_osteo__frail_",
                                                "osteo.png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:Osteoarthritis"))

diff_29_max <- diff_29_max[diff_29_max$adj.P.Val < 0.05, ]
nrow(diff_29_max)
# 13
rownames(diff_29_max)
# "A6KYJ4" "A9KJJ5" "P06702" "A6L792" "P05109" "A6L903" "Q5L8C5" "Q9AE24" "Q5L923" "A9KRZ4" "Q5LHW2" "A6L048" "Q8A9M2"

diff_30_max <- diff_30_max[diff_30_max$adj.P.Val < 0.05, ]
nrow(diff_30_max)
# 0
rownames(diff_30_max)

diff_31_max <- diff_31_max[diff_31_max$adj.P.Val < 0.05, ]
nrow(diff_31_max)
# 0
rownames(diff_31_max)

# ==================================
# B) Top-N intensities
# ==================================
# Build model matrix considering frailty, osteoarthritis and interaction
se_lfq_no_imp_osteo <- se_lfq_no_imp[, !is.na(colData(se_lfq_no_imp)$osteoarthritis)]
dim(se_lfq_no_imp_osteo)
# 442 196
design = model.matrix(~ colData(se_lfq_no_imp_osteo)$frailty*
                        colData(se_lfq_no_imp_osteo)$osteoarthritis)
colnames(design) = c("constant", "frailty", "osteoarthritis",
                     "frailty:osteoarthritis")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_lfq_no_imp_osteo), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_29_lfq <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_osteo)))
diff_30_lfq <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_osteo)))
diff_31_lfq <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_osteo)))

create_volcano_plot(diff_29_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_osteo__frail.",
                                                    "png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_30_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_osteo__osteo.",
                                                    "png"),
                    title = paste("Differentially expressed proteins",
                                  "- Osteoarthritis"))
create_volcano_plot(diff_31_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_osteo__frail_",
                                                    "osteo.png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:Osteoarthritis"))

diff_29_lfq <- diff_29_lfq[diff_29_lfq$adj.P.Val < 0.05, ]
nrow(diff_29_lfq)
# 1
rownames(diff_29_lfq)
# "Q8A015"

diff_30_lfq <- diff_30_lfq[diff_30_lfq$adj.P.Val < 0.05, ]
nrow(diff_30_lfq)
# 0
rownames(diff_30_lfq)

diff_31_lfq <- diff_31_lfq[diff_31_lfq$adj.P.Val < 0.05, ]
nrow(diff_31_lfq)
# 0
rownames(diff_31_lfq)

# ==================================
#   Frailty and hipfracture
# ==================================
#   A) MaxLFQ intensities
# ==================================
# Build model matrix considering frailty, hipfracture and interaction
se_maxlfq_no_imp_hip <- se_maxlfq_no_imp[, !is.na(colData(se_maxlfq_no_imp)$hipfracture)]
dim(se_maxlfq_no_imp_hip)
# 182 201
design = model.matrix(~ colData(se_maxlfq_no_imp_hip)$frailty*
                        colData(se_maxlfq_no_imp_hip)$hipfracture)
colnames(design) = c("constant", "frailty", "hipfracture",
                     "frailty:hipfracture")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_maxlfq_no_imp_hip), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_32_max <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_hip)))
diff_33_max <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_hip)))
diff_34_max <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_hip)))

create_volcano_plot(diff_32_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_hip__frail.",
                                                    "png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_33_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_hip__hip.",
                                                    "png"),
                    title = paste("Differentially expressed proteins",
                                  "- Hip fracture"))
create_volcano_plot(diff_34_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_hip__frail_",
                                                    "hip.png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:Hip fracture"))

diff_32_max <- diff_32_max[diff_32_max$adj.P.Val < 0.05, ]
nrow(diff_32_max)
# 7
rownames(diff_32_max)
# "P06702" "A6L792" "P05109" "A6KYJ0" "A9KRZ4" "A6KYJ4" "P94360"

diff_33_max <- diff_33_max[diff_33_max$adj.P.Val < 0.05, ]
nrow(diff_33_max)
# 0
rownames(diff_33_max)

diff_34_max <- diff_34_max[diff_34_max$adj.P.Val < 0.05, ]
nrow(diff_34_max)
# 0
rownames(diff_34_max)

# ==================================
#   B) Top-N intensities
# ==================================
# Build model matrix considering frailty, hipfracture and interaction
se_lfq_no_imp_hip <- se_lfq_no_imp[, !is.na(colData(se_lfq_no_imp)$hipfracture)]
dim(se_lfq_no_imp_hip)
# 442 201
design = model.matrix(~ colData(se_lfq_no_imp_hip)$frailty*
                        colData(se_lfq_no_imp_hip)$hipfracture)
colnames(design) = c("constant", "frailty", "hipfracture",
                     "frailty:hipfracture")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_lfq_no_imp_hip), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_32_lfq <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_hip)))
diff_33_lfq <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_hip)))
diff_34_lfq <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_hip)))

create_volcano_plot(diff_32_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_hip__frail.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_33_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_hip__hip.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Hip fracture"))
create_volcano_plot(diff_34_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_hip__frail_hip.",
                                                    "png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:Hip fracture"))

diff_32_lfq <- diff_32_lfq[diff_32_lfq$adj.P.Val < 0.05, ]
nrow(diff_32_lfq)
# 21
rownames(diff_32_lfq)
# "P06702" "A9KJL3" "P05109" "P33656" "B2UYT8" "A6L0U5" "Q8A015" "Q05650" "A9KMF6" "A9KNC4" "A9KJI4" "A6L4M1" "A6LPS3" "A6L7J5" "Q5L8B5" "A9KNK6" "Q5LHW2" "Q9Z9L6" "C4Z2R3" "A9KRZ1" "Q8A477"

diff_33_lfq <- diff_33_lfq[diff_33_lfq$adj.P.Val < 0.05, ]
nrow(diff_33_lfq)
# 2
rownames(diff_33_lfq)
# "Q3B6G3" "Q01523"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_33_lfq),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_33_lfq_ann <- merge(diff_33_lfq, prot_ann, by = "row.names")
row.names(diff_33_lfq_ann) <- diff_33_lfq_ann$Row.names
diff_33_lfq_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_33_lfq_ann <- diff_33_lfq_ann[order(diff_33_lfq_ann$adj.P.Val,
                                         -abs(diff_33_lfq_ann$logFC)), ]

# Save sig_prot_no_imp_lfq_hip.csv
write.csv(diff_33_lfq_ann, file = paste0(output_path, "lfq_hip.csv"),
          row.names = TRUE)

diff_34_lfq <- diff_34_lfq[diff_34_lfq$adj.P.Val < 0.05, ]
diff_34_lfq <- diff_34_lfq[!grepl("^NA(\\.|$)", rownames(diff_34_lfq)), ]
nrow(diff_34_lfq)
# 1
rownames(diff_34_lfq)
# "Q3B6G3"


# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_34_lfq),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_34_lfq_ann <- merge(diff_34_lfq, prot_ann, by = "row.names")
row.names(diff_34_lfq_ann) <- diff_34_lfq_ann$Row.names
diff_34_lfq_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_34_lfq_ann <- diff_34_lfq_ann[order(diff_34_lfq_ann$adj.P.Val,
                                         -abs(diff_34_lfq_ann$logFC)), ]

# Save sig_prot_no_imp_lfq_ft_hip.csv
write.csv(diff_34_lfq_ann, file = paste0(output_path, "lfq_ft_hip.csv"),
          row.names = TRUE)

# ==================================
#   Frailty and depression
# ==================================
#   A) MaxLFQ intensities
# ==================================
# Build model matrix considering frailty, depression and interaction
se_maxlfq_no_imp_depression <- se_maxlfq_no_imp[, !is.na(colData(se_maxlfq_no_imp)$depression)]
dim(se_maxlfq_no_imp_depression)
# 182 201
design = model.matrix(~ colData(se_maxlfq_no_imp_depression)$frailty*
                        colData(se_maxlfq_no_imp_depression)$depression)
colnames(design) = c("constant", "frailty", "depression", "frailty:depression")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_maxlfq_no_imp_depression), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_35_max <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_depression)))
diff_36_max <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_depression)))
diff_37_max <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_depression)))

create_volcano_plot(diff_35_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_depre__",
                                                    "frail.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_36_max, save_path = paste0(plot_path,"volcano_plot_",
                                                    "maxlfq_frail_depre__",
                                                    "depre.png"),
                    title = "Differentially expressed proteins - Depression")
create_volcano_plot(diff_37_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_depre__",
                                                    "frail_depre.png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:Depression"))

diff_35_max <- diff_35_max[diff_35_max$adj.P.Val < 0.05, ]
nrow(diff_35_max)
# 5
rownames(diff_35_max)
# "A6KYJ0" "P06702" "A9KRZ4" "P05109" "A6L792"

diff_36_max <- diff_36_max[diff_36_max$adj.P.Val < 0.05, ]
nrow(diff_36_max) 
# 0
rownames(diff_36_max)

diff_37_max <- diff_37_max[diff_37_max$adj.P.Val < 0.05, ]
nrow(diff_37_max)
# 0
rownames(diff_37_max)

# ==================================
#   B) Top-N intensities
# ==================================
# Build model matrix considering frailty, depression and interaction
se_lfq_no_imp_depression <- se_lfq_no_imp[, !is.na(colData(se_lfq_no_imp)$depression)]
dim(se_lfq_no_imp_depression)
# 442 201
design = model.matrix(~ colData(se_lfq_no_imp_depression)$frailty*
                        colData(se_lfq_no_imp_depression)$depression)
colnames(design) = c("constant", "frailty", "depression", "frailty:depression")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_lfq_no_imp_depression), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_35_lfq <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_depression)))
diff_36_lfq <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_depression)))
diff_37_lfq <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_depression)))

create_volcano_plot(diff_35_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_depre__frail.",
                                                    "png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_36_lfq, save_path = paste0(plot_path,"volcano_plot_",
                                                    "lfq_frail_depre__depre.",
                                                    "png"),
                    title = "Differentially expressed proteins - Depression")
create_volcano_plot(diff_37_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_depre__frail_",
                                                    "depre.png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:Depression"))

diff_35_lfq <- diff_35_lfq[diff_35_lfq$adj.P.Val < 0.05, ]
nrow(diff_35_lfq)
# 21
rownames(diff_35_lfq)
# "Q8A015" "A6L7J5" "P33656" "P05109" "Q9Z9L6" "A6L0U5" "A9KJL3" "A6LPS3" "P06702" "Q5LHW2" "Q5LH68" "B2UYT8" "A9KNK6" "A9KMF6" "A6LEJ1" "Q05650" "Q8A477" "A6L4M1" "A6KYJ0" "A9KNC4" "Q05203"

diff_36_lfq <- diff_36_lfq[diff_36_lfq$adj.P.Val < 0.05, ]
nrow(diff_36_lfq) 
# 0
rownames(diff_36_lfq)

diff_37_lfq <- diff_37_lfq[diff_37_lfq$adj.P.Val < 0.05, ]
diff_37_lfq <- diff_37_lfq[!grepl("^NA(\\.|$)", rownames(diff_37_lfq)), ]
nrow(diff_37_lfq)
# 0
rownames(diff_37_lfq)


# ==================================
#   Frailty and sarcopenia
# ==================================
#   A) MaxLFQ intensities
# ==================================
# Build model matrix considering frailty, sarcopenia and interaction
se_maxlfq_no_imp_sarco <- se_maxlfq_no_imp[, !is.na(colData(se_maxlfq_no_imp)$sarcopenia)]
dim(se_maxlfq_no_imp_sarco)
# 182 192
design = model.matrix(~ colData(se_maxlfq_no_imp_sarco)$frailty*
                        colData(se_maxlfq_no_imp_sarco)$sarcopenia)
colnames(design) = c("constant", "frailty", "sarcopenia", "frailty:sarcopenia")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_maxlfq_no_imp_sarco), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_38_max <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_sarco)))
diff_39_max <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_sarco)))
diff_40_max <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_sarco)))

create_volcano_plot(diff_38_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_sarco__",
                                                    "frail.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_39_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_sarco__",
                                                    "sarco.png"),
                    title = "Differentially expressed proteins - Sarcopenia")
create_volcano_plot(diff_40_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_sarco__frail",
                                                    "_sarco.png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:Sarcopenia"))

diff_38_max <- diff_38_max[diff_38_max$adj.P.Val < 0.05, ]
nrow(diff_38_max)
# 8
rownames(diff_38_max)
# "P06702" "A6KYJ0" "P05109" "A9KJL3" "O83023" "P94360" "A6L792" "C4Z2R3"

diff_39_max <- diff_39_max[diff_39_max$adj.P.Val < 0.05, ]
nrow(diff_39_max)
# 0
rownames(diff_39_max)

diff_40_max <- diff_40_max[diff_40_max$adj.P.Val < 0.05, ]
diff_40_max <- diff_40_max[!grepl("^NA(\\.|$)", rownames(diff_40_max)), ]
nrow(diff_40_max)
# 1
rownames(diff_40_max)
# "A6KYK7"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_40_max),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_40_max_ann <- merge(diff_40_max, prot_ann, by = "row.names")
row.names(diff_40_max_ann) <- diff_40_max_ann$Row.names
diff_40_max_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_40_max_ann <- diff_40_max_ann[order(diff_40_max_ann$adj.P.Val,
                                 -abs(diff_40_max_ann$logFC)), ]

# Save sig_prot_no_imp_maxlfq_ft_sarco.csv
write.csv(diff_40_max_ann, file = paste0(output_path, "maxlfq_ft_sarco.csv"),
          row.names = TRUE)

sel0 = which("A6KYK7"==rownames(assay(se_maxlfq_no_imp_sarco)))
df0 = data.frame(frailty=colData(se_maxlfq_no_imp_sarco)[,c("frailty")], 
                 sarco=colData(se_maxlfq_no_imp_sarco)[,c("sarcopenia")], 
                 expression=assay(se_maxlfq_no_imp_sarco)[sel0,])
df0$frailty_sarco <- paste(df0$frailty, df0$sarco, sep = "_")
ggplot(df0,aes(x=frailty_sarco,y=expression)) + geom_boxplot()

# ==================================
#   B) Top-N intensities
# ==================================
# Build model matrix considering frailty, sarcopenia and interaction
se_lfq_no_imp_sarco <- se_lfq_no_imp[, !is.na(colData(se_lfq_no_imp)$sarcopenia)]
dim(se_lfq_no_imp_sarco)
# 442 192
design = model.matrix(~ colData(se_lfq_no_imp_sarco)$frailty*
                        colData(se_lfq_no_imp_sarco)$sarcopenia)
colnames(design) = c("constant", "frailty", "sarcopenia", "frailty:sarcopenia")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_lfq_no_imp_sarco), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_38_lfq <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_sarco)))
diff_39_lfq <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_sarco)))
diff_40_lfq <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_sarco)))

create_volcano_plot(diff_38_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_sarco__frail.",
                                                    "png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_39_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_sarco__sarco.",
                                                    "png"),
                    title = "Differentially expressed proteins - Sarcopenia")
create_volcano_plot(diff_40_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_sarco__frail_",
                                                    "sarco.png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:Sarcopenia"))

diff_38_lfq <- diff_38_lfq[diff_38_lfq$adj.P.Val < 0.05, ]
nrow(diff_38_lfq)
# 10
rownames(diff_38_lfq)
# "P05109" "P06702" "A9KJL3" "Q05650" "A6L7J5" "C4Z2R3" "B2UYT8" "Q9Z9L6" "A9KMF6" "A9KNC4"

diff_39_lfq <- diff_39_lfq[diff_39_lfq$adj.P.Val < 0.05, ]
nrow(diff_39_lfq)
# 0
rownames(diff_39_lfq)

diff_40_lfq <- diff_40_lfq[diff_40_lfq$adj.P.Val < 0.05, ]
diff_40_lfq <- diff_40_lfq[!grepl("^NA(\\.|$)", rownames(diff_40_lfq)), ]
nrow(diff_40_lfq)
# 3
rownames(diff_40_lfq)
# "Q5L9E3" "A9KJI8" "A9KK94"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_40_lfq),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_40_lfq_ann <- merge(diff_40_lfq, prot_ann, by = "row.names")
row.names(diff_40_lfq_ann) <- diff_40_lfq_ann$Row.names
diff_40_lfq_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_40_lfq_ann <- diff_40_lfq_ann[order(diff_40_lfq_ann$adj.P.Val,
                                         -abs(diff_40_lfq_ann$logFC)), ]

# Save sig_prot_no_imp_lfq_ft_sarco.csv
write.csv(diff_40_lfq_ann, file = paste0(output_path, "lfq_ft_sarco.csv"),
          row.names = TRUE)

sel0 = which("A6KYK7"==rownames(assay(se_lfq_no_imp_sarco)))
df0 = data.frame(frailty=colData(se_lfq_no_imp_sarco)[,c("frailty")], 
                 sarco=colData(se_lfq_no_imp_sarco)[,c("sarcopenia")], 
                 expression=assay(se_lfq_no_imp_sarco)[sel0,])
df0$frailty_sarco <- paste(df0$frailty, df0$sarco, sep = "_")
ggplot(df0,aes(x=frailty_sarco,y=expression)) + geom_boxplot()

# ==================================
#   Frailty and ilef
# ==================================
#   A) MaxLFQ intensities
# ==================================
# Build model matrix considering frailty, ilef and interaction
se_maxlfq_no_imp_ilef <- se_maxlfq_no_imp[, !is.na(colData(se_maxlfq_no_imp)$ilef)]
dim(se_maxlfq_no_imp_ilef)
# 182 201
design = model.matrix(~ colData(se_maxlfq_no_imp_ilef)$frailty*
                        colData(se_maxlfq_no_imp_ilef)$ilef)
colnames(design) = c("constant", "frailty", "ilef", "frailty:ilef")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_maxlfq_no_imp_ilef), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_41_max <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_ilef)))
diff_42_max <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_ilef)))
diff_43_max <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_ilef)))

create_volcano_plot(diff_41_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_ilef__frail.",
                                                    "png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_42_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_ilef__ilef.",
                                                    "png"),
                    title = "Differentially expressed proteins - ILEF")
create_volcano_plot(diff_43_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_ilef__frail_",
                                                    "ilef.png"),
                    title = "Differentially expressed proteins - Frailty:ILEF")

diff_41_max <- diff_41_max[diff_41_max$adj.P.Val < 0.05, ]
nrow(diff_41_max)
# 2
rownames(diff_41_max)
# "A6KYJ0" "A6L792"

diff_42_max <- diff_42_max[diff_42_max$adj.P.Val < 0.05, ]
nrow(diff_42_max)
# 0
rownames(diff_42_max)

diff_43_max <- diff_43_max[diff_43_max$adj.P.Val < 0.05, ]
nrow(diff_43_max)
# 0
rownames(diff_43_max)

# ==================================
#   B) Top-N intensities
# ==================================
# Build model matrix considering frailty, ilef and interaction
se_lfq_no_imp_ilef <- se_lfq_no_imp[, !is.na(colData(se_lfq_no_imp)$ilef)]
dim(se_lfq_no_imp_ilef)
# 442 201
design = model.matrix(~ colData(se_lfq_no_imp_ilef)$frailty*
                        colData(se_lfq_no_imp_ilef)$ilef)
colnames(design) = c("constant", "frailty", "ilef", "frailty:ilef")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_lfq_no_imp_ilef), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_41_lfq <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_ilef)))
diff_42_lfq <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_ilef)))
diff_43_lfq <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_ilef)))

create_volcano_plot(diff_41_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_ilef__frail.",
                                                    "png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_42_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_ilef__ilef.png"),
                    title = "Differentially expressed proteins - ILEF")
create_volcano_plot(diff_43_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_ilef__frail_",
                                                    "ilef.png"),
                    title = "Differentially expressed proteins - Frailty:ILEF")

diff_41_lfq <- diff_41_lfq[diff_41_lfq$adj.P.Val < 0.05, ]
nrow(diff_41_lfq)
# 2
rownames(diff_41_lfq)
# "P05109" "P06702"

diff_42_lfq <- diff_42_lfq[diff_42_lfq$adj.P.Val < 0.05, ]
nrow(diff_42_lfq)
# 0
rownames(diff_42_lfq)

diff_43_lfq <- diff_43_lfq[diff_43_lfq$adj.P.Val < 0.05, ]
diff_43_lfq <- diff_43_lfq[!grepl("^NA(\\.|$)", rownames(diff_43_lfq)), ]
nrow(diff_43_lfq)
# 0
rownames(diff_43_lfq)


# ==================================
#   Frailty and age
# ==================================
#   A) MaxLFQ intensities
# ==================================
# Build model matrix considering frailty, age and interaction
se_maxlfq_no_imp_age <- se_maxlfq_no_imp[, !is.na(colData(se_maxlfq_no_imp)$age)]
dim(se_maxlfq_no_imp_age)
# 182 201
design = model.matrix(~ colData(se_maxlfq_no_imp_age)$frailty*
                        colData(se_maxlfq_no_imp_age)$age)
colnames(design) = c("constant", "frailty", "age", "frailty:age")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_maxlfq_no_imp_age), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_44_max <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_age)))
diff_45_max <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_age)))
diff_46_max <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_age)))

create_volcano_plot(diff_44_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_age__frail.",
                                                    "png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_45_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_age__age.",
                                                    "png"),
                    title = "Differentially expressed proteins- Age")
create_volcano_plot(diff_46_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_age__frail_",
                                                    "age.png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:age"))

diff_44_max <- diff_44_max[diff_44_max$adj.P.Val < 0.05, ]
nrow(diff_44_max)
# 0
rownames(diff_44_max)

diff_45_max <- diff_45_max[diff_45_max$adj.P.Val < 0.05, ]
nrow(diff_45_max)
# 0
rownames(diff_45_max)

diff_46_max <- diff_46_max[diff_46_max$adj.P.Val < 0.05, ]
nrow(diff_46_max)
# 0
rownames(diff_46_max)

# ==================================
#   B) Top-N intensities
# ==================================
# Build model matrix considering frailty, age and interaction
se_lfq_no_imp_age <- se_lfq_no_imp[, !is.na(colData(se_lfq_no_imp)$age)]
dim(se_lfq_no_imp_age)
# 442 201
design = model.matrix(~ colData(se_lfq_no_imp_age)$frailty*
                        colData(se_lfq_no_imp_age)$age)
colnames(design) = c("constant", "frailty", "age", "frailty:age")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_lfq_no_imp_age), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_44_lfq <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_age)))
diff_45_lfq <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_age)))
diff_46_lfq <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_age)))

create_volcano_plot(diff_44_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_age__frail.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_45_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_age__age.png"),
                    title = "Differentially expressed proteins- Age")
create_volcano_plot(diff_46_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_age__frail_age.",
                                                    "png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:age"))

diff_44_lfq <- diff_44_lfq[diff_44_lfq$adj.P.Val < 0.05, ]
nrow(diff_44_lfq)
# 0
rownames(diff_44_lfq)

diff_45_lfq <- diff_45_lfq[diff_45_lfq$adj.P.Val < 0.05, ]
nrow(diff_45_lfq)
# 0
rownames(diff_45_lfq)

diff_46_lfq <- diff_46_lfq[diff_46_lfq$adj.P.Val < 0.05, ]
nrow(diff_46_lfq)
# 0
rownames(diff_46_lfq)


# ==================================
#   Frailty and MEDAS
# ==================================
#   A) MaxLFQ intensities
# ==================================
# Build model matrix considering frailty, medas and interaction
se_maxlfq_no_imp_medas <- se_maxlfq_no_imp[, !is.na(colData(se_maxlfq_no_imp)$medas)]
dim(se_maxlfq_no_imp_medas)
# 182 189
design = model.matrix(~ colData(se_maxlfq_no_imp_medas)$frailty*
                        colData(se_maxlfq_no_imp_medas)$medas)
colnames(design) = c("constant", "frailty", "medas", "frailty:medas")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_maxlfq_no_imp_medas), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_47_max <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_medas)))
diff_48_max <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_medas)))
diff_49_max <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_medas)))

create_volcano_plot(diff_47_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_medas__",
                                                    "frail.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_48_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_medas__",
                                                    "medas.png"),
                    title = "Differentially expressed proteins- MEDAS")
create_volcano_plot(diff_49_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_medas__",
                                                    "frail_medas.png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:MEDAS"))

diff_47_max <- diff_47_max[diff_47_max$adj.P.Val < 0.05, ]
nrow(diff_47_max)
# 1
rownames(diff_47_max)
# "C4Z2V8"

diff_48_max <- diff_48_max[diff_48_max$adj.P.Val < 0.05, ]
nrow(diff_48_max)
# 0
rownames(diff_48_max)

diff_49_max <- diff_49_max[diff_49_max$adj.P.Val < 0.05, ]
nrow(diff_49_max)
# 1
rownames(diff_49_max)
# "C4Z2V8"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_49_max),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_49_max_ann <- merge(diff_49_max, prot_ann, by = "row.names")
row.names(diff_49_max_ann) <- diff_49_max_ann$Row.names
diff_49_max_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_49_max_ann <- diff_49_max_ann[order(diff_49_max_ann$adj.P.Val,
                                 -abs(diff_49_max_ann$logFC)), ]

# Save sig_prot_no_imp_maxlfq_ft_medas.csv
write.csv(diff_49_max_ann, file = paste0(output_path, "maxlfq_ft_medas.csv"),
          row.names = TRUE)

sel0 = which("C4Z2V8"==rownames(assay(se_maxlfq_no_imp_medas)))
df0 = data.frame(frailty=colData(se_maxlfq_no_imp_medas)[,c("frailty")], 
                 medas=colData(se_maxlfq_no_imp_medas)[,c("medas")], 
                 expression=assay(se_maxlfq_no_imp_medas)[sel0,])
ggplot(df0, aes(x = medas, y = expression, color = frailty)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", aes(linetype = frailty))

# ==================================
#   B) Top-N intensities
# ==================================
# Build model matrix considering frailty, medas and interaction
se_lfq_no_imp_medas <- se_lfq_no_imp[, !is.na(colData(se_lfq_no_imp)$medas)]
dim(se_lfq_no_imp_medas)
# 442 189
design = model.matrix(~ colData(se_lfq_no_imp_medas)$frailty*
                        colData(se_lfq_no_imp_medas)$medas)
colnames(design) = c("constant", "frailty", "medas", "frailty:medas")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_lfq_no_imp_medas), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_47_lfq <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_medas)))
diff_48_lfq <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_medas)))
diff_49_lfq <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_medas)))

create_volcano_plot(diff_47_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_medas__frail.",
                                                    "png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_48_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_medas__medas.",
                                                    "png"),
                    title = "Differentially expressed proteins- MEDAS")
create_volcano_plot(diff_49_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_medas__frail_",
                                                    "medas.png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:MEDAS"))

diff_47_lfq <- diff_47_lfq[diff_47_lfq$adj.P.Val < 0.05, ]
nrow(diff_47_lfq)
# 0
rownames(diff_47_lfq)

diff_48_lfq <- diff_48_lfq[diff_48_lfq$adj.P.Val < 0.05, ]
nrow(diff_48_lfq)
# 0
rownames(diff_48_lfq)

diff_49_lfq <- diff_49_lfq[diff_49_lfq$adj.P.Val < 0.05, ]
nrow(diff_49_lfq)
# 0
rownames(diff_49_lfq)

# ==================================
#   Frailty and energy
# ==================================
#   A) MaxLFQ intensities
# ==================================
# Build model matrix considering frailty, energy and interaction
se_maxlfq_no_imp_energy <- se_maxlfq_no_imp[, !is.na(colData(se_maxlfq_no_imp)$energy)]
dim(se_maxlfq_no_imp_energy)
# 182 189
design = model.matrix(~ colData(se_maxlfq_no_imp_energy)$frailty*
                        colData(se_maxlfq_no_imp_energy)$energy)
colnames(design) = c("constant", "frailty", "energy", "frailty:energy")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_maxlfq_no_imp_energy), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_50_max <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_energy)))
diff_51_max <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_energy)))
diff_52_max <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_energy)))

create_volcano_plot(diff_50_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_energy__",
                                                    "frail.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_51_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_energy__",
                                                    "energy.png"),
                    title = "Differentially expressed proteins - Energy")
create_volcano_plot(diff_52_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_energy__",
                                                    "frail_energy.png"),
                    title = paste("Differentially expressed proteins -",
                                  "Frailty:energy"))

diff_50_max <- diff_50_max[diff_50_max$adj.P.Val < 0.05, ]
nrow(diff_50_max)
# 0
rownames(diff_50_max)

diff_51_max <- diff_51_max[diff_51_max$adj.P.Val < 0.05, ]
nrow(diff_51_max)
# 0
rownames(diff_51_max)

diff_52_max <- diff_52_max[diff_52_max$adj.P.Val < 0.05, ]
nrow(diff_52_max)
# 0
rownames(diff_52_max)

# ==================================
#   B) Top-N intensities
# ==================================
# Build model matrix considering frailty, energy and interaction
se_lfq_no_imp_energy <- se_lfq_no_imp[, !is.na(colData(se_lfq_no_imp)$energy)]
dim(se_lfq_no_imp_energy)
# 442 189
design = model.matrix(~ colData(se_lfq_no_imp_energy)$frailty*
                        colData(se_lfq_no_imp_energy)$energy)
colnames(design) = c("constant", "frailty", "energy", "frailty:energy")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_lfq_no_imp_energy), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_50_lfq <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_energy)))
diff_51_lfq <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_energy)))
diff_52_lfq <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_energy)))

create_volcano_plot(diff_50_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_energy__frail.",
                                                    "png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_51_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_energy__energy.",
                                                    "png"),
                    title = "Differentially expressed proteins - Energy")
create_volcano_plot(diff_52_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_energy__frail_",
                                                    "energy.png"),
                    title = paste("Differentially expressed proteins -",
                                  "Frailty:energy"))

diff_50_lfq <- diff_50_lfq[diff_50_lfq$adj.P.Val < 0.05, ]
nrow(diff_50_lfq)
# 0
rownames(diff_50_lfq)

diff_51_lfq <- diff_51_lfq[diff_51_lfq$adj.P.Val < 0.05, ]
nrow(diff_51_lfq)
# 0
rownames(diff_51_lfq)

diff_52_lfq <- diff_52_lfq[diff_52_lfq$adj.P.Val < 0.05, ]
nrow(diff_52_lfq)
# 0
rownames(diff_52_lfq)


# ==================================
#   Frailty and bmi
# ==================================
#   A) MaxLFQ intensities
# ==================================
# Build model matrix considering frailty, bmi and interaction
se_maxlfq_no_imp_bmi <- se_maxlfq_no_imp[, !is.na(colData(se_maxlfq_no_imp)$bmi)]
dim(se_maxlfq_no_imp_bmi)
# 182 201
design = model.matrix(~ colData(se_maxlfq_no_imp_bmi)$frailty*
                        colData(se_maxlfq_no_imp_bmi)$bmi)
colnames(design) = c("constant", "frailty", "bmi", "frailty:bmi")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_maxlfq_no_imp_bmi), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_53_max <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_bmi)))
diff_54_max <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_bmi)))
diff_55_max <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_maxlfq_no_imp_bmi)))

create_volcano_plot(diff_53_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_bmi__frail.",
                                                    "png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_54_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_bmi__bmi.",
                                                    "png"),
                    title = "Differentially expressed proteins - BMI")
create_volcano_plot(diff_55_max, save_path = paste0(plot_path, "volcano_plot_",
                                                    "maxlfq_frail_bmi__frail_",
                                                    "bmi.png"),
                    title = paste("Differentially expressed proteins -",
                                  "Frailty:BMI"))

diff_53_max <- diff_53_max[diff_53_max$adj.P.Val < 0.05, ]
nrow(diff_53_max)
# 0
rownames(diff_53_max)

diff_54_max <- diff_54_max[diff_54_max$adj.P.Val < 0.05, ]
nrow(diff_54_max)
# 0
rownames(diff_54_max)

diff_55_max <- diff_55_max[diff_55_max$adj.P.Val < 0.05, ]
nrow(diff_55_max)
# 0
rownames(diff_55_max)

# ==================================
#   B) Top-N intensities
# ==================================
# Build model matrix considering frailty, bmi and interaction
se_lfq_no_imp_bmi <- se_lfq_no_imp[, !is.na(colData(se_lfq_no_imp)$bmi)]
dim(se_lfq_no_imp_bmi)
# 442 201
design = model.matrix(~ colData(se_lfq_no_imp_bmi)$frailty*
                        colData(se_lfq_no_imp_bmi)$bmi)
colnames(design) = c("constant", "frailty", "bmi", "frailty:bmi")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_lfq_no_imp_bmi), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_53_lfq <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_bmi)))
diff_54_lfq <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_bmi)))
diff_55_lfq <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                        number=nrow(assay(se_lfq_no_imp_bmi)))

create_volcano_plot(diff_53_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_bmi__frail.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_54_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_bmi__bmi.png"),
                    title = "Differentially expressed proteins - BMI")
create_volcano_plot(diff_55_lfq, save_path = paste0(plot_path, "volcano_plot_",
                                                    "lfq_frail_bmi__frail_",
                                                    "bmi.png"),
                    title = paste("Differentially expressed proteins -",
                                  "Frailty:BMI"))

diff_53_lfq <- diff_53_lfq[diff_53_lfq$adj.P.Val < 0.05, ]
nrow(diff_53_lfq)
# 4
rownames(diff_53_lfq)
# "P9WQH6" "E1WS50" "A6L048" "Q8WWA0"

diff_54_lfq <- diff_54_lfq[diff_54_lfq$adj.P.Val < 0.05, ]
nrow(diff_54_lfq)
# 0
rownames(diff_54_lfq)

diff_55_lfq <- diff_55_lfq[diff_55_lfq$adj.P.Val < 0.05, ]
nrow(diff_55_lfq)
# 3
rownames(diff_55_lfq)
# "E1WS50" "P9WQH6" "A6L048"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_55_lfq),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_55_lfq_ann <- merge(diff_55_lfq, prot_ann, by = "row.names")
row.names(diff_55_lfq_ann) <- diff_55_lfq_ann$Row.names
diff_55_lfq_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_55_lfq_ann <- diff_55_lfq_ann[order(diff_55_lfq_ann$adj.P.Val,
                                         -abs(diff_55_lfq_ann$logFC)), ]

# Save sig_prot_no_imp_lfq_ft_bmi.csv
write.csv(diff_55_lfq_ann, file = paste0(output_path, "lfq_ft_bmi.csv"),
          row.names = TRUE)
