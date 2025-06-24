# Load libraries
library(SummarizedExperiment)
library(ggplot2)
library(limma)
library(dplyr)
library(UniprotR)

# ==================================
# Running DE analysis with data from preprocessing
# ==================================
# A) Limma analysis from MaxLFQ
# ==================================
# - MaxLFQ,
# - log2transformed,
# - NA in at least 1 condition < 0.5,
# - only bacteria and human proteins,
# - no imputation
# ==================================
# B) Limma analysis from Top-N
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

# Choose option by introducing "A" or "B" in de_option variable:
de_option <- "A"

if (de_option == "A"){
  data_path <- paste0(work_path, "/prdea/data/maxlfq/")
  plot_path <- paste0(work_path, "/prdea/plots/limma_analysis/maxlfq/")
  output_path <- paste0(work_path, "/prdea/output/maxlfq/sig_prot_no_imp_")
  load(paste0(data_path, "/se_no_imp.RData"))
} else if (de_option == "B"){
  data_path <- paste0(work_path, "/data/lfq/")
  plot_path <- paste0(work_path, "/plots/limma_analysis/lfq/")
  output_file <- "sig_prot_no_imp_lfq_"
  load(paste0(data_path, "/se_no_imp.RData"))
} else {
  print("No valid option.")
}

# ==================================
# Frailty
# ==================================
# Build model matrix considering only frailty
design = model.matrix(~ colData(se_no_imp)$frailty)
colnames(design) = c("constant", "frailty")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp), design)

# Estimated coefficients for each protein
head(coef(fit))

# Moderate t-test of the differential expression by means of Bayes' empirical
# moderation of the standard errors to a global value
fit1 = eBayes(fit)

# p-values adjusted using the Benjamini-Hochberg method
diff_1 <- topTable(fit1, coef=2, adjust="BH", sort.by="P",
                   number=nrow(assay(se_no_imp)))

# Volcano plot
create_volcano_plot(diff_1, save_path = paste0(plot_path,
                                               "volcano_plot_frailty.png"),
                    title = "Differentially expressed proteins - Frailty")

diff_1 <- diff_1[diff_1$adj.P.Val < 0.05, ]
nrow(diff_1)
# A: 10
# B: 24
rownames(diff_1)
# A: "P06702" "A6KYJ0" "A6L792" "A9KRZ4" "P05109" "Q8A9M2" "P94360" "A6KYJ4" "C4KZP0" "A9KJL3"
# B: P06702" "A9KJL3" "P05109" "P33656" "Q8A015" "A9KMF6" "B2UYT8" "Q9Z9L6" "A6L4M1" "A6L0U5" "A6LPS3" "Q05650" "A9KNK6" "A9KNC4" "A6L7J5" "Q5LHW2" "Q5L8B5" "Q8A477" "Q5WLM8" "A9KRZ1" "A6L792" "C4Z2R3" "A9KJI4" "A6L2R5"

# Get annotation from Uniprot: gen, organism, protein, function, keywords
# and go terms
prot_ann <- GetProteinAnnontate(rownames(diff_1),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_1_ann <- merge(diff_1, prot_ann, by = "row.names")
row.names(diff_1_ann) <- diff_1_ann$Row.names
diff_1_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_1_ann <- diff_1_ann[order(diff_1_ann$adj.P.Val,
                                 -abs(diff_1_ann$logFC)), ]

# Save sig_prot_no_imp_(max)lfq_ft.csv
write.csv(diff_1_ann, file = paste0(output_path, "ft.csv"), row.names = TRUE)

sel0 = which("A9KJL3"==rownames(assay(se_no_imp)))
df0 = data.frame(frailty=colData(se_no_imp)[,c("frailty")],
                 expression=assay(se_no_imp)[sel0,])
ggplot(df0,aes(x=frailty,y=expression)) + geom_boxplot()

# ==================================
# Two predictor variables:
# ==================================
#   Fraily and sex
# ==================================
# Build model matrix considering frailty, sex and interaction
design = model.matrix(~ colData(se_no_imp)$frailty*colData(se_no_imp)$sex)
colnames(design) = c("constant", "frailty", "sex", "frailty:sex")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_2 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp)))
diff_3 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp)))
diff_4 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp)))

create_volcano_plot(diff_2, save_path = paste0(plot_path, "volcano_plot_",
                                               "frailty_sex__frailty.png"),
                     title = paste("Differentially expressed proteins",
                                   "- Frailty"))
create_volcano_plot(diff_3, save_path = paste0(plot_path, "volcano_plot_",
                                               "frailty_sex__sex.png"),
                    title = "Differentially expressed proteins - Sex")
create_volcano_plot(diff_4, save_path = paste0(plot_path, "volcano_plot_",
                                               "frailty_sex__fraity_sex.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Sex"))

diff_2 <- diff_2[diff_2$adj.P.Val < 0.05, ]
nrow(diff_2)
# A: 30
# B: 3
rownames(diff_2)
# A: "A6KYJ0" "C4KZP0" "A6L792" "A9KJJ5" "Q5L923" "A6KYJ7" "P05109" "A6KYJ4" "P95544" "A9KRZ4" "Q5L9B6" "C4ZBD3" "Q5L8C5" "Q8A9M2" "C4ZD46" "A8YXK9" "B9E9L7" "Q8A1A2" "A6KYH8" "A9KKU0" "P55990" "C4ZI85" "Q59199" "C4ZBG1" "P06702" "A6KYJ6" "A6KYH1" "C4ZF71" "Q18CF4" "C4Z0Q6"
# B: "A9KNC4" "Q05650" "A6LPS3"

diff_3 <- diff_3[diff_3$adj.P.Val < 0.05, ]
nrow(diff_3)
# A: 0
# B: 0
rownames(diff_3)

diff_4 <- diff_4[diff_4$adj.P.Val < 0.05, ]
nrow(diff_4)
# A: 1
# B: 0
rownames(diff_4)
# A: "C4Z0Q6"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_4),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_4_ann <- merge(diff_4, prot_ann, by = "row.names")
row.names(diff_4_ann) <- diff_4_ann$Row.names
diff_4_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_4_ann <- diff_4_ann[order(diff_4_ann$adj.P.Val,
                                 -abs(diff_4_ann$logFC)), ]

# Save sig_prot_no_imp_(max)lfq_ft_sex.csv
write.csv(diff_4_ann, file = paste0(output_path, "ft_sex.csv"), row.names = TRUE)

sel0 = which("C4Z0Q6"==rownames(assay(se_no_imp)))
df0 = data.frame(frailty=colData(se_no_imp)[,c("frailty")], 
                 sex=colData(se_no_imp)[,c("sex")], 
                 expression=assay(se_no_imp)[sel0,])
df0$frailty_sex <- paste(df0$frailty, df0$sex, sep = "_")
ggplot(df0,aes(x=frailty_sex,y=expression)) + geom_boxplot()

# ==================================
#   Frailty and education
# ==================================
# Build model matrix considering frailty, education and interaction
se_no_imp_educa <- se_no_imp[, !is.na(colData(se_no_imp)$education)]
dim(se_no_imp_educa)
# A: 182 201
# B: 442 201
design = model.matrix(~ colData(se_no_imp_educa)$frailty*
                        colData(se_no_imp_educa)$education)
colnames(design) = c("constant", "frailty", "secondary", "university",
                     "frailty:secondary", "frailty:university")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_educa), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 5 coefficients per each column of model matrix. Adjusted using BH
diff_5 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp_educa)))
diff_6 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp_educa)))
diff_7 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp_educa)))
diff_8 <- topTable(fit1,coef=5, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp_educa)))
diff_9 <- topTable(fit1,coef=6, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp_educa)))

create_volcano_plot(diff_5, save_path = paste0(plot_path, "volcano_plot_",
                                               "frailty_educa__frailty.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty"))
create_volcano_plot(diff_6, save_path = paste0(plot_path, "volcano_plot_",
                                               "frailty_educa__secondary.png"),
                    title = paste("Differentially expressed proteins -",
                                  "Secondary"))
create_volcano_plot(diff_7, save_path = paste0(plot_path, "volcano_plot_",
                                               "frailty_educa__university.png"),
                    title = paste("Differentially expressed proteins",
                                  "- University"))
create_volcano_plot(diff_8, save_path = paste0(plot_path, "volcano_plot_",
                                               "frailty_educa__frailty_",
                                               "secondary.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Secondary"))
create_volcano_plot(diff_9, save_path = paste0(plot_path, "volcano_plot_",
                                               "frailty_educa__fraity_",
                                               "university.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:University"))

diff_5 <- diff_5[diff_5$adj.P.Val < 0.05, ]
nrow(diff_5)
# A: 5
# B: 4
rownames(diff_5)
# A: "P06702" "A6KYJ4" "P05109" "A6L048" "A6L792"
# B: "P06702" "P05109" "Q8A015" "A2RC28"

diff_6 <- diff_6[diff_6$adj.P.Val < 0.05, ]
nrow(diff_6)
# A: 0
# B: 0
rownames(diff_6)

diff_7 <- diff_7[diff_7$adj.P.Val < 0.05, ]
nrow(diff_7)
# A: 0
# B: 0
rownames(diff_7)

diff_8 <- diff_8[diff_8$adj.P.Val < 0.05, ]
nrow(diff_8)
# A: 0
# B: 0
rownames(diff_8)

diff_9 <- diff_9[diff_9$adj.P.Val < 0.05, ]
nrow(diff_9)
# A: 0
# B: 0
rownames(diff_9)

# ==================================
#   Frailty and tobacco
# ==================================
# Build model matrix considering frailty, tobacco and interaction
se_no_imp_tobacco <- se_no_imp[, !is.na(colData(se_no_imp)$tobacco)]
dim(se_no_imp_tobacco)
# A: 182 201
# B: 442 201
design = model.matrix(~ colData(se_no_imp_tobacco)$frailty*
                        colData(se_no_imp_tobacco)$tobacco)
colnames(design) = c("constant", "frailty", "former", "current",
                     "frailty:former", "frailty:current")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_tobacco), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 6 coefficients per each column of the model matrix. Adjusted using BH
diff_10 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_tobacco)))
diff_11 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_tobacco)))
diff_12 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_tobacco)))
diff_13 <- topTable(fit1,coef=5, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_tobacco)))
diff_14 <- topTable(fit1,coef=6, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_tobacco)))

create_volcano_plot(diff_10, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_tobacco__frailty.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_11, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_tobacco__former.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Tobacco former"))
create_volcano_plot(diff_12, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_tobacco__current.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Tobacco current"))
create_volcano_plot(diff_13, save_path = paste0(plot_path, "/volcano_plot_",
                                                "frailty_tobacco__frailty_",
                                                "former.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Tobacco former"))
create_volcano_plot(diff_14, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_tobacco__frailty_",
                                                "current.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Tobacco current"))

diff_10 <- diff_10[diff_10$adj.P.Val < 0.05, ]
nrow(diff_10)
# A: 0
# B: 3
rownames(diff_10)
# B: "A2RC28" "Q8A015" "A9KJL3"

diff_11 <- diff_11[diff_11$adj.P.Val < 0.05, ]
nrow(diff_11)
# A: 0
# B: 0
rownames(diff_11)

diff_12 <- diff_12[diff_12$adj.P.Val < 0.05, ]
nrow(diff_12)
# A: 0
# B: 0
rownames(diff_12)

diff_13 <- diff_13[diff_13$adj.P.Val < 0.05, ]
nrow(diff_13)
# A: 0
# B: 0
rownames(diff_13)

diff_14 <- diff_14[diff_14$adj.P.Val < 0.05, ]
nrow(diff_14)
# A: 79
# B: 0
rownames(diff_14)
# A: "A6KYH6" "P0C2E7" "C4ZBT7" "C4ZBL1" "A6KYK2" "C4ZB90" "A6L1X1" "A6KYJ7" "A9VT65" "C4Z2R7" "A9KJL3" "A9KRZ2" "A6LFQ4" "P22983" "A6KYH8" "C4Z2R9" "C4ZBU3" "C4ZD46" "A8YXK9" "C4ZBS5" "C4ZBD3" "C4Z2T1" "C4Z3R4" "A9KK92" "Q59309" "A6KYK6" "A9KK94" "A6TH53" "A6KYI1" "C4ZBG1" "Q8A9M2" "P0C0G6" "A6KXA0" "C4ZBT1" "A6KYH0" "P94316" "A6KYJ8" "Q8RQP4" "A9KJI8" "C4Z2T8" "C4ZBS1" "A9KKU0" "Q1Q899" "C4Z2V8" "A6KYJ2" "Q88XX2" "A6L048" "Q8A1A2" "B9E9L7" "A6L1L8" "C0R090" "C4ZBD5" "A9KRZ3" "A6KXL2" "A6KYK3" "B8I7Y6" "C4Z1J4" "A6L792" "P95544" "C4KZP0" "A6TWI7" "Q5L8C5" "A9KJI0" "A6KYI8" "A6KYK9" "P19543" "A6KYJ6" "Q5L9B6" "A6L0V1" "A6L903" "P02768" "Q5L923" "C4ZBS3" "Q18CF4" "A4XI37" "C4ZF71" "A9KJJ5" "A5I7J4" "P24295"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_14),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_14_ann <- merge(diff_14, prot_ann, by = "row.names")
row.names(diff_14_ann) <- diff_14_ann$Row.names
diff_14_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_14_ann <- diff_14_ann[order(diff_14_ann$adj.P.Val,
                                 -abs(diff_14_ann$logFC)), ]

# Save sig_prot_no_imp_(max)lfq_ft_tobacco_current.csv
write.csv(diff_14_ann, file = paste0(output_path, "ft_tobacco_current.csv"),
          row.names = TRUE)

sel0 = which("P0C2E7"==rownames(assay(se_no_imp_tobacco)))
df0 = data.frame(frailty=colData(se_no_imp_tobacco)[,c("frailty")], 
                 tobacco=colData(se_no_imp_tobacco)[,c("tobacco")], 
                 expression=assay(se_no_imp_tobacco)[sel0,])
df0$frailty_tobacco <- paste(df0$frailty, df0$tobacco, sep = "_")
ggplot(df0,aes(x=frailty_tobacco,y=expression)) + geom_boxplot()

# ==================================
#   Frailty and alcohol
# ==================================
# Build model matrix considering frailty, alcohol and interaction
se_no_imp_alcohol <- se_no_imp[, !is.na(colData(se_no_imp)$alcohol)]
dim(se_no_imp_alcohol)
# A: 182 195
# B: 442 195
design = model.matrix(~ colData(se_no_imp_alcohol)$frailty*
                        colData(se_no_imp_alcohol)$alcohol)
colnames(design) = c("constant", "frailty", "alcohol_monthly", "alcohol_weekly",
                     "frailty:alcohol_monthly", "frailty:alcohol_weekly")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_alcohol), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 6 coefficients per each column of model matrix. Adjusted using BH
diff_15 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_alcohol)))
diff_16 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_alcohol)))
diff_17 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_alcohol)))
diff_18 <- topTable(fit1,coef=5, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_alcohol)))
diff_19 <- topTable(fit1,coef=6, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_alcohol)))

create_volcano_plot(diff_15, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_alcohol__frailty.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_16, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_alcohol__monthly.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Alcohol monthly"))
create_volcano_plot(diff_17, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_alcohol__weekly.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Alcohol weekly"))
create_volcano_plot(diff_18, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_alcohol__frailty_",
                                                "monthly.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Alcohol monthly"))
create_volcano_plot(diff_19, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_alcohol__frailty_",
                                                "weekly.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Frailty:Alcohol weekly"))

diff_15 <- diff_15[diff_15$adj.P.Val < 0.05, ]
nrow(diff_15) 
# A: 0
# B: 4
rownames(diff_15)
# B: "B2UYT8" "Q5L8B5" "A9KJL3" "Q8A477"


diff_16 <- diff_16[diff_16$adj.P.Val < 0.05, ]
nrow(diff_16)
# A: 0
# B: 0
rownames(diff_16)

diff_17 <- diff_17[diff_17$adj.P.Val < 0.05, ]
nrow(diff_17)
# A: 1
# B: 0
rownames(diff_17)
# A: P55259

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_17),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_17_ann <- merge(diff_17, prot_ann, by = "row.names")
row.names(diff_17_ann) <- diff_17_ann$Row.names
diff_17_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_17_ann <- diff_17_ann[order(diff_17_ann$adj.P.Val,
                                 -abs(diff_17_ann$logFC)), ]

# Save sig_prot_no_imp_(max)lfq_alcohol_weekly.csv
write.csv(diff_17_ann, file = paste0(output_path, "alcohol_weekly.csv"),
          row.names = TRUE)

diff_18 <- diff_18[diff_18$adj.P.Val < 0.05, ]
nrow(diff_18)
# A: 0
# B: 0
rownames(diff_18)

diff_19 <- diff_19[diff_19$adj.P.Val < 0.05, ]
nrow(diff_19) 
# A: 1
# B: 0
rownames(diff_19)
# A: "A6TH53"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_19),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_19_ann <- merge(diff_19, prot_ann, by = "row.names")
row.names(diff_19_ann) <- diff_19_ann$Row.names
diff_19_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_19_ann <- diff_19_ann[order(diff_19_ann$adj.P.Val,
                                 -abs(diff_19_ann$logFC)), ]

# Save sig_prot_no_imp_(max)lfq_ft_alcohol_weekly.csv
write.csv(diff_19_ann, file = paste0(output_path, "ft_alcohol_weekly.csv"),
          row.names = TRUE)


sel0 = which("A6TH53"==rownames(assay(se_no_imp_alcohol)))
df0 = data.frame(frailty=colData(se_no_imp_alcohol)[,c("frailty")], 
                 alcohol=colData(se_no_imp_alcohol)[,c("alcohol")], 
                 expression=assay(se_no_imp_alcohol)[sel0,])
df0$frailty_alcohol <- paste(df0$frailty, df0$alcohol, sep = "_")
ggplot(df0,aes(x=frailty_alcohol,y=expression)) + geom_boxplot()


# ==================================
#   Frailty and diabetes
# ==================================
# Build model matrix considering frailty, diabetes and interaction
se_no_imp_diabetes <- se_no_imp[, !is.na(colData(se_no_imp)$diabetes)]
dim(se_no_imp_diabetes)
# A: 182 201
# B: 442 201
design = model.matrix(~ colData(se_no_imp_diabetes)$frailty*
                        colData(se_no_imp_diabetes)$diabetes)
colnames(design) = c("constant", "frailty", "diabetes", "frailty:diabetes")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_diabetes), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_20 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp_diabetes)))
diff_21 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp_diabetes)))
diff_22 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp_diabetes)))

create_volcano_plot(diff_20, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_diabetes__frailty.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_21, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_diabetes__diabetes.png"),
                    title = "Differentially expressed proteins - Diabetes")
create_volcano_plot(diff_22, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_diabetes__frailty_",
                                                "diabetes.png"),
                    title = paste("Differentially expressed proteins -",
                                  "Frailty:Diabetes"))

diff_20 <- diff_20[diff_20$adj.P.Val < 0.05, ]
nrow(diff_20)
# A: 1
# B: 0
rownames(diff_20)
# A: "A6KYJ0"

diff_21 <- diff_21[diff_21$adj.P.Val < 0.05, ]
nrow(diff_21)
# A: 0
# B: 0
rownames(diff_21)

diff_22 <- diff_22[diff_22$adj.P.Val < 0.05, ]
nrow(diff_22)
# A: 0
# B: 0
rownames(diff_22)

# ==================================
#   Frailty and chf
# ==================================
# Build model matrix considering frailty, chf and interaction
se_no_imp_chf <- se_no_imp[, !is.na(colData(se_no_imp)$chf)]
dim(se_no_imp_chf)
# A: 182 199
# B: 442 199
design = model.matrix(~ colData(se_no_imp_chf)$frailty*
                        colData(se_no_imp_chf)$chf)
colnames(design) = c("constant", "frailty", "chf", "frailty:chf")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_chf), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_23 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_chf)))
diff_24 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_chf)))
diff_25 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_chf)))

create_volcano_plot(diff_23, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_chf__frailty.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_24, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_chf__chf.png"),
                    title = "Differentially expressed proteins - Chf")
create_volcano_plot(diff_25, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_chf__frailty_chf.png"),
                    title = paste("Differentially expressed proteins -",
                                  "Frailty:Chf"))

diff_23 <- diff_23[diff_23$adj.P.Val < 0.05, ]
nrow(diff_23)
# A: 5
# B: 22
rownames(diff_23)
# A: "P06702" "A6KYJ0" "P05109" "A9KJL3" "A6L792"
# B: "P06702" "P05109" "A9KJL3" "A6LPS3" "A6L4M1" "Q9Z9L6" "A6L0U5" "A9KMF6" "Q05650" "B2UYT8" "A9KNK6" "P33656" "Q5LHW2" "Q8A015" "A6L7J5" "Q5L8B5" "A6L792" "A9KRZ1" "A2RC28" "Q5LH68" "Q5WLM8" "Q8A477"

diff_24 <- diff_24[diff_24$adj.P.Val < 0.05, ]
nrow(diff_24)
# A: 0
# B: 1
rownames(diff_24)
# B: "A8F4S6"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_24),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_24_ann <- merge(diff_24, prot_ann, by = "row.names")
row.names(diff_24_ann) <- diff_24_ann$Row.names
diff_24_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_24_ann <- diff_24_ann[order(diff_24_ann$adj.P.Val,
                                 -abs(diff_24_ann$logFC)), ]

# Save sig_prot_no_imp_(max)lfq_chf.csv
write.csv(diff_24_ann, file = paste0(output_path, "chf.csv"),
          row.names = TRUE)

diff_25 <- diff_25[diff_25$adj.P.Val < 0.05, ]
nrow(diff_25)
# A: 0
# B: 0
rownames(diff_25)

# ==================================
#   Frailty and af
# ==================================
# Build model matrix considering frailty, af and interaction
se_no_imp_af <- se_no_imp[, !is.na(colData(se_no_imp)$af)]
dim(se_no_imp_af)
# A: 182 201
# B: 442 201
design = model.matrix(~ colData(se_no_imp_af)$frailty*
                        colData(se_no_imp_af)$af)
colnames(design) = c("constant", "frailty", "af", "frailty:af")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_af), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_26 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_af)))
diff_27 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_af)))
diff_28 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_af)))

create_volcano_plot(diff_26, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_af__frailty.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_27, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_af__af.png"),
                    title = "Differentially expressed proteins - AF")
create_volcano_plot(diff_28, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_af__frailty_af.png"),
                    title = paste("Differentially expressed proteins -",
                                  "Frailty:AF"))

diff_26 <- diff_26[diff_26$adj.P.Val < 0.05, ]
nrow(diff_26)
# A: 6
# B: 24
rownames(diff_26)
# A: "P06702" "A6L792" "A6KYJ0" "Q8A9M2" "A9KRZ4" "P05109"
# B: "A9KJL3" "P06702" "P05109" "Q8A015" "Q05650" "A6L4M1" "A6L7J5" "A9KMF6" "A6L0U5" "Q9Z9L6" "A6KYG9" "B2UYT8" "A6LPS3" "P33656" "A6L792" "Q5LHW2" "Q5L8B5" "A9KRZ1" "A9KNK6" "Q8A477" "Q5WLM8" "A6L2R5" "Q5LH68" "A9KJI4" "C4Z2R3"

diff_27 <- diff_27[diff_27$adj.P.Val < 0.05, ]
nrow(diff_27)
# A: 0
# B: 1
rownames(diff_27)
# B: "Q2T0I7"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_27),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_27_ann <- merge(diff_27, prot_ann, by = "row.names")
row.names(diff_27_ann) <- diff_27_ann$Row.names
diff_27_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_27_ann <- diff_27_ann[order(diff_27_ann$adj.P.Val,
                                 -abs(diff_27_ann$logFC)), ]

# Save sig_prot_no_imp_(max)lfq_af.csv
write.csv(diff_27_ann, file = paste0(output_path, "af.csv"),
          row.names = TRUE)

diff_28 <- diff_28[diff_28$adj.P.Val < 0.05, ]
nrow(diff_28)
# A: 0
# B: 2
rownames(diff_28)
# B: "Q2T0I7" "A6KYG9"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_28),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_28_ann <- merge(diff_28, prot_ann, by = "row.names")
row.names(diff_28_ann) <- diff_28_ann$Row.names
diff_28_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_28_ann <- diff_28_ann[order(diff_28_ann$adj.P.Val,
                                 -abs(diff_28_ann$logFC)), ]

# Save sig_prot_no_imp_(max)lfq_ft_af.csv
write.csv(diff_28_ann, file = paste0(output_path, "ft_af.csv"),
          row.names = TRUE)

# ==================================
#   Frailty and osteoarthritis
# ==================================
# Build model matrix considering frailty, osteoarthritis and interaction
se_no_imp_osteo <- se_no_imp[, !is.na(colData(se_no_imp)$osteoarthritis)]
dim(se_no_imp_osteo)
# A: 182 196
# B: 442 196
design = model.matrix(~ colData(se_no_imp_osteo)$frailty*
                        colData(se_no_imp_osteo)$osteoarthritis)
colnames(design) = c("constant", "frailty", "osteoarthritis",
                     "frailty:osteoarthritis")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_osteo), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_29 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_osteo)))
diff_30 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_osteo)))
diff_31 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_osteo)))

create_volcano_plot(diff_29, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_osteo__frailty.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_30, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_osteo__osteo.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Osteoarthritis"))
create_volcano_plot(diff_31, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_osteo__frailty_",
                                                "osteo.png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:Osteoarthritis"))

diff_29 <- diff_29[diff_29$adj.P.Val < 0.05, ]
nrow(diff_29)
# A: 13
# B: 1
rownames(diff_29)
# A: "A6KYJ4" "A9KJJ5" "P06702" "A6L792" "P05109" "A6L903" "Q5L8C5" "Q9AE24" "Q5L923" "A9KRZ4" "Q5LHW2" "A6L048" "Q8A9M2"
# B: "Q8A015"

diff_30 <- diff_30[diff_30$adj.P.Val < 0.05, ]
nrow(diff_30)
# A: 0
# B: 0
rownames(diff_30)

diff_31 <- diff_31[diff_31$adj.P.Val < 0.05, ]
nrow(diff_31)
# A: 0
# B: 0
rownames(diff_31)

# ==================================
#   Frailty and hipfracture
# ==================================
# Build model matrix considering frailty, hipfracture and interaction
se_no_imp_hip <- se_no_imp[, !is.na(colData(se_no_imp)$hipfracture)]
dim(se_no_imp_hip)
# A: 182 201
# B: 442 201
design = model.matrix(~ colData(se_no_imp_hip)$frailty*
                        colData(se_no_imp_hip)$hipfracture)
colnames(design) = c("constant", "frailty", "hipfracture",
                     "frailty:hipfracture")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_hip), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_32 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_hip)))
diff_33 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_hip)))
diff_34 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_hip)))

create_volcano_plot(diff_32, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_hip__frailty.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_33, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_hip__hip.png"),
                    title = paste("Differentially expressed proteins",
                                  "- Hip fracture"))
create_volcano_plot(diff_34, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_hip__frailty_",
                                                "hip.png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:Hip fracture"))

diff_32 <- diff_32[diff_32$adj.P.Val < 0.05, ]
nrow(diff_32)
# A: 7
# B: 21
rownames(diff_32)
# A: "P06702" "A6L792" "P05109" "A6KYJ0" "A9KRZ4" "A6KYJ4" "P94360"
# B: "P06702" "A9KJL3" "P05109" "P33656" "B2UYT8" "A6L0U5" "Q8A015" "Q05650" "A9KMF6" "A9KNC4" "A9KJI4" "A6L4M1" "A6LPS3" "A6L7J5" "Q5L8B5" "A9KNK6" "Q5LHW2" "Q9Z9L6" "C4Z2R3" "A9KRZ1" "Q8A477"

diff_33 <- diff_33[diff_33$adj.P.Val < 0.05, ]
nrow(diff_33)
# A: 0
# B: 2
rownames(diff_33)
# B: "Q3B6G3" "Q01523"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_33),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_33_ann <- merge(diff_33, prot_ann, by = "row.names")
row.names(diff_33_ann) <- diff_33_ann$Row.names
diff_33_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_33_ann <- diff_33_ann[order(diff_33_ann$adj.P.Val,
                                 -abs(diff_33_ann$logFC)), ]

# Save sig_prot_no_imp_(max)lfq_hip.csv
write.csv(diff_33_ann, file = paste0(output_path, "hip.csv"),
          row.names = TRUE)

diff_34 <- diff_34[diff_34$adj.P.Val < 0.05, ]
nrow(diff_34)
# A: 0
# B: 1
rownames(diff_34)
# B: "Q3B6G3"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_34),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_34_ann <- merge(diff_34, prot_ann, by = "row.names")
row.names(diff_34_ann) <- diff_34_ann$Row.names
diff_34_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_34_ann <- diff_34_ann[order(diff_34_ann$adj.P.Val,
                                 -abs(diff_34_ann$logFC)), ]

# Save sig_prot_no_imp_(max)lfq_ft_hip.csv
write.csv(diff_34_ann, file = paste0(output_path, "ft_hip.csv"),
          row.names = TRUE)

# ==================================
#   Frailty and depression
# ==================================
# Build model matrix considering frailty, depression and interaction
se_no_imp_depression <- se_no_imp[, !is.na(colData(se_no_imp)$depression)]
dim(se_no_imp_depression)
# A: 182 201
# B: 442 201
design = model.matrix(~ colData(se_no_imp_depression)$frailty*
                        colData(se_no_imp_depression)$depression)
colnames(design) = c("constant", "frailty", "depression", "frailty:depression")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_depression), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_35 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_depression)))
diff_36 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_depression)))
diff_37 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_depression)))

create_volcano_plot(diff_35, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_depression__",
                                                "frailty.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_36, save_path = paste0(plot_path,"volcano_plot_",
                                                "frailty_depression__",
                                                "depression.png"),
                    title = "Differentially expressed proteins - Depression")
create_volcano_plot(diff_37, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_depression__",
                                                "frailty_depression.png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:Depression"))

diff_35 <- diff_35[diff_35$adj.P.Val < 0.05, ]
nrow(diff_35)
# A: 5
# B: 21
rownames(diff_35)
# A: "A6KYJ0" "P06702" "A9KRZ4" "P05109" "A6L792"
# B: "Q8A015" "A6L7J5" "P33656" "P05109" "Q9Z9L6" "A6L0U5" "A9KJL3" "A6LPS3" "P06702" "Q5LHW2" "Q5LH68" "B2UYT8" "A9KNK6" "A9KMF6" "A6LEJ1" "Q05650" "Q8A477" "A6L4M1" "A6KYJ0" "A9KNC4" "Q05203"

diff_36 <- diff_36[diff_36$adj.P.Val < 0.05, ]
nrow(diff_36) 
# A: 0
# B: 0
rownames(diff_36)

diff_37 <- diff_37[diff_37$adj.P.Val < 0.05, ]
nrow(diff_37)
# A: 0
# B: 0
rownames(diff_37)

# ==================================
#   Frailty and sarcopenia
# ==================================
# Build model matrix considering frailty, sarcopenia and interaction
se_no_imp_sarco <- se_no_imp[, !is.na(colData(se_no_imp)$sarcopenia)]
dim(se_no_imp_sarco)
# A: 182 192
# B: 442 192
design = model.matrix(~ colData(se_no_imp_sarco)$frailty*
                        colData(se_no_imp_sarco)$sarcopenia)
colnames(design) = c("constant", "frailty", "sarcopenia", "frailty:sarcopenia")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_sarco), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_38 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_sarco)))
diff_39 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_sarco)))
diff_40 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_sarco)))

create_volcano_plot(diff_38, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_sarco__frailty.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_39, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_sarco__sarco.png"),
                    title = "Differentially expressed proteins - Sarcopenia")
create_volcano_plot(diff_40, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_sarco__",
                                                "frailty_sarco.png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:Sarcopenia"))

diff_38 <- diff_38[diff_38$adj.P.Val < 0.05, ]
nrow(diff_38)
# A: 8
# B: 10
rownames(diff_38)
# A: "P06702" "A6KYJ0" "P05109" "A9KJL3" "O83023" "P94360" "A6L792" "C4Z2R3"
# B: "P05109" "P06702" "A9KJL3" "Q05650" "A6L7J5" "C4Z2R3" "B2UYT8" "Q9Z9L6" "A9KMF6" "A9KNC4"

diff_39 <- diff_39[diff_39$adj.P.Val < 0.05, ]
nrow(diff_39)
# A: 0
# B: 0
rownames(diff_39)

diff_40 <- diff_40[diff_40$adj.P.Val < 0.05, ]
nrow(diff_40)
# A: 1
# B: 3
rownames(diff_40)
# A: "A6KYK7"
# B: "Q5L9E3" "A9KJI8" "A9KK94"

diff_40 <- diff_40[!grepl("^NA(\\.|$)", rownames(diff_40)), ]


# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_40),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_40_ann <- merge(diff_40, prot_ann, by = "row.names")
row.names(diff_40_ann) <- diff_40_ann$Row.names
diff_40_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_40_ann <- diff_40_ann[order(diff_40_ann$adj.P.Val,
                                 -abs(diff_40_ann$logFC)), ]

# Save sig_prot_no_imp_(max)lfq_ft_sarco.csv
write.csv(diff_40_ann, file = paste0(output_path, "ft_sarco.csv"),
          row.names = TRUE)

sel0 = which("A6KYK7"==rownames(assay(se_no_imp_sarco)))
df0 = data.frame(frailty=colData(se_no_imp_sarco)[,c("frailty")], 
                 sarco=colData(se_no_imp_sarco)[,c("sarcopenia")], 
                 expression=assay(se_no_imp_sarco)[sel0,])
df0$frailty_sarco <- paste(df0$frailty, df0$sarco, sep = "_")
ggplot(df0,aes(x=frailty_sarco,y=expression)) + geom_boxplot()

# ==================================
#   Frailty and ilef
# ==================================
# Build model matrix considering frailty, ilef and interaction
se_no_imp_ilef <- se_no_imp[, !is.na(colData(se_no_imp)$ilef)]
dim(se_no_imp_ilef)
# A: 182 201
# B: 442 201
design = model.matrix(~ colData(se_no_imp_ilef)$frailty*
                        colData(se_no_imp_ilef)$ilef)
colnames(design) = c("constant", "frailty", "ilef", "frailty:ilef")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_ilef), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_41 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_ilef)))
diff_42 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_ilef)))
diff_43 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_ilef)))

create_volcano_plot(diff_41, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_ilef__frailty.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_42, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_ilef__ilef.png"),
                    title = "Differentially expressed proteins - ILEF")
create_volcano_plot(diff_43, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_ilef__frailty_",
                                                "ilef.png"),
                    title = "Differentially expressed proteins - Frailty:ILEF")

diff_41 <- diff_41[diff_41$adj.P.Val < 0.05, ]
nrow(diff_41)
# A: 2
# B: 2
rownames(diff_41)
# A: "A6KYJ0" "A6L792"
# B: "P05109" "P06702"

diff_42 <- diff_42[diff_42$adj.P.Val < 0.05, ]
nrow(diff_42)
# A: 0
# B: 0
rownames(diff_42)

diff_43 <- diff_43[diff_43$adj.P.Val < 0.05, ]
nrow(diff_43)
# A: 0
# B: 0
rownames(diff_43)

# ==================================
#   Frailty and age
# ==================================
# Build model matrix considering frailty, age and interaction
se_no_imp_age <- se_no_imp[, !is.na(colData(se_no_imp)$age)]
dim(se_no_imp_age)
# A: 182 201
# B: 442 201
design = model.matrix(~ colData(se_no_imp_age)$frailty*
                        colData(se_no_imp_age)$age)
colnames(design) = c("constant", "frailty", "age", "frailty:age")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_age), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_44 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_age)))
diff_45 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_age)))
diff_46 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_age)))

create_volcano_plot(diff_44, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_age__frailty.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_45, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_age__age.png"),
                    title = "Differentially expressed proteins- Age")
create_volcano_plot(diff_46, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_age__frailty_age.png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:age"))

diff_44 <- diff_44[diff_44$adj.P.Val < 0.05, ]
nrow(diff_44)
# A: 2 -> 0
# B: 0
rownames(diff_44)
# A: "A6KYJ0" "A6L792"

diff_45 <- diff_45[diff_45$adj.P.Val < 0.05, ]
nrow(diff_45)
# A: 0
# B: 0
rownames(diff_45)

diff_46 <- diff_46[diff_46$adj.P.Val < 0.05, ]
nrow(diff_46)
# A: 0
# B: 0
rownames(diff_46)

# ==================================
#   Frailty and MEDAS
# ==================================
# Build model matrix considering frailty, medas and interaction
se_no_imp_medas <- se_no_imp[, !is.na(colData(se_no_imp)$medas)]
dim(se_no_imp_medas)
# A: 182 189
# B: 442 189
design = model.matrix(~ colData(se_no_imp_medas)$frailty*
                        colData(se_no_imp_medas)$medas)
colnames(design) = c("constant", "frailty", "medas", "frailty:medas")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_medas), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_47 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_medas)))
diff_48 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_medas)))
diff_49 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_medas)))

create_volcano_plot(diff_47, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_medas__frailty.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_48, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_medas__medas.png"),
                    title = "Differentially expressed proteins- MEDAS")
create_volcano_plot(diff_49, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_medas__frailty_",
                                                "medas.png"),
                    title = paste("Differentially expressed proteins - ",
                                  "Frailty:MEDAS"))

diff_47 <- diff_47[diff_47$adj.P.Val < 0.05, ]
nrow(diff_47)
# A: 1
# B: 0
rownames(diff_47)
# A: "C4Z2V8"

diff_48 <- diff_48[diff_48$adj.P.Val < 0.05, ]
nrow(diff_48)
# A: 0
# B: 0
rownames(diff_48)

diff_49 <- diff_49[diff_49$adj.P.Val < 0.05, ]
nrow(diff_49)
# A: 1
# B: 0
rownames(diff_49)
# A: "C4Z2V8"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_49),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_49_ann <- merge(diff_49, prot_ann, by = "row.names")
row.names(diff_49_ann) <- diff_49_ann$Row.names
diff_49_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_49_ann <- diff_49_ann[order(diff_49_ann$adj.P.Val,
                                 -abs(diff_49_ann$logFC)), ]

# Save sig_prot_no_imp_(max)lfq_ft_medas.csv
write.csv(diff_49_ann, file = paste0(output_path, "ft_medas.csv"),
          row.names = TRUE)

sel0 = which("C4Z2V8"==rownames(assay(se_no_imp_medas)))
df0 = data.frame(frailty=colData(se_no_imp_medas)[,c("frailty")], 
                 medas=colData(se_no_imp_medas)[,c("medas")], 
                 expression=assay(se_no_imp_medas)[sel0,])
ggplot(df0, aes(x = medas, y = expression, color = frailty)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", aes(linetype = frailty))

# ==================================
#  Frailty and energy
# ==================================
# Build model matrix considering frailty, energy and interaction
se_no_imp_energy <- se_no_imp[, !is.na(colData(se_no_imp)$energy)]
dim(se_no_imp_energy)
# A: 182 189
# B: 442 189
design = model.matrix(~ colData(se_no_imp_energy)$frailty*
                        colData(se_no_imp_energy)$energy)
colnames(design) = c("constant", "frailty", "energy", "frailty:energy")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_energy), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_50 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_energy)))
diff_51 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_energy)))
diff_52 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_energy)))

create_volcano_plot(diff_50, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_energy__frailty.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_51, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_energy__energy.png"),
                    title = "Differentially expressed proteins - Energy")
create_volcano_plot(diff_52, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_energy__frailty_",
                                                "energy.png"),
                    title = paste("Differentially expressed proteins -",
                                  "Frailty:energy"))

diff_50 <- diff_50[diff_50$adj.P.Val < 0.05, ]
nrow(diff_50)
# A: 0
# B: 0
rownames(diff_50)

diff_51 <- diff_51[diff_51$adj.P.Val < 0.05, ]
nrow(diff_51)
# A: 0
# B: 0
rownames(diff_51)

diff_52 <- diff_52[diff_52$adj.P.Val < 0.05, ]
nrow(diff_52)
# A: 0
# B: 0
rownames(diff_52)

# ==================================
#   Frailty and bmi
# ==================================
# Build model matrix considering frailty, bmi and interaction
se_no_imp_bmi <- se_no_imp[, !is.na(colData(se_no_imp)$bmi)]
dim(se_no_imp_bmi)
# A: 182 201
# B: 442 201
design = model.matrix(~ colData(se_no_imp_bmi)$frailty*
                        colData(se_no_imp_bmi)$bmi)
colnames(design) = c("constant", "frailty", "bmi", "frailty:bmi")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_bmi), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using BH
diff_53 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_bmi)))
diff_54 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_bmi)))
diff_55 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_bmi)))

create_volcano_plot(diff_53, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_bmi__frailty.png"),
                    title = "Differentially expressed proteins - Frailty")
create_volcano_plot(diff_54, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_bmi__bmi.png"),
                    title = "Differentially expressed proteins - BMI")
create_volcano_plot(diff_55, save_path = paste0(plot_path, "volcano_plot_",
                                                "frailty_bmi__frailty_bmi.png"),
                    title = paste("Differentially expressed proteins -",
                                  "Frailty:BMI"))

diff_53 <- diff_53[diff_53$adj.P.Val < 0.05, ]
nrow(diff_53)
# A: 0
# B: 4
rownames(diff_53)
# B: "P9WQH6" "E1WS50" "A6L048" "Q8WWA0"

diff_54 <- diff_54[diff_54$adj.P.Val < 0.05, ]
nrow(diff_54)
# A: 0
# B: 0
rownames(diff_54)

diff_55 <- diff_55[diff_55$adj.P.Val < 0.05, ]
nrow(diff_55)
# A: 0
# B: 3
rownames(diff_55)
# B: "E1WS50" "P9WQH6" "A6L048"

# Get annotation from Uniprot
prot_ann <- GetProteinAnnontate(rownames(diff_55),
                                columns = c("gene_primary", "organism_name",
                                            "protein_name", "cc_function",
                                            "keyword","go_p", "go_c", "go_f"))
colnames(prot_ann) <- c("gene_primary", "organism_name", "protein_name",
                        "cc_function", "keyword","go_p", "go_c", "go_f")
diff_55_ann <- merge(diff_55, prot_ann, by = "row.names")
row.names(diff_55_ann) <- diff_55_ann$Row.names
diff_55_ann$Row.names <- NULL

# Order by adj.P.Val and logFC
diff_55_ann <- diff_55_ann[order(diff_55_ann$adj.P.Val,
                                 -abs(diff_55_ann$logFC)), ]

# Save sig_prot_no_imp_(max)lfq_ft_bmi.csv
write.csv(diff_55_ann, file = paste0(output_path, "ft_bmi.csv"),
          row.names = TRUE)
