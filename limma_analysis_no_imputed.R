library(ggplot2)
library(limma)
library(dplyr)

# ==================================
# Limma analysis
# ==================================
# Running DE analysis with data from preprocessing no imputed: Norm MaxLFQ,
# log2transformed, NA in at least 1 condition < 0.5, no imputation
# ==================================
# Load filtered SummarizedExperiment
work_path <- getwd()
load(paste0(work_path, "/data/se_no_imp.RData"))

# ==================================
# 1. Frailty
# ==================================

# MODEL 1: FT vs NFT
##  Build model matrix considering only frailty
design = model.matrix(~ colData(se_no_imp)$frailty)
colnames(design) = c("constant", "frailty")
head(design)

## Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp), design)

## Estimated coefficients for each protein
head(coef(fit))

## Moderate t-test of the differential expression by means of Bayes' empirical
## moderation of the standard errors to a global value
fit1 = eBayes(fit)

## p-values adjusted using the Benjamini-Hochberg method
diff_1 <- topTable(fit1, coef=2, adjust="BH", sort.by="P",
                   number=nrow(assay(se_no_imp)))

create_volcano_plot(diff_1, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model1_frailty.png"),
                    title = "Model 1: Differentially expressed proteins - Frailty")

diff_1 <- diff_1[diff_1$adj.P.Val < 0.05, ]
nrow(diff_1) 
rownames(diff_1)
# 10 proteins (p.adj < 0.05): "P06702" "A6KYJ0" "A6L792" "A9KRZ4" "P05109"
# "Q8A9M2" "P94360" "A6KYJ4" "C4KZP0" "A9KJL3"

# Up vs down expressed
diff_1_up <- diff_1[diff_1$logFC > 0,]
nrow(diff_1_up) 
rownames(diff_1_up)
# 7: "P06702" "A6KYJ0" "A6L792" "P05109" "Q8A9M2" "A6KYJ4" "C4KZP0"
diff_1_down <- diff_1[diff_1$logFC < 0,]
nrow(diff_1_down) 
rownames(diff_1_down)
# 3: "A9KRZ4" "P94360" "A9KJL3"

sel0 = which("P06702"==rownames(assay(se_no_imp)))
df0 = data.frame(frailty=colData(se_no_imp)[,c("frailty")],
                 expression=assay(se_no_imp)[sel0,])
ggplot(df0,aes(x=frailty,y=expression)) + geom_boxplot()

# ==================================
# 2. Two predictor variables
# ==================================
#   2.1. Fraily and sex
# ==================================

# MODEL 2
# Build model matrix considering frailty and sex
design = model.matrix(~ colData(se_no_imp)$frailty + colData(se_no_imp)$sex)
colnames(design) = c("constant", "frailty", "sex")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 3 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_2 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp)))
diff_3 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp)))

create_volcano_plot(diff_2, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model2_frailty.png"),
                    title = "Model 2: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_3, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model2_sex.png"),
                    title = "Model 2: Differentially expressed proteins - Sex")

diff_2 <- diff_2[diff_2$adj.P.Val < 0.05, ]
nrow(diff_2) 
rownames(diff_2)
# 6: "P06702" "A6L792" "A6KYJ0" "A6KYJ4" "C4KZP0" "P05109"

diff_2_up <- diff_2[diff_2$logFC > 0,]
nrow(diff_2_up) 
rownames(diff_2_up)
# 6
diff_2_down <- diff_2[diff_2$logFC < 0,]
nrow(diff_2_down) 
rownames(diff_2_down)
# 0

diff_3 <- diff_3[diff_3$adj.P.Val < 0.05, ]
nrow(diff_3) 
rownames(diff_3)
# 0

# MODEL 3: Consider interaction between frailty and sex
# Build model matrix considering frailty, sex and interaction
design = model.matrix(~ colData(se_no_imp)$frailty*colData(se_no_imp)$sex)
colnames(design) = c("constant", "frailty", "sex", "frailty:sex")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_4 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp)))
diff_5 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp)))
diff_6 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp)))

create_volcano_plot(diff_4, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model3_frailty.png"),
                     title = "Model 3: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_5, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model3_sex.png"),
                    title = "Model 3: Differentially expressed proteins - Sex")
create_volcano_plot(diff_6, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model3_frailty_sex.png"),
                    title = "Model 3: Differentially expressed proteins - Frailty:Sex")

diff_4 <- diff_4[diff_4$adj.P.Val < 0.05, ]
nrow(diff_4) 
rownames(diff_4)
# 30: "A6KYJ0" "C4KZP0" "A6L792" "A9KJJ5" "Q5L923" "A6KYJ7" "P05109" "A6KYJ4"
# "P95544" "A9KRZ4" "Q5L9B6" "C4ZBD3" "Q5L8C5" "C4ZD46" "Q8A9M2" "A8YXK9"
# "B9E9L7" "Q8A1A2" "A6KYH8" "A9KKU0" "P55990" "C4ZI85" "C4ZBG1" "Q59199" 
# "P06702" "A6KYH1" "A6KYJ6" "C4ZF71" "Q18CF4" "C4Z0Q6

diff_4_up <- diff_4[diff_4$logFC > 0,]
nrow(diff_4_up) 
rownames(diff_4_up)
# 29: "A6KYJ0" "C4KZP0" "A6L792" "A9KJJ5" "Q5L923" "A6KYJ7" "P05109" "A6KYJ4"
# "P95544" "Q5L9B6" "C4ZBD3" "Q5L8C5" "C4ZD46" "Q8A9M2" "A8YXK9" "B9E9L7"
# "Q8A1A2" "A6KYH8" "A9KKU0" "P55990" "C4ZI85" "C4ZBG1" "Q59199" "P06702"
# "A6KYH1" "A6KYJ6" "C4ZF71" "Q18CF4" "C4Z0Q6"
diff_4_down <- diff_4[diff_4$logFC < 0,]
nrow(diff_4_down) 
rownames(diff_4_down)
# 1
# "A9KRZ4"

diff_5 <- diff_5[diff_5$adj.P.Val < 0.05, ]
nrow(diff_5) 
rownames(diff_5)
# 0

diff_6 <- diff_6[diff_6$adj.P.Val < 0.05, ]
nrow(diff_6) 
rownames(diff_6)
# 1: "C4Z0Q6"

diff_6_up <- diff_6[diff_6$logFC > 0,]
nrow(diff_6_up) 
rownames(diff_6_up)
# 0
diff_6_down <- diff_6[diff_6$logFC < 0,]
nrow(diff_6_down) 
rownames(diff_6_down)
# 1
# "C4Z0Q6"

sel0 = which("C4Z0Q6"==rownames(assay(se_no_imp)))
df0 = data.frame(frailty=colData(se_no_imp)[,c("frailty")], 
                 sex=colData(se_no_imp)[,c("sex")], 
                 expression=assay(se_no_imp)[sel0,])
df0$frailty_sex <- paste(df0$frailty, df0$sex, sep = "_")
ggplot(df0,aes(x=frailty_sex,y=expression)) + geom_boxplot()

# ==================================
#   2.2. Frailty and alcohol
# ==================================

# MODEL 4: 
# Build model matrix considering frailty and alcohol
se_no_imp_alcohol <- se_no_imp[, !is.na(colData(se_no_imp)$alcohol)]
dim(se_no_imp_alcohol)
design = model.matrix(~ colData(se_no_imp_alcohol)$frailty+
                        colData(se_no_imp_alcohol)$alcohol)
colnames(design) = c("constant", "frailty", "alcohol_monthly", "alcohol_weekly")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_alcohol), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_7 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp_alcohol)))
diff_8 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp_alcohol)))
diff_9 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp_alcohol)))

create_volcano_plot(diff_7, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model4_frailty.png"),
                    title = "Model 4: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_8, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model4_alcohol_monthly.png"),
                    title = "Model 4: Differentially expressed proteins - Alcohol monthly")
create_volcano_plot(diff_9, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model4_alcohol_weekly.png"),
                    title = "Model 4: Differentially expressed proteins - Alcohol weekly")


diff_7 <- diff_7[diff_7$adj.P.Val < 0.05, ]
nrow(diff_7) 
rownames(diff_7)
# 0

diff_8 <- diff_8[diff_8$adj.P.Val < 0.05, ]
nrow(diff_8) 
rownames(diff_8)
# 0

diff_9 <- diff_9[diff_9$adj.P.Val < 0.05, ]
nrow(diff_9) 
rownames(diff_9)
# 1: "P55259"

sel0 = which("P55259"==rownames(assay(se_no_imp)))
df0 = data.frame(alcohol=colData(se_no_imp)[,c("alcohol")],
                 expression=assay(se_no_imp)[sel0,])
ggplot(df0,aes(x=alcohol,y=expression)) + geom_boxplot()

diff_9_up <- diff_9[diff_9$logFC > 0,]
nrow(diff_9_up) 
rownames(diff_9_up)
# 1
# "P55259"
diff_9_down <- diff_9[diff_9$logFC < 0,]
nrow(diff_9_down) 
rownames(diff_9_down)
# 0

# MODEL 5: Consider interaction between frailty and alcohol
# Build model matrix considering frailty, alcohol and interaction
design = model.matrix(~ colData(se_no_imp_alcohol)$frailty*
                        colData(se_no_imp_alcohol)$alcohol)
colnames(design) = c("constant", "frailty", "alcohol_monthly", "alcohol_weekly",
                     "frailty:alcohol_monthly", "frailty:alcohol_weekly")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_alcohol), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 6 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_10 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_alcohol)))
diff_11 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_alcohol)))
diff_12 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_alcohol)))
diff_13 <- topTable(fit1,coef=5, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_alcohol)))
diff_14 <- topTable(fit1,coef=6, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_alcohol)))

create_volcano_plot(diff_10, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model5_frailty.png"),
                    title = "Model 5: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_11, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model5_alcohol_monthly.png"),
                    title = "Model 5: Differentially expressed proteins - Alcohol monthly")
create_volcano_plot(diff_12, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model5_alcohol_weekly.png"),
                    title = "Model 5: Differentially expressed proteins - Alcohol weekly")
create_volcano_plot(diff_13, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model5_frailty_alcohol_monthly.png"),
                    title = "Model 5: Differentially expressed proteins - Frailty:Alcohol monthly")
create_volcano_plot(diff_14, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model5_frailty_alcohol_weekly.png"),
                    title = "Model 5: Differentially expressed proteins - Frailty:Alcohol weekly")

diff_10 <- diff_10[diff_10$adj.P.Val < 0.05, ]
nrow(diff_10) 
rownames(diff_10)
# 0

diff_11 <- diff_11[diff_11$adj.P.Val < 0.05, ]
nrow(diff_11) 
rownames(diff_11)
# 0

diff_12 <- diff_12[diff_12$adj.P.Val < 0.05, ]
nrow(diff_12) 
rownames(diff_12)
# 1: "P55259"

diff_12_up <- diff_12[diff_12$logFC > 0,]
nrow(diff_12_up) 
rownames(diff_12_up)
# 1
# "P55259"
diff_12_down <- diff_12[diff_12$logFC < 0,]
nrow(diff_12_down) 
rownames(diff_12_down)
# 0

diff_13 <- diff_13[diff_13$adj.P.Val < 0.05, ]
nrow(diff_13) 
rownames(diff_13)
# 0

diff_14 <- diff_14[diff_14$adj.P.Val < 0.05, ]
nrow(diff_14) 
rownames(diff_14)
# 1: "A6TH53"

diff_14_up <- diff_14[diff_14$logFC > 0,]
nrow(diff_14_up) 
rownames(diff_14_up)
# 1
# "A6TH53"
diff_14_down <- diff_14[diff_14$logFC < 0,]
nrow(diff_14_down) 
rownames(diff_14_down)
# 0

sel0 = which("A6TH53"==rownames(assay(se_no_imp_alcohol)))
df0 = data.frame(frailty=colData(se_no_imp_alcohol)[,c("frailty")], 
                 alcohol=colData(se_no_imp_alcohol)[,c("alcohol")], 
                 expression=assay(se_no_imp_alcohol)[sel0,])
df0$frailty_alcohol <- paste(df0$frailty, df0$alcohol, sep = "_")
ggplot(df0,aes(x=frailty_alcohol,y=expression)) + geom_boxplot()


# ==================================
#   2.3. Frailty and tobacco
# ==================================

# MODEL 6
# Build model matrix considering frailty and tobacco
se_no_imp_tobacco <- se_no_imp[, !is.na(colData(se_no_imp)$tobacco)]
dim(se_no_imp_tobacco)
design = model.matrix(~ colData(se_no_imp_tobacco)$frailty+
                        colData(se_no_imp_tobacco)$tobacco)
colnames(design) = c("constant", "frailty", "former", "current")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_tobacco), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_15 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_tobacco)))
diff_16 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_tobacco)))
diff_17 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_tobacco)))

create_volcano_plot(diff_15, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model6_frailty.png"),
                    title = "Model 6: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_16, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model6_tabacco_former.png"),
                    title = "Model 6: Differentially expressed proteins - Tabacco former")
create_volcano_plot(diff_17, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model6_tabacco_current.png"),
                    title = "Model 6: Differentially expressed proteins - Tabacco current")

diff_15 <- diff_15[diff_15$adj.P.Val < 0.05, ]
nrow(diff_15) 
rownames(diff_15)
# 6: "P06702" "A6KYJ0" "A9KRZ4" "P94360" "A6L792" "P05109"

diff_15_up <- diff_15[diff_15$logFC > 0,]
nrow(diff_15_up)
rownames(diff_15_up)
# 4: "P06702" "A6KYJ0" "A6L792" "P05109"
diff_15_down <- diff_15[diff_15$logFC < 0,]
nrow(diff_15_down)
rownames(diff_15_down)
# 2: "A9KRZ4" "P94360"

diff_16 <- diff_16[diff_16$adj.P.Val < 0.05, ]
nrow(diff_16) 
rownames(diff_16)
# 0

diff_17 <- diff_17[diff_17$adj.P.Val < 0.05, ]
nrow(diff_17) 
rownames(diff_17)
# 1: "C4ZB90"

diff_17_up <- diff_17[diff_17$logFC > 0,]
nrow(diff_17_up) 
rownames(diff_17_up)
# 1
# "C4ZB90"
diff_17_down <- diff_17[diff_17$logFC < 0,]
nrow(diff_17_down) 
rownames(diff_17_down)
# 0

sel0 = which("C4ZB90"==rownames(assay(se_no_imp)))
df0 = data.frame(tobacco=colData(se_no_imp)[,c("tobacco")],
                 expression=assay(se_no_imp)[sel0,])
ggplot(df0,aes(x=tobacco,y=expression)) + geom_boxplot()

# MODEL 7: Consider interaction between frailty and tobacco
# Build model matrix considering frailty, tobacco and interaction
design = model.matrix(~ colData(se_no_imp_tobacco)$frailty*
                        colData(se_no_imp_tobacco)$tobacco)
colnames(design) = c("constant", "frailty", "former", "current",
                     "frailty:former", "frailty:current")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_tobacco), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 6 coefficients per each column of the model matrix. Adjusted using the Benjamini-Hochberg
diff_18 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_tobacco)))
diff_19 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_tobacco)))
diff_20 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_tobacco)))
diff_21 <- topTable(fit1,coef=5, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_tobacco)))
diff_22 <- topTable(fit1,coef=6, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_tobacco)))

create_volcano_plot(diff_18, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model7_frailty.png"),
                    title = "Model 7: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_19, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model7_tabacco_former.png"),
                    title = "Model 7: Differentially expressed proteins - Tabacco former")
create_volcano_plot(diff_20, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model7_tabacco_current.png"),
                    title = "Model 7: Differentially expressed proteins - Tabacco current")
create_volcano_plot(diff_21, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model7_frailty_tabacco_former.png"),
                    title = "Model 7: Differentially expressed proteins - Frailty:Tabacco former")
create_volcano_plot(diff_22, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model7_frailty_tabacco_current.png"),
                    title = "Model 7: Differentially expressed proteins - Frailty:Tabacco current")

diff_18 <- diff_18[diff_18$adj.P.Val < 0.05, ]
nrow(diff_18)
rownames(diff_18)
# 0

diff_19 <- diff_19[diff_19$adj.P.Val < 0.05, ]
nrow(diff_19)
rownames(diff_19)
# 0

diff_20 <- diff_20[diff_20$adj.P.Val < 0.05, ]
nrow(diff_20)
rownames(diff_20)
# 0

diff_21 <- diff_21[diff_21$adj.P.Val < 0.05, ]
nrow(diff_21)
rownames(diff_21)
# 0

diff_22 <- diff_22[diff_22$adj.P.Val < 0.05, ]
nrow(diff_22)
rownames(diff_22)
# 79: "A6KYH6" "P0C2E7" "C4ZBT7" "C4ZBL1" "A6KYK2" "C4ZB90" "A6L1X1" "A6KYJ7"
# "A9VT65" "C4Z2R7" "A9KJL3" "A9KRZ2" "A6LFQ4" "P22983" "A6KYH8" "C4Z2R9"
# "C4ZBU3" "C4ZD46" "A8YXK9" "C4ZBS5" "C4ZBD3" "C4Z2T1" "C4Z3R4" "A9KK92"
# "Q59309" "A6KYK6" "A9KK94" "A6TH53" "A6KYI1" "C4ZBG1" "Q8A9M2" "P0C0G6"
# "A6KXA0" "C4ZBT1" "A6KYH0" "P94316" "A6KYJ8" "Q8RQP4" "A9KJI8" "C4Z2T8"
# "C4ZBS1" "A9KKU0" "Q1Q899" "C4Z2V8" "A6KYJ2" "Q88XX2" "A6L048" "Q8A1A2"
# "B9E9L7" "A6L1L8" "C0R090" "C4ZBD5" "A9KRZ3" "A6KXL2" "A6KYK3" "B8I7Y6"
# "C4Z1J4" "A6L792" "P95544" "C4KZP0" "A6TWI7" "Q5L8C5" "A9KJI0" "A6KYI8"
# "P19543" "A6KYK9" "A6KYJ6" "Q5L9B6" "A6L0V1" "A6L903" "P02768" "Q5L923"
# "C4ZBS3" "Q18CF4" "A4XI37" "C4ZF71" "A9KJJ5" "A5I7J4" "P24295"

diff_22_up <- diff_22[diff_22$logFC > 0,]
nrow(diff_22_up)
rownames(diff_22_up)
# 79
diff_22_down <- diff_22[diff_22$logFC < 0,]
nrow(diff_22_up)
rownames(diff_22_down)
# 0

sel0 = which("A6KYH6"==rownames(assay(se_no_imp_tobacco)))
df0 = data.frame(frailty=colData(se_no_imp_tobacco)[,c("frailty")], 
                 tobacco=colData(se_no_imp_tobacco)[,c("tobacco")], 
                 expression=assay(se_no_imp_tobacco)[sel0,])
df0$frailty_tobacco <- paste(df0$frailty, df0$tobacco, sep = "_")
ggplot(df0,aes(x=frailty_tobacco,y=expression)) + geom_boxplot()

# ==================================
#   2.4. Frailty and diabetes
# ==================================
# MODEL 8
# Build model matrix considering frailty and diabetes
se_no_imp_diabetes <- se_no_imp[, !is.na(colData(se_no_imp)$diabetes)]
dim(se_no_imp_diabetes)
design = model.matrix(~ colData(se_no_imp_diabetes)$frailty+
                        colData(se_no_imp_diabetes)$diabetes)
colnames(design) = c("constant", "frailty", "diabetes")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_diabetes), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 3 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_23 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_diabetes)))
diff_24 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_diabetes)))

create_volcano_plot(diff_23, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model8_frailty.png"),
                    title = "Model 8: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_24, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model8_diabetes.png"),
                    title = "Model 8: Differentially expressed proteins - Diabetes")

diff_23 <- diff_23[diff_23$adj.P.Val < 0.05, ]
nrow(diff_23) 
rownames(diff_23)
# 5: "A6KYJ0" "P06702" "A6L792" "C4KZP0" "A9KRZ4"

diff_23_up <- diff_23[diff_23$logFC > 0,]
nrow(diff_23_up) 
rownames(diff_23_up)
# 4
# "A6KYJ0" "P06702" "A6L792" "C4KZP0"
diff_23_down <- diff_23[diff_23$logFC < 0,]
nrow(diff_23_down) 
rownames(diff_23_down)
# 1
# "A9KRZ4"

diff_24 <- diff_24[diff_24$adj.P.Val < 0.05, ]
nrow(diff_24) 
rownames(diff_24)
# 0

# MODEL 9
# Build model matrix considering frailty, diabetes and interaction
design = model.matrix(~ colData(se_no_imp_diabetes)$frailty*
                        colData(se_no_imp_diabetes)$diabetes)
colnames(design) = c("constant", "frailty", "diabetes", "frailty:diabetes")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_diabetes), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_25 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp_diabetes)))
diff_26 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp_diabetes)))
diff_27 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                   number=nrow(assay(se_no_imp_diabetes)))

create_volcano_plot(diff_25, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model9_frailty.png"),
                    title = "Model 9: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_26, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model9_diabetes.png"),
                    title = "Model 9: Differentially expressed proteins - Diabetes")
create_volcano_plot(diff_27, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model9_frailty_diabetes.png"),
                    title = "Model 9: Differentially expressed proteins - Frailty:Diabetes")

diff_25 <- diff_25[diff_25$adj.P.Val < 0.05, ]
nrow(diff_25) 
rownames(diff_25)
# 1: "A6KYJ0"

diff_25_up <- diff_25[diff_25$logFC > 0,]
nrow(diff_25_up)
rownames(diff_25_up)
# 1
# "A6KYJ0"
diff_25_down <- diff_25[diff_25$logFC < 0,]
nrow(diff_25_down) 
rownames(diff_25_down)
# 0

diff_26 <- diff_26[diff_26$adj.P.Val < 0.05, ]
nrow(diff_26) 
rownames(diff_26)
# 0

diff_27 <- diff_27[diff_27$adj.P.Val < 0.05, ]
nrow(diff_27)
rownames(diff_27)
# 0

# ==================================
#   2.5. Frailty and chf
# ==================================
# MODEL 10
# Build model matrix considering frailty and chf
se_no_imp_chf <- se_no_imp[, !is.na(colData(se_no_imp)$chf)]
dim(se_no_imp_chf)
design = model.matrix(~ colData(se_no_imp_chf)$frailty+
                        colData(se_no_imp_chf)$chf)
colnames(design) = c("constant", "frailty", "chf")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_chf), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 3 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_28 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_chf)))
diff_29 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_chf)))

create_volcano_plot(diff_28, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model10_frailty.png"),
                    title = "Model 10: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_29, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model10_chf.png"),
                    title = "Model 10: Differentially expressed proteins - Chf")

diff_28 <- diff_28[diff_28$adj.P.Val < 0.05, ]
nrow(diff_28) 
rownames(diff_28)
# 7: "P06702" "A6KYJ0" "A6L792" "A9KRZ4" "P05109" "Q8A9M2" "P94360"

diff_28_up <- diff_28[diff_28$logFC > 0,]
nrow(diff_28_up) 
rownames(diff_28_up)
# 5
# "P06702" "A6KYJ0" "A6L792" "P05109" "Q8A9M2"
diff_28_down <- diff_28[diff_28$logFC < 0,]
nrow(diff_28_down) 
rownames(diff_28_down)
# 2
# "A9KRZ4" "P94360"

diff_29 <- diff_29[diff_29$adj.P.Val < 0.05, ]
nrow(diff_29) 
rownames(diff_29)
# 0

# MODEL 11
# Build model matrix considering frailty, chf and interaction
design = model.matrix(~ colData(se_no_imp_chf)$frailty*
                        colData(se_no_imp_chf)$chf)
colnames(design) = c("constant", "frailty", "chf", "frailty:chf")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_chf), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_30 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_chf)))
diff_31 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_chf)))
diff_32 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_chf)))

create_volcano_plot(diff_30, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model11_frailty.png"),
                    title = "Model 11: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_31, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model11_chf.png"),
                    title = "Model 11: Differentially expressed proteins - Chf")
create_volcano_plot(diff_32, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model11_frailty_chf.png"),
                    title = "Model 11: Differentially expressed proteins - Frailty:Chf")

diff_30 <- diff_30[diff_30$adj.P.Val < 0.05, ]
nrow(diff_30) 
rownames(diff_30)
# 5: "P06702" "A6KYJ0" "P05109" "A9KJL3" "A6L792"

diff_30_up <- diff_30[diff_30$logFC > 0,]
nrow(diff_30_up)
rownames(diff_30_up)
# 4
# "P06702" "A6KYJ0" "P05109" "A6L792"
diff_30_down <- diff_30[diff_30$logFC < 0,]
nrow(diff_30_down) 
rownames(diff_30_down)
# 1
# "A9KJL3"

diff_31 <- diff_31[diff_31$adj.P.Val < 0.05, ]
nrow(diff_31) 
rownames(diff_31)
# 0

diff_32 <- diff_32[diff_32$adj.P.Val < 0.05, ]
nrow(diff_32)
rownames(diff_32)
# 0

# ==================================
#   2.6. Frailty and depression
# ==================================
# MODEL 12
# Build model matrix considering frailty and depression
se_no_imp_depression <- se_no_imp[, !is.na(colData(se_no_imp)$depression)]
dim(se_no_imp_depression)
design = model.matrix(~ colData(se_no_imp_depression)$frailty+
                        colData(se_no_imp_depression)$depression)
colnames(design) = c("constant", "frailty", "depression")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_depression), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 3 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_33 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_depression)))
diff_34 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_depression)))

create_volcano_plot(diff_33, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model12_frailty.png"),
                    title = "Model 12: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_34, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model12_depression.png"),
                    title = "Model 12: Differentially expressed proteins - Depression")

diff_33 <- diff_33[diff_33$adj.P.Val < 0.05, ]
nrow(diff_33) 
rownames(diff_33)
# 5: "A6KYJ0" "P06702" "A9KRZ4" "A6L792" "P05109"

diff_33_up <- diff_33[diff_33$logFC > 0,]
nrow(diff_33_up) 
rownames(diff_33_up)
# 4
# "A6KYJ0" "P06702" "A6L792" "P05109"
diff_33_down <- diff_33[diff_33$logFC < 0,]
nrow(diff_33_down) 
rownames(diff_33_down)
# 1
# "A9KRZ4"

diff_34 <- diff_34[diff_34$adj.P.Val < 0.05, ]
nrow(diff_34) 
rownames(diff_34)
# 0

# MODEL 13
# Build model matrix considering frailty, depression and interaction
design = model.matrix(~ colData(se_no_imp_depression)$frailty*
                        colData(se_no_imp_depression)$depression)
colnames(design) = c("constant", "frailty", "depression", "frailty:depression")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_depression), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_35 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_depression)))
diff_36 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_depression)))
diff_37 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_depression)))

create_volcano_plot(diff_35, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model13_frailty.png"),
                    title = "Model 13: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_36, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model13_depression.png"),
                    title = "Model 13: Differentially expressed proteins - Depression")
create_volcano_plot(diff_37, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model13_frailty_depression.png"),
                    title = "Model 13: Differentially expressed proteins - Frailty:Depression")

diff_35 <- diff_35[diff_35$adj.P.Val < 0.05, ]
nrow(diff_35) 
rownames(diff_35)
# 5: "A6KYJ0" "P06702" "A9KRZ4" "A6L792" "P05109"

diff_35_up <- diff_35[diff_35$logFC > 0,]
nrow(diff_35_up)
rownames(diff_35_up)
# 4
# A6KYJ0" "P06702" "A6L792" "P05109"
diff_35_down <- diff_35[diff_35$logFC < 0,]
nrow(diff_35_down) 
rownames(diff_35_down)
# 1
# "A9KRZ4"

diff_36 <- diff_36[diff_36$adj.P.Val < 0.05, ]
nrow(diff_36) 
rownames(diff_36)
# 0

diff_37 <- diff_37[diff_37$adj.P.Val < 0.05, ]
nrow(diff_37)
rownames(diff_37)
# 0

# ==================================
#   2.7. Frailty and osteoarthritis
# ==================================
# MODEL 14
# Build model matrix considering frailty and osteoarthritis
se_no_imp_osteo <- se_no_imp[, !is.na(colData(se_no_imp)$osteoarthritis)]
dim(se_no_imp_osteo)
design = model.matrix(~ colData(se_no_imp_osteo)$frailty+
                        colData(se_no_imp_osteo)$osteoarthritis)
colnames(design) = c("constant", "frailty", "osteoarthritis")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_osteo), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 3 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_38 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_osteo)))
diff_39 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_osteo)))

create_volcano_plot(diff_38, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model14_frailty.png"),
                    title = "Model 14: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_39, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model14_osteoarthritis.png"),
                    title = "Model 14: Differentially expressed proteins - Osteoarthritis")

diff_38 <- diff_38[diff_38$adj.P.Val < 0.05, ]
nrow(diff_38) 
rownames(diff_38)
# 9: "P06702" "A6KYJ4" "A6KYJ0" "A6L792" "P05109" "A9KRZ4" "A6L048" "A6L903"
# "Q5LHW2"

diff_38_up <- diff_38[diff_38$logFC > 0,]
nrow(diff_38_up) 
rownames(diff_38_up)
# 8
# "P06702" "A6KYJ4" "A6KYJ0" "A6L792" "P05109" "A6L048" "A6L903" "Q5LHW2"
diff_38_down <- diff_38[diff_38$logFC < 0,]
nrow(diff_38_down) 
rownames(diff_38_down)
# 1
# "A9KRZ4"

diff_39 <- diff_39[diff_39$adj.P.Val < 0.05, ]
nrow(diff_39) 
rownames(diff_39)
# 0

# MODEL 15
# Build model matrix considering frailty, osteoarthritis and interaction
design = model.matrix(~ colData(se_no_imp_osteo)$frailty*
                        colData(se_no_imp_osteo)$osteoarthritis)
colnames(design) = c("constant", "frailty", "osteoarthritis", "frailty:osteoarthritis")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_osteo), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_40 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_osteo)))
diff_41 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_osteo)))
diff_42 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_osteo)))

create_volcano_plot(diff_40, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model15_frailty.png"),
                    title = "Model 15: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_41, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model15_osteoarthritis.png"),
                    title = "Model 15: Differentially expressed proteins - Osteoarthritis")
create_volcano_plot(diff_42, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model15_frailty_osteoarthritis.png"),
                    title = "Model 15: Differentially expressed proteins - Frailty:Osteoarthritis")

diff_40 <- diff_40[diff_40$adj.P.Val < 0.05, ]
nrow(diff_40) 
rownames(diff_40)
# 13: "A6KYJ4" "A9KJJ5" "P06702" "A6L792" "P05109" "A6L903" "Q5L8C5" "Q9AE24"
# "Q5L923" "A9KRZ4" "Q5LHW2" "A6L048" "Q8A9M2"

diff_40_up <- diff_40[diff_40$logFC > 0,]
nrow(diff_40_up)
rownames(diff_40_up)
# 12
# "A6KYJ4" "A9KJJ5" "P06702" "A6L792" "P05109" "A6L903" "Q5L8C5" "Q9AE24"
# "Q5L923" "Q5LHW2" "A6L048" "Q8A9M2"
diff_40_down <- diff_40[diff_40$logFC < 0,]
nrow(diff_40_down) 
rownames(diff_40_down)
# 1
# "A9KRZ4"

diff_41 <- diff_41[diff_41$adj.P.Val < 0.05, ]
nrow(diff_41) 
rownames(diff_41)
# 0

diff_42 <- diff_42[diff_42$adj.P.Val < 0.05, ]
nrow(diff_42)
rownames(diff_42)
# 0

# ==================================
#   2.8. Frailty and sarcopenia
# ==================================
# MODEL 16
# Build model matrix considering frailty and sarcopenia
se_no_imp_sarco <- se_no_imp[, !is.na(colData(se_no_imp)$sarcopenia)]
dim(se_no_imp_sarco)
design = model.matrix(~ colData(se_no_imp_sarco)$frailty+
                        colData(se_no_imp_sarco)$sarcopenia)
colnames(design) = c("constant", "frailty", "sarcopenia")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_sarco), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 3 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_43 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_sarco)))
diff_44 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_sarco)))

create_volcano_plot(diff_43, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model16_frailty.png"),
                    title = "Model 16: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_44, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model16_sarcopenia.png"),
                    title = "Model 16: Differentially expressed proteins - Sarcopenia")

diff_43 <- diff_43[diff_43$adj.P.Val < 0.05, ]
nrow(diff_43) 
rownames(diff_43)
# 8: "A6KYJ0" "P06702" "A9KJL3" "O83023" "P94360" "A6L792" "C4Z2R3" "P05109"

diff_43_up <- diff_43[diff_43$logFC > 0,]
nrow(diff_43_up) 
rownames(diff_43_up)
# 4
# "A6KYJ0" "P06702" "A6L792" "P05109"
diff_43_down <- diff_43[diff_43$logFC < 0,]
nrow(diff_43_down) 
rownames(diff_43_down)
# 4
# "A9KJL3" "O83023" "P94360" "C4Z2R3"

diff_44 <- diff_44[diff_44$adj.P.Val < 0.05, ]
nrow(diff_44) 
rownames(diff_44)
# 0

# MODEL 17
# Build model matrix considering frailty, sarcopenia and interaction
design = model.matrix(~ colData(se_no_imp_sarco)$frailty*
                        colData(se_no_imp_sarco)$sarcopenia)
colnames(design) = c("constant", "frailty", "sarcopenia", "frailty:sarcopenia")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_sarco), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_45 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_sarco)))
diff_46 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_sarco)))
diff_47 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_sarco)))

create_volcano_plot(diff_45, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model17_frailty.png"),
                    title = "Model 17: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_46, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model17_sarcopenia.png"),
                    title = "Model 17: Differentially expressed proteins - Sarcopenia")
create_volcano_plot(diff_47, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model17_frailty_sarcopenia.png"),
                    title = "Model 17: Differentially expressed proteins - Frailty:Sarcopenia")

diff_45 <- diff_45[diff_45$adj.P.Val < 0.05, ]
nrow(diff_45) 
rownames(diff_45)
# 8
# "P06702" "A6KYJ0" "P05109" "A9KJL3" "O83023" "P94360" "A6L792" "C4Z2R3"

diff_45_up <- diff_45[diff_45$logFC > 0,]
nrow(diff_45_up)
rownames(diff_45_up)
# 4
# "P06702" "A6KYJ0" "P05109" "A6L792"
diff_45_down <- diff_45[diff_45$logFC < 0,]
nrow(diff_45_down) 
rownames(diff_45_down)
# 4
# "A9KJL3" "O83023" "P94360" "C4Z2R3"

diff_46 <- diff_46[diff_46$adj.P.Val < 0.05, ]
nrow(diff_46) 
rownames(diff_46)
# 0

diff_47 <- diff_47[diff_47$adj.P.Val < 0.05, ]
nrow(diff_47)
rownames(diff_47)
# 1: "A6KYK7"

diff_47_up <- diff_47[diff_47$logFC > 0,]
nrow(diff_47_up)
rownames(diff_47_up)
# 0

diff_47_down <- diff_47[diff_47$logFC < 0,]
nrow(diff_47_down) 
rownames(diff_47_down)
# 1
# "A6KYK7"

sel0 = which("A6KYK7"==rownames(assay(se_no_imp_sarco)))
df0 = data.frame(frailty=colData(se_no_imp_sarco)[,c("frailty")], 
                 sarco=colData(se_no_imp_sarco)[,c("sarcopenia")], 
                 expression=assay(se_no_imp_sarco)[sel0,])
df0$frailty_sarco <- paste(df0$frailty, df0$sarco, sep = "_")
ggplot(df0,aes(x=frailty_sarco,y=expression)) + geom_boxplot()

# ==================================
#   2.9. Frailty and bmi
# ==================================
# MODEL 18
# Build model matrix considering frailty and bmi
se_no_imp_bmi <- se_no_imp[, !is.na(colData(se_no_imp)$bmi)]
dim(se_no_imp_bmi)
design = model.matrix(~ colData(se_no_imp_bmi)$frailty+
                        colData(se_no_imp_bmi)$bmi)
colnames(design) = c("constant", "frailty", "bmi")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_bmi), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 3 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_48 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_bmi)))
diff_49 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_bmi)))

create_volcano_plot(diff_48, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model18_frailty.png"),
                    title = "Model 18: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_49, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model18_bmi.png"),
                    title = "Model 18: Differentially expressed proteins - BMI")

diff_48 <- diff_48[diff_48$adj.P.Val < 0.05, ]
nrow(diff_48) 
rownames(diff_48)
# 8: "P06702" "A6KYJ0" "C4KZP0" "A9KRZ4" "A6L792" "A6KYJ4" "P94360" "P05109"

diff_48_up <- diff_48[diff_48$logFC > 0,]
nrow(diff_48_up) 
rownames(diff_48_up)
# 6
# "P06702" "A6KYJ0" "C4KZP0" "A6L792" "A6KYJ4" "P05109"
diff_48_down <- diff_48[diff_48$logFC < 0,]
nrow(diff_48_down) 
rownames(diff_48_down)
# 2
# "A9KRZ4" "P94360"

diff_49 <- diff_49[diff_49$adj.P.Val < 0.05, ]
nrow(diff_49) 
rownames(diff_49)
# 0

# MODEL 19
# Build model matrix considering frailty, bmi and interaction
design = model.matrix(~ colData(se_no_imp_bmi)$frailty*
                        colData(se_no_imp_bmi)$bmi)
colnames(design) = c("constant", "frailty", "bmi", "frailty:bmi")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_bmi), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_50 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_bmi)))
diff_51 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_bmi)))
diff_52 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_bmi)))

create_volcano_plot(diff_50, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model19_frailty.png"),
                    title = "Model 19: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_51, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model19_bmi.png"),
                    title = "Model 19: Differentially expressed proteins - BMI")
create_volcano_plot(diff_52, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model19_frailty_bmi.png"),
                    title = "Model 19: Differentially expressed proteins - Frailty:BMI")

diff_50 <- diff_50[diff_50$adj.P.Val < 0.05, ]
nrow(diff_50) 
rownames(diff_50)
# 0

diff_51 <- diff_51[diff_51$adj.P.Val < 0.05, ]
nrow(diff_51) 
rownames(diff_51)
# 0

diff_52 <- diff_52[diff_52$adj.P.Val < 0.05, ]
nrow(diff_52)
rownames(diff_52)
# 0

# ==================================
#   2.10. Frailty and energy
# ==================================
# MODEL 20
# Build model matrix considering frailty and energy
se_no_imp_energy <- se_no_imp[, !is.na(colData(se_no_imp)$energy)]
dim(se_no_imp_energy)
design = model.matrix(~ colData(se_no_imp_energy)$frailty+
                        colData(se_no_imp_energy)$energy)
colnames(design) = c("constant", "frailty", "energy")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_energy), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 3 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_53 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_energy)))
diff_54 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_energy)))

create_volcano_plot(diff_53, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model20_frailty.png"),
                    title = "Model 20: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_54, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model20_energy.png"),
                    title = "Model 20: Differentially expressed proteins - Energy")

diff_53 <- diff_53[diff_53$adj.P.Val < 0.05, ]
nrow(diff_53) 
rownames(diff_53)
# 0

diff_54 <- diff_54[diff_54$adj.P.Val < 0.05, ]
nrow(diff_54) 
rownames(diff_54)
# 0

# MODEL 21
# Build model matrix considering frailty, energy and interaction
design = model.matrix(~ colData(se_no_imp_energy)$frailty*
                        colData(se_no_imp_energy)$energy)
colnames(design) = c("constant", "frailty", "energy", "frailty:energy")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_energy), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_55 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_energy)))
diff_56 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_energy)))
diff_57 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_energy)))

create_volcano_plot(diff_55, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model21_frailty.png"),
                    title = "Model 21: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_56, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model21_energy.png"),
                    title = "Model 21: Differentially expressed proteins - Energy")
create_volcano_plot(diff_57, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model21_frailty_energy.png"),
                    title = "Model 21: Differentially expressed proteins - Frailty:energy")

diff_55 <- diff_55[diff_55$adj.P.Val < 0.05, ]
nrow(diff_55) 
rownames(diff_55)
# 0

diff_56 <- diff_56[diff_56$adj.P.Val < 0.05, ]
nrow(diff_56) 
rownames(diff_56)
# 0

diff_57 <- diff_57[diff_57$adj.P.Val < 0.05, ]
nrow(diff_57)
rownames(diff_57)
# 0

# ==================================
#   2.11. Frailty and ilef
# ==================================
# MODEL 22
# Build model matrix considering frailty and ilef
se_no_imp_ilef <- se_no_imp[, !is.na(colData(se_no_imp)$ilef)]
dim(se_no_imp_ilef)
design = model.matrix(~ colData(se_no_imp_ilef)$frailty+
                        colData(se_no_imp_ilef)$ilef)
colnames(design) = c("constant", "frailty", "ilef")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_ilef), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 3 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_58 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_ilef)))
diff_59 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_ilef)))

create_volcano_plot(diff_58, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model22_frailty.png"),
                    title = "Model 22: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_59, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model22_ilef.png"),
                    title = "Model 22: Differentially expressed proteins - ILEF")

diff_58 <- diff_58[diff_58$adj.P.Val < 0.05, ]
nrow(diff_58) 
rownames(diff_58)
# 1: "A6L792"

diff_58_up <- diff_58[diff_58$logFC > 0,]
nrow(diff_58_up) 
rownames(diff_58_up)
# 1
# "A6L792"
diff_58_down <- diff_58[diff_58$logFC < 0,]
nrow(diff_58_down) 
rownames(diff_58_down)
# 0

diff_59 <- diff_59[diff_59$adj.P.Val < 0.05, ]
nrow(diff_59) 
rownames(diff_59)
# 0

# MODEL 23
# Build model matrix considering frailty, ilef and interaction
design = model.matrix(~ colData(se_no_imp_ilef)$frailty*
                        colData(se_no_imp_ilef)$ilef)
colnames(design) = c("constant", "frailty", "ilef", "frailty:ilef")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_ilef), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of model matrix. Adjusted using Benjamini-Hochberg
diff_60 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_ilef)))
diff_61 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_ilef)))
diff_62 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P",
                    number=nrow(assay(se_no_imp_ilef)))

create_volcano_plot(diff_60, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model23_frailty.png"),
                    title = "Model 23: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_61, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model23_ilef.png"),
                    title = "Model 23: Differentially expressed proteins - ILEF")
create_volcano_plot(diff_62, save_path = paste0(work_path, "/plots/limma_analysis/no_imputation/volcano_plot_model23_frailty_ilef.png"),
                    title = "Model 23: Differentially expressed proteins - Frailty:ILEF")

diff_60 <- diff_60[diff_60$adj.P.Val < 0.05, ]
nrow(diff_60) 
rownames(diff_60)
# 2: "A6KYJ0" "A6L792"

diff_60_up <- diff_60[diff_60$logFC > 0,]
nrow(diff_60_up) 
rownames(diff_60_up)
# 2
# "A6KYJ0" "A6L792"
diff_60_down <- diff_60[diff_60$logFC < 0,]
nrow(diff_60_down) 
rownames(diff_60_down)
# 0

diff_61 <- diff_61[diff_61$adj.P.Val < 0.05, ]
nrow(diff_61) 
rownames(diff_61)
# 0

diff_62 <- diff_62[diff_62$adj.P.Val < 0.05, ]
nrow(diff_62)
rownames(diff_62)
# 0

# ==================================
# 3. All predictor variables significantly related to frailty
# ==================================
# FINAL MODEL
# Keep cases with data for all predictor variables
sig_var <- c("alcohol", "tobacco", "diabetes", "chf", "depression",
             "osteoarthritis", "sarcopenia", "ilef", "bmi", "energy")
to_keep <- rowSums(!is.na(colData(se_no_imp)[, sig_var])) == length(sig_var)
se_no_imp_all <- se_no_imp[, to_keep]
dim(se_no_imp_all)

# Build model matrix considering all predictor variables and interaction
design = model.matrix(~ colData(se_no_imp_all)$frailty*colData(se_no_imp_all)$sex +
                        colData(se_no_imp_all)$frailty*colData(se_no_imp_all)$alcohol + 
                        colData(se_no_imp_all)$frailty*colData(se_no_imp_all)$tobacco +
                        colData(se_no_imp_all)$frailty*colData(se_no_imp_all)$diabetes +
                        colData(se_no_imp_all)$frailty*colData(se_no_imp_all)$chf +
                        colData(se_no_imp_all)$frailty*colData(se_no_imp_all)$depression +
                        colData(se_no_imp_all)$frailty*colData(se_no_imp_all)$osteoarthritis +
                        colData(se_no_imp_all)$frailty*colData(se_no_imp_all)$sarcopenia +
                        colData(se_no_imp_all)$frailty*colData(se_no_imp_all)$ilef +
                        colData(se_no_imp_all)$frailty*colData(se_no_imp_all)$bmi +
                        colData(se_no_imp_all)$frailty*colData(se_no_imp_all)$energy)
colnames(design) = c("constant", "frailty", "sex", "alcohol_monthly", "alcohol_weekly",
                     "former", "current", "diabetes", "chf", "depression", 
                     "osteoarthritis", "sarcopenia", "ilef", "bmi", "energy",
                     "frailty:sex","frailty:alcohol_monthly",
                     "frailty:alcohol_weekly", "frailty:former",
                     "frailty:current", "frailty:diabetes",
                     "frailty:chf", "frailty:depression",
                     "frailty:osteoarthritis", "frailty:sarcopenia",
                     "frailty:ilef", "frailty:bmi", "frailty:energy")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_no_imp_all), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 29 coefficients per each column of the model matrix. Adjusted using the Benjamini-Hochberg
diff_63 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_64 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_65 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_66 <- topTable(fit1,coef=5, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_67 <- topTable(fit1,coef=6, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_68 <- topTable(fit1,coef=7, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_69 <- topTable(fit1,coef=8, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_70 <- topTable(fit1,coef=9, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_71 <- topTable(fit1,coef=10, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_72 <- topTable(fit1,coef=11, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_73 <- topTable(fit1,coef=12, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_74 <- topTable(fit1,coef=13, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_75 <- topTable(fit1,coef=14, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_76 <- topTable(fit1,coef=15, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_77 <- topTable(fit1,coef=16, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_78 <- topTable(fit1,coef=17, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_79 <- topTable(fit1,coef=18, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_80 <- topTable(fit1,coef=19, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_81 <- topTable(fit1,coef=20, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_82 <- topTable(fit1,coef=21, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_83 <- topTable(fit1,coef=22, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_84 <- topTable(fit1,coef=23, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_85 <- topTable(fit1,coef=24, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_86 <- topTable(fit1,coef=25, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_87 <- topTable(fit1,coef=26, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_88 <- topTable(fit1,coef=27, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))
diff_89 <- topTable(fit1,coef=28, adjust ="BH", sort.by="P", number=nrow(assay(se_no_imp_all)))

# frailty
diff_63 <- diff_63[diff_63$adj.P.Val < 0.05, ]
nrow(diff_63) 
rownames(diff_63)
# 0

# sex
diff_64 <- diff_64[diff_64$adj.P.Val < 0.05, ]
nrow(diff_64)
rownames(diff_64)
# 0

# alcohol monthly
diff_65 <- diff_65[diff_65$adj.P.Val < 0.05, ]
nrow(diff_65) 
rownames(diff_65)
# 0

# alcohol weekly
diff_66 <- diff_66[diff_66$adj.P.Val < 0.05, ]
nrow(diff_66) 
rownames(diff_66)
# 1: "P55259"

diff_66_up <- diff_66[diff_66$logFC > 0,]
nrow(diff_66_up) 
rownames(diff_66_up)
# 1
# "P55259"
diff_66_down <- diff_66[diff_66$logFC < 0,]
nrow(diff_66_down) 
rownames(diff_66_down)
# 0

# former smoker
diff_67 <- diff_67[diff_67$adj.P.Val < 0.05, ]
nrow(diff_67) 
rownames(diff_67)
# 0

# current smoker
diff_68 <- diff_68[diff_68$adj.P.Val < 0.05, ]
nrow(diff_68)
rownames(diff_68)
# 0

# diabetes
diff_69 <- diff_69[diff_69$adj.P.Val < 0.05, ]
nrow(diff_69) 
rownames(diff_69)
# 0

# chf
diff_70 <- diff_70[diff_70$adj.P.Val < 0.05, ]
nrow(diff_70) 
rownames(diff_70)
# 0

# depression
diff_71 <- diff_71[diff_71$adj.P.Val < 0.05, ]
nrow(diff_71) 
rownames(diff_71)
# 1: "A9KNK6"

diff_71_up <- diff_71[diff_71$logFC > 0,]
nrow(diff_71_up) 
rownames(diff_71_up)
# 1
# "A9KNK6"
diff_71_down <- diff_71[diff_71$logFC < 0,]
nrow(diff_71_down) 
rownames(diff_71_down)
# 0

# osteoarthritis
diff_72 <- diff_72[diff_72$adj.P.Val < 0.05, ]
nrow(diff_72)
rownames(diff_72)
# 0

# sarcopenia
diff_73 <- diff_73[diff_73$adj.P.Val < 0.05, ]
nrow(diff_73)
rownames(diff_73)
# 0

# ilef
diff_74 <- diff_74[diff_74$adj.P.Val < 0.05, ]
nrow(diff_74)
rownames(diff_74)
# 0

# bmi
diff_75 <- diff_75[diff_75$adj.P.Val < 0.05, ]
nrow(diff_75)
rownames(diff_75)
# 0

# energy
diff_76 <- diff_76[diff_76$adj.P.Val < 0.05, ]
nrow(diff_76)
rownames(diff_76)
# 0

# frailty:sex
diff_77 <- diff_77[diff_77$adj.P.Val < 0.05, ]
nrow(diff_77) 
rownames(diff_77)
# 3: "P95544" "A6KYJ0" "Q189R2"

diff_77_up <- diff_77[diff_77$logFC > 0,]
nrow(diff_77_up) 
rownames(diff_77_up)
# 0
diff_77_down <- diff_77[diff_77$logFC < 0,]
nrow(diff_77_down) 
rownames(diff_77_down)
# 3
# "P95544" "A6KYJ0" "Q189R2"

# frailty:alcohol_monthly
diff_78 <- diff_78[diff_78$adj.P.Val < 0.05, ]
nrow(diff_78)
rownames(diff_78)
# 0

# frailty:alcohol_weekly
diff_79 <- diff_79[diff_79$adj.P.Val < 0.05, ]
nrow(diff_79)
rownames(diff_79)
# 0

# frailty:former
diff_80 <- diff_80[diff_80$adj.P.Val < 0.05, ]
nrow(diff_80) 
rownames(diff_80)
# 1: "P95544"

diff_80_up <- diff_80[diff_80$logFC > 0,]
nrow(diff_80_up) 
rownames(diff_80_up)
# 0
diff_80_down <- diff_80[diff_80$logFC < 0,]
nrow(diff_80_down) 
rownames(diff_80_down)
# 1
# "P95544"

# frailty:current
diff_81 <- diff_81[diff_81$adj.P.Val < 0.05, ]
nrow(diff_81)
rownames(diff_81)
# 1: "A6KXA0"

diff_81_up <- diff_81[diff_81$logFC > 0,]
nrow(diff_81_up) 
rownames(diff_81_up)
# 1
# "A6KXA0"
diff_81_down <- diff_81[diff_81$logFC < 0,]
nrow(diff_81_down) 
rownames(diff_81_down)
# 0

# frailty:diabetes
diff_82 <- diff_82[diff_82$adj.P.Val < 0.05, ]
nrow(diff_82)
rownames(diff_82)
# 0

# frailty:chf
diff_83 <- diff_83[diff_83$adj.P.Val < 0.05, ]
nrow(diff_83) 
rownames(diff_83)
# 0

# frailty:depression
diff_84 <- diff_84[diff_84$adj.P.Val < 0.05, ]
nrow(diff_84) 
rownames(diff_84)
# 0

# frailty:osteoarthritis
diff_85 <- diff_85[diff_85$adj.P.Val < 0.05, ]
nrow(diff_85) 
rownames(diff_85)
# 0

# frailty:sarcopenia
diff_86 <- diff_86[diff_86$adj.P.Val < 0.05, ]
nrow(diff_86)
rownames(diff_86)
# 0

# frailty:ilef
diff_87 <- diff_87[diff_87$adj.P.Val < 0.05, ]
nrow(diff_87) 
rownames(diff_87)
# 0

# frailty:bmi
diff_88 <- diff_88[diff_88$adj.P.Val < 0.05, ]
nrow(diff_88) 
rownames(diff_88)
# 2: "Q59309" "P0C2E7"

diff_88_up <- diff_88[diff_88$logFC > 0,]
nrow(diff_88_up) 
rownames(diff_88_up)
# 2
# "Q59309" "P0C2E7"
diff_88_down <- diff_88[diff_88$logFC < 0,]
nrow(diff_88_down) 
rownames(diff_88_down)
# 0

sel0 = which("P0C2E7"==rownames(assay(se_no_imp_all)))
df0 = data.frame(frailty=colData(se_no_imp_all)[,c("frailty")], 
                 bmi=colData(se_no_imp_all)[,c("bmi")], 
                 expression=assay(se_no_imp_all)[sel0,])
ggplot(df0, aes(x = bmi, y = expression, color = frailty)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", aes(linetype = frailty))

# frailty:energy
diff_89 <- diff_89[diff_89$adj.P.Val < 0.05, ]
nrow(diff_89)
rownames(diff_89)
# 0

