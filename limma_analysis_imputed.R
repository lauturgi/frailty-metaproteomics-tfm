library(ggplot2)
library(limma)
library(dplyr)

# ==================================
# Limma analysis
# ==================================
# Running DE analysis with data from preprocessing imputed: Norm MaxLFQ,
# log2transformed, NA in at least 1 condition < 0.5, perseus-type imputation
# ==================================
# Load filtered and imputed SummarizedExperiment
work_path <- getwd()
load(paste0(work_path, "/data/se_perseus.RData"))
# ==================================

# ==================================
# 1. Frailty
# ==================================

# MODEL 1: FT vs NFT
##  Build model matrix considering only frailty
design = model.matrix(~ colData(se_perseus)$frailty)
colnames(design) = c("constant", "frailty")
head(design)

## Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_perseus), design)

## Estimated coefficients for each protein
head(coef(fit))

## Moderate t-test of the differential expression by means of Bayes' empirical
## moderation of the standard errors to a global value
fit1 = eBayes(fit)

## p-values adjusted using the Benjamini-Hochberg method
diff_1 <- topTable(fit1, coef=2, adjust="BH", sort.by="P", number=nrow(assay(se_perseus)))

create_volcano_plot(diff_1, save_path = paste0(work_path, "/plots/volcano_plot_model1_frailty.png"),
                    title = "Model 1: Differentially expressed proteins - Frailty")

diff_1 <- diff_1[diff_1$adj.P.Val < 0.05, ]
nrow(diff_1) 
rownames(diff_1)
# 3 proteins (p.adj < 0.05): "P05109" "P06702" "A9KRZ1"

# Up vs down expressed
diff_1_up <- diff_1[diff_1$logFC > 0,]
# 2
diff_1_down <- diff_1[diff_1$logFC < 0,]
# 1

sel0 = which("P05109"==rownames(assay(se_perseus)))
df0 = data.frame(frailty=colData(se_perseus)[,c("frailty")], expression=assay(se_perseus)[sel0,])
ggplot(df0,aes(x=frailty,y=expression)) + geom_boxplot()

# ==================================
# 2. Two predictor variables
# ==================================
#   2.1. Fraily and sex
# ==================================

# MODEL 2
# Build model matrix considering frailty and sex
design = model.matrix(~ colData(se_perseus)$frailty + colData(se_perseus)$sex)
colnames(design) = c("constant", "frailty", "sex")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_perseus), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 3 coefficients per each column of the model matrix. Adjusted using the Benjamini-Hochberg
diff_2 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus)))
diff_3 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus)))

create_volcano_plot(diff_2, save_path = paste0(work_path, "/plots/volcano_plot_model2_frailty.png"),
                    title = "Model 2: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_3, save_path = paste0(work_path, "/plots/volcano_plot_model2_sex.png"),
                    title = "Model 2: Differentially expressed proteins - Sex")

diff_2 <- diff_2[diff_2$adj.P.Val < 0.05, ]
nrow(diff_2) 
rownames(diff_2)
# 2: "P06702" "P05109"

diff_2_up <- diff_2[diff_2$logFC > 0,]
# 2
diff_2_down <- diff_2[diff_2$logFC < 0,]
# 0

diff_3 <- diff_3[diff_3$adj.P.Val < 0.05, ]
nrow(diff_3) 
rownames(diff_3)
# 3: "A9KJL3" "A4W5A0" "A0Q2T1"

sel0 = which("A9KJL3"==rownames(assay(se_perseus)))
df0 = data.frame(sex=colData(se_perseus)[,c("sex")], expression=assay(se_perseus)[sel0,])
ggplot(df0,aes(x=sex,y=expression)) + geom_boxplot()

diff_3_up <- diff_3[diff_3$logFC > 0,]
# 0
diff_3_down <- diff_3[diff_3$logFC < 0,]
# 3

# MODEL 3: Consider interaction between frailty and sex
# Build model matrix considering frailty, sex and interaction
design = model.matrix(~ colData(se_perseus)$frailty*colData(se_perseus)$sex)
colnames(design) = c("constant", "frailty", "sex", "frailty:sex")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_perseus), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of the model matrix. Adjusted using the Benjamini-Hochberg
diff_4 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus)))
diff_5 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus)))
diff_6 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus)))

create_volcano_plot(diff_4, save_path = paste0(work_path, "/plots/volcano_plot_model3_frailty.png"),
                    title = "Model 3: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_5, save_path = paste0(work_path, "/plots/volcano_plot_model3_sex.png"),
                    title = "Model 3: Differentially expressed proteins - Sex")
create_volcano_plot(diff_6, save_path = paste0(work_path, "/plots/volcano_plot_model3_frailty_sex.png"),
                    title = "Model 3: Differentially expressed proteins - Frailty:Sex")

diff_4 <- diff_4[diff_4$adj.P.Val < 0.05, ]
nrow(diff_4) 
rownames(diff_4)
# 0

diff_5 <- diff_5[diff_5$adj.P.Val < 0.05, ]
nrow(diff_5) 
rownames(diff_5)
# 0

diff_6 <- diff_6[diff_6$adj.P.Val < 0.05, ]
nrow(diff_6) 
rownames(diff_6)
# 0

sel0 = which("C4Z0Q6"==rownames(assay(se_perseus)))
df0 = data.frame(frailty=colData(se_perseus)[,c("frailty")], 
                 sex=colData(se_perseus)[,c("sex")], 
                 expression=assay(se_perseus)[sel0,])
df0$frailty_sex <- paste(df0$frailty, df0$sex, sep = "_")
ggplot(df0,aes(x=frailty_sex,y=expression)) + geom_boxplot()

# ==================================
#   2.2. Frailty and alcohol
# ==================================

# MODEL 4: 
# Build model matrix considering frailty and alcohol
se_perseus_alcohol <- se_perseus[, !is.na(colData(se_perseus)$alcohol)]
dim(se_perseus_alcohol)
design = model.matrix(~ colData(se_perseus_alcohol)$frailty + colData(se_perseus_alcohol)$alcohol)
colnames(design) = c("constant", "frailty", "alcohol_monthly", "alcohol_weekly")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_perseus_alcohol), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of the model matrix. Adjusted using the Benjamini-Hochberg
diff_7 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_alcohol)))
diff_8 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_alcohol)))
diff_9 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_alcohol)))

create_volcano_plot(diff_7, save_path = paste0(work_path, "/plots/volcano_plot_model4_frailty.png"),
                    title = "Model 4: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_8, save_path = paste0(work_path, "/plots/volcano_plot_model4_alcohol_monthly.png"),
                    title = "Model 4: Differentially expressed proteins - Alcohol monthly")
create_volcano_plot(diff_9, save_path = paste0(work_path, "/plots/volcano_plot_model4_alcohol_weekly.png"),
                    title = "Model 4: Differentially expressed proteins - Alcohol weekly")


diff_7 <- diff_7[diff_7$adj.P.Val < 0.05, ]
nrow(diff_7) 
rownames(diff_7)
# 2: "P05109" "P06702"

diff_8 <- diff_8[diff_8$adj.P.Val < 0.05, ]
nrow(diff_8) 
rownames(diff_8)
# 0

diff_9 <- diff_9[diff_9$adj.P.Val < 0.05, ]
nrow(diff_9) 
rownames(diff_9)
# 0


# MODEL 5: Consider interaction between frailty and alcohol
# Build model matrix considering frailty, alcohol and interaction
design = model.matrix(~ colData(se_perseus_alcohol)$frailty*colData(se_perseus_alcohol)$alcohol)
colnames(design) = c("constant", "frailty", "alcohol_monthly", "alcohol_weekly", "frailty:alcohol_monthly", "frailty:alcohol_weekly")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_perseus_alcohol), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 6 coefficients per each column of the model matrix. Adjusted using the Benjamini-Hochberg
diff_10 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_alcohol)))
diff_11 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_alcohol)))
diff_12 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_alcohol)))
diff_13 <- topTable(fit1,coef=5, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_alcohol)))
diff_14 <- topTable(fit1,coef=6, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_alcohol)))

create_volcano_plot(diff_10, save_path = paste0(work_path, "/plots/volcano_plot_model5_frailty.png"),
                    title = "Model 5: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_11, save_path = paste0(work_path, "/plots/volcano_plot_model5_alcohol_monthly.png"),
                    title = "Model 5: Differentially expressed proteins - Alcohol monthly")
create_volcano_plot(diff_12, save_path = paste0(work_path, "/plots/volcano_plot_model5_alcohol_weekly.png"),
                    title = "Model 5: Differentially expressed proteins - Alcohol weekly")
create_volcano_plot(diff_13, save_path = paste0(work_path, "/plots/volcano_plot_model5_frailty_alcohol_monthly.png"),
                    title = "Model 5: Differentially expressed proteins - Frailty:Alcohol monthly")
create_volcano_plot(diff_14, save_path = paste0(work_path, "/plots/volcano_plot_model5_frailty_alcohol_weekly.png"),
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
# 0

diff_13 <- diff_13[diff_13$adj.P.Val < 0.05, ]
nrow(diff_13) 
rownames(diff_13)
# 0

diff_14 <- diff_14[diff_14$adj.P.Val < 0.05, ]
nrow(diff_14) 
rownames(diff_14)
# 0

sel0 = which("P22983"==rownames(assay(se_perseus_alcohol)))
df0 = data.frame(frailty=colData(se_perseus_alcohol)[,c("frailty")], 
                 alcohol=colData(se_perseus_alcohol)[,c("alcohol")], 
                 expression=assay(se_perseus_alcohol)[sel0,])
df0$frailty_alcohol <- paste(df0$frailty, df0$alcohol, sep = "_")
ggplot(df0,aes(x=frailty_alcohol,y=expression)) + geom_boxplot()


# ==================================
#   2.3. Frailty and tobacco
# ==================================

# MODEL 6
# Build model matrix considering frailty and tobacco
se_perseus_tobacco <- se_perseus[, !is.na(colData(se_perseus)$tobacco)]
dim(se_perseus_tobacco)
design = model.matrix(~ colData(se_perseus_tobacco)$frailty + colData(se_perseus_tobacco)$tobacco)
colnames(design) = c("constant", "frailty", "former", "current")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_perseus_tobacco), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 4 coefficients per each column of the model matrix. Adjusted using the Benjamini-Hochberg
diff_15 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_tobacco)))
diff_16 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_tobacco)))
diff_17 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_tobacco)))

create_volcano_plot(diff_15, save_path = paste0(work_path, "/plots/volcano_plot_model6_frailty.png"),
                    title = "Model 6: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_16, save_path = paste0(work_path, "/plots/volcano_plot_model6_tabacco_former.png"),
                    title = "Model 6: Differentially expressed proteins - Tabacco former")
create_volcano_plot(diff_17, save_path = paste0(work_path, "/plots/volcano_plot_model6_tabacco_current.png"),
                    title = "Model 6: Differentially expressed proteins - Tabacco current")

diff_15 <- diff_15[diff_15$adj.P.Val < 0.05, ]
nrow(diff_15) 
rownames(diff_15)
# 2: "P05109" "P06702"

diff_15_up <- diff_15[diff_15$logFC > 0,]
# 2
diff_15_down <- diff_15[diff_15$logFC < 0,]
# 0

diff_16 <- diff_16[diff_16$adj.P.Val < 0.05, ]
nrow(diff_16) 
rownames(diff_16)
# 0

diff_17 <- diff_17[diff_17$adj.P.Val < 0.05, ]
nrow(diff_17) 
rownames(diff_17)
# 0

# MODEL 7: Consider interaction between frailty and tobacco
# Build model matrix considering frailty, tobacco and interaction
design = model.matrix(~ colData(se_perseus_tobacco)$frailty*colData(se_perseus_tobacco)$tobacco)
colnames(design) = c("constant", "frailty", "former", "current", "frailty:former", "frailty:current")
head(design)

# Fit a linear model for each row of the expression matrix
fit = lmFit(assay(se_perseus_tobacco), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 6 coefficients per each column of the model matrix. Adjusted using the Benjamini-Hochberg
diff_18 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_tobacco)))
diff_19 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_tobacco)))
diff_20 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_tobacco)))
diff_21 <- topTable(fit1,coef=5, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_tobacco)))
diff_22 <- topTable(fit1,coef=6, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_tobacco)))

create_volcano_plot(diff_18, save_path = paste0(work_path, "/plots/volcano_plot_model7_frailty.png"),
                    title = "Model 7: Differentially expressed proteins - Frailty")
create_volcano_plot(diff_19, save_path = paste0(work_path, "/plots/volcano_plot_model7_tabacco_former.png"),
                    title = "Model 7: Differentially expressed proteins - Tabacco former")
create_volcano_plot(diff_20, save_path = paste0(work_path, "/plots/volcano_plot_model7_tabacco_current.png"),
                    title = "Model 7: Differentially expressed proteins - Tabacco current")
create_volcano_plot(diff_21, save_path = paste0(work_path, "/plots/volcano_plot_model7_frailty_tabacco_former.png"),
                    title = "Model 7: Differentially expressed proteins - Frailty:Tabacco former")
create_volcano_plot(diff_22, save_path = paste0(work_path, "/plots/volcano_plot_model7_frailty_tabacco_current.png"),
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
# 106: 
# "A6KYH6" "Q1Q899" "C4ZBT7" "A6KYJ6" "P0C2E7" "P19543" "Q5L9B6" "A6L1X1" "A6KYJ4" "A9KK94" "C4ZBS1" "A9KSJ1"
# "A6L0V1" "A9KNK6" "C4ZB90" "A9KKU0" "A9KK92" "C4Z2V8" "C4Z3R4" "B8I7Y6" "A9KJJ1" "P00698" "A9KRZ2" "C4ZD46"
# "C4ZBS3" "C4Z2T1" "A9KJI4" "C4Z2T9" "A6TH53" "C4ZBL1" "C4ZBG1" "A9KJL3" "A6LPS3" "C4Z2R7" "Q0SQC8" "Q5L8C5"
# "A9KNC4" "G3KIM7" "C4Z9C9" "Q59309" "C0R090" "P95544" "Q5L923" "Q8A294" "O83023" "A9VT65" "Q57079" "A6KYK3"
# "C4Z0D7" "A6KYH8" "P01012" "Q88XX2" "Q8R7U7" "A4W5A0" "A6KYH0" "C4Z1J4" "Q8RQP4" "P22983" "A6KYJ8" "C4ZBD3"
# "Q5M0M5" "A4XI37" "P0C0G6" "Q59199" "A6LFQ4" "C4ZF71" "A6L1L8" "P02768" "A7I3U7" "O09460" "A6KYI1" "A6KYJ2"
# "A6L048" "A3DJI2" "A6KYJ7" "A6KYK2" "A5I7J4" "Q04FQ4" "B2UYT8" "A9VP74" "C4ZBD5" "A0Q2T1" "A8YXK9" "A6KXA0"
# "C4Z2U0" "C4ZBT1" "C4Z2T8" "P69776" "C4Z5P8" "P26823" "P24295" "C4KZP0" "C4ZBS5" "A6KXL2" "Q8A753" "C4Z2T0"
# "Q8A1A2" "Q5LHW2" "A9KJI0" "P55990" "Q8A463" "A6KYK6" "A7HBL7" "A6KYH1" "B3DPK4" "H9L478"

sel0 = which("Q1Q899"==rownames(assay(se_perseus_tobacco)))
df0 = data.frame(frailty=colData(se_perseus_tobacco)[,c("frailty")], 
                 tobacco=colData(se_perseus_tobacco)[,c("tobacco")], 
                 expression=assay(se_perseus_tobacco)[sel0,])
df0$frailty_tobacco <- paste(df0$frailty, df0$tobacco, sep = "_")
ggplot(df0,aes(x=frailty_tobacco,y=expression)) + geom_boxplot()

diff_22_up <- diff_22[diff_22$logFC > 0,]
# 106
diff_22_down <- diff_22[diff_22$logFC < 0,]
# 0

# ==================================
#   2.4. Frailty and diabetes
# ==================================
# MODEL 8
# Build model matrix considering frailty and diabetes
# MODEL 9
# Build model matrix considering frailty, diabetes and interaction
# ==================================
#   2.5. Frailty and chf
# ==================================
# MODEL 10
# Build model matrix considering frailty and chf
# MODEL 11
# Build model matrix considering frailty, chf and interaction
# ==================================
#   2.6. Frailty and depression
# ==================================
# MODEL 12
# Build model matrix considering frailty and depression
# MODEL 13
# Build model matrix considering frailty, depression and interaction
# ==================================
#   2.7. Frailty and osteoarthritis
# ==================================
# MODEL 14
# Build model matrix considering frailty and osteoarthritis
# MODEL 15
# Build model matrix considering frailty, osteoarthritis and interaction
# ==================================
#   2.8. Frailty and sarcopenia
# ==================================
# MODEL 16
# Build model matrix considering frailty and sarcopenia
# MODEL 17
# Build model matrix considering frailty, sarcopenia and interaction
# ==================================
#   2.9. Frailty and bmi
# ==================================
# MODEL 18
# Build model matrix considering frailty and bmi
# MODEL 19
# Build model matrix considering frailty, bmi and interaction
# ==================================
#   2.10. Frailty and energy
# ==================================
# MODEL 20
# Build model matrix considering frailty and energy
# MODEL 21
# Build model matrix considering frailty, energy and interaction
# ==================================
#   2.11. Frailty and ilef
# ==================================
# MODEL 22
# Build model matrix considering frailty and ilef
# MODEL 23
# Build model matrix considering frailty, ilef and interaction

# ==================================
# 3. All predictor variables significantly related to frailty
# ==================================
# FINAL MODEL
# Keep cases with data for all predictor variables
sig_var <- c("alcohol", "tobacco", "diabetes", "chf", "depression",
             "osteoarthritis", "sarcopenia", "ilef", "bmi", "energy")
to_keep <- rowSums(!is.na(colData(se_perseus)[, sig_var])) == length(sig_var)
se_perseus_all <- se_perseus[, to_keep]
dim(se_perseus_all)

# Build model matrix considering all predictor variables and interaction
design = model.matrix(~ colData(se_perseus_all)$frailty*colData(se_perseus_all)$sex +
                        colData(se_perseus_all)$frailty*colData(se_perseus_all)$alcohol + 
                        colData(se_perseus_all)$frailty*colData(se_perseus_all)$tobacco +
                        colData(se_perseus_all)$frailty*colData(se_perseus_all)$diabetes +
                        colData(se_perseus_all)$frailty*colData(se_perseus_all)$chf +
                        colData(se_perseus_all)$frailty*colData(se_perseus_all)$depression +
                        colData(se_perseus_all)$frailty*colData(se_perseus_all)$osteoarthritis +
                        colData(se_perseus_all)$frailty*colData(se_perseus_all)$sarcopenia +
                        colData(se_perseus_all)$frailty*colData(se_perseus_all)$ilef +
                        colData(se_perseus_all)$frailty*colData(se_perseus_all)$bmi +
                        colData(se_perseus_all)$frailty*colData(se_perseus_all)$energy)
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
fit = lmFit(assay(se_perseus_all), design)
fit1 = eBayes(fit)
head(coef(fit1))

# 29 coefficients per each column of the model matrix. Adjusted using the Benjamini-Hochberg
diff_23 <- topTable(fit1,coef=2, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_24 <- topTable(fit1,coef=3, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_25 <- topTable(fit1,coef=4, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_26 <- topTable(fit1,coef=5, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_27 <- topTable(fit1,coef=6, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_28 <- topTable(fit1,coef=7, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_29 <- topTable(fit1,coef=8, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_30 <- topTable(fit1,coef=9, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_31 <- topTable(fit1,coef=10, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_32 <- topTable(fit1,coef=11, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_33 <- topTable(fit1,coef=12, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_34 <- topTable(fit1,coef=13, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_35 <- topTable(fit1,coef=14, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_36 <- topTable(fit1,coef=15, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_37 <- topTable(fit1,coef=16, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_38 <- topTable(fit1,coef=17, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_39 <- topTable(fit1,coef=18, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_40 <- topTable(fit1,coef=19, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_41 <- topTable(fit1,coef=20, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_42 <- topTable(fit1,coef=21, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_43 <- topTable(fit1,coef=22, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_44 <- topTable(fit1,coef=23, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_45 <- topTable(fit1,coef=24, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_46 <- topTable(fit1,coef=25, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_47 <- topTable(fit1,coef=26, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_48 <- topTable(fit1,coef=27, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))
diff_49 <- topTable(fit1,coef=28, adjust ="BH", sort.by="P", number=nrow(assay(se_perseus_all)))

# frailty
diff_23 <- diff_23[diff_23$adj.P.Val < 0.05, ]
nrow(diff_23) 
rownames(diff_23)
# 0

# sex
diff_24 <- diff_24[diff_24$adj.P.Val < 0.05, ]
nrow(diff_24) 
rownames(diff_24)
# 0

# alcohol monthly
diff_25 <- diff_25[diff_25$adj.P.Val < 0.05, ]
nrow(diff_25) 
rownames(diff_25)
# 0

# alcohol weekly
diff_26 <- diff_26[diff_26$adj.P.Val < 0.05, ]
nrow(diff_26) 
rownames(diff_26)
# 0

# former smoker
diff_27 <- diff_27[diff_27$adj.P.Val < 0.05, ]
nrow(diff_27) 
rownames(diff_27)
# 0

# current smoker
diff_28 <- diff_28[diff_28$adj.P.Val < 0.05, ]
nrow(diff_28) 
rownames(diff_28)
# 0

# diabetes
diff_29 <- diff_29[diff_29$adj.P.Val < 0.05, ]
nrow(diff_29) 
rownames(diff_29)
# 0

# chf
diff_30 <- diff_30[diff_30$adj.P.Val < 0.05, ]
nrow(diff_30) 
rownames(diff_30)
# 0

# depression
diff_31 <- diff_31[diff_31$adj.P.Val < 0.05, ]
nrow(diff_31) 
rownames(diff_31)
# 0

# osteoarthritis
diff_32 <- diff_32[diff_32$adj.P.Val < 0.05, ]
nrow(diff_32) 
rownames(diff_32)
# 0

# sarcopenia
diff_33 <- diff_33[diff_33$adj.P.Val < 0.05, ]
nrow(diff_33) 
rownames(diff_33)
# 0

# ilef
diff_34 <- diff_34[diff_34$adj.P.Val < 0.05, ]
nrow(diff_34) 
rownames(diff_34)
# 0

# bmi
diff_35 <- diff_35[diff_35$adj.P.Val < 0.05, ]
nrow(diff_35) 
rownames(diff_35)
# 0

# energy
diff_36 <- diff_36[diff_36$adj.P.Val < 0.05, ]
nrow(diff_36) 
rownames(diff_36)
# 1: "A9KRZ1"

# frailty:sex
diff_37 <- diff_37[diff_37$adj.P.Val < 0.05, ]
nrow(diff_37) 
rownames(diff_37)
# 0

sel0 = which("Q189R2"==rownames(assay(se_perseus_all)))
df0 = data.frame(frailty=colData(se_perseus_all)[,c("frailty")], 
                 sex=colData(se_perseus_all)[,c("sex")], 
                 expression=assay(se_perseus_all)[sel0,])
df0$frailty_sex <- paste(df0$frailty, df0$sex, sep = "_")
ggplot(df0,aes(x=frailty_sex,y=expression)) + geom_boxplot()

# frailty:alcohol_monthly
diff_38 <- diff_38[diff_38$adj.P.Val < 0.05, ]
nrow(diff_38) 
rownames(diff_38)
# 0

# frailty:alcohol_weekly
diff_39 <- diff_39[diff_39$adj.P.Val < 0.05, ]
nrow(diff_39) 
rownames(diff_39)
# 0

# frailty:former
diff_40 <- diff_40[diff_40$adj.P.Val < 0.05, ]
nrow(diff_40) 
rownames(diff_40)
# 0

# frailty:current
diff_41 <- diff_41[diff_41$adj.P.Val < 0.05, ]
nrow(diff_41) 
rownames(diff_41)
# 0

sel0 = which("A6KXA0"==rownames(assay(se_perseus_all)))
df0 = data.frame(frailty=colData(se_perseus_all)[,c("frailty")], 
                 tobacco=colData(se_perseus_all)[,c("tobacco")], 
                 expression=assay(se_perseus_all)[sel0,])
df0$frailty_tobacco <- paste(df0$frailty, df0$tobacco, sep = "_")
ggplot(df0,aes(x=frailty_tobacco,y=expression)) + geom_boxplot()


# frailty:diabetes
diff_42 <- diff_42[diff_42$adj.P.Val < 0.05, ]
nrow(diff_42) 
rownames(diff_42)
# 0

# frailty:chf
diff_43 <- diff_43[diff_43$adj.P.Val < 0.05, ]
nrow(diff_43) 
rownames(diff_43)
# 0

# frailty:depression
diff_44 <- diff_44[diff_44$adj.P.Val < 0.05, ]
nrow(diff_44) 
rownames(diff_44)
# 0

# frailty:osteoarthritis
diff_45 <- diff_45[diff_45$adj.P.Val < 0.05, ]
nrow(diff_45) 
rownames(diff_45)
# 0

# frailty:sarcopenia
diff_46 <- diff_46[diff_46$adj.P.Val < 0.05, ]
nrow(diff_46) 
rownames(diff_46)
# 0

# frailty:ilef
diff_47 <- diff_47[diff_47$adj.P.Val < 0.05, ]
nrow(diff_47) 
rownames(diff_47)
# 0

# frailty:bmi
diff_48 <- diff_48[diff_48$adj.P.Val < 0.05, ]
nrow(diff_48) 
rownames(diff_48)
# 0

sel0 = which("Q59309"==rownames(assay(se_perseus_all)))
df0 = data.frame(frailty=colData(se_perseus_all)[,c("frailty")], 
                 bmi=colData(se_perseus_all)[,c("bmi")], 
                 expression=assay(se_perseus_all)[sel0,])
ggplot(df0, aes(x = bmi, y = expression, color = frailty)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", aes(linetype = frailty))


# frailty:energy
diff_49 <- diff_49[diff_49$adj.P.Val < 0.05, ]
nrow(diff_49) 
rownames(diff_49)
# 0

