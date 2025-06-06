#install.packages("ggplot2")
# Load libraries
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# Working directory
work_path <- getwd()

# Load functions
source(paste0(work_path, "/functions/create_bar_plot.R"))
source(paste0(work_path, "/functions/create_density_plot.R"))
source(paste0(work_path, "/functions/check_normal_and_test.R"))

# Read metadata
metadata <- read.csv(paste0(work_path, "/data/Tabla metadatos ENRICA.csv"),
                     sep=";")

# Get frailty cohort cases and columns of interest
metadata_ft <- metadata[metadata$Fragilidad %in% c("FT", "NFT"),
                        -ncol(metadata)]

# Class of each variable
str(metadata_ft)

# 203 in total, but 2 cases are missing because clinical info not available:
# EE266 	MALE 	NFT 	
# EE324 	MALE 	FT 

# Include these two cases to metadata
missing_cases <- data.frame(
  ID.FISABIO = c("EE266", "EE324"),
  w23sex = c(0, 0),
  Fragilidad = c("NFT", "FT")
)

# Ensure that missing_cases has the same columns as metadata_ft
missing_columns <- setdiff(names(metadata_ft), names(missing_cases))
# Add any missing columns to missing_cases with NA values
missing_cases[missing_columns] <- NA

# Add the new rows
metadata_ft <- rbind(metadata_ft, missing_cases)

# Replicate column
metadata_ft$replicate <- substr(metadata_ft$ID.FISABIO, 3, 5)

# Reorganise metadata
metadata_ft <- metadata_ft %>%
  rename(
    frailty = Fragilidad,  # Rename variables
    sex = w23sex,
    education = w23education,
    tobacco = w23tobacco,
    alcohol = w23alcohol,
    diabetes = w23diabetes,
    chronicrespdisease = w23chronicrespdisease,
    ami = w23ami, 
    stroke = w23stroke, 
    chf = w23chf,
    af = w23af, 
    osteoarthritis = w23osteoarthritis,
    reumatarthitis = w23reumatarthitis,
    hipfracture = w23hipfracture,
    depression = w23depression,
    cancer = w23cancer,
    sarcopenia = w23sarcopenia,
    age = w23age,
    medas = w23medas,
    energy = w23energy,
    weight = w23weight,
    height = w23height,
    bmi = w23bmi,
    sppbscore = w23sppbscore,
    sppb_balance = w23sppb_balance,
    sppb_gaitspeed = w23sppb_gaitspeed,
    sppb_sitstand = w23sppb_sitstand,
    fried_exhaustion = w23fried_exhaustion,
    fried_weightloss = w23fried_weightloss,
    fried_lowactivity = w23fried_lowactivity,
    fried_lowactivity_mod = w23fried_lowactivity_mod,
    fried_gaitspeed_chs = w23fried_gaitspeed_chs,
    fried_gaitspeed_mod = w23fried_gaitspeed_mod,
    fried_gripstrength_chs = w23fried_gripstrength_chs,
    fried_gripstrength_mod = w23fried_gripstrength_mod,
    heces = w23heces,
    frailty_fried_chs_2c = w23frailty_fried_chs_2c, 
    frailty_fried_chs_3c = w23frailty_fried_chs_3c,
    frailty_fried_mod_2c = w23frailty_fried_mod_2c, 
    frailty_fried_mod_3c = w23frailty_fried_mod_3c
  ) %>%
  mutate(
    frailty = factor(frailty, levels = c("NFT", "FT")),
    # Group alcohol consumption in 3 levels
    # 0: Rare | 1: Monthly | 2: Weekly
    # 1: A diario o casi a diario --> 2
    # 2: 5-6 días por semana --> 2
    # 3: 3-4 días por semana --> 2
    # 4: 1-2 días por semana --> 2
    # 5: 2-3 días en un mes --> 1
    # 6: 1 vez al mes --> 1
    # 7: Menos de una vez al mes --> 1
    # 8: No en los últimos 12 meses, he dejado de tomar alcohol --> 0
    # 9: Nunca o solamente unos sorbos para probar a lo largo de toda la vida --> 0
    alcohol = factor(ifelse(alcohol %in% c("1", "2", "3", "4"), "2",
                            ifelse(alcohol %in% c("5", "6", "7"), "1",
                                   ifelse(alcohol %in% c("8", "9"), "0", NA)
                            )
    ), levels = c("0", "1", "2")),
    
    ilef = factor(ifelse(sppbscore > 9, "0",  # ILEF calculated from SPPB score
                         ifelse(sppbscore <= 9, "1", NA)
    ))
  )


# Change variables with discrete values to factor
discrete_var <- c("sex",
                  "education",
                  "tobacco",
                  "alcohol",
                  "diabetes",
                  "chronicrespdisease",
                  "ami",
                  "stroke",
                  "chf",
                  "af",
                  "osteoarthritis",
                  "reumatarthitis",
                  "hipfracture",
                  "depression",
                  "cancer",
                  "sppb_balance",
                  "sppb_gaitspeed",
                  "sppb_sitstand",
                  "fried_exhaustion",
                  "fried_weightloss",
                  "fried_lowactivity",
                  "fried_lowactivity_mod",
                  "fried_gaitspeed_chs",
                  "fried_gaitspeed_mod",
                  "fried_gripstrength_chs",
                  "fried_gripstrength_mod",
                  "frailty_fried_chs_2c",
                  "frailty_fried_chs_3c",
                  "frailty_fried_mod_2c",
                  "frailty_fried_mod_3c",
                  "sarcopenia",
                  "heces",
                  "frailty",
                  "ilef")

for (var in discrete_var) {
  metadata_ft[[var]] <- as.factor(metadata_ft[[var]])
}

# Check classes
str(metadata_ft)

# Number of cases
dim(metadata_ft)
# 203  45

# Get summary
summary(metadata_ft)

# Plot discrete variables
for (var in discrete_var) {
  p <- create_bar_plot(metadata_ft, var, title = paste("Bar plot for", var))
  save_path <- paste0(work_path,"/plots/metadata/bar_plot_",var,".png")
  ggsave(filename = save_path, plot = p, width = 6, height = 4)
}

#Plot continuous variables
continuous_var <- c("age",
                   "medas",
                   "energy",
                   "weight",
                   "height",
                   "bmi",
                   "sppbscore")

for (var in continuous_var) {
  p <- create_density_plot(metadata_ft, var, title = paste("Density plot for", var))
  save_path <- paste0(work_path,"/plots/metadata/density_plot_",var,".png")
  ggsave(filename = save_path, plot = p, width = 6, height = 4)
}

# Save metadatos
save_path <- paste0(work_path,"/data/metadata.RData")
save(metadata_ft, file = save_path)

# Number of men (0) and women (1)
table(metadata_ft$sex)
# 0   1 
# 101 102 

# Mean age and sd
mean(metadata_ft[metadata_ft$sex == "0", "age"], na.rm = TRUE)
# 75.82
sd(metadata_ft[metadata_ft$sex == "0", "age"], na.rm = TRUE)
# 3.77
mean(metadata_ft[metadata_ft$sex == "1", "age"])
# 76.98
sd(metadata_ft[metadata_ft$sex == "1", "age"])
# 4.65

# ==================================
# Test correlation
# ==================================
colnames(metadata_ft)
vars <- metadata_ft[, c("sex", "education", "tobacco", "alcohol","diabetes",
                        "chronicrespdisease", "ami", "stroke", "chf", "af",
                        "osteoarthritis", "reumatarthitis", "hipfracture", 
                        "depression", "cancer", "sarcopenia", "ilef", "age",
                        "medas", "energy", #"weight", "height",
                        "bmi", #"sppbscore", "frailty_fried_chs_2c",
                        #"frailty_fried_chs_3c",
                        "frailty_fried_mod_2c"#, "frailty_fried_mod_3c"
                        )]

vars[] <- lapply(vars, as.numeric)

# Correlation and p-values matrix
n <- ncol(vars)
cor_matrix <- matrix(NA, n, n)
p_matrix <- matrix(NA, n, n)
colnames(cor_matrix) <- rownames(cor_matrix) <- colnames(vars)
colnames(p_matrix) <- rownames(p_matrix) <- colnames(vars)

# Test association between paired samples by Spearman correlation coefficient
for (i in 1:n) {
  for (j in 1:n) {
    test <- cor.test(vars[[i]], vars[[j]], method = "spearman", exact = FALSE)
    cor_matrix[i, j] <- test$estimate
    p_matrix[i, j] <- test$p.value
  }
}

# Create labels for p-values
p_text <- ifelse(p_matrix < 0.001, "***",
                 ifelse(p_matrix < 0.01, "**",
                        ifelse(p_matrix < 0.05, "*", "")))

# Plot heatmap
save_path <- paste0(work_path, "/plots/metadata/plot_heatmap_spearman_corr.png")
png(save_path, width = 8, height = 6, units = "in", res = 300)
Heatmap(cor_matrix,
       name = "Spearman correlation",
       col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
       cell_fun = function(j, i, x, y, ...)
         grid.text(p_text[i, j], x, y, gp = gpar(fontsize = 6))
       ,
       column_names_gp = gpar(fontsize = 8), 
       row_names_gp = gpar(fontsize = 8),
       show_row_dend = FALSE,
       clustering_distance_rows = "euclidean",
       clustering_distance_columns = "euclidean",
       use_raster = FALSE)
dev.off()

# Keep p-values < 0.05 from p_matrix and "frailty_fried_mod_2c"
p_fried <- p_matrix["frailty_fried_mod_2c",]
p_fried_sig <- p_fried[p_fried < 0.05]
# Keep significant correlations for "frailty_fried_mod_2c" from cor_matrix
cor_fried <- cor_matrix["frailty_fried_mod_2c",]
cor_fried_sig <- cor_fried[names(p_fried_sig)]
# Order by significance
ordered_vars <- names(sort(p_fried_sig))
cor_fried_sig <- cor_fried_sig[ordered_vars]
p_fried_sig   <- p_fried_sig[ordered_vars]
# Dataframe from cor_fried_sig and p_fried_sig
cor_fried_df <- data.frame(
  variable = names(cor_fried_sig),
  rho = cor_fried_sig,
  p_value = p_fried_sig[names(cor_fried_sig)]
)

# ==================================
# Test continuous variables 
# ==================================
# Age vs frailty
ggplot(metadata_ft, aes(x = age, fill = frailty)) +
  geom_density(alpha = 0.4) + 
  theme_minimal()
check_normal_and_test(metadata_ft, "age", "frailty")
# Wilcoxon rank sum test with continuity correction
# data:  df[[var_col]] by df[[group_col]]
# W = 2290.5, p-value = 4.618e-08
# alternative hypothesis: true location shift is not equal to 0

## Mean and sd
mean(metadata_ft[metadata_ft$frailty == "FT", "age"], na.rm = TRUE)
# 79.01562
sd(metadata_ft[metadata_ft$frailty == "FT", "age"], na.rm = TRUE)
# 4.547779
mean(metadata_ft[metadata_ft$frailty == "NFT", "age"], na.rm = TRUE)
# 75.18978
sd(metadata_ft[metadata_ft$frailty == "NFT", "age"], na.rm = TRUE)
# 3.532481

# MEDAS vs frailty
ggplot(metadata_ft, aes(x = medas, fill = frailty)) +
  geom_density(alpha = 0.4) + 
  theme_minimal()
check_normal_and_test(metadata_ft, "medas", "frailty")
# Wilcoxon rank sum test with continuity correction
# data:  df[[var_col]] by df[[group_col]]
# W = 4693, p-value = 0.002641
# alternative hypothesis: true location shift is not equal to 0

## Mean and sd
mean(metadata_ft[metadata_ft$frailty == "FT", "medas"], na.rm = TRUE)
# 7.309091
sd(metadata_ft[metadata_ft$frailty == "FT", "medas"], na.rm = TRUE)
# 1.399375
mean(metadata_ft[metadata_ft$frailty == "NFT", "medas"], na.rm = TRUE)
# 8.044776
sd(metadata_ft[metadata_ft$frailty == "NFT", "medas"], na.rm = TRUE)
# 1.662798

# Energy vs frailty
ggplot(metadata_ft, aes(x = energy, fill = frailty)) +
  geom_density(alpha = 0.4) + 
  theme_minimal()
check_normal_and_test(metadata_ft, "energy", "frailty")
# Welch Two Sample t-test
# data:  df[[var_col]] by df[[group_col]]
# t = 3.5316, df = 96.345, p-value = 0.0006357
# alternative hypothesis: true difference in means between group NFT and group FT is not equal to 0
# 95 percent confidence interval:
#   63.34243 225.91428
# sample estimates:
#   mean in group NFT  mean in group FT 
# 1934.828          1790.200

wilcox.test(metadata_ft[["energy"]] ~ metadata_ft[["frailty"]], data = metadata_ft)
#Wilcoxon rank sum test with continuity correction
#data:  metadata_ft[["energy"]] by metadata_ft[["frailty"]]
#W = 4890, p-value = 0.0004218
#alternative hypothesis: true location shift is not equal to 0

## Mean and sd
mean(metadata_ft[metadata_ft$frailty == "FT", "energy"], na.rm = TRUE)
# 1790.2
sd(metadata_ft[metadata_ft$frailty == "FT", "energy"], na.rm = TRUE)
# 259.2014
mean(metadata_ft[metadata_ft$frailty == "NFT", "energy"], na.rm = TRUE)
# 1934.828
sd(metadata_ft[metadata_ft$frailty == "NFT", "energy"], na.rm = TRUE)
# 247.0667

# Weight vs frailty
ggplot(metadata_ft, aes(x = weight, fill = frailty)) +
  geom_density(alpha = 0.4) + 
  theme_minimal()
check_normal_and_test(metadata_ft, "weight", "frailty")
# Wilcoxon rank sum test with continuity correction
# data:  df[[var_col]] by df[[group_col]]
# W = 4890, p-value = 0.188
# alternative hypothesis: true location shift is not equal to 0

# Height vs frailty
ggplot(metadata_ft, aes(x = height, fill = frailty)) +
  geom_density(alpha = 0.4) + 
  theme_minimal()
check_normal_and_test(metadata_ft, "height", "frailty")
# Wilcoxon rank sum test with continuity correction
# data:  df[[var_col]] by df[[group_col]]
# W = 6552, p-value = 1.65e-08
# alternative hypothesis: true location shift is not equal to 0

# BMI vs frailty
ggplot(metadata_ft, aes(x = bmi, fill = frailty)) +
  geom_density(alpha = 0.4) + 
  theme_minimal()
check_normal_and_test(metadata_ft, "bmi", "frailty")
# Wilcoxon rank sum test with continuity correction
# data:  df[[var_col]] by df[[group_col]]
# W = 3095.5, p-value = 0.0007527
# alternative hypothesis: true location shift is not equal to 0

## Mean and sd
mean(metadata_ft[metadata_ft$frailty == "FT", "bmi"], na.rm = TRUE)
# 28.76562
sd(metadata_ft[metadata_ft$frailty == "FT", "bmi"], na.rm = TRUE)
# 4.891659
mean(metadata_ft[metadata_ft$frailty == "NFT", "bmi"], na.rm = TRUE)
# 26.43796
sd(metadata_ft[metadata_ft$frailty == "NFT", "bmi"], na.rm = TRUE)
# 3.131579

# SPPB score vs frailty
ggplot(metadata_ft, aes(x = sppbscore, fill = frailty)) +
  geom_density(alpha = 0.4) + 
  theme_minimal()
check_normal_and_test(metadata_ft, "sppbscore", "frailty")
# Wilcoxon rank sum test with continuity correction
# data:  df[[var_col]] by df[[group_col]]
# W = 8352.5, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

# ==================================
# Test categorical variables 
# ==================================
# Contingency table: sex vs frailty
cont <- table(metadata_ft$frailty, metadata_ft$sex)
cont
#     0  1
# NFT 84 54
# FT  17 48
#chisq.test(cont)
fisher.test(cont)
# Fisher's Exact Test for Count Data
# data:  cont
# p-value = 4.886e-06
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.199041 8.971141
# sample estimates:
# odds ratio 
#   4.358703 

# Contingency table: education vs frailty
cont <- table(metadata_ft$frailty, metadata_ft$education)
cont
# 0  1  2
# NFT 73 24 40
# FT  43 14  7
#chisq.test(cont)
fisher.test(cont)
# Fisher's Exact Test for Count Data
# data:  cont
# p-value = 0.01287
# alternative hypothesis: two.sided

# Contingency table: tobacco vs frailty
cont <- table(metadata_ft$frailty, metadata_ft$tobacco)
cont
# 0  1  2
# NFT 68 56 13
# FT  44 15  5
fisher.test(cont)
# Fisher's Exact Test for Count Data
# data:  cont
# p-value = 0.03465
# alternative hypothesis: two.sided

# Contingency table: alcohol vs frailty
cont <- table(metadata_ft$frailty, metadata_ft$alcohol)
cont
# 0  1  2
# NFT 27 26 80
# FT  31 13 18
fisher.test(cont)
# Fisher's Exact Test for Count Data
# data:  cont
# p-value = 3.132e-05
# alternative hypothesis: two.sided

# Contingency table: diabetes vs frailty
cont <- table(metadata_ft$frailty, metadata_ft$diabetes)
cont
# 0   1
# NFT 124  13
# FT   45  19
fisher.test(cont)
# Fisher's Exact Test for Count Data
# data:  cont
# p-value = 0.0006715
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.714061 9.596599
# sample estimates:
# odds ratio 
#  3.994854
  
# Contingency table: chronicrespdisease vs frailty
cont <- table(metadata_ft$frailty, metadata_ft$chronicrespdisease)
cont
# 0   1
# NFT 124  13
# FT   55   8
fisher.test(cont)
# Fisher's Exact Test for Count Data
# data:  cont
# p-value = 0.4695
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 0.469226 3.850142
# sample estimates:
# odds ratio 
#  1.385076 
  
# Contingency table: ami vs frailty
cont <- table(metadata_ft$frailty, metadata_ft$ami)
cont
# 0   1
# NFT 132   5
# FT   60   4
fisher.test(cont)
# Fisher's Exact Test for Count Data
# data:  cont
# p-value = 0.4698
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.3358585 8.4731555
# sample estimates:
# odds ratio 
#   1.754583 

# Contingency table: stroke vs frailty
cont <- table(metadata_ft$frailty, metadata_ft$stroke)
cont
# 0   1
# NFT 135   2
# FT   60   3
fisher.test(cont)
# Fisher's Exact Test for Count Data
# data:  cont
# p-value = 0.1809
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.3739862 41.1401062
# sample estimates:
# odds ratio 
#   3.351773

# Contingency table: chf vs frailty
cont <- table(metadata_ft$frailty, metadata_ft$chf)
cont
# 0   1
# NFT 134   3
# FT   55   7
fisher.test(cont)
# Fisher's Exact Test for Count Data
# data:  cont
# p-value = 0.01128
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.230308 34.937916
# sample estimates:
# odds ratio 
#    5.62866 

# Contingency table: af vs frailty
cont <- table(metadata_ft$frailty, metadata_ft$af)
cont
# 0   1
# NFT 134   3
# FT   57   7
fisher.test(cont)
# Fisher's Exact Test for Count Data
# data:  cont
# p-value = 0.01293
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.18886 33.68034
# sample estimates:
# odds ratio 
#  5.433394 

# Contingency table: osteoarthritis vs frailty
cont <- table(metadata_ft$frailty, metadata_ft$osteoarthritis)
cont
# 0  1
# NFT 98 37
# FT  19 42
fisher.test(cont)
# Fisher's Exact Test for Count Data
# data:  cont
# p-value = 5.894e-08
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   2.88285 12.02329
# sample estimates:
# odds ratio 
#   5.794481 

# Contingency table: reumatarthitis vs frailty
cont <- table(metadata_ft$frailty, metadata_ft$reumatarthitis)
cont
# 0   1
# NFT 125   9
# FT   53   7
fisher.test(cont)
# Fisher's Exact Test for Count Data
# data:  cont
# p-value = 0.266
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.5479213 5.8483006
# sample estimates:
# odds ratio 
#   1.828251 

# Contingency table: depression vs frailty
cont <- table(metadata_ft$frailty, metadata_ft$depression)
cont
# 0   1
# NFT 132   5
# FT   52  12
fisher.test(cont)
# Fisher's Exact Test for Count Data
# data:  cont
# p-value = 0.0007165
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.864796 22.959462
# sample estimates:
# odds ratio 
#    6.02959 

# Contingency table: hipfracture vs frailty
cont <- table(metadata_ft$frailty, metadata_ft$hipfracture)
cont
# 0   1
# NFT 136   1
# FT   60   4
fisher.test(cont)
# Fisher's Exact Test for Count Data
# data:  cont
# p-value = 0.03641
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#    0.8634421 448.5673789
# sample estimates:
# odds ratio 
#   8.962024 

# Contingency table: cancer vs frailty
cont <- table(metadata_ft$frailty, metadata_ft$cancer)
cont
# 0   1
# NFT 127  10
# FT   54  10
fisher.test(cont)
# Fisher's Exact Test for Count Data
# data:  cont
# p-value = 0.07883
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.8221271 6.6735986
# sample estimates:
# odds ratio 
#  2.340832 

# Contingency table: sarcopenia vs frailty
cont <- table(metadata_ft$frailty, metadata_ft$sarcopenia)
cont
# 0   1
# NFT 133   1
# FT   44  14
fisher.test(cont)
# Fisher's Exact Test for Count Data
# data:  cont
# p-value = 1.792e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#     5.9947 1784.9038
# sample estimates:
# odds ratio 
#   41.48235 

# Contingency table: ILEF vs frailty
cont <- table(metadata_ft$frailty, metadata_ft$ilef)
cont
# 0   1
# NFT 132   5
# FT   17  47
fisher.test(cont)
# Fisher's Exact Test for Count Data
# data:  cont
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   23.81665 257.44851
# sample estimates:
# odds ratio 
#   70.00736 