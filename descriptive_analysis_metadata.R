#install.packages("ggplot2")
library(ggplot2)
library(ComplexHeatmap)

# Working directory
work_path <- "/home/lauturgi/workspace/tfm"

# Metadata directory
metadata_path <- "/mnt/d/proteomica/fragilidad/datos/"
metadata_name <- "Tabla metadatos ENRICA.csv"

# Read metadata
metadata <- read.csv(paste0(metadata_path,metadata_name),sep=";")

# Get frailty cohort cases and columns of interest
metadata_ft <- metadata[metadata$Fragilidad %in% c("FT", "NFT"),
                        -ncol(metadata)]

# Class of each variable
str(metadata_ft)

# Change variables with discrete values to factor
discrete_var <- c("w23sex",
                  "w23education",
                  "w23tobacco",
                  "w23alcohol",
                  "w23diabetes",
                  "w23chronicrespdisease",
                  "w23ami",
                  "w23stroke",
                  "w23chf",
                  "w23af",
                  "w23osteoarthritis",
                  "w23reumatarthitis",
                  "w23hipfracture",
                  "w23depression",
                  "w23cancer",
                  "w23sppb_balance",
                  "w23sppb_gaitspeed",
                  "w23sppb_sitstand",
                  "w23fried_exhaustion",
                  "w23fried_weightloss",
                  "w23fried_lowactivity",
                  "w23fried_lowactivity_mod",
                  "w23fried_gaitspeed_chs",
                  "w23fried_gaitspeed_mod",
                  "w23fried_gripstrength_chs",
                  "w23fried_gripstrength_mod",
                  "w23frailty_fried_chs_2c",
                  "w23frailty_fried_chs_3c",
                  "w23frailty_fried_mod_2c",
                  "w23frailty_fried_mod_3c",
                  "w23sarcopenia",
                  "w23heces",
                  "Fragilidad")

for (var in discrete_var) {
  metadata_ft[[var]] <- as.factor(metadata_ft[[var]])
}

# Check classes
str(metadata_ft)

# Number of cases
dim(metadata_ft)

# Get summary
summary(metadata_ft)

# Plot discrete variables
for (var in discrete_var) {
  create_bar_plot(metadata_ft, var, title = paste("Bar plot for", var),
                  paste0(work_path,"/plots/bar_plot_",var,".png"))
}

#Plot continuous variables
continuous_var <- c("w23age",
                   "w23medas",
                   "w23energy",
                   "w23weight",
                   "w23height",
                   "w23bmi",
                   "w23sppbscore")

for (var in continuous_var) {
  create_density_plot(metadata_ft, var, title = paste("Density plot for", var),
                      paste0(work_path,"/plots/density_plot_",var,".png"))
}

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

# Save metadatos
save_path <- paste0(work_path,"/data/metadatos.RData")
save(metadata_ft, file = save_path)

load(paste0(work_path, "/data/metadatos.RData"))

# ==================================
# Test continuos variables: FT vs NFT 
# ==================================
cor_matrix <- cor(metadata_ft[, c("w23age", "w23medas", "w23energy",
                                  "w23weight", "w23height", "w23bmi",
                                  "w23sppbscore")], 
                  method = "spearman", use = "complete.obs")
Heatmap(cor_matrix, 
        name = "Spearman correlation", 
        col = colorRampPalette(c("white", "blue", "red"))(100),
        show_row_dend = FALSE, 
        column_names_gp = gpar(fontsize = 6), 
        row_names_gp = gpar(fontsize = 6), 
        # cluster_columns = FALSE,
        # cluster_rows = FALSE,
        clustering_distance_rows = "euclidean", 
        clustering_distance_columns = "euclidean",
        use_raster = FALSE)

# Age vs frailty
ggplot(metadata_ft, aes(x = w23age, fill = Fragilidad)) +
  geom_density(alpha = 0.4) + 
  theme_minimal()

summary(metadata_ft[metadata_ft$Fragilidad == "FT", "w23age"])
summary(metadata_ft[metadata_ft$Fragilidad == "NFT", "w23age"])

check_normal_and_test(metadata_ft, "w23age", "Fragilidad")
#Wilcoxon rank sum test with continuity correction
#data:  df[[var_col]] by df[[group_col]]
#W = 6477.5, p-value = 4.618e-08
#alternative hypothesis: true location shift is not equal to 0

# medas vs frailty
ggplot(metadata_ft, aes(x = w23medas, fill = Fragilidad)) +
  geom_density(alpha = 0.4) + 
  theme_minimal()
check_normal_and_test(metadata_ft, "w23medas", "Fragilidad")
#Wilcoxon rank sum test with continuity correction
#W = 2677, p-value = 0.002641
#alternative hypothesis: true location shift is not equal to 0

# energy vs frailty
ggplot(metadata_ft, aes(x = w23energy, fill = Fragilidad)) +
  geom_density(alpha = 0.4) + 
  theme_minimal()
check_normal_and_test(metadata_ft, "w23energy", "Fragilidad")
#Welch Two Sample t-test
#t = -3.5316, df = 96.345, p-value = 0.0006357
#alternative hypothesis: true difference in means between group FT and group NFT is not equal to 0
#95 percent confidence interval:
#-225.91428  -63.34243
#sample estimates:
#mean in group FT mean in group NFT 
#1790.200          1934.828

# Weight vs frailty
ggplot(metadata_ft, aes(x = w23weight, fill = Fragilidad)) +
  geom_density(alpha = 0.4) + 
  theme_minimal()
check_normal_and_test(metadata_ft, "w23weight", "Fragilidad")
#Wilcoxon rank sum test with continuity correction
#W = 3878, p-value = 0.188
#alternative hypothesis: true location shift is not equal to 0


# Height vs frailty
ggplot(metadata_ft, aes(x = w23height, fill = Fragilidad)) +
  geom_density(alpha = 0.4) + 
  theme_minimal()
check_normal_and_test(metadata_ft, "w23height", "Fragilidad")
#Wilcoxon rank sum test with continuity correction
#W = 2216, p-value = 1.65e-08
#alternative hypothesis: true location shift is not equal to 0

# BMI vs frailty
ggplot(metadata_ft, aes(x = w23bmi, fill = Fragilidad)) +
  geom_density(alpha = 0.4) + 
  theme_minimal()
check_normal_and_test(metadata_ft, "w23bmi", "Fragilidad")
#Wilcoxon rank sum test with continuity correction
#W = 5672.5, p-value = 0.0007527
#alternative hypothesis: true location shift is not equal to 0

#w23sppbscore vs frailty
ggplot(metadata_ft, aes(x = w23sppbscore, fill = Fragilidad)) +
  geom_density(alpha = 0.4) + 
  theme_minimal()
check_normal_and_test(metadata_ft, "w23sppbscore", "Fragilidad")
#Wilcoxon rank sum test with continuity correction
#W = 415.5, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0
