#install.packages("ggplot2")
library(ggplot2)

# Working directory
work_path <- "/home/lauturgi/workspace/tfm"

# Metadata directory
metadata_path <- "/mnt/d/proteomica/fragilidad/datos/"
metadata_name <- "Tabla metadatos ENRICA.csv"

# Read metadata
metadata <- read.csv(paste0(metadata_path,metadata_name),sep=";")

# Get frailty cohort cases and columns of interest
metadata_ft <- metadata[metadata$Fragilidad %in% c("FT", "NFT"), -ncol(metadata)]

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
  create_bar_plot(metadata_ft, var, title = paste("Bar plot for", var), paste0(work_path,"/plots/bar_plot_",var,".png"))
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
  create_density_plot(metadata_ft, var, title = paste("Density plot for", var), paste0(work_path,"/plots/density_plot_",var,".png"))
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

