#renv::install("DESeq2")
# Load libraries
library(DESeq2)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)

# Working directory
work_path <- getwd()

# Load functions
source(paste0(work_path, "/functions/filter_valids.R"))

# Load metadata
load(paste0(work_path, "/data/metadata.RData"))

# Path to combined_peptide.tsv
fragpipe_path <- paste0("/mnt/d/proteomica/fragilidad/datos/ProteinIdent/",
                        "MSfraggerSPbacteria/DatosProtIdentCombinados/")
# Load combined_peptide.tsv
peptides  <- read.delim(paste0(fragpipe_path,"combined_peptide.tsv"),
                        row.names = 1) %>%
  as.data.frame() %>% 
  select(ends_with("Intensity")) %>%
  select(-contains("MaxLFQ")) %>%  # Unselect MaxLFQ
  select(contains(metadata_ft$replicate))

# Number of peptides
nrow(peptides)
# 29259

# Number of intensities
sum(peptides != 0)
# 312065

# Change colnames to include the group each replicate belongs to
replicates_ft <- metadata_ft[metadata_ft$frailty == "FT", "replicate"]
replicates_nft <- metadata_ft[metadata_ft$frailty == "NFT", "replicate"]
for (i in seq_along(colnames(peptides))) {
  colname <- colnames(peptides)[i]  # Get column name
  
  if (substr(colname, 5, 7) %in% replicates_ft) {
    colnames(peptides)[i] <- paste0("FT", substr(colname, 4, nchar(colname)))
  } else {
    colnames(peptides)[i] <- paste0("NFT", substr(colname, 4, nchar(colname)))
  }
}

# Number of peptides with all NA
peptides[peptides==0] <- NA
peptides_to_keep <- apply(peptides, 1, function(x) !all(is.na(x)))
sum(!peptides_to_keep)
# 

# Filter peptides with all NA
peptides <- peptides[peptides_to_keep, ]

# Number of peptides after filtering those with all NA
nrow(peptides)
# 

# Number of intensities after filtering those with all NA
table(!is.na(peptides))
# 

# Number of peptides detected per sample
num_pep_sample <- apply(peptides, 2, function(x) sum(!is.na(x)))
num_pep_sample <- as.data.frame(num_pep_sample)
colnames(num_pep_sample) <- "Peptides"
num_pep_sample$Sample <- sapply(strsplit(rownames(num_pep_sample), "\\."), `[`, 1)
rownames(num_pep_sample) <- NULL
# Add frailty group
num_pep_sample$Frailty <- sapply(strsplit(num_pep_sample$Sample, "_"), `[`, 1)

# Summary total peptides detected per sample
summary(num_pep_sample$Peptides)
# 

# Density plot number of peptides per sample in FT vs NFT
ggplot(num_pep_sample, aes(x = Peptides, fill = Frailty)) +
  geom_density(alpha = 0.4, linewidth = 0.2) +
  labs(title = "Number of peptides per sample", x = "Number of peptides", 
       y = "Density") +
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue")) +
  theme_minimal()

# Number of missing per sample
num_na_sample <- apply(peptides, 2, function(x) sum(is.na(x)))
num_na_sample <- as.data.frame(num_na_sample)
colnames(num_na_sample) <- "NAs"
num_na_sample$Sample <- sapply(strsplit(rownames(num_na_sample), "\\."), `[`, 1)
rownames(num_na_sample) <- NULL
# Add frailty group
num_na_sample$Frailty <- sapply(strsplit(num_na_sample$Sample, "_"), `[`, 1)

# Density plot number of missing per sample in FT vs NFT
ggplot(num_na_sample, aes(x = NAs, fill = Frailty)) +
  geom_density(alpha = 0.4, linewidth = 0.2) +
  labs(title = "Number of NAs per sample", x = "Number of NAs", y = "Density") +
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue")) +
  theme_minimal()

# Make binary assay (1 = valid value, 0 = missing value)
missval <- ifelse(is.na(peptides), 0, 1)
colnames(missval) <- sapply(strsplit(colnames(missval), "\\."), `[`, 1)

# Convert matrix to long format dataframe
df <- melt(missval)

# Rename columns
colnames(df) <- c("peptide", "sample", "intensity")
df$sample <- as.character(df$sample)


# Add frailty group
df$frailty <- sapply(strsplit(df$sample, "_"), `[`, 1)

# Plot heatmap
ggplot(df, aes(sample, peptide, fill = intensity)) + 
  geom_tile() +
  facet_grid(. ~ frailty, scales = "free_x", space = "free_x") +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 2, angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_gradient(low = "white", high = "black")

# Complex Heatmap
# Clustering rows and columns by binary distance
# Column split by frailty groups
Heatmap(missval, name = "NA", col = c("white", "black"),
             show_row_dend = FALSE, 
             column_names_gp = gpar(fontsize = 4),
             show_row_names = FALSE, 
             # cluster_columns = FALSE,
             # cluster_rows = FALSE,
             column_split = metadata_ft$Fragilidad, 
             clustering_distance_rows = "binary", 
             clustering_distance_columns = "binary", 
             use_raster = FALSE)


# Annotation of columns by frailty group
column_ann <- HeatmapAnnotation(frailty = metadata_ft$Fragilidad,
                                col = list(frailty = c("FT" = "red", 
                                                       "NFT" = "blue")))
Heatmap(missval, name = "NA", col = c("white", "black"), show_row_dend = FALSE, 
        column_names_gp = gpar(fontsize = 4), show_row_names = FALSE, 
        # cluster_columns = FALSE,
        # cluster_rows = FALSE,
        clustering_distance_rows = "binary", 
        clustering_distance_columns = "binary",
        top_annotation = column_ann, 
        use_raster = FALSE)


# ==================================
# Filter by missing
# ==================================
# Frailty column as grouping
cond_options <- table(metadata_ft[,"frailty"]) %>% as.data.frame()

# Filtering by minimum count in at least 1 condition
cond_opts <- cond_options$Var1
cond_count <- cond_options$Freq * 0.5

pep_exp <- filter_valids(peptides,
                         conditions = cond_opts,
                         min_count = cond_count,
                         at_least_one = T)

save(pep_exp, file = paste0(work_path, "/psva/data/lfq/pep_exp.RData"))

# Number of peptides after filtering by sample proportion
nrow(pep_exp)
# 699

# Number of intensities
pep_exp[pep_exp==0] <- NA
table(!is.na(pep_exp))
# FALSE    TRUE 
# 48147   93750 

# Number of peptides detected per sample
num_pep_sample <- apply(pep_exp, 2, function(x) sum(!is.na(x)))
num_pep_sample <- as.data.frame(num_pep_sample)
colnames(num_pep_sample) <- "Peptides"
num_pep_sample$Sample <- sapply(strsplit(rownames(num_pep_sample), "\\."), `[`, 1)
rownames(num_pep_sample) <- NULL
# Add frailty group
num_pep_sample$Frailty <- sapply(strsplit(num_pep_sample$Sample, "_"), `[`, 1)

# Summary total peptides detected per sample
summary(num_pep_sample$Peptides)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 55.0   380.5   500.0   461.8   569.0   678.0  

# Density plot number of peptides per sample in FT vs NFT
ggplot(num_pep_sample, aes(x = Peptides, fill = Frailty)) +
  geom_density(alpha = 0.4, linewidth = 0.2) +
  labs(title = "Number of peptides per sample", x = "Number of peptides", 
       y = "Density") +
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue")) +
  theme_minimal()

# Number of missing per sample
num_na_sample <- apply(pep_exp, 2, function(x) sum(is.na(x)))
num_na_sample <- as.data.frame(num_na_sample)
colnames(num_na_sample) <- "NAs"
num_na_sample$Sample <- sapply(strsplit(rownames(num_na_sample), "\\."), `[`, 1)
rownames(num_na_sample) <- NULL
# Add frailty group
num_na_sample$Frailty <- sapply(strsplit(num_na_sample$Sample, "_"), `[`, 1)

# Density plot number of missing per sample in FT vs NFT
ggplot(num_na_sample, aes(x = NAs, fill = Frailty)) +
  geom_density(alpha = 0.4, linewidth = 0.2) +
  labs(title = "Number of NAs per sample", x = "Number of NAs", y = "Density") +
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue")) +
  theme_minimal()

# Make binary assay (1 = valid value, 0 = missing value)
missval <- ifelse(is.na(pep_exp), 0, 1)
colnames(missval) <- sapply(strsplit(colnames(missval), "\\."), `[`, 1)

# Convert matrix to long format dataframe
df <- melt(missval)

# Rename columns
colnames(df) <- c("peptide", "sample", "intensity")
df$sample <- as.character(df$sample)

# Add frailty group
df$frailty <- sapply(strsplit(df$sample, "_"), `[`, 1)

# Plot heatmap
ggplot(df, aes(sample, peptide, fill = intensity)) + 
  geom_tile() +
  facet_grid(. ~ frailty, scales = "free_x", space = "free_x") +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 2, angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_gradient(low = "white", high = "black")

# Complex Heatmap
# Clustering rows and columns by binary distance
# Column split by frailty groups
Heatmap(missval, name = "NA", col = c("white", "black"),
        show_row_dend = FALSE, 
        column_names_gp = gpar(fontsize = 4),
        show_row_names = FALSE, 
        # cluster_columns = FALSE,
        # cluster_rows = FALSE,
        column_split = metadata_ft$Fragilidad, 
        clustering_distance_rows = "binary", 
        clustering_distance_columns = "binary", 
        use_raster = FALSE)


# Annotation of columns by frailty group
column_ann <- HeatmapAnnotation(frailty = metadata_ft$Fragilidad,
                                col = list(frailty = c("FT" = "red", 
                                                       "NFT" = "blue")))
Heatmap(missval, name = "NA", col = c("white", "black"), show_row_dend = FALSE, 
        column_names_gp = gpar(fontsize = 4), show_row_names = FALSE, 
        # cluster_columns = FALSE,
        # cluster_rows = FALSE,
        clustering_distance_rows = "binary", 
        clustering_distance_columns = "binary",
        top_annotation = column_ann, 
        use_raster = FALSE)

# ==================================
# Normalization
# ==================================
# Transform pep_exp to long
pep_exp_long <- pep_exp
pep_exp_long$Peptides <- rownames(pep_exp_long)
rownames(pep_exp_long) <- NULL
pep_exp_long <- melt(pep_exp_long)
colnames(pep_exp_long) <- c("Peptides", "Samples", "Intensities")
pep_exp_long$Samples <- as.character(pep_exp_long$Samples)
pep_exp_long$Samples <- sapply(strsplit(pep_exp_long$Samples, "\\."), `[`, 1)
# Add frailty group
pep_exp_long$Frailty <- sapply(strsplit(pep_exp_long$Samples, "_"), `[`, 1)

summary(pep_exp_long$Intensities)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# 486    34113    60737   147707   120817 76101048    48147 

# Plot intensities vs samples
ggplot(pep_exp_long, aes(x = Samples, y = Intensities, fill = Frailty)) + 
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    legend.position = c(0.95, 0.95),
    legend.background = element_rect(fill = alpha("white", 0.5))) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple")) +
  ylim(0, 250000) +
  facet_wrap(~Frailty, scales = "free_x", ncol = 1)

# Plot intensities vs frailty group
ggplot(pep_exp_long, aes(x = Frailty, y = Intensities, fill = Frailty)) + 
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    legend.position = c(0.95, 0.95),
    legend.background = element_rect(fill = alpha("white", 0.5))) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple")) +
  ylim(0, 250000)

# Intensities density plot
ggplot(pep_exp_long, aes(x = Intensities)) +
  geom_density(fill = "blue", alpha = 0.4) +
  theme_minimal() + 
  xlim(0, 700000)

# Estimate size factor for each sample -> norm_pep (vector)
norm_pep <- estimateSizeFactorsForMatrix(as.matrix(pep_exp)#,
                                         #type = "poscounts"
                                        )
# Normalise dividing each column by its size factor
norm_pep_exp <- sweep(as.matrix(pep_exp), 2, norm_pep, "/")
norm_pep_exp <- as.data.frame(norm_pep_exp)

save(norm_pep_exp, file = paste0(work_path, "/psva/data/lfq/norm_pep_exp.RData"))


# Log2 transform
log_pep_exp <- norm_pep_exp
log_pep_exp[, grep("Intensity", names(norm_pep_exp))] <- lapply(
  norm_pep_exp[, grep("Intensity", names(norm_pep_exp))], 
  function(x) log2(1 + x)
)

# Transform log_pep_exp to long
log_pep_exp_long <- log_pep_exp
log_pep_exp_long$Peptides <- rownames(log_pep_exp_long)
rownames(log_pep_exp_long) <- NULL
log_pep_exp_long <- melt(log_pep_exp_long)
colnames(log_pep_exp_long) <- c("Peptides", "Samples", "Intensities")
log_pep_exp_long$Samples <- as.character(log_pep_exp_long$Samples)
log_pep_exp_long$Samples <- sapply(strsplit(log_pep_exp_long$Samples, "\\."), `[`, 1)
# Add frailty group
log_pep_exp_long$Frailty <- sapply(strsplit(log_pep_exp_long$Samples, "_"), `[`, 1)

summary(log_pep_exp_long$Intensities)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  6.90   14.83   16.07   16.14   17.38   26.52   48147 

# Plot normalised intensities vs samples
ggplot(log_pep_exp_long, aes(x = Samples, y = Intensities, fill = Frailty)) + 
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    legend.position = c(0.95, 0.95),
    legend.background = element_rect(fill = alpha("white", 0.5))) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple")) +
  facet_wrap(~Frailty, scales = "free_x", ncol = 1)

# Plot normalised intensities vs frailty group
ggplot(log_pep_exp_long, aes(x = Frailty, y = Intensities, fill = Frailty)) + 
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    legend.position = c(0.95, 0.95),
    legend.background = element_rect(fill = alpha("white", 0.5))) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple"))

# Normalised intensities density plot
ggplot(log_pep_exp_long, aes(x = Intensities)) +
  geom_density(fill = "blue", alpha = 0.4) +
  theme_minimal()

#================================
# PCA
#================================
# Principal component analysis (PCA) biplot visualizing PC1 and PC2, the two
# components that explain the most variation in the dataset. 
# Data exploration technique that can be used for quality control and pattern
# discover.
#================================
numcond <- length(cond_opts)
seqcond <- 1:numcond
# colours_to_plot <- lacroix_palette("Pamplemousse", n = numcond,
#                                   type = "continuous")

colnames(metadata) <- c("Samples", "Condition")
conditions <- metadata$Condition

pca<- prcomp(t(log_pep_exp), center=T, scale=F)
sampleVals<-data.frame(pca$x)
exprVals<-data.frame(pca$rotation)
PoV <- (pca$sdev^2/sum(pca$sdev^2))*100


coords<-data.frame(sampleVals, Condition = conditions,
                   samplename = rownames(sampleVals))
numPCs <- 1:length(PoV)

for (i in 1:length(PoV)) {
  percent <- paste0("(", round(PoV[i],2), "%)")
  percentNoBrack <- paste0(round(PoV[i],2), "%")
  name <- paste0("PC", i, "per")
  name2 <- paste0("PC",i,"per_short")
  assign(name, percent)
  assign(name2, percentNoBrack)
}

(pcaplot <- ggplot(coords, aes(x = PC1, y = PC2)) +
    stat_ellipse(geom = "polygon", alpha=.2, aes(color=Condition,
                                                 fill=Condition)) +
    geom_point(size=5, aes(colour=Condition, shape=Condition)) + 
    # scale_color_manual(values=c(colours_to_plot)) +
    # scale_fill_manual(values=c(colours_to_plot)) +
    scale_x_continuous(name= paste0("PC1", " ", PC1per))+ # labels depend on selected PCs
    scale_y_continuous(name= paste0("PC2", " ", PC2per))+ theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank(), 
          text=element_text(size=12)) )

#================================
# Hierachical clustering
#================================
# A hierarchical clustering dendrogram visualizing the relationships between
# samples as explained by peptide intensity similarity.
# Euclidian distance and the complete agglomeration hierachical clustering method.
#================================
dd <- dist(log_pep_exp %>% t(), method = "euclidian") # be able to chose distance
hc <- hclust(dd, method = "complete") # be able to chose method
condcolours <- data.frame(Condition = cond_opts, Colour = as.character(colours_to_plot))
condcolours <- merge(metadata, condcolours, by = "Condition")

colnames(log_pep_exp) <- colnames(log_pep_exp) %>% substr(., 1, 7)
dend <- log_pep_exp %>% t() %>% dist(method = "euclidian") %>% 
  hclust(method = "complete") %>% as.dendrogram(hang=0.1) %>%
  set("leaves_pch", 19) %>% 
  set("leaves_col", rainbow(length(labels(.))), order_value = T) %>% 
  set('branches_lwd', 0.7) %>%
  set('labels_cex', 0.7)

#plot(dend)
dend_ggplot <- as.ggdend(dend)
ggplot(dend_ggplot, horiz=T, theme=theme_dendro(),  offset_labels = -10)

dend <- log_pep_exp %>% t() %>% dist(method = "euclidian") %>% 
  hclust(method = "complete") %>% as.dendrogram(hang=0.1) %>%
  set("leaves_pch", 19) %>% 
  set("leaves_col", rainbow(length(labels(.))), order_value = T) %>% 
  set('branches_lwd', 0.7) %>%
  set('labels_cex', 0.7)

