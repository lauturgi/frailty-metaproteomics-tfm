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
  select(contains("MaxLFQ")) %>%  # MaxLFQ intensities
  select(contains(metadata_ft$replicate))


# Change 0 values to NA
peptides[peptides==0] <- NA

# Number of peptides
nrow(peptides)
# 29259

# Number of MaxLFQ intensities
sum(!is.na(peptides))
# 305082

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
peptides_to_keep <- apply(peptides, 1, function(x) !all(is.na(x)))
sum(!peptides_to_keep)
# 10089

# Filter peptides with all NA
peptides <- peptides[peptides_to_keep, ]

# Number of peptides after filtering those with all NA
nrow(peptides)
# 19170

# Number of MaxLFQ intensities after filtering those with all NA
table(!is.na(peptides))
# TRUE    FALSE 
# 3586428  305082 

# Number of peptides per sample
num_pep_sample <- apply(peptides, 2, function(x) sum(!is.na(x)))
num_pep_sample <- as.data.frame(num_pep_sample)
colnames(num_pep_sample) <- "peptides"
num_pep_sample$samples <- sapply(strsplit(rownames(num_pep_sample), "\\."), `[`, 1)
rownames(num_pep_sample) <- NULL
# Add frailty group
num_pep_sample$frailty <- sapply(strsplit(num_pep_sample$samples, "_"), `[`, 1)

# Summary total peptides per sample
summary(num_pep_sample$peptides)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 64     922    1341    1503    1813    5971

# Density plot number of peptides per sample FT vs NFT
ggplot(num_pep_sample, aes(x = peptides, fill = frailty)) +
  geom_density(alpha = 0.4, linewidth = 0.2) +
  labs(title = "Number of peptides per sample", x = "Number of peptides", 
       y = "Density") +
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue")) +
  theme_minimal()

# Number of missing per sample
num_na_sample <- apply(peptides, 2, function(x) sum(is.na(x)))
num_na_sample <- as.data.frame(num_na_sample)
colnames(num_na_sample) <- "NAs"
num_na_sample$samples <- sapply(strsplit(rownames(num_na_sample), "\\."), `[`, 1)
rownames(num_na_sample) <- NULL
# Add frailty group
num_na_sample$frailty <- sapply(strsplit(num_na_sample$samples, "_"), `[`, 1)

# Density plot number of missing per sample FT vs NFT
ggplot(num_na_sample, aes(x = NAs, fill = frailty)) +
  geom_density(alpha = 0.4, linewidth = 0.2) +
  labs(title = "Number of NAs per sample", x = "Number of NAs", y = "Density") +
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue")) +
  theme_minimal()

# Make binary (1 = valid value, 0 = missing value)
missval <- ifelse(is.na(peptides), 0, 1)
colnames(missval) <- sapply(strsplit(colnames(missval), "\\."), `[`, 1)

# Convert dataframe to long format
df <- melt(missval)

# Rename columns
colnames(df) <- c("peptides", "samples", "maxlfq")
df$samples <- as.character(df$samples)


# Add frailty group
df$frailty <- sapply(strsplit(df$samples, "_"), `[`, 1)

# Plot heatmap
ggplot(df, aes(samples, peptides, fill = maxlfq)) + 
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
# Fragilidad column as grouping
cond_options <- table(metadata_ft[,"frailty"]) %>% as.data.frame()

# Filtering by minimum count in at least 1 condition
cond_opts <- cond_options$Var1
cond_count <- cond_options$Freq * 0.5

pep_exp <- filter_valids(peptides,
                         conditions = cond_opts,
                         min_count = cond_count,
                         at_least_one = T)

save(pep_exp, file = paste0(work_path, "/psva/data/maxlfq/pep_exp.RData"))

# Number of peptides after filtering by sample proportion
pep_exp[pep_exp==0] <- NA
nrow(pep_exp)
# 699

# Number of MaxLFQ intensities
table(!is.na(pep_exp))
# FALSE    TRUE 
# 48148   93749 

# Number of peptides per sample
num_pep_sample <- apply(pep_exp, 2, function(x) sum(!is.na(x)))
num_pep_sample <- as.data.frame(num_pep_sample)
colnames(num_pep_sample) <- "peptides"
num_pep_sample$samples <- sapply(strsplit(rownames(num_pep_sample), "\\."), `[`, 1)
rownames(num_pep_sample) <- NULL
# Add frailty group
num_pep_sample$frailty <- sapply(strsplit(num_pep_sample$samples, "_"), `[`, 1)

# Summary total peptides per sample
summary(num_pep_sample$peptides)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 55.0   380.5   500.0   461.8   569.0   678.0  

# Density plot number of peptides per sample FT vs NFT
ggplot(num_pep_sample, aes(x = peptides, fill = frailty)) +
  geom_density(alpha = 0.4, linewidth = 0.2) +
  labs(title = "Number of peptides per sample", x = "Number of peptides", 
       y = "Density") +
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue")) +
  theme_minimal()

# Number of missing per sample
num_na_sample <- apply(pep_exp, 2, function(x) sum(is.na(x)))
num_na_sample <- as.data.frame(num_na_sample)
colnames(num_na_sample) <- "NAs"
num_na_sample$samples <- sapply(strsplit(rownames(num_na_sample), "\\."), `[`, 1)
rownames(num_na_sample) <- NULL
# Add frailty group
num_na_sample$frailty <- sapply(strsplit(num_na_sample$samples, "_"), `[`, 1)

# Density plot number of missing per sample in FT vs NFT
ggplot(num_na_sample, aes(x = NAs, fill = frailty)) +
  geom_density(alpha = 0.4, linewidth = 0.2) +
  labs(title = "Number of NAs per sample", x = "Number of NAs", y = "Density") +
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue")) +
  theme_minimal()

# Make binary (1 = valid value, 0 = missing value)
missval <- ifelse(is.na(pep_exp), 0, 1)
colnames(missval) <- sapply(strsplit(colnames(missval), "\\."), `[`, 1)

# Convert dataframe to long format
df <- melt(missval)

# Rename columns
colnames(df) <- c("peptides", "samples", "maxlfq")
df$samples <- as.character(df$samples)

# Add frailty group
df$frailty <- sapply(strsplit(df$samples, "_"), `[`, 1)

# Plot heatmap
ggplot(df, aes(samples, peptides, fill = maxlfq)) + 
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
pep_exp_long$peptides <- rownames(pep_exp_long)
rownames(pep_exp_long) <- NULL
pep_exp_long <- melt(pep_exp_long)
colnames(pep_exp_long) <- c("peptides", "samples", "maxlfq")
pep_exp_long$samples <- as.character(pep_exp_long$samples)
pep_exp_long$samples <- sapply(strsplit(pep_exp_long$samples, "\\."), `[`, 1)
# Add frailty group
pep_exp_long$frailty <- sapply(strsplit(pep_exp_long$samples, "_"), `[`, 1)

summary(pep_exp_long$maxlfq)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# 486    32564    54325   126635   106158 58777692    48148 

# Plot MaxLFQ vs samples
ggplot(pep_exp_long, aes(x = samples, y = maxlfq, fill = frailty)) + 
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    legend.position = c(0.95, 0.95),
    legend.background = element_rect(fill = alpha("white", 0.5))) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple")) +
  ylim(0, 250000) +
  facet_wrap(~frailty, scales = "free_x", ncol = 1)

# Plot MaxLFQ vs frailty group
ggplot(pep_exp_long, aes(x = frailty, y = maxlfq, fill = frailty)) + 
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    legend.position = c(0.95, 0.95),
    legend.background = element_rect(fill = alpha("white", 0.5))) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple")) +
  ylim(0, 250000)

# MaxLFQ density plot
ggplot(pep_exp_long, aes(x = maxlfq)) +
  geom_density(fill = "blue", alpha = 0.4) +
  theme_minimal() + 
  xlim(0, 700000)

# Log2 transform
log_pep_exp <- pep_exp
log_pep_exp[, grep("Intensity", names(log_pep_exp))] <- lapply(
  log_pep_exp[, grep("Intensity", names(log_pep_exp))], 
  function(x) log2(x)
)

# Transform log_pep_exp to long
log_pep_exp_long <- log_pep_exp
log_pep_exp_long$peptides <- rownames(log_pep_exp_long)
rownames(log_pep_exp_long) <- NULL
log_pep_exp_long <- melt(log_pep_exp_long)
colnames(log_pep_exp_long) <- c("peptides", "samples", "maxlfq")
log_pep_exp_long$samples <- as.character(log_pep_exp_long$samples)
log_pep_exp_long$samples <- sapply(strsplit(log_pep_exp_long$samples, "\\."), `[`, 1)
# Add frailty group
log_pep_exp_long$frailty <- sapply(strsplit(log_pep_exp_long$samples, "_"), `[`, 1)

summary(log_pep_exp_long$maxlfq)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  8.93   14.99   15.73   15.94   16.70   25.81   48148 

# Plot log2MaxLFQ vs samples
ggplot(log_pep_exp_long, aes(x = samples, y = maxlfq, fill = frailty)) + 
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    legend.position = c(0.95, 0.95),
    legend.background = element_rect(fill = alpha("white", 0.5))) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple")) +
  facet_wrap(~frailty, scales = "free_x", ncol = 1)

# Plot log2MaxLFQ vs frailty group
ggplot(log_pep_exp_long, aes(x = frailty, y = maxlfq, fill = frailty)) + 
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    legend.position = c(0.95, 0.95),
    legend.background = element_rect(fill = alpha("white", 0.5))) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple"))

# Log2MaxLFQ density plot
ggplot(log_pep_exp_long, aes(x = maxlfq)) +
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

