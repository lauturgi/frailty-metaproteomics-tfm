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
source(paste0(work_path, "/functions/build_summarized_experiment.R"))

# Data, plot and output paths
data_path <- paste0(work_path, "/psva/data/")
plot_path <- paste0(work_path, "/psva/plots/")
output_path <- paste0(work_path, "/psva/output/")

# Path to combined_protein.tsv and experiment_annotation.tsv
combined_peptide <- paste0(data_path, "combined_peptide.tsv")
experiment_ann <- paste0(data_path, "experiment_annotation.tsv")

# ==================================
# 1. Preprocessing MaxLFQ intensities
# ==================================

# Make SummarizedExperiment (SE)
se_pep_maxlfq <- build_summarized_experiment(combined_file = combined_peptide,
                                             experiment_ann = experiment_ann,
                                             level = "peptide")

# Check class of se_pep_maxlfq object
class(se_pep_maxlfq)

# Check log2, exp, lfq_type and level
metadata(se_pep_maxlfq)
#$log2transform
# FALSE
#$lfq_type
# "MaxLFQ"
#$level
# "peptide"

# Check number of rows and columns
dim(se_pep_maxlfq)
# 29259  274

# ==================================
# 1.1. Filtering frailty cohort
# ==================================

# Load metadata
load(paste0(work_path, "/data/metadata.RData"))

# Align metadata with same order of `replicate` in colData(se_pep_maxlfq)
replicate_se <- colData(se_pep_maxlfq)$replicate
replicate_metadata <- metadata_ft$replicate
metadata_ft <- metadata_ft[match(replicate_se, replicate_metadata), ]
metadata_ft <- metadata_ft[!is.na(metadata_ft$frailty), ]
# Keep only data from FT cohort
se_pep_maxlfq <- se_pep_maxlfq[, colData(se_pep_maxlfq)$replicate %in%
                                 metadata_ft$replicate]
unique(colData(se_pep_maxlfq)$replicate)

# Check number of rows and columns
dim(se_pep_maxlfq)
# 29259  203

# Get intensities matrix from SE and convert to dataframe
pep_maxlfq  <- as.data.frame(assay(se_pep_maxlfq))

# Change colnames to include the group each replicate belongs to
replicates_ft <- metadata_ft[metadata_ft$frailty == "FT", "replicate"]
replicates_nft <- metadata_ft[metadata_ft$frailty == "NFT", "replicate"]
for (i in seq_along(colnames(pep_maxlfq))) {
  colname <- colnames(pep_maxlfq)[i]  # Get column name
  
  if (substr(colname, 5, 7) %in% replicates_ft) {
    colnames(pep_maxlfq)[i] <- paste0("FT", substr(colname, 4, nchar(colname)))
  } else {
    colnames(pep_maxlfq)[i] <- paste0("NFT", substr(colname, 4, nchar(colname)))
  }
}

# ==================================
# 1.2. Filtering by missing
# ==================================

# Number of peptides
nrow(pep_maxlfq)
# 29259

# Number of NA values in assay
na_table <- table(is.na(pep_maxlfq))
na_table
# FALSE    TRUE 
# 305082 5634495 

# NA proportion
prop_na <- na_table["TRUE"] / sum(na_table)
prop_na
# 0.9486357

# Peptides with at least one intensity value
peptides_to_keep <- apply(pep_maxlfq, 1, function(x) !all(is.na(x)))

# Number of peptides to keep
sum(peptides_to_keep)
# 19170

# Number of peptides with all NA
sum(!peptides_to_keep)
# 10089

# Filter peptides with all NA
pep_maxlfq <- pep_maxlfq[peptides_to_keep, ]

# Check number of peptides after filtering those wil all NA
nrow(pep_maxlfq)
# 19170

# Number of NA values in assay after filtering those with all NA
na_table <- table(is.na(pep_maxlfq))
na_table
# FALSE    TRUE 
# 305082 3586428 

# NA proportion
prop_na <- na_table["TRUE"] / sum(na_table)
prop_na
# 0.9216032

# Number of peptides detected per sample
num_pep_sample <- apply(pep_maxlfq, 2, function(x) sum(!is.na(x)))
num_pep_sample <- as.data.frame(num_pep_sample)
colnames(num_pep_sample) <- "peptides"
num_pep_sample$sample <- rownames(num_pep_sample)
rownames(num_pep_sample) <- NULL
# Add frailty group
num_pep_sample$frailty <- sapply(strsplit(num_pep_sample$sample, "_"), `[`, 1)

# Summary total peptides detected per sample
summary(num_pep_sample$peptides)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 64     922    1341    1503    1813    5971

# Density plot number of peptides per sample FT vs NFT
p <- ggplot(num_pep_sample, aes(x = peptides, fill = frailty)) +
  geom_density(linewidth = 0.4, alpha = 0.3) +
  labs(title = "Distribution of peptides per sample",
       x = "Number of peptides", 
       y = "Density") +
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

save_path <- paste0(plot_path, "preprocessing/density_plot_maxlfq_num_pep_",
                    "sample.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- ggplot(num_pep_sample, aes(x = frailty, y = peptides, fill = frailty)) + 
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Frailty group", 
       y = "Number of peptides") +
  coord_cartesian(y = c(0, 7200)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail"))

save_path <- paste0(plot_path, "preprocessing/boxplot_maxlfq_num_pep_",
                    "sample.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Mean number of peptides in frail
mean(num_pep_sample[num_pep_sample$frailty == "FT", "peptides"])
# 1362.569
# Standard deviation number of peptides in frail
sd(num_pep_sample[num_pep_sample$frailty == "FT", "peptides"])
# 872.471

# Mean number of peptides in non-frail
mean(num_pep_sample[num_pep_sample$frailty == "NFT", "peptides"])
# 1568.949
# Standard deviation number of peptides in non-frail
sd(num_pep_sample[num_pep_sample$frailty == "NFT", "peptides"])
# 975.23

# Wilcox test number of peptides FT vs NFT
wilcox.test(peptides ~ frailty, data = num_pep_sample)
# Wilcoxon rank sum test with continuity correction
# data:  peptides by frailty
# W = 3844, p-value = 0.101
# alternative hypothesis: true location shift is not equal to 0

# Make binary (1 = valid value, 0 = missing value)
missval <- ifelse(is.na(pep_maxlfq), 0, 1)

# Convert dataframe to long format
df <- melt(missval)

# Rename columns
colnames(df) <- c("peptides", "sample", "maxlfq")
df$sample <- as.character(df$sample)

# Add frailty group
df$frailty <- sapply(strsplit(df$sample, "_"), `[`, 1)

# Plot heatmap
p <- ggplot(df, aes(sample, peptides, fill = maxlfq)) + 
  geom_tile() +
  facet_grid(. ~ frailty, scales = "free_x", space = "free_x") +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 2, angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_gradient(low = "white", high = "black")

save_path <- paste0(plot_path,"preprocessing/heatmap_maxlfq_na.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Complex Heatmap
# Clustering rows and columns by binary distance
# Column split by frailty groups
p <- Heatmap(missval, name = "NA", col = c("white", "black"),
             show_row_dend = FALSE, 
             column_names_gp = gpar(fontsize = 4),
             show_row_names = FALSE, 
             # cluster_columns = FALSE,
             # cluster_rows = FALSE,
             column_split = metadata_ft$frailty, 
             clustering_distance_rows = "binary", 
             clustering_distance_columns = "binary", 
             use_raster = FALSE)

save_path <- paste0(plot_path,"preprocessing/heatmap_maxlfq_clust_na.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
draw(p)
dev.off()

# Annotation of columns by frailty group
column_ann <- HeatmapAnnotation(
  frailty = metadata_ft$frailty,
  col = list(frailty = c("FT" = "red", "NFT" = "blue")),
  show_annotation_name = FALSE,
  show_legend = FALSE
  )

p <- Heatmap(missval, col = c("white", "black"),
            show_row_dend = FALSE, 
            column_names_gp = gpar(fontsize = 4),
            show_row_names = FALSE, 
            show_column_names = FALSE,
            # cluster_columns = FALSE,
            # cluster_rows = FALSE,
            clustering_distance_rows = "binary", 
            clustering_distance_columns = "binary",
            top_annotation = column_ann, 
            use_raster = FALSE,
            show_heatmap_legend = FALSE
            )

# Heatmap legend
legend_1 <- Legend(at = c(0, 1),
                   labels = c("Missing value", "Valid value"),
                   legend_gp = gpar(fill = c("white", "black")),
                   title = NULL)
# Annotation legend
legend_2 <- Legend(at = c("FT", "NFT"),
                   labels = c("Frail", "Non-Frail"),
                   legend_gp = gpar(fill = c("red", "blue")),
                   title = NULL)

# Combine legends
legends <- packLegend(legend_2, legend_1, direction = "horizontal")

# Plot heatmap
save_path <- paste0(plot_path,"preprocessing/heatmap_maxlfq_clust_na_ann.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
draw(p, heatmap_legend_list = legends, heatmap_legend_side = "top")
dev.off()

# Frailty column as grouping
cond_options <- table(metadata_ft[,"frailty"]) %>% as.data.frame()

# Filtering by minimum count in at least 1 condition
cond_opts <- cond_options$Var1
cond_count <- cond_options$Freq * 0.5

pep_exp <- filter_valids(pep_maxlfq,
                         conditions = cond_opts,
                         min_count = cond_count,
                         at_least_one = T)

save(pep_exp, file = paste0(work_path, "/psva/data/pep_exp.RData"))

# Number of peptides after filtering by sample proportion
nrow(pep_exp)
# 699

# Number of NA values in assay after filtering by sample proportion
na_table <- table(is.na(pep_exp))
na_table
# FALSE  TRUE 
# 93749 48148  

# NA proportion
prop_na <- na_table["TRUE"] / sum(na_table)
prop_na
# 0.3393165

# Number of peptides detected per sample
num_pep_sample <- apply(pep_exp, 2, function(x) sum(!is.na(x)))
num_pep_sample <- as.data.frame(num_pep_sample)
colnames(num_pep_sample) <- "peptides"
num_pep_sample$sample <- rownames(num_pep_sample)
rownames(num_pep_sample) <- NULL
# Add frailty group
num_pep_sample$frailty <- sapply(strsplit(num_pep_sample$sample, "_"), `[`, 1)

# Summary total peptides detected per sample
summary(num_pep_sample$peptides)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 55.0   380.5   500.0   461.8   569.0   678.0  

# Density plot number of peptides per sample FT vs NFT
p <- ggplot(num_pep_sample, aes(x = peptides, fill = frailty)) +
  geom_density(linewidth = 0.4, alpha = 0.3) +
  labs(title = "Distribution of peptides per sample",
       x = "Number of peptides", 
       y = "Density") +
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

save_path <- paste0(plot_path, "preprocessing/density_plot_maxlfq_num_pep_",
                    "sample_filt.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- ggplot(num_pep_sample, aes(x = frailty, y = peptides, fill = frailty)) + 
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Frailty group", 
       y = "Number of peptides") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail"))

save_path <- paste0(plot_path, "preprocessing/boxplot_maxlfq_num_pep_sample_",
                    "filt.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Mean number of peptides in frail
mean(num_pep_sample[num_pep_sample$frailty == "FT", "peptides"])
# 437.3385
# Standard deviation number of peptides in frail
sd(num_pep_sample[num_pep_sample$frailty == "FT", "peptides"])
# 131.0829

# Mean number of peptides in non-frail
mean(num_pep_sample[num_pep_sample$frailty == "NFT", "peptides"])
# 473.3478
# Standard deviation number of peptides in non-frail
sd(num_pep_sample[num_pep_sample$frailty == "NFT", "peptides"])
# 147.0389

# Wilcox test number of peptides FT vs NFT
wilcox.test(peptides ~ frailty, data = num_pep_sample)
# Wilcoxon rank sum test with continuity correction
# data:  peptides by frailty
# W = 3571, p-value = 0.01932
# alternative hypothesis: true location shift is not equal to 0

# Make binary
missval <- ifelse(is.na(pep_exp), 0, 1)

# Convert dataframe to long format
df <- melt(missval)

# Rename columns
colnames(df) <- c("peptides", "sample", "maxlfq")
df$sample <- as.character(df$sample)

# Add frailty group
df$frailty <- sapply(strsplit(df$sample, "_"), `[`, 1)

# Plot heatmap
p <- ggplot(df, aes(sample, peptides, fill = maxlfq)) + 
  geom_tile() +
  facet_grid(. ~ frailty, scales = "free_x", space = "free_x") +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 2, angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_gradient(low = "white", high = "black")

save_path <- paste0(plot_path, "preprocessing/heatmap_maxlfq_na_filt.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Complex Heatmap
# Clustering rows and columns by binary distance
# Column split by frailty groups
p <- Heatmap(missval, name = "NA", col = c("white", "black"),
        show_row_dend = FALSE, 
        column_names_gp = gpar(fontsize = 4),
        show_row_names = FALSE, 
        # cluster_columns = FALSE,
        # cluster_rows = FALSE,
        column_split = metadata_ft$Fragilidad, 
        clustering_distance_rows = "binary", 
        clustering_distance_columns = "binary", 
        use_raster = FALSE)
save_path <- paste0(plot_path, "preprocessing/heatmap_maxlfq_clust_na_filt_",
                    "min.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
draw(p)
dev.off()

# Annotation of columns by frailty group
column_ann <- HeatmapAnnotation(
  frailty = metadata_ft$frailty,
  col = list(frailty = c("FT" = "red", "NFT" = "blue")),
  show_annotation_name = FALSE,
  show_legend = FALSE
  )

p <- Heatmap(missval, col = c("white", "black"),
             show_row_dend = FALSE, 
             column_names_gp = gpar(fontsize = 4),
             show_row_names = FALSE, 
             show_column_names = FALSE,
             # cluster_columns = FALSE,
             # cluster_rows = FALSE,
             clustering_distance_rows = "binary", 
             clustering_distance_columns = "binary",
             top_annotation = column_ann, 
             use_raster = FALSE,
             show_heatmap_legend = FALSE)

# Heatmap legend
legend_1 <- Legend(at = c(0, 1),
                   labels = c("Missing value", "Valid value"),
                   legend_gp = gpar(fill = c("white", "black")),
                   title = NULL)
# Annotation legend
legend_2 <- Legend(at = c("FT", "NFT"),
                   labels = c("Frail", "Non-Frail"),
                   legend_gp = gpar(fill = c("red", "blue")),
                   title = NULL)

# Combine legends
legends <- packLegend(legend_2, legend_1, direction = "horizontal")

# Plot heatmap
save_path <- paste0(plot_path,"preprocessing/heatmap_maxlfq_clust_na_ann_",
                    "filt.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
draw(p, heatmap_legend_list = legends, heatmap_legend_side = "top")
dev.off()

# ==================================
# 1.3. Normalization
# ==================================
# Transform pep_exp to long
pep_exp_long <- pep_exp
pep_exp_long$peptides <- rownames(pep_exp_long)
rownames(pep_exp_long) <- NULL
pep_exp_long <- melt(pep_exp_long)
colnames(pep_exp_long) <- c("peptides", "sample", "maxlfq")
pep_exp_long$sample <- as.character(pep_exp_long$sample)
# Add frailty group
pep_exp_long$frailty <- sapply(strsplit(pep_exp_long$sample, "_"), `[`, 1)

summary(pep_exp_long$maxlfq)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# 486    32564    54325   126635   106158 58777692    48148 

sd(pep_exp_long$maxlfq, na.rm = TRUE)
# 573228.6

# Plot maxlfq vs samples
p <- ggplot(pep_exp_long, aes(x = sample, y = maxlfq, fill = frailty)) + 
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Samples", 
       y = "Peptide MaxLFQ intensities") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),  #element_text(angle = 90, size = 6),
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail")) +
  #ylim(0, 250000) +
  facet_wrap("frailty", scales = "free_x", ncol = 1)

save_path <- paste0(plot_path,"preprocessing/boxplot_maxlfq_samples.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Plot maxlfq vs frailty group
p <- ggplot(pep_exp_long, aes(x = frailty, y = maxlfq, fill = frailty)) + 
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Frailty group", 
       y = "Peptide MaxLFQ intensities") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 8),
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail"))

save_path <- paste0(plot_path,"preprocessing/boxplot_maxlfq_groups.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Intensities density plot
p <- ggplot(pep_exp_long, aes(x = maxlfq)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  theme_minimal(base_size = 12) + 
  labs(title = "Distribution of peptide MaxLFQ intensities",
       x = "Peptide MaxLFQ intensities", 
       y = "Density")

save_path <- paste0(plot_path,"preprocessing/density_plot_maxlfq.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Log2 transform
log_pep_exp <- pep_exp
log_pep_exp[sapply(log_pep_exp, is.numeric)] <- lapply(
  log_pep_exp[sapply(log_pep_exp, is.numeric)],
  function(x) log2(x)
)

# Transform log_pep_exp to long
log_pep_exp_long <- log_pep_exp
log_pep_exp_long$peptides <- rownames(log_pep_exp_long)
rownames(log_pep_exp_long) <- NULL
log_pep_exp_long <- melt(log_pep_exp_long)
colnames(log_pep_exp_long) <- c("peptides", "sample", "maxlfq")
log_pep_exp_long$sample <- as.character(log_pep_exp_long$sample)
# Add frailty group
log_pep_exp_long$frailty <- sapply(strsplit(log_pep_exp_long$sample, "_"), `[`, 1)

summary(log_pep_exp_long$maxlfq)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  8.92   14.99   15.73   15.94   16.70   25.81   48148 

sd(log_pep_exp_long$maxlfq, na.rm = TRUE)
# 1.379822

# Plot log2maxlfq vs samples
p <- ggplot(log_pep_exp_long, aes(x = sample, y = maxlfq, fill = frailty)) + 
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Samples", 
       y = "Peptide log2 MaxLFQ intensity") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x =  element_blank(), #element_text(angle = 90, size = 6),
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail")) +
  facet_wrap("frailty", scales = "free_x", ncol = 1)

save_path <- paste0(plot_path,"preprocessing/boxplot_maxlfq_samples_log2.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Plot log2maxlfq vs frailty group
p <- ggplot(log_pep_exp_long, aes(x = frailty, y = maxlfq, fill = frailty)) + 
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Frailty group", 
       y = "Peptide log2 MaxLFQ intensity") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 8),
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail"))

save_path <- paste0(plot_path,"preprocessing/boxplot_maxlfq_groups_log2.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Log2maxlfq density plot
p <- ggplot(log_pep_exp_long, aes(x = maxlfq)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  labs(title = "Distribution of peptide log2 MaxLFQ intensities",
      x = "Peptide log2 MaxLFQ intensity", 
      y = "Density") +
  theme_minimal(base_size = 12) +
  coord_cartesian(xlim = c(9, 26), ylim = c(0, 0.35))

save_path <- paste0(plot_path,"preprocessing/density_plot_maxlfq_log2.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# ==================================
# 2. Preprocessing Top-N intensities
# ==================================

# Make SummarizedExperiment (SE)
se_pep_lfq <- build_summarized_experiment(combined_file = combined_peptide,
                                          experiment_ann = experiment_ann,
                                          lfq_type = "LFQ",
                                          level = "peptide")

# Check class of se_pep_lfq object
class(se_pep_lfq)

# Check log2, exp, lfq_type and level
metadata(se_pep_lfq)
#$log2transform
# FALSE
#$lfq_type
# "LFQ"
#$level
# "peptide"

# Check number of rows and columns
dim(se_pep_lfq)
# 29259  274

# ==================================
# 2.1. Filtering frailty cohort
# ==================================

# Keep only data from FT cohort
se_pep_lfq <- se_pep_lfq[, colData(se_pep_lfq)$replicate %in%
                           metadata_ft$replicate]
unique(colData(se_pep_lfq)$replicate)

# Check number of rows and columns
dim(se_pep_lfq)
# 29259  203

# Get intensities matrix from SE and convert to dataframe
pep_lfq  <- as.data.frame(assay(se_pep_lfq))

# Change colnames to include the group each replicate belongs to
for (i in seq_along(colnames(pep_lfq))) {
  colname <- colnames(pep_lfq)[i]  # Get column name
  
  if (substr(colname, 5, 7) %in% replicates_ft) {
    colnames(pep_lfq)[i] <- paste0("FT", substr(colname, 4, nchar(colname)))
  } else {
    colnames(pep_lfq)[i] <- paste0("NFT", substr(colname, 4, nchar(colname)))
  }
}

# ==================================
# 2.2. Filtering by missing
# ==================================

# Number of peptides
nrow(pep_lfq)
# 29259

# Number of NA values in assay
na_table <- table(is.na(pep_lfq))
na_table
# FALSE    TRUE 
# 312065 5627512 

# NA proportion
prop_na <- na_table["TRUE"] / sum(na_table)
prop_na
# 0.9474601

# Peptides with at least one intensity value
peptides_to_keep <- apply(pep_lfq, 1, function(x) !all(is.na(x)))

# Number of peptides to keep
sum(peptides_to_keep)
# 25742

# Number of peptides with all NA
sum(!peptides_to_keep)
# 3517

# Filter peptides with all NA
pep_lfq <- pep_lfq[peptides_to_keep, ]

# Number of peptides after filtering those with all NA
nrow(pep_lfq)
# 25742

# Number of NA values in assay after filtering those with all NA
na_table <- table(is.na(pep_lfq))
na_table
#  FALSE    TRUE 
# 312065 4913561 

# NA proportion
prop_na <- na_table["TRUE"] / sum(na_table)
prop_na
# 0.9402818

# Number of peptides detected per sample
num_pep_sample <- apply(pep_lfq, 2, function(x) sum(!is.na(x)))
num_pep_sample <- as.data.frame(num_pep_sample)
colnames(num_pep_sample) <- "peptides"
num_pep_sample$sample <- rownames(num_pep_sample)
rownames(num_pep_sample) <- NULL
# Add frailty group
num_pep_sample$frailty <- sapply(strsplit(num_pep_sample$sample, "_"), `[`, 1)

# Summary total peptides detected per sample
summary(num_pep_sample$peptides)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 64.0   934.5  1360.0  1537.3  1872.0  7152.0

# Density plot number of peptides per sample FT vs NFT
p <- ggplot(num_pep_sample, aes(x = peptides, fill = frailty)) +
  geom_density(linewidth = 0.4, alpha = 0.3) +
  labs(title = "Distribution of peptides per sample",
       x = "Number of peptides", 
       y = "Density") +
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

save_path <- paste0(plot_path, "preprocessing/density_plot_lfq_num_pep_",
                    "sample.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- ggplot(num_pep_sample, aes(x = frailty, y = peptides, fill = frailty)) + 
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Frailty group", 
       y = "Number of peptides") +
  coord_cartesian(y = c(0, 7200)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail"))

save_path <- paste0(plot_path, "preprocessing/boxplot_lfq_num_pep_",
                    "sample.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Mean number of peptides in frail
mean(num_pep_sample[num_pep_sample$frailty == "FT", "peptides"])
# 1400.969
# Standard deviation number of peptides in frail
sd(num_pep_sample[num_pep_sample$frailty == "FT", "peptides"])
# 959.4652

# Mean number of peptides in non-frail
mean(num_pep_sample[num_pep_sample$frailty == "NFT", "peptides"])
# 1601.464
# Standard deviation number of peptides in non-frail
sd(num_pep_sample[num_pep_sample$frailty == "NFT", "peptides"])
# 1049.815

# Wilcox test number of peptides FT vs NFT
wilcox.test(peptides ~ frailty, data = num_pep_sample)
#Wilcoxon rank sum test with continuity correction
#data:  peptides by frailty
#W = 3864.5, p-value = 0.1124
#alternative hypothesis: true location shift is not equal to 0

# Make binary assay
missval <- ifelse(is.na(pep_lfq), 0, 1)

# Convert matrix to long format dataframe
df <- melt(missval)

# Rename columns
colnames(df) <- c("peptides", "sample", "intensity")
df$sample <- as.character(df$sample)

# Add frailty group
df$frailty <- sapply(strsplit(df$sample, "_"), `[`, 1)

# Plot heatmap
p <- ggplot(df, aes(sample, peptides, fill = intensity)) + 
  geom_tile() +
  facet_grid(. ~ frailty, scales = "free_x", space = "free_x") +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 2, angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_gradient(low = "white", high = "black")

save_path <- paste0(plot_path,"preprocessing/heatmap_lfq_na.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Complex Heatmap
# Clustering rows and columns by binary distance
# Column split by frailty groups
p <- Heatmap(missval, name = "NA", col = c("white", "black"),
             show_row_dend = FALSE, 
             column_names_gp = gpar(fontsize = 4),
             show_row_names = FALSE, 
             # cluster_columns = FALSE,
             # cluster_rows = FALSE,
             column_split = metadata_ft$frailty, 
             clustering_distance_rows = "binary", 
             clustering_distance_columns = "binary", 
             use_raster = FALSE)

save_path <- paste0(plot_path,"preprocessing/heatmap_lfq_clust_na.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
draw(p)
dev.off()

# Annotation of columns by frailty group
column_ann <- HeatmapAnnotation(
  frailty = metadata_ft$frailty,
  col = list(frailty = c("FT" = "red", "NFT" = "blue")),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

p <- Heatmap(missval, col = c("white", "black"),
             show_row_dend = FALSE, 
             column_names_gp = gpar(fontsize = 4),
             show_row_names = FALSE, 
             show_column_names = FALSE,
             # cluster_columns = FALSE,
             # cluster_rows = FALSE,
             clustering_distance_rows = "binary", 
             clustering_distance_columns = "binary",
             top_annotation = column_ann, 
             use_raster = FALSE,
             show_heatmap_legend = FALSE
)

# Heatmap legend
legend_1 <- Legend(at = c(0, 1),
                   labels = c("Missing value", "Valid value"),
                   legend_gp = gpar(fill = c("white", "black")),
                   title = NULL)
# Annotation legend
legend_2 <- Legend(at = c("FT", "NFT"),
                   labels = c("Frail", "Non-Frail"),
                   legend_gp = gpar(fill = c("red", "blue")),
                   title = NULL)

# Combine legends
legends <- packLegend(legend_2, legend_1, direction = "horizontal")

# Plot heatmap
save_path <- paste0(plot_path,"preprocessing/heatmap_lfq_clust_na_ann.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
draw(p, heatmap_legend_list = legends, heatmap_legend_side = "top")
dev.off()

# Frailty column as grouping
cond_options <- table(metadata_ft[,"frailty"]) %>% as.data.frame()

# Filtering by minimum count in at least 1 condition
cond_opts <- cond_options$Var1
cond_count <- cond_options$Freq * 0.5

pep_exp <- filter_valids(pep_lfq,
                         conditions = cond_opts,
                         min_count = cond_count,
                         at_least_one = T)

# Number of peptides after filtering by sample proportion
nrow(pep_exp)
# 699

# Number of NA values in assay after filtering by sample proportion
na_table <- table(is.na(pep_exp))
na_table
# FALSE  TRUE 
# 93750 48147 

# NA proportion
prop_na <- na_table["TRUE"] / sum(na_table)
prop_na
# 0.3393095

# Number of peptides detected per sample
num_pep_sample <- apply(pep_exp, 2, function(x) sum(!is.na(x)))
num_pep_sample <- as.data.frame(num_pep_sample)
colnames(num_pep_sample) <- "peptides"
num_pep_sample$sample <- rownames(num_pep_sample)
rownames(num_pep_sample) <- NULL
# Add frailty group
num_pep_sample$frailty <- sapply(strsplit(num_pep_sample$sample, "_"), `[`, 1)

# Summary total peptides detected per sample
summary(num_pep_sample$peptides)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 55.0   380.5   500.0   461.8   569.0   678.0  

# Density plot number of peptides per sample FT vs NFT
p <- ggplot(num_pep_sample, aes(x = peptides, fill = frailty)) +
  geom_density(linewidth = 0.4, alpha = 0.3) +
  labs(title = "Distribution of peptides per sample",
       x = "Number of peptides", 
       y = "Density") +
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

save_path <- paste0(plot_path, "preprocessing/density_plot_lfq_num_pep_",
                    "sample_filt.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

p <- ggplot(num_pep_sample, aes(x = frailty, y = peptides, fill = frailty)) + 
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Frailty group", 
       y = "Number of peptides") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail"))

save_path <- paste0(plot_path, "preprocessing/boxplot_lfq_num_pep_sample_",
                    "filt.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Mean number of peptides in frail
mean(num_pep_sample[num_pep_sample$frailty == "FT", "peptides"])
# 437.3385
# Standard deviation number of peptides in frail
sd(num_pep_sample[num_pep_sample$frailty == "FT", "peptides"])
# 131.0829

# Mean number of peptides in non-frail
mean(num_pep_sample[num_pep_sample$frailty == "NFT", "peptides"])
# 473.3551
# Standard deviation number of peptides in non-frail
sd(num_pep_sample[num_pep_sample$frailty == "NFT", "peptides"])
# 147.0446

# Wilcox test number of peptides FT vs NFT
wilcox.test(peptides ~ frailty, data = num_pep_sample)
# Wilcoxon rank sum test with continuity correction
# data:  peptides by frailty
# W = 3570.5, p-value = 0.01925
# alternative hypothesis: true location shift is not equal to 0

# Make binary assay
missval <- ifelse(is.na(pep_exp), 0, 1)

# Convert matrix to long format dataframe
df <- melt(missval)

# Rename columns
colnames(df) <- c("peptides", "sample", "intensity")
df$sample <- as.character(df$sample)

# Add frailty group
df$frailty <- sapply(strsplit(df$sample, "_"), `[`, 1)

# Plot heatmap
p <- ggplot(df, aes(sample, peptides, fill = intensity)) + 
  geom_tile() +
  facet_grid(. ~ frailty, scales = "free_x", space = "free_x") +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 2, angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_gradient(low = "white", high = "black")

save_path <- paste0(plot_path, "preprocessing/heatmap_lfq_na_filt.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Complex Heatmap
# Clustering rows and columns by binary distance
# Column split by frailty groups
p <- Heatmap(missval, name = "NA", col = c("white", "black"),
             show_row_dend = FALSE, 
             column_names_gp = gpar(fontsize = 4),
             show_row_names = FALSE, 
             # cluster_columns = FALSE,
             # cluster_rows = FALSE,
             column_split = metadata_ft$Fragilidad, 
             clustering_distance_rows = "binary", 
             clustering_distance_columns = "binary", 
             use_raster = FALSE)
save_path <- paste0(plot_path, "preprocessing/heatmap_lfq_clust_na_filt_",
                    "min.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
draw(p)
dev.off()

# Annotation of columns by frailty group
column_ann <- HeatmapAnnotation(
  frailty = metadata_ft$frailty,
  col = list(frailty = c("FT" = "red", "NFT" = "blue")),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

p <- Heatmap(missval, col = c("white", "black"),
             show_row_dend = FALSE, 
             column_names_gp = gpar(fontsize = 4),
             show_row_names = FALSE, 
             show_column_names = FALSE,
             # cluster_columns = FALSE,
             # cluster_rows = FALSE,
             clustering_distance_rows = "binary", 
             clustering_distance_columns = "binary",
             top_annotation = column_ann, 
             use_raster = FALSE,
             show_heatmap_legend = FALSE)

# Heatmap legend
legend_1 <- Legend(at = c(0, 1),
                   labels = c("Missing value", "Valid value"),
                   legend_gp = gpar(fill = c("white", "black")),
                   title = NULL)
# Annotation legend
legend_2 <- Legend(at = c("FT", "NFT"),
                   labels = c("Frail", "Non-Frail"),
                   legend_gp = gpar(fill = c("red", "blue")),
                   title = NULL)

# Combine legends
legends <- packLegend(legend_2, legend_1, direction = "horizontal")

# Plot heatmap
save_path <- paste0(plot_path,"preprocessing/heatmap_lfq_clust_na_ann_",
                    "filt.png")
png(save_path, width = 8, height = 4, units = "in", res = 300)
draw(p, heatmap_legend_list = legends, heatmap_legend_side = "top")
dev.off()


# ==================================
# 2.3. Normalization
# ==================================
# Transform pep_exp to long
pep_exp_long <- pep_exp
pep_exp_long$peptides <- rownames(pep_exp_long)
rownames(pep_exp_long) <- NULL
pep_exp_long <- melt(pep_exp_long)
colnames(pep_exp_long) <- c("peptides", "sample", "intensity")
pep_exp_long$sample <- as.character(pep_exp_long$sample)
# Add frailty group
pep_exp_long$frailty <- sapply(strsplit(pep_exp_long$sample, "_"), `[`, 1)

summary(pep_exp_long$intensity)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# 486    34113    60737   147707   120817 76101048    48147 
sd(pep_exp_long$intensity, na.rm = TRUE)
# 733799.6

# Plot intensities vs samples
p <- ggplot(pep_exp_long, aes(x = sample, y = intensity, fill = frailty)) + 
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Samples", 
       y = "Peptide Top-N intensity") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),  #element_text(angle = 90, size = 6),
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail")) +
  #ylim(0, 250000) +
  facet_wrap("frailty", scales = "free_x", ncol = 1)

save_path <- paste0(plot_path,"preprocessing/boxplot_lfq_samples.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Plot intensities vs frailty group
p <- ggplot(pep_exp_long, aes(x = frailty, y = intensity, fill = frailty)) + 
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Frailty group", 
       y = "Peptide Top-N intensity") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 8),
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail"))

save_path <- paste0(plot_path,"preprocessing/boxplot_lfq_groups.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Intensities density plot
p <- ggplot(pep_exp_long, aes(x = intensity)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  theme_minimal(base_size = 12) + 
  labs(title = "Distribution of peptide Top-N intensities",
       x = "Peptide Top-N intensity", 
       y = "Density")

save_path <- paste0(plot_path,"preprocessing/density_plot_lfq.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)


# Estimate size factor for each sample -> norm_pep (vector)
norm_pep <- estimateSizeFactorsForMatrix(as.matrix(pep_exp))

# Normalise dividing each column by its size factor
norm_pep_exp <- sweep(as.matrix(pep_exp), 2, norm_pep, "/")
norm_pep_exp <- as.data.frame(norm_pep_exp)

# Log2 transform
log_pep_exp <- norm_pep_exp
log_pep_exp[sapply(log_pep_exp, is.numeric)] <- lapply(
  log_pep_exp[sapply(log_pep_exp, is.numeric)],
  function(x) log2(x)
)

# Transform log_pep_exp to long
log_pep_exp_long <- log_pep_exp
log_pep_exp_long$peptides <- rownames(log_pep_exp_long)
rownames(log_pep_exp_long) <- NULL
log_pep_exp_long <- melt(log_pep_exp_long)
colnames(log_pep_exp_long) <- c("peptides", "sample", "intensity")
log_pep_exp_long$sample <- as.character(log_pep_exp_long$sample)
# Add frailty group
log_pep_exp_long$frailty <- sapply(strsplit(log_pep_exp_long$sample, "_"), `[`, 1)

summary(log_pep_exp_long$intensity)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  6.89   14.83   16.07   16.14   17.38   26.52   48147 

sd(log_pep_exp_long$intensity, na.rm = TRUE)
# 1.91668

# Plot log2 normalized intensities vs samples
p <- ggplot(log_pep_exp_long, aes(x = sample, y = intensity, fill = frailty)) + 
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Samples", 
       y = "Peptide log2 SF normalized Top-N intensity") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x =  element_blank(), #element_text(angle = 90, size = 6),
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail")) +
  facet_wrap("frailty", scales = "free_x", ncol = 1)

save_path <- paste0(plot_path,"preprocessing/boxplot_lfq_samples_norm_log2.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Plot log2 normalized intensities vs frailty group
p <- ggplot(log_pep_exp_long, aes(x = frailty, y = intensity, fill = frailty)) + 
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Frailty group", 
       y = "Peptide log2 SF normalized Top-N intensity") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 8),
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail"))

save_path <- paste0(plot_path,"preprocessing/boxplot_lfq_groups_norm_log2.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Log2 normalized intensities density plot
p <- ggplot(log_pep_exp_long, aes(x = intensity)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  labs(title = "Distribution of peptide log2 normalized Top-N intensities",
       x = "Peptide log2 SF normalized Top-N intensity", 
       y = "Density") +
  theme_minimal(base_size = 11) +
  coord_cartesian(xlim = c(8.5, 27), ylim = c(0, 0.35))

save_path <- paste0(plot_path,"preprocessing/density_plot_lfq_norm_log2.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Normalise by vsn
library(vsn)
fit <- vsnMatrix(as.matrix(pep_exp))
vsn_pep_exp <- vsn::predict(fit, as.matrix(pep_exp))

# Transform vsn_pep_exp to long
vsn_pep_exp_long <- as.data.frame(vsn_pep_exp)
vsn_pep_exp_long$peptides <- rownames(vsn_pep_exp_long)
rownames(vsn_pep_exp_long) <- NULL
vsn_pep_exp_long <- melt(vsn_pep_exp_long)
colnames(vsn_pep_exp_long) <- c("peptides", "sample", "intensity")
vsn_pep_exp_long$sample <- as.character(vsn_pep_exp_long$sample)
# Add frailty group
vsn_pep_exp_long$frailty <- sapply(strsplit(vsn_pep_exp_long$sample, "_"), `[`, 1)

summary(vsn_pep_exp_long$intensity)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 9.09   15.13   15.92   16.10   16.89   25.48   48147 

sd(vsn_pep_exp_long$intensity, na.rm = TRUE)
# 1.382352

# Plot VSN intensities vs samples
p <- ggplot(vsn_pep_exp_long, aes(x = sample, y = intensity, fill = frailty)) + 
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Samples", 
       y = "Peptide VSN Top-N intensity") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x =  element_blank(), #element_text(angle = 90, size = 6),
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail")) +
  facet_wrap("frailty", scales = "free_x", ncol = 1)

save_path <- paste0(plot_path,"preprocessing/boxplot_lfq_samples_vsn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Plot log2 normalized intensities vs frailty group
p <- ggplot(vsn_pep_exp_long, aes(x = frailty, y = intensity, fill = frailty)) + 
  geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.5) +
  labs(x = "Frailty group", 
       y = "Peptide VSN Top-N intensity") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 8),
    legend.position = "top",
    legend.title = element_blank()) + 
  scale_fill_manual(values = c("FT" = "red", "NFT" = "purple"),
                    labels = c("FT" = "Frail", "NFT" = "Non-Frail"))

save_path <- paste0(plot_path,"preprocessing/boxplot_lfq_groups_vsn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Log2 normalized intensities density plot
p <- ggplot(vsn_pep_exp_long, aes(x = intensity)) +
  geom_density(fill = "blue", linewidth = 0.4, alpha = 0.4) +
  labs(title = "Distribution of peptide VSN Top-N intensities",
       x = "Peptide VSN Top-N intensity", 
       y = "Density") +
  theme_minimal(base_size = 12) +
  coord_cartesian(xlim = c(9, 26), ylim = c(0, 0.35))

save_path <- paste0(plot_path,"preprocessing/density_plot_lfq_vsn.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)
