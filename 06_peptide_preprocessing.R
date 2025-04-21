#renv::install("DESeq2")
# Load libraries
library(DESeq2)
library(dplyr)
#library(tibble)

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
  select(-contains("MaxLFQ")) %>%  # Exclude columns that contain MaxLFQ
  select(contains(metadata_ft$replicate))

# Change colnames to include the group each replicate belongs to
replicates_ft <- metadata_ft[metadata_ft$Fragilidad == "FT", "replicate"]
replicates_nft <- metadata_ft[metadata_ft$Fragilidad == "NFT", "replicate"]
for (i in seq_along(colnames(peptides))) {
  colname <- colnames(peptides)[i]  # Get column name
  
  if (substr(colname, 5, 7) %in% replicates_ft) {
    colnames(peptides)[i] <- paste0("FT", substr(colname, 4, nchar(colname)))
  } else {
    colnames(peptides)[i] <- paste0("NFT", substr(colname, 4, nchar(colname)))
  }
}

# Fragilidad column as grouping
cond_options <- table(metadata_ft[,"Fragilidad"]) %>% as.data.frame()

# Filtering by minimum count in at least 1 condition
cond_opts <- cond_options$Var1
cond_count <- cond_options$Freq * 0.5
pep_exp <- filter_valids(peptides,
                         conditions = cond_opts,
                         min_count = cond_count,
                         at_least_one = T)

# ==================================
# Normalization
# ==================================
# Estimate size factor for each sample -> norm_pep (vector)
norm_pep <- estimateSizeFactorsForMatrix(as.matrix(pep_exp),
                                         type = "poscounts")
# Normalise dividing each column by its size factor
norm_pep_exp <- sweep(as.matrix(pep_exp), 2, norm_pep, "/")
norm_pep_exp <- as.data.frame(norm_pep_exp)

save(norm_pep_exp, file = paste0(work_path, "/data/norm_pep_exp.RData"))

# Log2 transform ?????
log_pep_exp <- norm_pep_exp
log_pep_exp[, grep("Intensity", names(norm_pep_exp))] <- lapply(
  norm_pep_exp[, grep("Intensity", names(norm_pep_exp))], 
  function(x) log2(1 + x)
)

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

