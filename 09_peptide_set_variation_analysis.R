#renv::install("GSVA")
# Load libraries
library(GSVA)
library(dplyr)
library(limma)
library(tidyr)
library(tibble)
library(ggplot2)


# Working directory
work_path <- getwd()


# Load KEGG peptide sets
load(paste0(work_path, "/data/kegg_pepsets.RData"))

# Load peptide kegg intensities
load(paste0(work_path, "/data/pep_kegg_exp.RData"))


# GSVA
# Parameters for running GSVA
params <- gsvaParam(
  as.matrix(pep_kegg_exp),
  kegg_pepsets,
  minSize = 10,
  kcdf = "Gaussian"
)

# Calculate GSVA enrichment scores
gsva_kegg <- gsva(params)

# Save kegg_pepsets
# save(gsva_kegg, file = paste0(work_path, "/data/gsva_kegg.RData"))
# load(paste0(work_path, "/data/gsva_kegg.RData"))

# Number of peptide sets
nrow(gsva_kegg)
# 49

#================================
# Functional peptide enrichment analysis
#================================
# Load metadata
load(paste0(work_path, "/data/metadata.RData"))

# Build model matrix
control <- "NFT"
cond <- factor(metadata_ft$Fragilidad) %>%
  relevel(. , ref=control) # NFT control
design <- model.matrix(~  cond)
colnames(design) <- c("NFT", "FT")
head(design)

# Fit a linear model for each peptide set of GSVA scores matrix
fit <- lmFit(gsva_kegg, design)
# Estimated coefficients for each peptide set
head(coef(fit))

# Moderate t-test of differential KEGG pathways by means of Bayes' empirical
# moderation of the standard errors to a global value
fit <- eBayes(fit)
res <- decideTests(fit, p.value=0.05, adjust="BH")
#res <- topTable(fit, coef=2, adjust="BH", sort.by="P", number=nrow(gsva_kegg))
summary(res)

#        NFT FT
# Down     0  4
# NotSig  48 44
# Up       1  1

res <- res %>% as.data.frame()

# Get significant peptide sets from GSVA scores matrix
sig_tests <- res[abs(res[,2]) > 0,]
sig_gsva <- gsva_kegg[rownames(gsva_kegg) %in% rownames(sig_tests),]

# Rownames to column "Pathway", reshape to long format and add "Condition"
sig_gsvaplot_data <- data.frame(sig_gsva) %>%
  rownames_to_column(., var="Pathway") %>%
  pivot_longer(!Pathway, names_to="Samples", values_to = "GSVA")
sig_gsvaplot_data$Condition <- sapply(strsplit(sig_gsvaplot_data$Samples, "_"),
                                      `[`, 1)

# Compute distance (euclidean) between pathways (rows) of GSVA matrix and cluster
cluster_path <- rownames(gsva_kegg)[hclust(dist(gsva_kegg))$order]
sig_gsvaplot_data$Pathway<- factor(sig_gsvaplot_data$Pathway,
                                   levels = cluster_path)

# Compute distance between samples (columns) of GSVA matrix and cluster
cluster_sample <- colnames(gsva_kegg)[hclust(dist(gsva_kegg %>% t()))$order]
sig_gsvaplot_data$Samples<- factor(sig_gsvaplot_data$Samples,
                                   levels = cluster_sample)

# Rownames to column "Pathway" and reshape to long format
sig_tests <- as.data.frame(sig_tests) %>%
  rownames_to_column(., var = "Pathway") %>%  
  pivot_longer(!Pathway, names_to= "Condition", values_to = "Significance")


# Heatmap to plot significantly different KEGG pathways
ggplot(data = sig_gsvaplot_data, mapping = aes(x = Samples, y = Pathway,
                                               fill = GSVA)) +
  facet_grid(~ Condition, switch = "x", scales = "free_x", space = "free_x") +
  scale_fill_gradientn(colours=c("#67A7C1","white","#FF6F59"),
                       space = "Lab", name="GSVA enrichment score",
                       na.value = "lightgrey") + 
  geom_tile(na.rm = TRUE, color="white") +
  xlab(label = "Sample") +
  ylab(label="") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


# Peptides and KOs from significant KEGG pathways
sig_path <- rownames(sig_gsva)
sig_kegg_pepsets <- kegg_pepsets[sig_path]
