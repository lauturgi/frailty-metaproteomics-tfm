#renv::install("DESeq2")
#renv::install("GSVA")
# Load libraries
library(DESeq2)
library(GSVA)
library(limma)
library(dplyr)
library(tibble)


# Working directory
work_path <- getwd()

# Load metadata
load(paste0(work_path, "/metadata.RData"))

# Load KEGG peptide sets
load(paste0(work_path, "/kegg_pepsets.RData"))

# Load peptide kegg intensities
load(paste0(work_path, "/pep_kegg_exp.RData"))


# GSVA
# Parameters for running GSVA
params <- gsvaParam(
  as.matrix(pep_kegg_exp),
  kegg_pepsets,
  minSize = 10,
  kcdf = "Gaussian"
)

# Execute GSVA with gsvaParam
gsva_kegg <- gsva(params)

num_peptidesets <- nrow(gsva_kegg)

#================================
# Functional peptide enrichment analysis
#================================

# KEGG pathways
control_condition <- "NFT"

cond <- factor(metadata$frailty) %>% relevel(. , ref=control_condition) # NFT is the control
design <- model.matrix(~  cond)
colnames(design)[2:ncol(design)] <- substr(colnames(design)[2:ncol(design)], 5, 
                                           nchar(colnames(design)[2:ncol(design)])) #just removing "cond"
colnames(design)[1] <- control_condition
fit <- lmFit(gsva_kegg, design)
fit <- eBayes(fit, trend=T)
res <- decideTests(fit, p.value=0.05, adjust="BH")
summary(res)
#NFT FT
#Down    10  5
#NotSig  17 36
#Up      15  1

res <- res %>% as.data.frame()

# A simple heatmap describing significantly different KEGG pathways when compared to a reference condition

sig_tests <- res[abs(res[,2]) > 0,]

sig_gsva <- gsva_kegg[rownames(gsva_kegg) %in% rownames(sig_tests),]

sig_gsvaplot_data <- data.frame(sig_gsva) %>% rownames_to_column(., var="Pathway") %>%
  pivot_longer(!Pathway, names_to="Samples", values_to = "GSVA")
sig_gsvaplot_data$Condition <- sapply(strsplit(sig_gsvaplot_data$Samples, "_"), `[`, 1)

# clustering pathways
clusterdata <- rownames(gsva_kegg)[hclust(dist(gsva_kegg))$order]
sig_gsvaplot_data$Pathway<- factor(sig_gsvaplot_data$Pathway, levels = clusterdata)

# clustering samples
clustersamples <- colnames(gsva_kegg)[hclust(dist(gsva_kegg %>% t()))$order]
sig_gsvaplot_data$Samples<- factor(sig_gsvaplot_data$Samples, levels = clustersamples)

sig_tests <- as.data.frame(sig_tests) %>% rownames_to_column(., var = "Pathway") %>%  
  pivot_longer(!Pathway, names_to= "Condition", values_to = "Significance")

sig_gsvaplot_data <- sig_gsvaplot_data %>% merge(., sig_tests, by=c("Pathway", "Condition"))  %>%
  dplyr::mutate(GSVA = ifelse(Significance == 0 & Condition != control_condition, NA, GSVA))

(sigplot <- ggplot(data = sig_gsvaplot_data, mapping = aes(x = Samples, y = Pathway, fill = GSVA)) + 
    facet_grid(~ Condition, switch = "x", scales = "free_x", space = "free_x") +
    scale_fill_gradientn(colours=c("#67A7C1","white","#FF6F59"),
                         space = "Lab", name="GSVA enrichment score", na.value = "lightgrey") + 
    geom_tile(na.rm = TRUE, color="white") +
    xlab(label = "Sample") +
    ylab(label="") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))

# Although less informative, all identified peptide gene sets and their
# accompanying GSVA scores are plotted below.
# A simple heatmap describing all KEGG pathways identified in the dataset
gsvaplot_data <- data.frame(gsva_kegg) %>% rownames_to_column(., var="Pathway") %>%
  pivot_longer(!Pathway, names_to="Samples", values_to = "GSVA")

gsvaplot_data$Condition <- sapply(strsplit(gsvaplot_data$Samples, "_"), `[`, 1)


gsvaplot_data$Pathway<- factor(gsvaplot_data$Pathway, levels = clusterdata)
gsvaplot_data$Samples<- factor(gsvaplot_data$Samples, levels = clustersamples)


(gsvaplot <- ggplot(data = gsvaplot_data, mapping = aes(x = Samples, y = Pathway, fill = GSVA)) + 
    facet_grid(~ Condition, switch = "x", scales = "free_x", space = "free_x") +
    scale_fill_gradientn(colours=c("#67A7C1","white","#FF6F59"),
                         space = "Lab", name="GSVA enrichment score", na.value = "#64686F") + 
    geom_tile(na.rm = TRUE, color="white") +
    xlab(label = "Sample") +
    ylab(label="") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))


