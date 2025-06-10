#renv::install("GSVA")
# Load libraries
library(GSVA)
library(dplyr)
library(limma)
library(tidyr)
library(tibble)
library(ggplot2)
library(KEGGREST)


# Working directory
work_path <- getwd()


# Load KEGG peptide sets
load(paste0(work_path, "/psva/data/kegg_pepsets.RData"))

# Load peptide kegg intensities
load(paste0(work_path, "/psva/data/maxlfq/pep_kegg_exp.RData"))

# Convert 0 to NA
pep_kegg_exp[pep_kegg_exp==0] <- NA

# GSVA
# Parameters for running GSVA
params <- gsvaParam(
  as.matrix(pep_kegg_exp),
  kegg_pepsets,
  minSize = 10,
  kcdf = "Gaussian", 
  use = "na.rm"  # NA values data will be removed from calculations
)

# Calculate GSVA enrichment scores
gsva_kegg <- gsva(params)

# Save kegg_pepsets
# save(gsva_kegg, file = paste0(work_path, "/psva/data/maxlfq/gsva_kegg.RData"))
# load(paste0(work_path, "/psva/data/maxlfq/gsva_kegg.RData"))

# Number of peptide sets
nrow(gsva_kegg)
# 71

#================================
# Functional peptide enrichment analysis
#================================
# Load metadata
load(paste0(work_path, "/data/metadata.RData"))

# Build model matrix
control <- "NFT"
cond <- factor(metadata_ft$frailty) %>%
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
res_2 <- topTable(fit, coef=2, adjust="BH", sort.by="P", number=nrow(gsva_kegg))
write.csv(res_2, file = paste0(work_path, "/psva/data/maxlfq/",
                                "res_es_gsva.csv"))
summary(res)

#        NFT FT
# Down     15  3
# NotSig  47 63
# Up       9  5

res <- res %>% as.data.frame()

# Get significant peptide sets from GSVA scores matrix
sig_tests <- res[abs(res[,2]) > 0,]
sig_gsva <- gsva_kegg[rownames(gsva_kegg) %in% rownames(sig_tests),]

save(sig_gsva, file = paste0(work_path, "/psva/data/maxlfq/sig_gsva.RData"))

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


# Peptides from significant KEGG pathways
sig_path <- rownames(sig_gsva)
sig_kegg_pepsets <- kegg_pepsets[sig_path]

# Keep only peptides present in pep_kegg_exp
pep_to_keep <- rownames(pep_kegg_exp)
sig_kegg_pepsets <- lapply(sig_kegg_pepsets, function(pepset) {
  intersect(pepset, pep_to_keep)
})
#load(paste0(work_path, "/psva/data/maxlfq/sig_kegg_pepsets.RData"))

# Load matched_core_pep_ko
load(paste0(work_path, "/psva/data/maxlfq/matched_core_pep_ko.RData"))

# Limma analysis for peptides of each path

## ABC transporters
pep_kegg_exp_abc <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                   sig_kegg_pepsets[[1]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_abc, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_abc))

res$pep <- rownames(res)
rownames(res) <- NULL

# Join res and matched_core_pep_ko to get ko
res_ko <- inner_join(res, matched_core_pep_ko, by = c("pep" = "newpep_name"))
res_ko <- unique(res_ko[, c("logFC", "adj.P.Val", "pep", "ko")])

# Get ko description
ko_ids <- unique(res_ko$ko)
ko_descr <- sapply(ko_ids, function(ko) {
  keggList(ko)
  })
ko_df <- data.frame(
  ko = sapply(names(ko_descr), function(x) strsplit(x, "\\.")[[1]][1]),
  descr = as.character(ko_descr)
)
res_ko <- merge(res_ko, ko_df, by = "ko")

# Order by adj.P.Val and logFC
res_ko <- res_ko[order(res_ko$adj.P.Val, -abs(res_ko$logFC)), ]

# Save res_ko csv file
write.csv(res_ko, file = paste0(work_path, "/psva/data/maxlfq/",
                                    "res_dea_abc.csv"))

## Fructose and mannose metabolism
pep_kegg_exp_fm <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                  sig_kegg_pepsets[[2]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_fm, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_fm))

res$pep <- rownames(res)
rownames(res) <- NULL

# Join res and matched_core_pep_ko to get ko
res_ko <- inner_join(res, matched_core_pep_ko, by = c("pep" = "newpep_name"))
res_ko <- unique(res_ko[, c("logFC", "adj.P.Val", "pep", "ko")])

# Get ko description
ko_ids <- unique(res_ko$ko)
ko_descr <- sapply(ko_ids, function(ko) {
  keggList(ko)
})
ko_df <- data.frame(
  ko = sapply(names(ko_descr), function(x) strsplit(x, "\\.")[[1]][1]),
  descr = as.character(ko_descr)
)
res_ko <- merge(res_ko, ko_df, by = "ko")

# Order by adj.P.Val and logFC
res_ko <- res_ko[order(res_ko$adj.P.Val, -abs(res_ko$logFC)), ]

# Save res_ko csv file
write.csv(res_ko, file = paste0(work_path, "/psva/data/maxlfq/",
                                    "res_dea_fruc_mano.csv"))

## Galactose metabolism
pep_kegg_exp_gal <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                  sig_kegg_pepsets[[3]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_gal, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_gal))

res$pep <- rownames(res)
rownames(res) <- NULL

# Join res and matched_core_pep_ko to get ko
res_ko <- inner_join(res, matched_core_pep_ko, by = c("pep" = "newpep_name"))
res_ko <- unique(res_ko[, c("logFC", "adj.P.Val", "pep", "ko")])

# Get ko description
ko_ids <- unique(res_ko$ko)
ko_descr <- sapply(ko_ids, function(ko) {
  keggList(ko)
})
ko_df <- data.frame(
  ko = sapply(names(ko_descr), function(x) strsplit(x, "\\.")[[1]][1]),
  descr = as.character(ko_descr)
)
res_ko <- merge(res_ko, ko_df, by = "ko")

# Order by adj.P.Val and logFC
res_ko <- res_ko[order(res_ko$adj.P.Val, -abs(res_ko$logFC)), ]

# Save res_ko csv file
write.csv(res_ko, file = paste0(work_path, "/psva/data/maxlfq/",
                                    "res_dea_gal.csv"))

## Necroptosis
pep_kegg_exp_necro <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                   sig_kegg_pepsets[[4]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_necro, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_necro))

res$pep <- rownames(res)
rownames(res) <- NULL

# Join res and matched_core_pep_ko to get ko
res_ko <- inner_join(res, matched_core_pep_ko, by = c("pep" = "newpep_name"))
res_ko <- unique(res_ko[, c("logFC", "adj.P.Val", "pep", "ko")])

# Get ko description
ko_ids <- unique(res_ko$ko)
ko_descr <- sapply(ko_ids, function(ko) {
  keggList(ko)
})
ko_df <- data.frame(
  ko = sapply(names(ko_descr), function(x) strsplit(x, "\\.")[[1]][1]),
  descr = as.character(ko_descr)
)
res_ko <- merge(res_ko, ko_df, by = "ko")

# Order by adj.P.Val and logFC
res_ko <- res_ko[order(res_ko$adj.P.Val, -abs(res_ko$logFC)), ]

# Save res_ko csv file
write.csv(res_ko, file = paste0(work_path, "/psva/data/maxlfq/",
                                    "res_dea_necro.csv"))

## Pentose and glucuronate interconversions
pep_kegg_exp_pg_inter <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                     sig_kegg_pepsets[[5]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_pg_inter, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_pg_inter))

res$pep <- rownames(res)
rownames(res) <- NULL

# Join res and matched_core_pep_ko to get ko
res_ko <- inner_join(res, matched_core_pep_ko, by = c("pep" = "newpep_name"))
res_ko <- unique(res_ko[, c("logFC", "adj.P.Val", "pep", "ko")])

# Get ko description
ko_ids <- unique(res_ko$ko)
ko_descr <- sapply(ko_ids, function(ko) {
  keggList(ko)
})
ko_df <- data.frame(
  ko = sapply(names(ko_descr), function(x) strsplit(x, "\\.")[[1]][1]),
  descr = as.character(ko_descr)
)
res_ko <- merge(res_ko, ko_df, by = "ko")

# Order by adj.P.Val and logFC
res_ko <- res_ko[order(res_ko$adj.P.Val, -abs(res_ko$logFC)), ]

# Save res_ko csv file
write.csv(res_ko, file = paste0(work_path, "/psva/data/maxlfq/",
                                    "res_dea_pg_inter.csv"))

## Proximal tubule bicarbonate reclamation
pep_kegg_exp_prox_tub <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                        sig_kegg_pepsets[[6]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_prox_tub, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_prox_tub))

res$pep <- rownames(res)
rownames(res) <- NULL

# Join res and matched_core_pep_ko to get ko
res_ko <- inner_join(res, matched_core_pep_ko, by = c("pep" = "newpep_name"))
res_ko <- unique(res_ko[, c("logFC", "adj.P.Val", "pep", "ko")])

# Get ko description
ko_ids <- unique(res_ko$ko)
ko_descr <- sapply(ko_ids, function(ko) {
  keggList(ko)
})
ko_df <- data.frame(
  ko = sapply(names(ko_descr), function(x) strsplit(x, "\\.")[[1]][1]),
  descr = as.character(ko_descr)
)
res_ko <- merge(res_ko, ko_df, by = "ko")

# Order by adj.P.Val and logFC
res_ko <- res_ko[order(res_ko$adj.P.Val, -abs(res_ko$logFC)), ]


# Save res_ko csv file
write.csv(res_ko, file = paste0(work_path, "/psva/data/maxlfq/",
                                    "res_dea_prox_tub.csv"))

## Secretion system
pep_kegg_exp_secre <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                        sig_kegg_pepsets[[7]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_secre, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_secre))

res$pep <- rownames(res)
rownames(res) <- NULL

# Join res and matched_core_pep_ko to get ko
res_ko <- inner_join(res, matched_core_pep_ko, by = c("pep" = "newpep_name"))
res_ko <- unique(res_ko[, c("logFC", "adj.P.Val", "pep", "ko")])

# Get ko description
ko_ids <- unique(res_ko$ko)
ko_descr <- sapply(ko_ids, function(ko) {
  keggList(ko)
})
ko_df <- data.frame(
  ko = sapply(names(ko_descr), function(x) strsplit(x, "\\.")[[1]][1]),
  descr = as.character(ko_descr)
)
res_ko <- merge(res_ko, ko_df, by = "ko")

# Order by adj.P.Val and logFC
res_ko <- res_ko[order(res_ko$adj.P.Val, -abs(res_ko$logFC)), ]

# Save res_ko csv file
write.csv(res_ko, file = paste0(work_path, "/psva/data/maxlfq/",
                                    "res_dea_secre.csv"))

## Transporters
pep_kegg_exp_trans <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                     sig_kegg_pepsets[[8]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_trans, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_trans))

res$pep <- rownames(res)
rownames(res) <- NULL

# Join res and matched_core_pep_ko to get ko
res_ko <- inner_join(res, matched_core_pep_ko, by = c("pep" = "newpep_name"))
res_ko <- unique(res_ko[, c("logFC", "adj.P.Val", "pep", "ko")])

# Get ko description
ko_ids <- unique(res_ko$ko)
ko_descr <- sapply(ko_ids, function(ko) {
  keggList(ko)
})
ko_df <- data.frame(
  ko = sapply(names(ko_descr), function(x) strsplit(x, "\\.")[[1]][1]),
  descr = as.character(ko_descr)
)
res_ko <- merge(res_ko, ko_df, by = "ko")

# Order by adj.P.Val and logFC
res_ko <- res_ko[order(res_ko$adj.P.Val, -abs(res_ko$logFC)), ]

# Save res_ko csv file
write.csv(res_ko, file = paste0(work_path, "/psva/data/maxlfq/",
                                    "res_dea_trans.csv"))
