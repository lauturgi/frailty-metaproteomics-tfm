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

# Data, output and plot paths
data_path <- paste0(work_path, "/psva/data/")
output_path <- paste0(work_path, "/psva/output/")
plot_path <- paste0(work_path, "/psva/plots/enrichment/")

# Load KEGG peptide sets
load(paste0(data_path, "kegg_pepsets.RData"))

# Load peptide kegg intensities
load(paste0(data_path, "pep_kegg_exp.RData"))

# GSVA
# Parameters for running GSVA
params <- gsvaParam(
  as.matrix(pep_kegg_exp),
  kegg_pepsets,
  minSize = 10,
  kcdf = "Gaussian", 
  use = "na.rm"  # NA values data will be removed from calculations
)

# Calculate GSVA ranks
gsva_ranks <- gsva_ranks(params)

# Save GSVA ranks
save(gsva_ranks, file = paste0(data_path, "gsva_ranks.RData"))
gsva_ranks_expr <- gsva_ranks@exprData
save(gsva_ranks_expr, file = paste0(data_path, "gsva_ranks_expr.RData"))


# Calculate GSVA enrichment scores (ES) from ranks
gsva_scores <- gsva_scores(gsva_ranks)

# Save GSVA scores
save(gsva_scores, file = paste0(data_path, "gsva_scores.RData"))

# Calculate directly GSVA ES
gsva_kegg <- gsva(params)

# NA ES in gene sets with less than 10 genes after removing missing values

# Save GSVA ES
save(gsva_kegg, file = paste0(data_path, "gsva_kegg.RData"))

# Number of peptide sets
nrow(gsva_kegg)
# 71

#================================
# Functional peptide enrichment analysis
#================================
# Load metadata
load(paste0(work_path, "/data/metadata.RData"))

# Align metadata with same order of `replicate` in pep_kegg_exp
replicate_pep <- sapply(strsplit(colnames(pep_kegg_exp), "_"), function(x) x[2])
replicate_metadata <- metadata_ft$replicate
metadata_ft <- metadata_ft[match(replicate_pep, replicate_metadata), ]

# Build model matrix
control <- "NFT"
cond <- factor(metadata_ft$frailty) %>%
  relevel(. , ref=control) # NFT control
design <- model.matrix(~cond)
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
summary(res)
# NFT FT
# Down     0  0
# NotSig   6 56
# Up      65 15

res_2 <- topTable(fit, coef=2, adjust="BH", sort.by="P", number=nrow(gsva_kegg))

# Order by adj.P.Val and logFC
res_2 <- res_2[order(res_2$adj.P.Val, -abs(res_2$logFC)), ]

# Write csv with top table results
write.csv(res_2, file = paste0(output_path, "res_es_gsva.csv"))

res <- res %>% as.data.frame()

# Get significant peptide sets from GSVA scores matrix
sig_tests <- res[abs(res[,2]) > 0,]
sig_gsva <- gsva_kegg[rownames(gsva_kegg) %in% rownames(sig_tests),]

save(sig_gsva, file = paste0(data_path, "sig_gsva.RData"))

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

# Annotation
ann_plot <- ggplot(sig_gsvaplot_data, aes(x = Samples, y = 1, fill = Condition)) +
  geom_tile() +
  scale_fill_manual(values = c("FT" = "red", "NFT" = "blue"),
                    labels = c("FT" = "Frail", "NFT" = "Non-frail")
                    ) +
  theme_void() +
  theme(legend.position = "top")

# Heatmap to plot significantly different KEGG pathways
p <- ggplot(data = sig_gsvaplot_data, aes(x = Samples, y = Pathway, 
                                          fill = GSVA)) +
  #facet_grid(~ Condition, switch = "x", scales = "free_x", space = "free_x") +
  scale_fill_gradientn(colours=c("blue","white","red"),
                       space = "Lab", name="GSVA enrichment score",
                       na.value = "lightgrey") + 
  theme_minimal(base_size = 11) +
  geom_tile(na.rm = TRUE, color="white") +
  xlab(label = "Sample") +
  ylab(label="") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ann_plot / p + plot_layout(heights = c(0.05, 1))

save_path <- paste0(plot_path, "heatmap_plot_gsva_es.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Peptides from significant KEGG pathways
sig_path <- rownames(sig_gsva)
sig_kegg_pepsets <- kegg_pepsets[sig_path]

# Keep only peptides present in pep_kegg_exp
pep_to_keep <- rownames(pep_kegg_exp)
sig_kegg_pepsets <- lapply(sig_kegg_pepsets, function(pepset) {
  intersect(pepset, pep_to_keep)
})

# Load matched_core_pep_ko
load(paste0(data_path, "/matched_core_pep_ko.RData"))

# Limma analysis for peptides of each path

## Alanine, aspartate and glutamate metabolism
pep_kegg_exp_ade <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                   sig_kegg_pepsets[[1]], ]

gsva_ranks_ade <- gsva_ranks_expr[rownames(gsva_ranks_expr) %in%
                                    sig_kegg_pepsets[[1]], ]

test <- cor.test(vars[[i]], vars[[j]], method = "spearman", exact = FALSE)
cor_matrix[i, j] <- test$estimate
p_matrix[i, j] <- test$p.value

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_ade, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_ade))

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
write.csv(res_ko, file = paste0(output_path, "res_dea_ade.csv"))

## Arginine biosynthesis
pep_kegg_exp_r <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                  sig_kegg_pepsets[[2]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_r, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_r))

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
write.csv(res_ko, file = paste0(output_path, "res_dea_r.csv"))

## Chromosome and associated proteins
pep_kegg_exp_chr <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                 sig_kegg_pepsets[[3]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_chr, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_chr))

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
write.csv(res_ko, file = paste0(output_path, "res_dea_chr.csv"))

## Cysteine and methionine metabolism
pep_kegg_exp_cm <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                   sig_kegg_pepsets[[4]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_cm, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_cm))

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
write.csv(res_ko, file = paste0(output_path, "res_dea_cm.csv"))

## Cytoskeleton in muscle cells
pep_kegg_exp_cyto <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                  sig_kegg_pepsets[[5]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_cyto, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_cyto))

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
write.csv(res_ko, file = paste0(output_path, "res_dea_cyto.csv"))

## Enzymes with EC numbers
pep_kegg_exp_ec <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                    sig_kegg_pepsets[[6]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_ec, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_ec))

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
write.csv(res_ko, file = paste0(output_path, "res_dea_ec.csv"))

## Fructose and mannose metabolism
pep_kegg_exp_fm <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                  sig_kegg_pepsets[[7]], ]

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
write.csv(res_ko, file = paste0(output_path, "res_dea_fruc_mano.csv"))

## Glycine, serine and threonine metabolism 
pep_kegg_exp_gst <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                   sig_kegg_pepsets[[8]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_gst, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_gst))

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
write.csv(res_ko, file = paste0(output_path, "res_dea_gst.csv"))

## HIF-1 signaling pathway 
pep_kegg_exp_hif1 <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                   sig_kegg_pepsets[[9]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_hif1, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_hif1))

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
write.csv(res_ko, file = paste0(output_path, "res_dea_hif1.csv"))

## Methane metabolism 
pep_kegg_exp_met <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                    sig_kegg_pepsets[[10]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_met, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_met))

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
write.csv(res_ko, file = paste0(output_path, "res_dea_met.csv"))

## Necroptosis 
pep_kegg_exp_necro <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                     sig_kegg_pepsets[[11]], ]

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
write.csv(res_ko, file = paste0(output_path, "res_dea_necro.csv"))

## Nitrogen metabolism
pep_kegg_exp_nitro <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                     sig_kegg_pepsets[[12]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_nitro, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_nitro))

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
write.csv(res_ko, file = paste0(output_path, "res_dea_nitro.csv"))


## Pentose and glucuronate interconversions
pep_kegg_exp_pg_inter <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                     sig_kegg_pepsets[[13]], ]

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
write.csv(res_ko, file = paste0(output_path, "res_dea_pg_inter.csv"))

## Proximal tubule bicarbonate reclamation
pep_kegg_exp_prox_tub <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                        sig_kegg_pepsets[[14]], ]

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
write.csv(res_ko, file = paste0(output_path, "res_dea_prox_tub.csv"))

## Taurine and hypotaurine metabolism
pep_kegg_exp_tau <- pep_kegg_exp[rownames(pep_kegg_exp) %in%
                                        sig_kegg_pepsets[[15]], ]

# Fit a linear model for each peptide
fit <- lmFit(pep_kegg_exp_tau, design)
# Estimated coefficients for each peptide
head(coef(fit))

# Moderated t-test with empirical Bayes shrinkage of standard errors
fit <- eBayes(fit)
res <- topTable(fit, coef=2, adjust="BH", sort.by="P",
                number=nrow(pep_kegg_exp_tau))

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
write.csv(res_ko, file = paste0(output_path, "res_dea_tau.csv"))
