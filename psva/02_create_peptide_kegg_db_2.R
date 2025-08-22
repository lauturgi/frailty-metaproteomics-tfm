# Load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)


# Working directory
work_path <- getwd()

# Plot path
plot_path <- paste0(work_path, "/psva/plots/")

# Read UHGP annotation file
uhgp_ann <- fread(file.path(paste0(work_path, "/psva/data/uhgp-90_eggNOG.tsv")),
                  sep = "\t", header = TRUE) %>%
  select(1, 12) %>%
  filter(KEGG_ko != "-")  # Filter UHGP proteins without KEGG KO
colnames(uhgp_ann) <- c("Protein", "KO")
nrow(uhgp_ann)
# 4851159

# Separate rows with more than one KO in several rows
uhgp_ann <- uhgp_ann %>%
  separate_rows(KO, sep = ",")
nrow(uhgp_ann)
# 5379202

# Count number of KO terms per protein
num_ko_prot <- uhgp_ann %>%
  group_by(Protein) %>%
  summarise(KO_count = n())

# Count number of proteins with 1, 2, 3 and 4 or more KO terms
df <- num_ko_prot %>%
  mutate(num_ko_group = ifelse(KO_count >= 4, ">=4", 
                               as.character(KO_count))) %>%
  count(num_ko_group) %>%
  mutate(percent = (n / sum(n)) * 100) 
df$num_ko_group <- factor(df$num_ko_group, levels = df$num_ko_group)

# Get label positions for ggrepel
df2 <- df %>% 
  mutate(csum = rev(cumsum(rev(percent))), 
         pos = percent/2 + lead(csum, 1),
         pos = if_else(is.na(pos), percent/2, pos))

# Pie plot
p <- ggplot(df, aes(x = "" , y = percent, fill = num_ko_group)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_label_repel(data = df2, aes(y = pos,
                                   label = paste0(round(percent, 1), "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Number KO")) +
  theme_void() +
  scale_fill_manual(values = c("1" = "#7570b3",
                               "2" = "#b3b3e6", 
                               "3" = "#f4a582",
                               ">=4" = "#92c5de"))

save_path <- paste0(plot_path, "/pie_plot_num_ko_protein.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)


# Read diamond_uhgp_out
Prot2UR <- fread(file.path(work_path,
                           "/psva/data/diamond_uhgp_best_hits_out.tsv"),
                 header = FALSE) %>%
  select(1, 3) %>%
  rename(Protein = V1, ID = V3)
nrow(Prot2UR)
# 12317685


# Dataframe with total number of UHGP proteins and subsets annotated with UniRef
# and KO
n_prot_df <- data.frame(
  variable = c("Total UHGP proteins", "Annotated with UniRef", "Annotated with KO"),
  count = c(13811247, 12317685, 4851159)
)
n_prot_df$variable <- factor(n_prot_df$variable, levels = n_prot_df$variable)

# Calculate percentage of proteins with UniRef and KO annotation
n_prot_df$percent <- n_prot_df$count/
  n_prot_df$count[n_prot_df$variable == "Total UHGP proteins"]*100

# Bar plot
p <- ggplot(n_prot_df, aes(x = variable, y = count, fill=variable)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(
    scales::comma(count), " (", round(percent, 1), "%)")), vjust = -0.5) +
  scale_y_continuous(
    labels = scales::comma
  ) +
  labs(y = "Number of proteins") +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = c("Total UHGP proteins" = "lightgrey", 
                               "Annotated with UniRef" = "lightgrey", 
                               "Annotated with KO" = "#7570b3")) + 
  theme(
    axis.title.x = element_blank(), 
    legend.position = "none"
  )

save_path <- paste0(plot_path, "/bar_plot_num_proteins_uhgp_uniref_ko.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Read uhgp_pep
Pep2Prot <- fread(file.path(work_path, "/psva/data/uhgp_pep.tsv"),
                  header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)
nrow(Pep2Prot)
# 1027319075

nrow(unique(Pep2Prot[, "Protein"]))
# 13797153

# Join by "Protein" to get UniRef90
Pep2Prot2UR <- left_join(Pep2Prot, Prot2UR, by = "Protein",
                         relationship = "many-to-many")
nrow(Pep2Prot2UR)
# 1027319075

nrow(unique(Pep2Prot2UR[, "Protein"]))
# 13797153

nrow(unique(Pep2Prot2UR[!is.na(Pep2Prot2UR$ID), "Protein"]))
# 12310128

# Join by "Protein" to get KO
Pep2Prot2UR2KEGG <- left_join(Pep2Prot2UR, uhgp_ann, by = "Protein",
                              relationship = "many-to-many")
nrow(Pep2Prot2UR2KEGG)
# 1079321301

# Number missing KO
sum(is.na(Pep2Prot2UR2KEGG$KO))
# 601822433

# Number proteins with KO
nrow(unique(Pep2Prot2UR2KEGG[!is.na(Pep2Prot2UR2KEGG$KO), "Protein"]))
# 4850365

# Number proteins without KO
nrow(unique(Pep2Prot2UR2KEGG[is.na(Pep2Prot2UR2KEGG$KO), "Protein"]))
# 8946788

# Number peptides with KO
nrow(unique(Pep2Prot2UR2KEGG[!is.na(Pep2Prot2UR2KEGG$KO), "Peptide"]))
# 283000810

# Number peptides without KO
nrow(unique(Pep2Prot2UR2KEGG[is.na(Pep2Prot2UR2KEGG$KO), "Peptide"]))
# 435375925

# Number missing UniRef90
sum(is.na(Pep2Prot2UR2KEGG$ID))
# 54138006

# Number proteins with UniRef90
nrow(unique(Pep2Prot2UR2KEGG[!is.na(Pep2Prot2UR2KEGG$ID), "Protein"]))
# 12310128

# Number proteins without UniRef90
nrow(unique(Pep2Prot2UR2KEGG[is.na(Pep2Prot2UR2KEGG$ID), "Protein"]))
# 1487025

# Number peptides with UniRef90
nrow(unique(Pep2Prot2UR2KEGG[!is.na(Pep2Prot2UR2KEGG$ID), "Peptide"]))
# 662797145

# Number peptides without UniRef90
nrow(unique(Pep2Prot2UR2KEGG[is.na(Pep2Prot2UR2KEGG$ID), "Peptide"]))
# 48064982

# Number proteins with KO but not UniRef90
nrow(unique(Pep2Prot2UR2KEGG[!is.na(Pep2Prot2UR2KEGG$KO) &
                               is.na(Pep2Prot2UR2KEGG$ID), "Protein"]))
# 344

# Number peptides with KO but not UniRef90
nrow(unique(Pep2Prot2UR2KEGG[!is.na(Pep2Prot2UR2KEGG$KO) &
                               is.na(Pep2Prot2UR2KEGG$ID), "Peptide"]))
# 11323

# Number of proteins with KO and UniRef90
nrow(unique(Pep2Prot2UR2KEGG[!is.na(Pep2Prot2UR2KEGG$KO) &
                               !is.na(Pep2Prot2UR2KEGG$ID), "Protein"]))
# 4850021

# Number of peptides with KO and UniRef90
nrow(unique(Pep2Prot2UR2KEGG[!is.na(Pep2Prot2UR2KEGG$KO) &
                               !is.na(Pep2Prot2UR2KEGG$ID), "Peptide"]))
# 282991120

# Filter peptides without KO annotation
Pep2Prot2UR2KEGG <- Pep2Prot2UR2KEGG %>% filter(!is.na(KO))
nrow(Pep2Prot2UR2KEGG)
# 477498868
nrow(unique(Pep2Prot2UR2KEGG[, "Protein"]))
# 4850365

# Dataframe with total number of UHGP peptides and subsets annotated with UniRef
# and KO
n_pep_df <- data.frame(
  variable = c("Total UHGP peptides", "Annotated with UniRef", "Annotated with KO"),
  count = c(1027319075, 662797145, 283000810)
)
n_pep_df$variable <- factor(n_pep_df$variable, levels = n_pep_df$variable)

# Calculate percentage of peptides with UniRef and KO annotation
n_pep_df$percent <- n_pep_df$count/
  n_pep_df$count[n_pep_df$variable == "Total UHGP peptides"]*100

# Bar plot
p <- ggplot(n_pep_df, aes(x = variable, y = count, fill=variable)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(
    scales::comma(count), " (", round(percent, 1), "%)")), vjust = -0.5) +
  scale_y_continuous(
    labels = scales::comma
  ) +
  labs(y = "Number of peptides") +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = c("Total UHGP peptides" = "lightgrey", 
                               "Annotated with UniRef" = "lightgrey", 
                               "Annotated with KO" = "#a4133c")) + 
  theme(
    axis.title.x = element_blank(), 
    legend.position = "none"
  )

save_path <- paste0(plot_path, "/bar_plot_num_peptides_uhgp_uniref_ko.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# Count occurrences of each unique peptide–KO term combination
count_ko_pep <- Pep2Prot2UR2KEGG %>%
  group_by(Peptide, KO) %>%
  summarise(KO_count = n())

# Total KO inherited by suming KO_count per peptide
total_ko_pep <- count_ko_pep %>%
  group_by(Peptide) %>%
  summarise(Total_KO_count = sum(KO_count))

# Percentage of peptides with 1, 2, 3, 4, 5, 6 or more KO terms inherited
total_ko_df <- total_ko_pep %>%
  mutate(total_ko_group = ifelse(Total_KO_count >= 4, ">=4", 
                               as.character(Total_KO_count))) %>%
  count(total_ko_group) %>%
  mutate(percent = (n / sum(n)) * 100) 
total_ko_df$total_ko_group <- factor(total_ko_df$total_ko_group, 
                                     levels = total_ko_df$total_ko_group)

# Get label positions for ggrepel
df <- total_ko_df %>% 
  mutate(csum = rev(cumsum(rev(percent))), 
         pos = percent/2 + lead(csum, 1),
         pos = if_else(is.na(pos), percent/2, pos))

# Pie plot
p <- ggplot(total_ko_df, aes(x = "" , y = percent, fill = total_ko_group)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_label_repel(data = df, aes(y = pos,
                                   label = paste0(round(percent, 1), "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Total KO")) +
  theme_void() +
  scale_fill_manual(values = c("1" = "#a4133c",
                               "2" = "#e63946", 
                               "3" = "#f1a208",
                               ">=4" = "#ffb5a7"))

save_path <- paste0(plot_path, "/pie_plot_total_ko_peptide.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# max 9264
# Number of distint KO assigned per peptide
dist_ko_pep <- count_ko_pep %>%
  group_by(Peptide) %>%
  summarise(KO_dist = n_distinct(KO))

# Percentage of peptides with 1, 2, 3, 4 or more distinct KO terms
dist_ko_df <- dist_ko_pep %>%
  mutate(dist_ko_group = ifelse(KO_dist >= 4, ">=4", 
                                 as.character(KO_dist))) %>%
  count(dist_ko_group) %>%
  mutate(percent = (n / sum(n)) * 100) 
dist_ko_df$dist_ko_group <- factor(dist_ko_df$dist_ko_group, 
                                     levels = dist_ko_df$dist_ko_group)

# Get label positions for ggrepel
df <- dist_ko_df %>% 
  mutate(csum = rev(cumsum(rev(percent))), 
         pos = percent/2 + lead(csum, 1),
         pos = if_else(is.na(pos), percent/2, pos))

# Pie plot
p <- ggplot(dist_ko_df, aes(x = "" , y = percent, fill = dist_ko_group)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_label_repel(data = df, aes(y = pos,
                                  label = paste0(round(percent, 1), "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Distinct KO")) +
  theme_void() +
  scale_fill_manual(values = c("1" = "#a4133c",
                               "2" = "#e63946", 
                               "3" = "#f1a208",
                               ">=4" = "#ffb5a7"))

save_path <- paste0(plot_path, "/pie_plot_dist_ko_peptide.png")
ggsave(filename = save_path, plot = p, width = 8, height = 4, dpi = 300)

# max 1548

# Save peptide-protein-UniRef90-KEGG dataset as csv file
write.csv(Pep2Prot2UR2KEGG, file = paste0(work_path,
                                          "/psva/data/core_pep_kegg_db_2.csv"))