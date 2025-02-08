library(plotly)
library(dplyr)

# Working directory
work_path <- "/home/lauturgi/workspace/tfm"

# Unipept LCA analysis output directory
unipept_path <- "/mnt/d/proteomica/fragilidad/datos/PeptideFilesUnipept/lca_output/"

# Read Unipept LCA output files from frailty cohort
nrow_samples <- c()  # Vector to store number of peptides of each sample
for (file in list.files(unipept_path)) {
  replicate <- substring(file, 4, 6)
  if (replicate %in% metadata_ft$replicate) {
    df_name <- paste0("unipept_lca_", replicate)
    df <- read.csv(paste0(unipept_path,file),sep=",")
    nrow_sample <- nrow(df)
    print(paste0("Number of unique peptides between 5 and 50 aa of sample ",
                 replicate, ": ", nrow_sample))
    nrow_samples <- c(nrow_samples, nrow_sample)
    assign(df_name, df)
    # Save LCA output as RData
    save_path <- paste0(work_path,paste0("/data/",df_name,".RData"))
    save(df, file = save_path)
  }
}

# Summary number of peptides
summary(nrow_samples)

# Group by superkingdom, phylum, class, order, family, genus and species and
# count number of peptides that belong to each LCA
taxonomy_118 <- unipept_lca_118 %>%
  group_by(superkingdom_name, phylum_name, class_name, order_name, family_name,
           genus_name, species_name) %>%
  summarize(count = n())

# Create taxonomy for sunburst plot
# Number of total peptides that constitute the root
root_118 <- unipept_lca_118 %>%
  mutate(parent = "", taxon = "root") %>%
  group_by(taxon, parent) %>%
  summarize(count = n())
# Number of peptides that belong to each superkingdom and set root as parent
superkingdom_118 <- unipept_lca_118 %>%
  mutate(parent = "root") %>%
  group_by(superkingdom_name, parent) %>%
  summarize(count = n()) %>%
  filter(superkingdom_name != "") %>%
  rename(taxon = superkingdom_name)

# Number of peptides that belong to each phylum and set superkingdom as parent
phylum_118 <- unipept_lca_118 %>%
  mutate(parent = superkingdom_name) %>%
  group_by(phylum_name, parent) %>%
  summarize(count = n()) %>%
  filter(phylum_name != "") %>%
  rename(taxon = phylum_name)

# Number of peptides that belong to each class and set phylum as parent
class_118 <- unipept_lca_118 %>%
  mutate(parent = phylum_name) %>%
  group_by(class_name, parent) %>%
  summarize(count = n()) %>%
  filter(class_name != "") %>%
  rename(taxon = class_name)

# Number of peptides that belong to each order and set class as parent
order_118 <- unipept_lca_118 %>%
  mutate(parent = class_name) %>%
  group_by(order_name, parent) %>%
  summarize(count = n()) %>%
  filter(order_name != "") %>%
  rename(taxon = order_name)

# Number of peptides that belong to each family and set order as parent
family_118 <- unipept_lca_118 %>%
  mutate(parent = order_name) %>%
  group_by(family_name, parent) %>%
  summarize(count = n()) %>%
  filter(family_name != "") %>%
  rename(taxon = family_name)

# Number of peptides that belong to each genus and set family as parent
genus_118 <- unipept_lca_118 %>%
  mutate(parent = family_name) %>%
  group_by(genus_name, parent) %>%
  summarize(count = n()) %>%
  filter(genus_name != "") %>%
  rename(taxon = genus_name)

# Number of peptides that belong to each species and set genus as parent
species_118 <- unipept_lca_118 %>%
  mutate(parent = genus_name) %>%
  group_by(species_name, parent) %>%
  summarize(count = n()) %>%
  filter(species_name != "") %>%
  rename(taxon = species_name)

# Bind rows for each level
sunburst_118 <- bind_rows(
  root_118,
  superkingdom_118,
  phylum_118,
  class_118,
  order_118,
  family_118,
  genus_118,
  species_118
)

# Plot sunburst
fig <- plot_ly(
  data = sunburst_118,
  labels = ~taxon,
  parents = ~parent,
  values = ~count,
  type = 'sunburst',
  branchvalues = 'total'
)

fig
