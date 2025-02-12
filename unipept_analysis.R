library(plotly)
library(dplyr)
library(htmlwidgets)

# Working directory
work_path <- "/home/lauturgi/workspace/tfm"

# Unipept LCA analysis output directory
unipept_path <- "/mnt/d/proteomica/fragilidad/datos/PeptideFilesUnipept/lca_output_no_dup/"

# Read Unipept LCA output files from frailty cohort
nrow_samples <- c()  # Vector to store number of peptides of each sample
for (file in list.files(unipept_path)) {
  replicate <- substring(file, 4, 6)
  if (replicate %in% metadata_ft$replicate) {
    df_name <- paste0("unipept_lca_no_dup_", replicate)
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

# List of LCA dataframes without duplicates
unipept_dfs <- ls(pattern = "unipept_lca_no_dup")

# Loop over the dataframes
for (i in 1:length(unipept_dfs)) {
  # Get replicate and dataframe from sample i
  replicate <- substr(unipept_dfs[i], 20, 22)
  df <- get(unipept_dfs[i])

  # Group by superkingdom, phylum, class, order, family, genus and species and
  # count number of peptides that belong to each LCA
  taxonomy <- df %>%
    group_by(superkingdom_name, phylum_name, class_name, order_name, family_name,
            genus_name, species_name) %>%
    summarize(count = n())

  # Create taxonomy for sunburst plot
  # Number of total peptides that constitute the root
  root <- df %>%
    mutate(parent = "", taxon = "root") %>%
    group_by(taxon, parent) %>%
    summarize(count = n())

  # Number of peptides that belong to each superkingdom and set root as parent
  superkingdom <- df %>%
    mutate(parent = "root") %>%
    group_by(superkingdom_name, parent) %>%
    summarize(count = n()) %>%
    filter(superkingdom_name != "") %>%
    rename(taxon = superkingdom_name)

  # Number of peptides that belong to each phylum and set superkingdom as parent
  phylum <- df %>%
    mutate(parent = superkingdom_name) %>%
    group_by(phylum_name, parent) %>%
    summarize(count = n()) %>%
    filter(phylum_name != "") %>%
    rename(taxon = phylum_name)

  # Number of peptides that belong to each class and set phylum as parent
  class <- df %>%
    mutate(parent = phylum_name) %>%
    group_by(class_name, parent) %>%
    summarize(count = n()) %>%
    filter(class_name != "") %>%
    rename(taxon = class_name)

  # Number of peptides that belong to each order and set class as parent
  order <- df %>%
    mutate(parent = class_name) %>%
    group_by(order_name, parent) %>%
    summarize(count = n()) %>%
    filter(order_name != "") %>%
    rename(taxon = order_name)

  # Number of peptides that belong to each family and set order as parent
  family <- df %>%
    mutate(parent = order_name) %>%
    group_by(family_name, parent) %>%
    summarize(count = n()) %>%
    filter(family_name != "") %>%
    rename(taxon = family_name)

  # Number of peptides that belong to each genus and set family as parent
  genus <- df %>%
    mutate(parent = family_name) %>%
    group_by(genus_name, parent) %>%
    summarize(count = n()) %>%
    filter(genus_name != "") %>%
    rename(taxon = genus_name)

  # Number of peptides that belong to each species and set genus as parent
  species <- df %>%
    mutate(parent = genus_name) %>%
    group_by(species_name, parent) %>%
    summarize(count = n()) %>%
    filter(species_name != "") %>%
    rename(taxon = species_name)

  # Bind rows for each level
  sunburst <- bind_rows(
    root,
    superkingdom,
    phylum,
    class,
    order,
    family,
    genus,
    species
  )

  # Plot sunburst
  fig <- plot_ly(
    data = sunburst,
    labels = ~taxon,
    parents = ~parent,
    values = ~count,
    type = 'sunburst',
    branchvalues = 'total'
  ) %>%
    layout(
      title = paste("Sunburst Plot - Replicate", replicate)
    )

  # Save the plot
  file_name <- paste0(work_path,"/plots/unipept_analysis/sunburst_replicate_", replicate, ".png")
  saveWidget(
    widget = fig,
    file = paste0(work_path,"/plots/unipept_analysis/sunburst_",replicate,".html"),
    selfcontained = TRUE
  )
}
