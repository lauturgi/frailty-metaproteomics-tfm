# Load libraries
library(dplyr)
library(stringr)
library(plotly)

# Get the current working directory
work_path <- getwd()

# Create dataframe to store KO levels
ko_df <- data.frame()

# KO levels downloaded from: https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001&format=htext&filedir=
# Read ko00001.keg as a vector of lines
lines <- readLines(paste0(work_path, "/psva/data/ko00001.keg"))

for (line in lines) {
  # If line starts with "A", KO Level 1 (L1)
  if (startsWith(line, "A")) {
    # Remove "A" and leading spaces
    line <- str_trim(sub("^A\\s*", "", line))
    # Split line to extract L1 code and description
    line_s <- str_split(line, "\\s+", n = 2)
    L1 <- line_s[[1]][1]
    L1_DESC <- line_s[[1]][2]
  }
  # If line starts with "B", KO Level 2 (L2)
  else if (startsWith(line, "B")){
    # Remove "B" and leading spaces
    line <- str_trim(sub("^B\\s*", "", line))
    # Split line to extract L1 code and description
    line_s <- str_split(line, "\\s+", n = 2)
    L2 <- line_s[[1]][1]
    L2_DESC <- line_s[[1]][2]
  }
  # If line starts with "C", KO Level 3 (L3)
  else if (startsWith(line, "C")) {
    # Remove "C" and leading spaces
    line <- str_trim(sub("^C\\s*", "", line))
    # Split line to extract L3 code and description
    line_s <- str_split(line, "\\s+", n = 2)
    L3 <- line_s[[1]][1]
    L3_DESC <- line_s[[1]][2]
  }
  # If line starts with "D", KO Level 4 (L4) within the current L3
  else if (startsWith(line, "D")) {
    # Remove "D" and leading spaces
    line <- str_trim(sub("^D\\s*", "", line))
    # Split line to extract L4 code and description
    line_s <- str_split(line, "\\s+", n = 2)
    # Bind row to ko_df
    ko_df <- bind_rows(ko_df, data.frame(
      L1 = L1,
      L1_DESC = L1_DESC,
      L2 = L2,
      L2_DESC = L2_DESC,
      L3 = L3,
      L3_DESC = L3_DESC,
      L4 = line_s[[1]][1],
      L4_DESC = line_s[[1]][2]
    ))
  }
}


# Filter ko_df to keep KO hierarchy for L3 = 00430 
ko_sub <- ko_df %>%
  filter(L3 == "00430") %>%
  mutate(
    L1_DESC = paste(L1, L1_DESC, sep = "; "),
    L2_DESC = paste(L2, L2_DESC, sep = "; "),
    L3_DESC = paste(L3, L3_DESC, sep = "; "),
    L4_DESC = paste(L4, L4_DESC, sep = "; ")
  )

# List of unique nodes
nodes <- unique(c(ko_sub$L1_DESC, ko_sub$L2_DESC, ko_sub$L3_DESC, ko_sub$L4_DESC))

# Index for each node
node_index <- setNames(seq_along(nodes) - 1, nodes)

# Links between levels
l1_l2 <- ko_sub %>%
  distinct(L1_DESC, L2_DESC) %>%
  mutate(source = node_index[L1_DESC],
         target = node_index[L2_DESC])

l2_l3 <- ko_sub %>%
  distinct(L2_DESC, L3_DESC) %>%
  mutate(source = node_index[L2_DESC],
         target = node_index[L3_DESC])

l3_l4 <- ko_sub %>%
  distinct(L3_DESC, L4_DESC) %>%
  mutate(source = node_index[L3_DESC],
         target = node_index[L4_DESC])

# Bind all links
links <- bind_rows(l1_l2, l2_l3, l3_l4) %>%
  mutate(value = 1)


fig <- plot_ly(
  type = "sankey",
  orientation = "h",
  node = list(
    label = nodes
  ),
  link = list(
    source = links$source,
    target = links$target,
    value = links$value
  )
)

fig <- fig %>% layout(
  font = list(
    size = 14,
    color = 'black'
  )
)

fig
