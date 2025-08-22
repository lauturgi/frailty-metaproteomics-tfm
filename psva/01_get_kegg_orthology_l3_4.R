# Load libraries
library(dplyr)
library(stringr)
library(plotly)

# Get the current working directory
work_path <- getwd()

# Create dataframe to store KO level 4 and 3
ko_df <- data.frame()

# KO levels downloaded from: https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001&format=htext&filedir=
# Read ko00001.keg as a vector of lines
lines <- readLines(paste0(work_path, "/psva/data/ko00001.keg"))

for (line in lines) {
  # If line starts with "C", KO Level 3 (L3)
  if (startsWith(line, "C")) {
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
      L3 = L3,
      L3_DESC = L3_DESC,
      L4 = line_s[[1]][1],
      L4_DESC = line_s[[1]][2]
    ))
  }
}

save(ko_df, file = paste0(work_path, "/psva/data/kegg_ko_l3_4.RData"))
