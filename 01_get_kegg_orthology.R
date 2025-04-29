# Load libraries
library(dplyr)
library(stringr)

# Get the current working directory
work_path <- getwd()

# KO levels downloaded from: https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001&format=htext&filedir=
# Read .keg file as a vector of lines of text
lines <- readLines(paste0("/mnt/d/proteomica/fragilidad/Pepfunk/kegg/",
                          "ko00001.keg"))

# Create dataframe to store KO level 4 and 3
ko_df <- data.frame()

# Read ko00001.keg as a vector of lines
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
      L4 = line_s[[1]][1],
      L4_DESC = line_s[[1]][2],
      L3 = L3,
      L3_DESC = L3_DESC
    ))
  }
}

save(ko_df, file = paste0(work_path, "/data/kegg_ko_l3_4.RData"))

# Compare with kegg_L3.txt (old) from pepfunk
pepfunk_path <- "/mnt/d/proteomica/fragilidad/Pepfunk/database/"
kegg_old <- read.csv2(paste0(pepfunk_path,"kegg_L3.txt"), sep = "\t", 
                      header = FALSE, colClasses = "character")
colnames(kegg_old) <- c("L4", "L4_DESC", "L3", "L3_DESC")

kegg_check <- merge(ko_df, kegg_old, by = c("L4", "L3"), all = TRUE)

kegg_check <- kegg_check[is.na(kegg_check$L3_DESC.x),]
unique(kegg_check$L3)

# "00072" --> Deleted; merged into 00650
# "00281" --> Deleted; merged into 00907
# "00260" --> K00050 missing
# "00630" --> K01450 and K00050 missing
# "99980" --> KO moved to other pathway or bride
# "99982" --> K00183, K18500, K18501 missing, K00184, K00185, K04755 moved to pathway, K12262, K16183 moved to 99980, K23302, K23303 moved to 02000
# "00400" --> K15652 moved to 01053, K00210 missing
# "00401" --> K00210 missing
# "00410" "00640" "00471" "00472" "00520" "00190"
# "99996" "00590" "00360" "00071" "00380" "00627" "00965" "09113" "00670" "02000" "00473" "00230" "00540" "00140"
# "02042" "00940" "00760" "99992" "00270" "00430" "00710" "00903" "00999" "00052" "99981" "00860" "03400" "02022"
# "02035" "00500" "00010" "00051" "03029" "03021" "99986" "99977" "99973" "99979" "00240" "00310" "99995" "99985"
# "99975" "03036" "99993" "99988" "99978" "04515" "99997" "99994" "99974" "04031" "03110" "04091" "03019" "03032"
# "02044" "00231" "01058" "00330" "00440" "03009" "03012" "00073" "00562" "01054" "04131" "99987" "00280" "03200"
# "99983" "04040" "04052"

# Compare with ko group by pathways and bride from KEGGREST
load(paste0(work_path, "/data/kegg_ko_path_old.RData"))

kegg_check_2 <- merge(ko_df, kegg_ko_path, by.x = c("L4", "L3"), 
                    by.y = c("KEGG", "PATHWAY"), all = TRUE)

kegg_check_2 <- kegg_check_2[is.na(kegg_check_2$L3_DESC),]
unique(kegg_check_2$L3)
# "01000(L2 --> NO)" "01100" "01110" "01120" "01220" "01240" "01230" "01250" "01200" "01212" "01210" "01232" "01310" "01320"