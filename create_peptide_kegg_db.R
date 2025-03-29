# Load libraries
library(dplyr)
library(KEGGREST)

# ==================================
# KEGG annotation
# ==================================
# Select proteins with KEGG annotation
# ==================================
# Working directory
work_path <- getwd()

# Load KEGG_UniRef90.idmapping.dat fractions
#for (file in list.files(work_path)) {
#  if (substr(file, 1, 7) == "UR2KEGG"){
#    df_name <- file
#    df <- read.delim(paste0(work_path, "/", df_name),
#                     header = FALSE)
#    df <- df %>% rename(ID = V1,
#                        AnType = V2,
#                        AnID = V3) #rename columns
#    df <- df %>% 
#      pivot_wider(names_from = AnType,
#                  values_from = AnID) #matches kegg to UniRef90 by UniProtKB ID
#    df <- df %>% filter(KEGG != "NULL") #get rid of rows w/out kegg
#    df <- unnest(df, UniRef90) #UniRef90 from "list" to "character"
#    assign(df_name, df)
#  }
#}
kegg_ko_gene <- keggLink("ko", "cmy:102947773")
keggLink("cmy:102947773")


# Load libraries
#library(dplyr)
#library(tidyr)
#
#UR2KEGGaa <- read.delim(paste0(work_path, "/", "UR2KEGGaa"),
#                        header = FALSE)
#colnames(UR2KEGGaa) <- c("ID", "AnType", "AnID")  #rename columns
#
#UR2KEGGaa <- UR2KEGGaa %>% 
#  pivot_wider(names_from = AnType,
#              values_from = AnID) #matches kegg terms to UniRef90 terms using UniProtKB ID
#UR2KEGGaa <- UR2KEGGaa %>% filter(KEGG != "NULL") #get rid of rows w/out kegg
#UR2KEGGaa <- unnest(UR2KEGGaa, UniRef90) #UniRef90 from "list" to "character"
#UR2KEGGaa$KO <- NULL
#for (i in 1:nrow(UR2KEGGaa)){
#  print(i)
#  UR2KEGGaa[[i, "KO"]] <- keggLink("ko", UR2KEGGaa[i, "KEGG"])
#}


UR2KEGGaa <- read.delim(file.path(work_path, "UR2KEGGaa"), header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  filter(AnType == "KEGG" & AnID != "NULL") %>%
  group_by(ID) %>%
  summarize(KEGG = first(AnID), UniRef90 = first(ID), .groups = "drop") %>%
  mutate(KO = sapply(KEGG, function(x) {
    Sys.sleep(1)
    keggLink("ko", x)
  }))

