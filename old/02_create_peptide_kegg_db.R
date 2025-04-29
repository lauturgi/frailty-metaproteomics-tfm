# Load libraries
library(dplyr)
library(tidyr)
library(stringr)


# ==================================
# KEGG annotation
# ==================================
# Select proteins with KEGG annotation
# ==================================
# Working directory
work_path <- getwd()

# Load KEGG_UniRef90.idmapping.dat fractions
UR2KEGGaa <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGaa"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGab <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGab"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGac <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGac"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGad <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGad"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGae <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGae"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGaf <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGaf"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGag <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGag"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGah <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGah"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGai <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGai"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGaj <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGaj"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGak <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGak"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGal <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGal"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGam <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGam"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGan <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGan"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGao <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGao"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGap <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGap"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGaq <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGaq"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGar <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGar"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGas <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGas"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGat <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGat"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGau <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGau"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGav <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGav"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGaw <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGaw"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGax <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGax"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGay <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGay"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGaz <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGaz"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

print("UR2KEGGa ingested.")

UR2KEGGba <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGba"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbb <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbb"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbc <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbc"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbd <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbd"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbe <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbe"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbf <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbf"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbg <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbg"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbh <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbh"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbi <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbi"), header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbj <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbj"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbk <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbk"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbl <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbl"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbm <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbm"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbn <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbn"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbo <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbo"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbp <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbp"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbq <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbq"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbr <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbr"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbs <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbs"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbt <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbt"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbu <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbu"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbv <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbv"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbw <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbw"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbx <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbx"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGby <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGby"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGbz <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGbz"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

print("UR2KEGGb ingested.")

UR2KEGGca <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGca"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcb <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcb"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcc <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcc"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcd <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcd"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGce <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGce"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcf <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcf"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcg <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcg"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGch <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGch"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGci <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGci"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcj <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcj"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGck <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGck"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcl <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcl"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcm <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcm"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcn <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcn"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGco <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGco"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcp <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcp"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcq <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcq"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcr <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcr"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcs <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcs"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGct <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGct"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcu <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcu"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcv <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcv"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcw <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcw"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcx <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcx"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcy <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcy"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGcz <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGcz"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

print("UR2KEGGc ingested.")

UR2KEGGda <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGda"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdb <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdb"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdc <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdc"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdd <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdd"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGde <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGde"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdf <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdf"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdg <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdg"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdh <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdh"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdi <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdi"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdj <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdj"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdk <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdk"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdl <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdl"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdm <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdm"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdn <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdn"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdo <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdo"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdp <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdp"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdq <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdq"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdr <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdr"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGds <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGds"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdt <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdt"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdu <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdu"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdv <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdv"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdw <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdw"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdx <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdx"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdy <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdy"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGdz <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGdz"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

print("UR2KEGGd ingested.")

UR2KEGGea <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGea"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGeb <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGeb"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGec <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGec"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGed <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGed"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGee <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGee"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGef <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGef"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGeg <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGeg"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGeh <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGeh"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGei <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGei"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGej <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGej"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGek <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGek"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGel <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGel"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGem <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGem"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGen <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGen"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGeo <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGeo"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGep <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGep"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGeq <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGeq"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGer <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGer"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGes <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGes"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGet <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGet"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGeu <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGeu"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGev <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGev"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGew <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGew"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGex <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGex"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGey <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGey"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGez <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGez"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

print("UR2KEGGe ingested.")

UR2KEGGfa <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfa"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfb <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfb"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfc <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfc"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfd <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfd"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfe <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfe"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGff <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGff"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfg <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfg"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfh <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfh"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfi <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfi"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfj <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfj"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfk <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfk"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfl <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfl"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfm <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfm"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfn <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfn"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfo <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfo"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfp <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfp"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfq <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfq"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfr <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfr"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfs <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfs"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGft <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGft"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfu <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfu"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfv <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfv"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfw <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfw"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfx <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfx"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfy <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfy"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGfz <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGfz"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

print("UR2KEGGf ingested.")

UR2KEGGga <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGga"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgb <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgb"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgc <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgc"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgd <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgd"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGge <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGge"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgf <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgf"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgg <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgg"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgh <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgh"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgi <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgi"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgj <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgj"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgk <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgk"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgl <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgl"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgm <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgm"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgn <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgn"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgo <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgo"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgp <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgp"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgq <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgq"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgr <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgr"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgs <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgs"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgt <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgt"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgu <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgu"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgv <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgv"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgw <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgw"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgx <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgx"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgy <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgy"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGgz <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGgz"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

print("UR2KEGGg ingested.")

UR2KEGGha <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGha"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhb <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhb"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhc <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhc"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhd <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhd"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhe <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhe"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhf <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhf"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhg <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhg"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhh <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhh"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhi <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhi"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhj <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhj"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhk <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhk"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhl <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhl"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhm <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhm"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhn <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhn"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGho <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGho"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhp <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhp"), 
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhq <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhq"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhr <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhr"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhs <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhs"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGht <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGht"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhu <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhu"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhv <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhv"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhw <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhw"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhx <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhx"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhy <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhy"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGhz <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGhz"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

print("UR2KEGGh ingested.")

UR2KEGGia <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGia"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGib <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGib"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGic <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGic"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGid <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGid"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGie <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGie"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGif <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGif"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGig <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGig"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGih <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGih"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGii <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGii"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGij <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGij"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGik <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGik"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGil <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGil"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGim <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGim"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGin <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGin"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGio <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGio"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGip <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGip"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGiq <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGiq"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGir <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGir"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGis <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGis"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGit <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGit"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGiu <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGiu"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGiv <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGiv"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGiw <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGiw"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGix <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGix"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGiy <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGiy"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGiz <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGiz"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

print("UR2KEGGi ingested.")

UR2KEGGja <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGja"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjb <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjb"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjc <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjc"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjd <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjd"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGje <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGje"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjf <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjf"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjg <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjg"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjh <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjh"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGji <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGji"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjj <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjj"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjk <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjk"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjl <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjl"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjm <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjm"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjn <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjn"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjo <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjo"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjp <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjp"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjq <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjq"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjr <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjr"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjs <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjs"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjt <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjt"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGju <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGju"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjv <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjv"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

#UR2KEGGjw <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjw"),
#                        header = FALSE) %>%
#  rename(ID = V1, AnType = V2, AnID = V3) %>%
#  pivot_wider(names_from = AnType, values_from = AnID)  %>%
#  filter(KEGG != "NULL")

UR2KEGGjx <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjx"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjy <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjy"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGjz <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGjz"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

print("UR2KEGGj ingested.")

UR2KEGGka <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGka"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGkb <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGkb"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGkc <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGkc"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGkd <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGkd"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGke <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGke"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGkf <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGkf"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGkg <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGkg"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGkh <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGkh"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGki <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGki"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGkj <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGkj"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGkk <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGkk"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGkl <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGkl"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGkm <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGkm"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGkn <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGkn"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGko <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGko"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGkp <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGkp"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGkq <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGkq"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

UR2KEGGkr <- read.delim(file.path(work_path, "/UR2KEGG/UR2KEGGkr"),
                        header = FALSE) %>%
  rename(ID = V1, AnType = V2, AnID = V3) %>%
  pivot_wider(names_from = AnType, values_from = AnID, values_fn = list) %>%
  filter(KEGG != "NULL")

print("UR2KEGGk ingested.")

# Unnest KEGG
UR2KEGGaa <- unnest(UR2KEGGaa, c(UniRef90, KEGG))
UR2KEGGab <- unnest(UR2KEGGab, c(UniRef90, KEGG))
UR2KEGGac <- unnest(UR2KEGGac, c(UniRef90, KEGG))
UR2KEGGad <- unnest(UR2KEGGad, c(UniRef90, KEGG))
UR2KEGGae <- unnest(UR2KEGGae, c(UniRef90, KEGG))
UR2KEGGaf <- unnest(UR2KEGGaf, c(UniRef90, KEGG))
UR2KEGGag <- unnest(UR2KEGGag, c(UniRef90, KEGG))
UR2KEGGah <- unnest(UR2KEGGah, c(UniRef90, KEGG))
UR2KEGGai <- unnest(UR2KEGGai, c(UniRef90, KEGG))
UR2KEGGaj <- unnest(UR2KEGGaj, c(UniRef90, KEGG))
UR2KEGGak <- unnest(UR2KEGGak, c(UniRef90, KEGG))
UR2KEGGal <- unnest(UR2KEGGal, c(UniRef90, KEGG))
UR2KEGGam <- unnest(UR2KEGGam, c(UniRef90, KEGG))
UR2KEGGan <- unnest(UR2KEGGan, c(UniRef90, KEGG))
UR2KEGGao <- unnest(UR2KEGGao, c(UniRef90, KEGG))
UR2KEGGap <- unnest(UR2KEGGap, c(UniRef90, KEGG))
UR2KEGGaq <- unnest(UR2KEGGaq, c(UniRef90, KEGG))
UR2KEGGar <- unnest(UR2KEGGar, c(UniRef90, KEGG))
UR2KEGGas <- unnest(UR2KEGGas, c(UniRef90, KEGG))
UR2KEGGat <- unnest(UR2KEGGat, c(UniRef90, KEGG))
UR2KEGGau <- unnest(UR2KEGGau, c(UniRef90, KEGG))
UR2KEGGav <- unnest(UR2KEGGav, c(UniRef90, KEGG))
UR2KEGGaw <- unnest(UR2KEGGaw, c(UniRef90, KEGG))
UR2KEGGax <- unnest(UR2KEGGax, c(UniRef90, KEGG))
UR2KEGGay <- unnest(UR2KEGGay, c(UniRef90, KEGG))
UR2KEGGaz <- unnest(UR2KEGGaz, c(UniRef90, KEGG))
UR2KEGGba <- unnest(UR2KEGGba, c(UniRef90, KEGG))
UR2KEGGbb <- unnest(UR2KEGGbb, c(UniRef90, KEGG))
UR2KEGGbc <- unnest(UR2KEGGbc, c(UniRef90, KEGG))
UR2KEGGbd <- unnest(UR2KEGGbd, c(UniRef90, KEGG))
UR2KEGGbe <- unnest(UR2KEGGbe, c(UniRef90, KEGG))
UR2KEGGbf <- unnest(UR2KEGGbf, c(UniRef90, KEGG))
UR2KEGGbg <- unnest(UR2KEGGbg, c(UniRef90, KEGG))
UR2KEGGbh <- unnest(UR2KEGGbh, c(UniRef90, KEGG))
UR2KEGGbi <- unnest(UR2KEGGbi, c(UniRef90, KEGG))
UR2KEGGbj <- unnest(UR2KEGGbj, c(UniRef90, KEGG))
UR2KEGGbk <- unnest(UR2KEGGbk, c(UniRef90, KEGG))
UR2KEGGbl <- unnest(UR2KEGGbl, c(UniRef90, KEGG))
UR2KEGGbm <- unnest(UR2KEGGbm, c(UniRef90, KEGG))
UR2KEGGbn <- unnest(UR2KEGGbn, c(UniRef90, KEGG))
UR2KEGGbo <- unnest(UR2KEGGbo, c(UniRef90, KEGG))
UR2KEGGbp <- unnest(UR2KEGGbp, c(UniRef90, KEGG))
UR2KEGGbq <- unnest(UR2KEGGbq, c(UniRef90, KEGG))
UR2KEGGbr <- unnest(UR2KEGGbr, c(UniRef90, KEGG))
UR2KEGGbs <- unnest(UR2KEGGbs, c(UniRef90, KEGG))
UR2KEGGbt <- unnest(UR2KEGGbt, c(UniRef90, KEGG))
UR2KEGGbu <- unnest(UR2KEGGbu, c(UniRef90, KEGG))
UR2KEGGbv <- unnest(UR2KEGGbv, c(UniRef90, KEGG))
UR2KEGGbw <- unnest(UR2KEGGbw, c(UniRef90, KEGG))
UR2KEGGbx <- unnest(UR2KEGGbx, c(UniRef90, KEGG))
UR2KEGGby <- unnest(UR2KEGGby, c(UniRef90, KEGG))
UR2KEGGbz <- unnest(UR2KEGGbz, c(UniRef90, KEGG))
UR2KEGGca <- unnest(UR2KEGGca, c(UniRef90, KEGG))
UR2KEGGcb <- unnest(UR2KEGGcb, c(UniRef90, KEGG))
UR2KEGGcc <- unnest(UR2KEGGcc, c(UniRef90, KEGG))
UR2KEGGcd <- unnest(UR2KEGGcd, c(UniRef90, KEGG))
UR2KEGGce <- unnest(UR2KEGGce, c(UniRef90, KEGG))
UR2KEGGcf <- unnest(UR2KEGGcf, c(UniRef90, KEGG))
UR2KEGGcg <- unnest(UR2KEGGcg, c(UniRef90, KEGG))
UR2KEGGch <- unnest(UR2KEGGch, c(UniRef90, KEGG))
UR2KEGGci <- unnest(UR2KEGGci, c(UniRef90, KEGG))
UR2KEGGcj <- unnest(UR2KEGGcj, c(UniRef90, KEGG))
UR2KEGGck <- unnest(UR2KEGGck, c(UniRef90, KEGG))
UR2KEGGcl <- unnest(UR2KEGGcl, c(UniRef90, KEGG))
UR2KEGGcm <- unnest(UR2KEGGcm, c(UniRef90, KEGG))
UR2KEGGcn <- unnest(UR2KEGGcn, c(UniRef90, KEGG))
UR2KEGGco <- unnest(UR2KEGGco, c(UniRef90, KEGG))
UR2KEGGcp <- unnest(UR2KEGGcp, c(UniRef90, KEGG))
UR2KEGGcq <- unnest(UR2KEGGcq, c(UniRef90, KEGG))
UR2KEGGcr <- unnest(UR2KEGGcr, c(UniRef90, KEGG))
UR2KEGGcs <- unnest(UR2KEGGcs, c(UniRef90, KEGG))
UR2KEGGct <- unnest(UR2KEGGct, c(UniRef90, KEGG))
UR2KEGGcu <- unnest(UR2KEGGcu, c(UniRef90, KEGG))
UR2KEGGcv <- unnest(UR2KEGGcv, c(UniRef90, KEGG))
UR2KEGGcw <- unnest(UR2KEGGcw, c(UniRef90, KEGG))
UR2KEGGcx <- unnest(UR2KEGGcx, c(UniRef90, KEGG))
UR2KEGGcy <- unnest(UR2KEGGcy, c(UniRef90, KEGG))
UR2KEGGcz <- unnest(UR2KEGGcz, c(UniRef90, KEGG))
UR2KEGGda <- unnest(UR2KEGGda, c(UniRef90, KEGG))
UR2KEGGdb <- unnest(UR2KEGGdb, c(UniRef90, KEGG))
UR2KEGGdc <- unnest(UR2KEGGdc, c(UniRef90, KEGG))
UR2KEGGdd <- unnest(UR2KEGGdd, c(UniRef90, KEGG))
UR2KEGGde <- unnest(UR2KEGGde, c(UniRef90, KEGG))
UR2KEGGdf <- unnest(UR2KEGGdf, c(UniRef90, KEGG))
UR2KEGGdg <- unnest(UR2KEGGdg, c(UniRef90, KEGG))
UR2KEGGdh <- unnest(UR2KEGGdh, c(UniRef90, KEGG))
UR2KEGGdi <- unnest(UR2KEGGdi, c(UniRef90, KEGG))
UR2KEGGdj <- unnest(UR2KEGGdj, c(UniRef90, KEGG))
UR2KEGGdk <- unnest(UR2KEGGdk, c(UniRef90, KEGG))
UR2KEGGdl <- unnest(UR2KEGGdl, c(UniRef90, KEGG))
UR2KEGGdm <- unnest(UR2KEGGdm, c(UniRef90, KEGG))
UR2KEGGdn <- unnest(UR2KEGGdn, c(UniRef90, KEGG))
UR2KEGGdo <- unnest(UR2KEGGdo, c(UniRef90, KEGG))
UR2KEGGdp <- unnest(UR2KEGGdp, c(UniRef90, KEGG))
UR2KEGGdq <- unnest(UR2KEGGdq, c(UniRef90, KEGG))
UR2KEGGdr <- unnest(UR2KEGGdr, c(UniRef90, KEGG))
UR2KEGGds <- unnest(UR2KEGGds, c(UniRef90, KEGG))
UR2KEGGdt <- unnest(UR2KEGGdt, c(UniRef90, KEGG))
UR2KEGGdu <- unnest(UR2KEGGdu, c(UniRef90, KEGG))
UR2KEGGdv <- unnest(UR2KEGGdv, c(UniRef90, KEGG))
UR2KEGGdw <- unnest(UR2KEGGdw, c(UniRef90, KEGG))
UR2KEGGdx <- unnest(UR2KEGGdx, c(UniRef90, KEGG))
UR2KEGGdy <- unnest(UR2KEGGdy, c(UniRef90, KEGG))
UR2KEGGdz <- unnest(UR2KEGGdz, c(UniRef90, KEGG))
UR2KEGGea <- unnest(UR2KEGGea, c(UniRef90, KEGG))
UR2KEGGeb <- unnest(UR2KEGGeb, c(UniRef90, KEGG))
UR2KEGGec <- unnest(UR2KEGGec, c(UniRef90, KEGG))
UR2KEGGed <- unnest(UR2KEGGed, c(UniRef90, KEGG))
UR2KEGGee <- unnest(UR2KEGGee, c(UniRef90, KEGG))
UR2KEGGef <- unnest(UR2KEGGef, c(UniRef90, KEGG))
UR2KEGGeg <- unnest(UR2KEGGeg, c(UniRef90, KEGG))
UR2KEGGeh <- unnest(UR2KEGGeh, c(UniRef90, KEGG))
UR2KEGGei <- unnest(UR2KEGGei, c(UniRef90, KEGG))
UR2KEGGej <- unnest(UR2KEGGej, c(UniRef90, KEGG))
UR2KEGGek <- unnest(UR2KEGGek, c(UniRef90, KEGG))
UR2KEGGel <- unnest(UR2KEGGel, c(UniRef90, KEGG))
UR2KEGGem <- unnest(UR2KEGGem, c(UniRef90, KEGG))
UR2KEGGen <- unnest(UR2KEGGen, c(UniRef90, KEGG))
UR2KEGGeo <- unnest(UR2KEGGeo, c(UniRef90, KEGG))
UR2KEGGep <- unnest(UR2KEGGep, c(UniRef90, KEGG))
UR2KEGGeq <- unnest(UR2KEGGeq, c(UniRef90, KEGG))
UR2KEGGer <- unnest(UR2KEGGer, c(UniRef90, KEGG))
UR2KEGGes <- unnest(UR2KEGGes, c(UniRef90, KEGG))
UR2KEGGet <- unnest(UR2KEGGet, c(UniRef90, KEGG))
UR2KEGGeu <- unnest(UR2KEGGeu, c(UniRef90, KEGG))
UR2KEGGev <- unnest(UR2KEGGev, c(UniRef90, KEGG))
UR2KEGGew <- unnest(UR2KEGGew, c(UniRef90, KEGG))
UR2KEGGex <- unnest(UR2KEGGex, c(UniRef90, KEGG))
UR2KEGGey <- unnest(UR2KEGGey, c(UniRef90, KEGG))
UR2KEGGez <- unnest(UR2KEGGez, c(UniRef90, KEGG))
UR2KEGGfa <- unnest(UR2KEGGfa, c(UniRef90, KEGG))
UR2KEGGfb <- unnest(UR2KEGGfb, c(UniRef90, KEGG))
UR2KEGGfc <- unnest(UR2KEGGfc, c(UniRef90, KEGG))
UR2KEGGfd <- unnest(UR2KEGGfd, c(UniRef90, KEGG))
UR2KEGGfe <- unnest(UR2KEGGfe, c(UniRef90, KEGG))
UR2KEGGff <- unnest(UR2KEGGff, c(UniRef90, KEGG))
UR2KEGGfg <- unnest(UR2KEGGfg, c(UniRef90, KEGG))
UR2KEGGfh <- unnest(UR2KEGGfh, c(UniRef90, KEGG))
UR2KEGGfi <- unnest(UR2KEGGfi, c(UniRef90, KEGG))
UR2KEGGfj <- unnest(UR2KEGGfj, c(UniRef90, KEGG))
UR2KEGGfk <- unnest(UR2KEGGfk, c(UniRef90, KEGG))
UR2KEGGfl <- unnest(UR2KEGGfl, c(UniRef90, KEGG))
UR2KEGGfm <- unnest(UR2KEGGfm, c(UniRef90, KEGG))
UR2KEGGfn <- unnest(UR2KEGGfn, c(UniRef90, KEGG))
UR2KEGGfo <- unnest(UR2KEGGfo, c(UniRef90, KEGG))
UR2KEGGfp <- unnest(UR2KEGGfp, c(UniRef90, KEGG))
UR2KEGGfq <- unnest(UR2KEGGfq, c(UniRef90, KEGG))
UR2KEGGfr <- unnest(UR2KEGGfr, c(UniRef90, KEGG))
UR2KEGGfs <- unnest(UR2KEGGfs, c(UniRef90, KEGG))
UR2KEGGft <- unnest(UR2KEGGft, c(UniRef90, KEGG))
UR2KEGGfu <- unnest(UR2KEGGfu, c(UniRef90, KEGG))
UR2KEGGfv <- unnest(UR2KEGGfv, c(UniRef90, KEGG))
UR2KEGGfw <- unnest(UR2KEGGfw, c(UniRef90, KEGG))
UR2KEGGfx <- unnest(UR2KEGGfx, c(UniRef90, KEGG))
UR2KEGGfy <- unnest(UR2KEGGfy, c(UniRef90, KEGG))
UR2KEGGfz <- unnest(UR2KEGGfz, c(UniRef90, KEGG))
UR2KEGGga <- unnest(UR2KEGGga, c(UniRef90, KEGG))
UR2KEGGgb <- unnest(UR2KEGGgb, c(UniRef90, KEGG))
UR2KEGGgc <- unnest(UR2KEGGgc, c(UniRef90, KEGG))
UR2KEGGgd <- unnest(UR2KEGGgd, c(UniRef90, KEGG))
UR2KEGGge <- unnest(UR2KEGGge, c(UniRef90, KEGG))
UR2KEGGgf <- unnest(UR2KEGGgf, c(UniRef90, KEGG))
UR2KEGGgg <- unnest(UR2KEGGgg, c(UniRef90, KEGG))
UR2KEGGgh <- unnest(UR2KEGGgh, c(UniRef90, KEGG))
UR2KEGGgi <- unnest(UR2KEGGgi, c(UniRef90, KEGG))
UR2KEGGgj <- unnest(UR2KEGGgj, c(UniRef90, KEGG))
UR2KEGGgk <- unnest(UR2KEGGgk, c(UniRef90, KEGG))
UR2KEGGgl <- unnest(UR2KEGGgl, c(UniRef90, KEGG))
UR2KEGGgm <- unnest(UR2KEGGgm, c(UniRef90, KEGG))
UR2KEGGgn <- unnest(UR2KEGGgn, c(UniRef90, KEGG))
UR2KEGGgo <- unnest(UR2KEGGgo, c(UniRef90, KEGG))
UR2KEGGgp <- unnest(UR2KEGGgp, c(UniRef90, KEGG))
UR2KEGGgq <- unnest(UR2KEGGgq, c(UniRef90, KEGG))
UR2KEGGgr <- unnest(UR2KEGGgr, c(UniRef90, KEGG))
UR2KEGGgs <- unnest(UR2KEGGgs, c(UniRef90, KEGG))
UR2KEGGgt <- unnest(UR2KEGGgt, c(UniRef90, KEGG))
UR2KEGGgu <- unnest(UR2KEGGgu, c(UniRef90, KEGG))
UR2KEGGgv <- unnest(UR2KEGGgv, c(UniRef90, KEGG))
UR2KEGGgw <- unnest(UR2KEGGgw, c(UniRef90, KEGG))
UR2KEGGgx <- unnest(UR2KEGGgx, c(UniRef90, KEGG))
UR2KEGGgy <- unnest(UR2KEGGgy, c(UniRef90, KEGG))
UR2KEGGgz <- unnest(UR2KEGGgz, c(UniRef90, KEGG))
UR2KEGGha <- unnest(UR2KEGGha, c(UniRef90, KEGG))
UR2KEGGhb <- unnest(UR2KEGGhb, c(UniRef90, KEGG))
UR2KEGGhc <- unnest(UR2KEGGhc, c(UniRef90, KEGG))
UR2KEGGhd <- unnest(UR2KEGGhd, c(UniRef90, KEGG))
UR2KEGGhe <- unnest(UR2KEGGhe, c(UniRef90, KEGG))
UR2KEGGhf <- unnest(UR2KEGGhf, c(UniRef90, KEGG))
UR2KEGGhg <- unnest(UR2KEGGhg, c(UniRef90, KEGG))
UR2KEGGhh <- unnest(UR2KEGGhh, c(UniRef90, KEGG))
UR2KEGGhi <- unnest(UR2KEGGhi, c(UniRef90, KEGG))
UR2KEGGhj <- unnest(UR2KEGGhj, c(UniRef90, KEGG))
UR2KEGGhk <- unnest(UR2KEGGhk, c(UniRef90, KEGG))
UR2KEGGhl <- unnest(UR2KEGGhl, c(UniRef90, KEGG))
UR2KEGGhm <- unnest(UR2KEGGhm, c(UniRef90, KEGG))
UR2KEGGhn <- unnest(UR2KEGGhn, c(UniRef90, KEGG))
UR2KEGGho <- unnest(UR2KEGGho, c(UniRef90, KEGG))
UR2KEGGhp <- unnest(UR2KEGGhp, c(UniRef90, KEGG))
UR2KEGGhq <- unnest(UR2KEGGhq, c(UniRef90, KEGG))
UR2KEGGhr <- unnest(UR2KEGGhr, c(UniRef90, KEGG))
UR2KEGGhs <- unnest(UR2KEGGhs, c(UniRef90, KEGG))
UR2KEGGht <- unnest(UR2KEGGht, c(UniRef90, KEGG))
UR2KEGGhu <- unnest(UR2KEGGhu, c(UniRef90, KEGG))
UR2KEGGhv <- unnest(UR2KEGGhv, c(UniRef90, KEGG))
UR2KEGGhw <- unnest(UR2KEGGhw, c(UniRef90, KEGG))
UR2KEGGhx <- unnest(UR2KEGGhx, c(UniRef90, KEGG))
UR2KEGGhy <- unnest(UR2KEGGhy, c(UniRef90, KEGG))
UR2KEGGhz <- unnest(UR2KEGGhz, c(UniRef90, KEGG))
UR2KEGGia <- unnest(UR2KEGGia, c(UniRef90, KEGG))
UR2KEGGib <- unnest(UR2KEGGib, c(UniRef90, KEGG))
UR2KEGGic <- unnest(UR2KEGGic, c(UniRef90, KEGG))
UR2KEGGid <- unnest(UR2KEGGid, c(UniRef90, KEGG))
UR2KEGGie <- unnest(UR2KEGGie, c(UniRef90, KEGG))
UR2KEGGif <- unnest(UR2KEGGif, c(UniRef90, KEGG))
UR2KEGGig <- unnest(UR2KEGGig, c(UniRef90, KEGG))
UR2KEGGih <- unnest(UR2KEGGih, c(UniRef90, KEGG))
UR2KEGGii <- unnest(UR2KEGGii, c(UniRef90, KEGG))
UR2KEGGij <- unnest(UR2KEGGij, c(UniRef90, KEGG))
UR2KEGGik <- unnest(UR2KEGGik, c(UniRef90, KEGG))
UR2KEGGil <- unnest(UR2KEGGil, c(UniRef90, KEGG))
UR2KEGGim <- unnest(UR2KEGGim, c(UniRef90, KEGG))
UR2KEGGin <- unnest(UR2KEGGin, c(UniRef90, KEGG))
UR2KEGGio <- unnest(UR2KEGGio, c(UniRef90, KEGG))
UR2KEGGip <- unnest(UR2KEGGip, c(UniRef90, KEGG))
UR2KEGGiq <- unnest(UR2KEGGiq, c(UniRef90, KEGG))
UR2KEGGir <- unnest(UR2KEGGir, c(UniRef90, KEGG))
UR2KEGGis <- unnest(UR2KEGGis, c(UniRef90, KEGG))
UR2KEGGit <- unnest(UR2KEGGit, c(UniRef90, KEGG))
UR2KEGGiu <- unnest(UR2KEGGiu, c(UniRef90, KEGG))
UR2KEGGiv <- unnest(UR2KEGGiv, c(UniRef90, KEGG))
UR2KEGGiw <- unnest(UR2KEGGiw, c(UniRef90, KEGG))
UR2KEGGix <- unnest(UR2KEGGix, c(UniRef90, KEGG))
UR2KEGGiy <- unnest(UR2KEGGiy, c(UniRef90, KEGG))
UR2KEGGiz <- unnest(UR2KEGGiz, c(UniRef90, KEGG))
UR2KEGGja <- unnest(UR2KEGGja, c(UniRef90, KEGG))
UR2KEGGjb <- unnest(UR2KEGGjb, c(UniRef90, KEGG))
UR2KEGGjc <- unnest(UR2KEGGjc, c(UniRef90, KEGG))
UR2KEGGjd <- unnest(UR2KEGGjd, c(UniRef90, KEGG))
UR2KEGGje <- unnest(UR2KEGGje, c(UniRef90, KEGG))
UR2KEGGjf <- unnest(UR2KEGGjf, c(UniRef90, KEGG))
UR2KEGGjg <- unnest(UR2KEGGjg, c(UniRef90, KEGG))
UR2KEGGjh <- unnest(UR2KEGGjh, c(UniRef90, KEGG))
UR2KEGGji <- unnest(UR2KEGGji, c(UniRef90, KEGG))
UR2KEGGjj <- unnest(UR2KEGGjj, c(UniRef90, KEGG))
UR2KEGGjk <- unnest(UR2KEGGjk, c(UniRef90, KEGG))
UR2KEGGjl <- unnest(UR2KEGGjl, c(UniRef90, KEGG))
UR2KEGGjm <- unnest(UR2KEGGjm, c(UniRef90, KEGG))
UR2KEGGjn <- unnest(UR2KEGGjn, c(UniRef90, KEGG))
UR2KEGGjo <- unnest(UR2KEGGjo, c(UniRef90, KEGG))
UR2KEGGjp <- unnest(UR2KEGGjp, c(UniRef90, KEGG))
UR2KEGGjq <- unnest(UR2KEGGjq, c(UniRef90, KEGG))
UR2KEGGjr <- unnest(UR2KEGGjr, c(UniRef90, KEGG))
UR2KEGGjs <- unnest(UR2KEGGjs, c(UniRef90, KEGG))
UR2KEGGjt <- unnest(UR2KEGGjt, c(UniRef90, KEGG))
UR2KEGGju <- unnest(UR2KEGGju, c(UniRef90, KEGG))
UR2KEGGjv <- unnest(UR2KEGGjv, c(UniRef90, KEGG))
#UR2KEGGjw <- unnest(UR2KEGGjw, c(UniRef90, KEGG))
UR2KEGGjx <- unnest(UR2KEGGjx, c(UniRef90, KEGG))
UR2KEGGjy <- unnest(UR2KEGGjy, c(UniRef90, KEGG))
UR2KEGGjz <- unnest(UR2KEGGjz, c(UniRef90, KEGG))
UR2KEGGka <- unnest(UR2KEGGka, c(UniRef90, KEGG))
UR2KEGGkb <- unnest(UR2KEGGkb, c(UniRef90, KEGG))
UR2KEGGkc <- unnest(UR2KEGGkc, c(UniRef90, KEGG))
UR2KEGGkd <- unnest(UR2KEGGkd, c(UniRef90, KEGG))
UR2KEGGke <- unnest(UR2KEGGke, c(UniRef90, KEGG))
UR2KEGGkf <- unnest(UR2KEGGkf, c(UniRef90, KEGG))
UR2KEGGkg <- unnest(UR2KEGGkg, c(UniRef90, KEGG))
UR2KEGGkh <- unnest(UR2KEGGkh, c(UniRef90, KEGG))
UR2KEGGki <- unnest(UR2KEGGki, c(UniRef90, KEGG))
UR2KEGGkj <- unnest(UR2KEGGkj, c(UniRef90, KEGG))
UR2KEGGkk <- unnest(UR2KEGGkk, c(UniRef90, KEGG))
UR2KEGGkl <- unnest(UR2KEGGkl, c(UniRef90, KEGG))
UR2KEGGkm <- unnest(UR2KEGGkm, c(UniRef90, KEGG))
UR2KEGGkn <- unnest(UR2KEGGkn, c(UniRef90, KEGG))
UR2KEGGko <- unnest(UR2KEGGko, c(UniRef90, KEGG))
UR2KEGGkp <- unnest(UR2KEGGkp, c(UniRef90, KEGG))
UR2KEGGkq <- unnest(UR2KEGGkq, c(UniRef90, KEGG))
UR2KEGGkr <- unnest(UR2KEGGkr, c(UniRef90, KEGG))

# Combine all the processed fractions into one dataset
UR2KEGG <- bind_rows(UR2KEGGaa, UR2KEGGab, UR2KEGGac, UR2KEGGad, UR2KEGGae, 
                     UR2KEGGaf, UR2KEGGag, UR2KEGGah, UR2KEGGai, UR2KEGGaj, 
                     UR2KEGGak, UR2KEGGal, UR2KEGGam, UR2KEGGan, UR2KEGGao,
                     UR2KEGGap, UR2KEGGaq, UR2KEGGar, UR2KEGGas, UR2KEGGat, 
                     UR2KEGGau, UR2KEGGav, UR2KEGGaw, UR2KEGGax, UR2KEGGay, 
                     UR2KEGGaz, UR2KEGGba, UR2KEGGbb, UR2KEGGbc, UR2KEGGbd,
                     UR2KEGGbe, UR2KEGGbf, UR2KEGGbg, UR2KEGGbh, UR2KEGGbi,
                     UR2KEGGbj, UR2KEGGbk, UR2KEGGbl, UR2KEGGbm, UR2KEGGbn,
                     UR2KEGGbo, UR2KEGGbp, UR2KEGGbq, UR2KEGGbr, UR2KEGGbs,
                     UR2KEGGbt, UR2KEGGbu, UR2KEGGbv, UR2KEGGbw, UR2KEGGbx,
                     UR2KEGGby, UR2KEGGbz, UR2KEGGca, UR2KEGGcb, UR2KEGGcc,
                     UR2KEGGcd, UR2KEGGce, UR2KEGGcf, UR2KEGGcg, UR2KEGGch,
                     UR2KEGGci, UR2KEGGcj, UR2KEGGck, UR2KEGGcl, UR2KEGGcm,
                     UR2KEGGcn, UR2KEGGco, UR2KEGGcp, UR2KEGGcq, UR2KEGGcr,
                     UR2KEGGcs, UR2KEGGct, UR2KEGGcu, UR2KEGGcv, UR2KEGGcw,
                     UR2KEGGcx, UR2KEGGcy, UR2KEGGcz, UR2KEGGda, UR2KEGGdb,
                     UR2KEGGdc, UR2KEGGdd, UR2KEGGde, UR2KEGGdf, UR2KEGGdg,
                     UR2KEGGdh, UR2KEGGdi, UR2KEGGdj, UR2KEGGdk, UR2KEGGdl,
                     UR2KEGGdm, UR2KEGGdn, UR2KEGGdo, UR2KEGGdp, UR2KEGGdq,
                     UR2KEGGdr, UR2KEGGds, UR2KEGGdt, UR2KEGGdu, UR2KEGGdv,
                     UR2KEGGdw, UR2KEGGdx, UR2KEGGdy, UR2KEGGdz, UR2KEGGea,
                     UR2KEGGeb, UR2KEGGec, UR2KEGGed, UR2KEGGee, UR2KEGGef,
                     UR2KEGGeg, UR2KEGGeh, UR2KEGGei, UR2KEGGej, UR2KEGGek,
                     UR2KEGGel, UR2KEGGem, UR2KEGGen, UR2KEGGeo, UR2KEGGep,
                     UR2KEGGeq, UR2KEGGer, UR2KEGGes, UR2KEGGet, UR2KEGGeu,
                     UR2KEGGev, UR2KEGGew, UR2KEGGex, UR2KEGGey, UR2KEGGez, 
                     UR2KEGGfa, UR2KEGGfb, UR2KEGGfc, UR2KEGGfd, UR2KEGGfe,
                     UR2KEGGff, UR2KEGGfg, UR2KEGGfh, UR2KEGGfi, UR2KEGGfj,
                     UR2KEGGfk, UR2KEGGfl, UR2KEGGfm, UR2KEGGfn, UR2KEGGfo,
                     UR2KEGGfp, UR2KEGGfq, UR2KEGGfr, UR2KEGGfs, UR2KEGGft,
                     UR2KEGGfu, UR2KEGGfv, UR2KEGGfw, UR2KEGGfx, UR2KEGGfy,
                     UR2KEGGfz, UR2KEGGga, UR2KEGGgb, UR2KEGGgc, UR2KEGGgd, 
                     UR2KEGGge, UR2KEGGgf, UR2KEGGgg, UR2KEGGgh, UR2KEGGgi,
                     UR2KEGGgj, UR2KEGGgk, UR2KEGGgl, UR2KEGGgm, UR2KEGGgn,
                     UR2KEGGgo, UR2KEGGgp, UR2KEGGgq, UR2KEGGgr, UR2KEGGgs,
                     UR2KEGGgt, UR2KEGGgu, UR2KEGGgv, UR2KEGGgw, UR2KEGGgx,
                     UR2KEGGgy, UR2KEGGgz, UR2KEGGha, UR2KEGGhb, UR2KEGGhc,
                     UR2KEGGhd, UR2KEGGhe, UR2KEGGhf, UR2KEGGhg, UR2KEGGhh, 
                     UR2KEGGhi, UR2KEGGhj, UR2KEGGhk, UR2KEGGhl, UR2KEGGhm,
                     UR2KEGGhn, UR2KEGGho, UR2KEGGhp, UR2KEGGhq, UR2KEGGhr,
                     UR2KEGGhs, UR2KEGGht, UR2KEGGhu, UR2KEGGhv, UR2KEGGhw,
                     UR2KEGGhx, UR2KEGGhy, UR2KEGGhz, UR2KEGGia, UR2KEGGib,
                     UR2KEGGic, UR2KEGGid, UR2KEGGie, UR2KEGGif, UR2KEGGig,
                     UR2KEGGih, UR2KEGGii, UR2KEGGij, UR2KEGGik, UR2KEGGil,
                     UR2KEGGim, UR2KEGGin, UR2KEGGio, UR2KEGGip, UR2KEGGiq,
                     UR2KEGGir, UR2KEGGis, UR2KEGGit, UR2KEGGiu, UR2KEGGiv,
                     UR2KEGGiw, UR2KEGGix, UR2KEGGiy, UR2KEGGiz, UR2KEGGja,
                     UR2KEGGjb, UR2KEGGjc, UR2KEGGjd, UR2KEGGje, UR2KEGGjf,
                     UR2KEGGjg, UR2KEGGjh, UR2KEGGji, UR2KEGGjj, UR2KEGGjk,
                     UR2KEGGjl, UR2KEGGjm, UR2KEGGjn, UR2KEGGjo, UR2KEGGjp, 
                     UR2KEGGjq, UR2KEGGjr, UR2KEGGjs, UR2KEGGjt, UR2KEGGju,
                     UR2KEGGjv, #UR2KEGGjw,
                     UR2KEGGjx, UR2KEGGjy, UR2KEGGjz,
                     UR2KEGGka, UR2KEGGkb, UR2KEGGkc, UR2KEGGkd, UR2KEGGke,
                     UR2KEGGkf, UR2KEGGkg, UR2KEGGkh, UR2KEGGki, UR2KEGGkj,
                     UR2KEGGkk, UR2KEGGkl, UR2KEGGkm, UR2KEGGkn, UR2KEGGko,
                     UR2KEGGkp, UR2KEGGkq, UR2KEGGkr
)

print("UR2KEGGa_k combined.")

remove(UR2KEGGaa, UR2KEGGab, UR2KEGGac, UR2KEGGad, UR2KEGGae, UR2KEGGaf,
       UR2KEGGag, UR2KEGGah, UR2KEGGai, UR2KEGGaj, UR2KEGGak, UR2KEGGal,
       UR2KEGGam, UR2KEGGan, UR2KEGGao, UR2KEGGap, UR2KEGGaq, UR2KEGGar,
       UR2KEGGas, UR2KEGGat, UR2KEGGau, UR2KEGGav, UR2KEGGaw, UR2KEGGax,
       UR2KEGGay, UR2KEGGaz, UR2KEGGba, UR2KEGGbb, UR2KEGGbc, UR2KEGGbd,
       UR2KEGGbe, UR2KEGGbf, UR2KEGGbg, UR2KEGGbh, UR2KEGGbi, UR2KEGGbj,
       UR2KEGGbk, UR2KEGGbl, UR2KEGGbm, UR2KEGGbn, UR2KEGGbo, UR2KEGGbp,
       UR2KEGGbq, UR2KEGGbr, UR2KEGGbs, UR2KEGGbt, UR2KEGGbu, UR2KEGGbv,
       UR2KEGGbw, UR2KEGGbx, UR2KEGGby, UR2KEGGbz, UR2KEGGca, UR2KEGGcb,
       UR2KEGGcc, UR2KEGGcd, UR2KEGGce, UR2KEGGcf, UR2KEGGcg, UR2KEGGch, 
       UR2KEGGci, UR2KEGGcj, UR2KEGGck, UR2KEGGcl, UR2KEGGcm, UR2KEGGcn,
       UR2KEGGco, UR2KEGGcp, UR2KEGGcq, UR2KEGGcr, UR2KEGGcs, UR2KEGGct,
       UR2KEGGcu, UR2KEGGcv, UR2KEGGcw, UR2KEGGcx, UR2KEGGcy, UR2KEGGcz,
       UR2KEGGda, UR2KEGGdb, UR2KEGGdc, UR2KEGGdd, UR2KEGGde, UR2KEGGdf,
       UR2KEGGdg, UR2KEGGdh, UR2KEGGdi, UR2KEGGdj, UR2KEGGdk, UR2KEGGdl,
       UR2KEGGdm, UR2KEGGdn, UR2KEGGdo, UR2KEGGdp, UR2KEGGdq, UR2KEGGdr,
       UR2KEGGds, UR2KEGGdt, UR2KEGGdu, UR2KEGGdv, UR2KEGGdw, UR2KEGGdx,
       UR2KEGGdy, UR2KEGGdz, UR2KEGGea, UR2KEGGeb, UR2KEGGec, UR2KEGGed,
       UR2KEGGee, UR2KEGGef, UR2KEGGeg, UR2KEGGeh, UR2KEGGei, UR2KEGGej,
       UR2KEGGek, UR2KEGGel, UR2KEGGem, UR2KEGGen, UR2KEGGeo, UR2KEGGep, 
       UR2KEGGeq, UR2KEGGer, UR2KEGGes, UR2KEGGet, UR2KEGGeu, UR2KEGGev,
       UR2KEGGew, UR2KEGGex, UR2KEGGey, UR2KEGGez, UR2KEGGfa, UR2KEGGfb,
       UR2KEGGfc, UR2KEGGfd, UR2KEGGfe, UR2KEGGff, UR2KEGGfg, UR2KEGGfh,
       UR2KEGGfi, UR2KEGGfj, UR2KEGGfk, UR2KEGGfl, UR2KEGGfm, UR2KEGGfn,
       UR2KEGGfo, UR2KEGGfp, UR2KEGGfq, UR2KEGGfr, UR2KEGGfs, UR2KEGGft,
       UR2KEGGfu, UR2KEGGfv, UR2KEGGfw, UR2KEGGfx, UR2KEGGfy, UR2KEGGfz,
       UR2KEGGga, UR2KEGGgb, UR2KEGGgc, UR2KEGGgd, UR2KEGGge, UR2KEGGgf,
       UR2KEGGgg, UR2KEGGgh, UR2KEGGgi, UR2KEGGgj, UR2KEGGgk, UR2KEGGgl,
       UR2KEGGgm, UR2KEGGgn, UR2KEGGgo, UR2KEGGgp, UR2KEGGgq, UR2KEGGgr,
       UR2KEGGgs, UR2KEGGgt, UR2KEGGgu, UR2KEGGgv, UR2KEGGgw, UR2KEGGgx,
       UR2KEGGgy, UR2KEGGgz, UR2KEGGha, UR2KEGGhb, UR2KEGGhc, UR2KEGGhd,
       UR2KEGGhe, UR2KEGGhf, UR2KEGGhg, UR2KEGGhh, UR2KEGGhi, UR2KEGGhj,
       UR2KEGGhk, UR2KEGGhl, UR2KEGGhm, UR2KEGGhn, UR2KEGGho, UR2KEGGhp,
       UR2KEGGhq, UR2KEGGhr, UR2KEGGhs, UR2KEGGht, UR2KEGGhu, UR2KEGGhv,
       UR2KEGGhw, UR2KEGGhx, UR2KEGGhy, UR2KEGGhz, UR2KEGGia, UR2KEGGib,
       UR2KEGGic, UR2KEGGid, UR2KEGGie, UR2KEGGif, UR2KEGGig, UR2KEGGih,
       UR2KEGGii, UR2KEGGij, UR2KEGGik, UR2KEGGil, UR2KEGGim, UR2KEGGin,
       UR2KEGGio, UR2KEGGip, UR2KEGGiq, UR2KEGGir, UR2KEGGis, UR2KEGGit,
       UR2KEGGiu, UR2KEGGiv, UR2KEGGiw, UR2KEGGix, UR2KEGGiy, UR2KEGGiz,
       UR2KEGGja, UR2KEGGjb, UR2KEGGjc, UR2KEGGjd, UR2KEGGje, UR2KEGGjf,
       UR2KEGGjg, UR2KEGGjh, UR2KEGGji, UR2KEGGjj, UR2KEGGjk, UR2KEGGjl,
       UR2KEGGjm, UR2KEGGjn, UR2KEGGjo, UR2KEGGjp, UR2KEGGjq, UR2KEGGjr,
       UR2KEGGjs, UR2KEGGjt, UR2KEGGju, UR2KEGGjv, #UR2KEGGjw,
       UR2KEGGjx,
       UR2KEGGjy, UR2KEGGjz, UR2KEGGka, UR2KEGGkb, UR2KEGGkc, UR2KEGGkd,
       UR2KEGGke, UR2KEGGkf, UR2KEGGkg, UR2KEGGkh, UR2KEGGki, UR2KEGGkj,
       UR2KEGGkk, UR2KEGGkl, UR2KEGGkm, UR2KEGGkn, UR2KEGGko, UR2KEGGkp,
       UR2KEGGkq, UR2KEGGkr
       )

# Read uhgp_pep fractions
Pep2Protaa <- read.delim(file.path(work_path, "uhgp_pepaa"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)

print("uhgp_pepaa.tsv ingested (Pep2Protaa).")

Pep2Protab <- read.delim(file.path(work_path, "uhgp_pepab"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)

print("uhgp_pepab.tsv ingested (Pep2Protab).")

Pep2Protac <- read.delim(file.path(work_path, "uhgp_pepac"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)

print("uhgp_pepac.tsv ingested (Pep2Protac).")

Pep2Protad <- read.delim(file.path(work_path, "uhgp_pepad"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)

print("uhgp_pepad.tsv ingested (Pep2Protad).")

Pep2Protae <- read.delim(file.path(work_path, "uhgp_pepae"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)

print("uhgp_pepae.tsv ingested (Pep2Protae).")

Pep2Protaf <- read.delim(file.path(work_path, "uhgp_pepaf"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)

print("uhgp_pepaf.tsv ingested (Pep2Protaf).")

Pep2Protag <- read.delim(file.path(work_path, "uhgp_pepag"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)

print("uhgp_pepag.tsv ingested (Pep2Protag).")

Pep2Protah <- read.delim(file.path(work_path, "uhgp_pepah"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)

print("uhgp_pepah.tsv ingested (Pep2Protah).")

Pep2Protai <- read.delim(file.path(work_path, "uhgp_pepai"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)

print("uhgp_pepai.tsv ingested (Pep2Protai).")

Pep2Protaj <- read.delim(file.path(work_path, "uhgp_pepaj"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)

print("uhgp_pepaj.tsv ingested (Pep2Protaj).")

Pep2Protak <- read.delim(file.path(work_path, "uhgp_pepak"), header = FALSE) %>%
  rename(Peptide = V1, Protein = V2)

print("uhgp_pepak.tsv ingested (Pep2Protak).")

# Read diamond_uhgp_out
Prot2UR <- read.delim(file.path(work_path, "diamond_uhgp_best_hits_out.tsv"),
                      header = FALSE) %>%
  select(1, 3) %>%
  rename(Protein = V1, ID = V3)

print("diamond_uhgp_best_hits_out.tsv ingested (Prot2UR).")

# Join by values of the "Protein"
print("Join Pep2Protaa and Prot2UR by Protein...")
Pep2Prot2URaa <- left_join(Pep2Protaa, Prot2UR, by = "Protein",
                         relationship = "many-to-many")
print("Done.")

print("Join Pep2Protab and Prot2UR by Protein...")
Pep2Prot2URab <- left_join(Pep2Protab, Prot2UR, by = "Protein",
                         relationship = "many-to-many")
print("Done.")

print("Join Pep2Protac and Prot2UR by Protein...")
Pep2Prot2URac <- left_join(Pep2Protac, Prot2UR, by = "Protein",
                          relationship = "many-to-many")
print("Done.")

print("Join Pep2Protad and Prot2UR by Protein...")
Pep2Prot2URad <- left_join(Pep2Protad, Prot2UR, by = "Protein",
                          relationship = "many-to-many")
print("Done.")

print("Join Pep2Protae and Prot2UR by Protein...")
Pep2Prot2URae <- left_join(Pep2Protae, Prot2UR, by = "Protein",
                          relationship = "many-to-many")
print("Done.")

print("Join Pep2Protaf and Prot2UR by Protein...")
Pep2Prot2URaf <- left_join(Pep2Protaf, Prot2UR, by = "Protein",
                          relationship = "many-to-many")
print("Done.")

print("Join Pep2Protag and Prot2UR by Protein...")
Pep2Prot2URag <- left_join(Pep2Protag, Prot2UR, by = "Protein",
                          relationship = "many-to-many")
print("Done.")

print("Join Pep2Protah and Prot2UR by Protein...")
Pep2Prot2URah <- left_join(Pep2Protah, Prot2UR, by = "Protein",
                          relationship = "many-to-many")
print("Done.")

print("Join Pep2Protai and Prot2UR by Protein...")
Pep2Prot2URai <- left_join(Pep2Protai, Prot2UR, by = "Protein",
                          relationship = "many-to-many")
print("Done.")

print("Join Pep2Protaj and Prot2UR by Protein...")
Pep2Prot2URaj <- left_join(Pep2Protaj, Prot2UR, by = "Protein",
                          relationship = "many-to-many")
print("Done.")

print("Join Pep2Protak and Prot2UR by Protein...")
Pep2Prot2URak <- left_join(Pep2Protak, Prot2UR, by = "Protein",
                          relationship = "many-to-many")
print("Done.")

# Remove "UniRef90_" prefix
Pep2Prot2URaa$ID <- substr(Pep2Prot2URaa$ID, 10, length(Pep2Prot2URaa$ID))
Pep2Prot2URab$ID <- substr(Pep2Prot2URab$ID, 10, length(Pep2Prot2URab$ID))
Pep2Prot2URac$ID <- substr(Pep2Prot2URac$ID, 10, length(Pep2Prot2URac$ID))
Pep2Prot2URad$ID <- substr(Pep2Prot2URad$ID, 10, length(Pep2Prot2URad$ID))
Pep2Prot2URae$ID <- substr(Pep2Prot2URae$ID, 10, length(Pep2Prot2URae$ID))
Pep2Prot2URaf$ID <- substr(Pep2Prot2URaf$ID, 10, length(Pep2Prot2URaf$ID))
Pep2Prot2URag$ID <- substr(Pep2Prot2URag$ID, 10, length(Pep2Prot2URag$ID))
Pep2Prot2URah$ID <- substr(Pep2Prot2URah$ID, 10, length(Pep2Prot2URah$ID))
Pep2Prot2URai$ID <- substr(Pep2Prot2URai$ID, 10, length(Pep2Prot2URai$ID))
Pep2Prot2URaj$ID <- substr(Pep2Prot2URaj$ID, 10, length(Pep2Prot2URaj$ID))
Pep2Prot2URak$ID <- substr(Pep2Prot2URak$ID, 10, length(Pep2Prot2URak$ID))

# Join by "UniRef90"
print("Join Pep2Prot2URaa and UR2KEGG by UniRef90_ID...")
Pep2Prot2UR2KEGGaa <- left_join(Pep2Prot2URaa, UR2KEGG, by = "ID",
                              relationship = "many-to-many")
print("Done.")

print("Join Pep2Prot2URab and UR2KEGG by UniRef90_ID...")
Pep2Prot2UR2KEGGab <- left_join(Pep2Prot2URab, UR2KEGG, by = "ID",
                              relationship = "many-to-many")
print("Done.")
print("Join Pep2Prot2URac and UR2KEGG by UniRef90_ID...")
Pep2Prot2UR2KEGGac <- left_join(Pep2Prot2URac, UR2KEGG, by = "ID",
                              relationship = "many-to-many")
print("Done.")

print("Join Pep2Prot2URad and UR2KEGG by UniRef90_ID...")
Pep2Prot2UR2KEGGad <- left_join(Pep2Prot2URad, UR2KEGG, by = "ID",
                              relationship = "many-to-many")
print("Done.")

print("Join Pep2Prot2URae and UR2KEGG by UniRef90_ID...")
Pep2Prot2UR2KEGGae <- left_join(Pep2Prot2URae, UR2KEGG, by = "ID",
                              relationship = "many-to-many")
print("Done.")

print("Join Pep2Prot2URaf and UR2KEGG by UniRef90_ID...")
Pep2Prot2UR2KEGGaf <- left_join(Pep2Prot2URaf, UR2KEGG, by = "ID",
                              relationship = "many-to-many")
print("Done.")

print("Join Pep2Prot2URag and UR2KEGG by UniRef90_ID...")
Pep2Prot2UR2KEGGag <- left_join(Pep2Prot2URag, UR2KEGG, by = "ID",
                              relationship = "many-to-many")
print("Done.")

print("Join Pep2Prot2URah and UR2KEGG by UniRef90_ID...")
Pep2Prot2UR2KEGGah <- left_join(Pep2Prot2URah, UR2KEGG, by = "ID",
                              relationship = "many-to-many")
print("Done.")

print("Join Pep2Prot2URai and UR2KEGG by UniRef90_ID...")
Pep2Prot2UR2KEGGai <- left_join(Pep2Prot2URai, UR2KEGG, by = "ID",
                              relationship = "many-to-many")
print("Done.")

print("Join Pep2Prot2URaj and UR2KEGG by UniRef90_ID...")
Pep2Prot2UR2KEGGaj <- left_join(Pep2Prot2URaj, UR2KEGG, by = "ID",
                              relationship = "many-to-many")
print("Done.")

print("Join Pep2Prot2URak and UR2KEGG by UniRef90_ID...")
Pep2Prot2UR2KEGGak <- left_join(Pep2Prot2URak, UR2KEGG, by = "ID",
                              relationship = "many-to-many")
print("Done.")

# Filter peptides without KEGG annotation
Pep2Prot2UR2KEGGaa <- Pep2Prot2UR2KEGGaa %>% filter(KEGG != "NULL")
Pep2Prot2UR2KEGGab <- Pep2Prot2UR2KEGGab %>% filter(KEGG != "NULL")
Pep2Prot2UR2KEGGac <- Pep2Prot2UR2KEGGac %>% filter(KEGG != "NULL")
Pep2Prot2UR2KEGGad <- Pep2Prot2UR2KEGGad %>% filter(KEGG != "NULL")
Pep2Prot2UR2KEGGae <- Pep2Prot2UR2KEGGae %>% filter(KEGG != "NULL")
Pep2Prot2UR2KEGGaf <- Pep2Prot2UR2KEGGaf %>% filter(KEGG != "NULL")
Pep2Prot2UR2KEGGag <- Pep2Prot2UR2KEGGag %>% filter(KEGG != "NULL")
Pep2Prot2UR2KEGGah <- Pep2Prot2UR2KEGGah %>% filter(KEGG != "NULL")
Pep2Prot2UR2KEGGai <- Pep2Prot2UR2KEGGai %>% filter(KEGG != "NULL")
Pep2Prot2UR2KEGGaj <- Pep2Prot2UR2KEGGaj %>% filter(KEGG != "NULL")
Pep2Prot2UR2KEGGak <- Pep2Prot2UR2KEGGak %>% filter(KEGG != "NULL")

print("Filtered peptides without KEGG annotation.")

# Bind uhgp_pep fractions
Pep2Prot2UR2KEGG <- bind_rows(Pep2Prot2UR2KEGGaa, Pep2Prot2UR2KEGGab,
                              Pep2Prot2UR2KEGGac, Pep2Prot2UR2KEGGad, 
                              Pep2Prot2UR2KEGGae, Pep2Prot2UR2KEGGaf, 
                              Pep2Prot2UR2KEGGag, Pep2Prot2UR2KEGGah, 
                              Pep2Prot2UR2KEGGai, Pep2Prot2UR2KEGGaj, 
                              Pep2Prot2UR2KEGGak)

Pep2Prot2UR2KEGG <- select(Pep2Prot2UR2KEGG, -ID) #get rid of ID column

# Expand peptides with >1 KEGG terms, which are in a list, to give each KEGG
# term its own row
Pep2Prot2UR2KEGG <- unnest(Pep2Prot2UR2KEGG, KEGG)

# Save peptide-protein-UniRef90-KEGG dataset as csv file
write.csv(Pep2Prot2UR2KEGG, file = paste0(work_path, "/core_pep_kegg_db.csv"))
print("core_pep_kegg_db.csv saved.")
