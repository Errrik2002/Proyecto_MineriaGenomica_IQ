library(Biostrings)
library(dplyr)
library(tidyverse)

metadatos_nosiheptide_cluster <- read.csv("01_RawData/macrociclacion_BLAST_02_Oct_2026/metadata_macrociclacion/metadata_1_nosiheptide_ACR48344.1_CLUSTER.csv")
View(metadatos_nosiheptide_cluster)
metadatos_nosiheptide_cluster <- metadatos_nosiheptide_cluster[1:26,]



metadatos_nosiheptide_hit <- read.csv("01_RawData/macrociclacion_BLAST_02_Oct_2026/metadata_macrociclacion/metadata_2_nosiheptide_ACR48344.1_HIT.csv")
View(metadatos_nosiheptide_hit)

metadatos_nosiheptide <- bind_cols(metadatos_nosiheptide_cluster, metadatos_nosiheptide_hit)


