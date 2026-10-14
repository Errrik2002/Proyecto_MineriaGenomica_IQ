library(tidyverse)
library(Biostrings)

metadata_comparar <- readRDS("03_Results/conjunto_noactinomicetos_YcaO/metadata_homologos.rds")
View(metadata_comparar)
metadata_bacillus <- subset(metadata_comparar,metadata_comparar$query_homologos == "Bacillus cereus")
View(metadata_bacillus)
print(metadata_bacillus$accession)

secuencias_bacillus <- readRDS("03_Results/conjunto_noactinomicetos_YcaO/secuencias_unicas_YcaO.rds")

secuencias_bacillus <- as.data.frame(secuencias_bacillus)
dim(secuencias_bacillus)

secuencias_bacillus1 <- subset(secuencias_bacillus, secuencias_bacillus$x)

