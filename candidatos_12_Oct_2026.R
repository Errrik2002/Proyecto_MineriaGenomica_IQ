BiocManager::install("biomaRt")
install.packages("rentrez")

library(biomaRt)
library(Biostrings)
library(rentrez)

metadata_Query <- readRDS("03_Results/conjunto_noactinomicetos_YcaO/metadata_homologos.rds")
View(metadata_Query)

clado2a <- subset(metadata_Query, metadata_Query$accession == "MBP3040987.1")
clado2a$assembly
clado2b <- subset(metadata_Query, metadata_Query$accession == "WP_239432773.1")
clado2b$assembly
clado2c <- subset(metadata_Query, metadata_Query$accession ==  "WP_101569076.1")
clado2c$assembly

