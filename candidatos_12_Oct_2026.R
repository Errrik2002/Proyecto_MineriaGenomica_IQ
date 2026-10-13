install.packages("phytools")
install.packages("Preorder")
BiocManager::install("TreeTools")
BiocManager::install("heattree")
BiocManager::install("ggtree")
BiocManager::install("metacoder")
BiocManager::install("phylotate")

library(phytools)
library(Biostrings)
library(heattree)
library(TreeTools)
library(ggtree)
#library(metacoder)
#library(phylotate)

iqtree_conjunto_noactinomicetos_Ycao_ML <- read.tree("03_Results/conjunto_noactinomicetos_YcaO/conjunto_MSA_ajustado_clipkit.fasta.treefile")

metadata_Query <- readRDS("03_Results/conjunto_noactinomicetos_YcaO/metadata_homologos.rds")
View(metadata_Query)

#clado2a <- subset(metadata_Query, metadata_Query$accession == "MBP3040987.1")
#clado2b <- subset(metadata_Query, metadata_Query$accession == "WP_239432773.1")
#clado2c <- subset(metadata_Query, metadata_Query$accession ==  "WP_101569076.1")

plot(iqtree_conjunto_noactinomicetos_Ycao_ML)

arbol_ordenado <- Preorder(iqtree_conjunto_noactinomicetos_Ycao_ML)


heat_tree(iqtree_conjunto_noactinomicetos_Ycao_ML)


tips_arbol <- c("MBP3040987.1", "WP_239432773.1", "WP_101569076.1")
buscar_nodo <- getMRCA(iqtree_conjunto_noactinomicetos_Ycao_ML, tips_arbol)
?getMRCA()
#indice del nodo
buscar_nodo
