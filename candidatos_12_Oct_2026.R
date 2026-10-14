install.packages("phytools")
install.packages("Preorder")
BiocManager::install("TreeTools")
#BiocManager::install("heattree")
#BiocManager::install("ggtree")
BiocManager::install("metacoder")
BiocManager::install("phylotate")

library(phytools)
library(Biostrings)
library(heattree)
library(TreeTools)
library(ggtree)
#library(metacoder)
#library(phylotate)
library(tidyverse)

iqtree_conjunto_noactinomicetos_Ycao_ML <- read.tree("03_Results/conjunto_noactinomicetos_YcaO/conjunto_MSA_ajustado_clipkit.fasta.treefile")

metadata_Query <- readRDS("03_Results/conjunto_noactinomicetos_YcaO/metadata_homologos.rds")
View(metadata_Query)



#arbol_ordenado <- Preorder(iqtree_conjunto_noactinomicetos_Ycao_ML)
#arbol_ordenado$tip.label
#arbol_ordenado2 <- Cladewise(iqtree_conjunto_noactinomicetos_Ycao_ML)
#class(arbol_ordenado2)

#Ver los clados
ggtree(iqtree_conjunto_noactinomicetos_Ycao_ML, layout = "rectangular", size= 0.2) + 
  geom_text2(aes(subset = !isTip, label = node), size = 0.5, color = "red", hjust = -0.3) +
  geom_tiplab(size= 0.5) +
  theme_tree2()


#lista_clados <- subtrees(iqtree_conjunto_noactinomicetos_Ycao_ML)
#iqtree_conjunto_noactinomicetos_Ycao_ML$Nnode
#extraer los clados una vez visto su numero del clado/nodo

###Clado A
cladoA <- extract.clade(iqtree_conjunto_noactinomicetos_Ycao_ML, node = 573)$tip.label
cladoA[1]
length(cladoA)#49

### Clado B
cladoB1 <- extract.clade(iqtree_conjunto_noactinomicetos_Ycao_ML, node = 782)$tip.label
cladoB2 <- extract.clade(iqtree_conjunto_noactinomicetos_Ycao_ML, node = 800)$tip.label
cladoB3 <- extract.clade(iqtree_conjunto_noactinomicetos_Ycao_ML, node = 814)$tip.label
cladoB4 <- extract.clade(iqtree_conjunto_noactinomicetos_Ycao_ML, node = 816)$tip.label
cladoB_total <- c(cladoB1,cladoB2, cladoB3, cladoB4)
cladoB_total
cladoB_total[1]
length(cladoB_total)#42

############## obtener los datos de los clados
#pasarlos a un solo csv

####CICLO FOR PARA EL CLADo A
a <- 1
clado_A_total_FOR <- list()
clado_A_total_FOR <- vector("list", length = length(cladoA))


for (i in 1:length(cladoA)){
  clado <- subset(metadata_Query, metadata_Query$accession == cladoA[a])
  clado <- as.data.frame(clado)
  print(i)
  clado_A_total_FOR[[a]] <- clado
  a <- a + 1
}
clado_df_A <- do.call(rbind, clado_A_total_FOR)
View(clado_df_A)
class(clado_df_A)

#escribe los accession numbers sin header en el csv
write.table(clado_df_A$assembly,sep=",", "01_RawData/candidatos_genomas/candidatos_genomas_cladoA.csv", row.names=FALSE, col.names= FALSE)
