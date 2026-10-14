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


#Ver los clados
ggtree(iqtree_conjunto_noactinomicetos_Ycao_ML, layout = "rectangular", size= 0.2) + 
  geom_text2(aes(subset = !isTip, label = node), size = 0.5, color = "red", hjust = -0.3) +
  geom_tiplab(size= 0.5) +
  theme_tree2()

##############################################################
#extraer los clados una vez visto su numero del clado/nodo####
##############################################################

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

##########################################################################
############## obtener los datos de los clados ###########################
##############pasarlos a un solo csv######################################
##########################################################################

#### CICLO FOR PARA EL CLADo A ####


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
saveRDS(clado_df_A, "03_Results/candidatos_antismash/clado_A.rds")

#escribe los accession numbers sin header en el csv
write.table(clado_df_A$assembly,sep=",", "01_RawData/candidatos_genomas/candidatos_genomas_cladoA.csv", row.names=FALSE, col.names= FALSE)


####################################################
################ Proceso clado B ###################
####################################################

b <- 1
clado_B_FOR <- list()
clado_B_FOR <- vector("list", length= length(cladoB_total))

for (i in 1:length(cladoB_total)){
  clado_dos <- subset(metadata_Query, metadata_Query$accession == cladoB_total[b])
  clado_dos <- as.data.frame(clado_dos)
  print(i)

  clado_B_FOR[[b]] <- clado_dos

  b <- b + 1
}

clado_B_df <- do.call(rbind, clado_B_FOR)
View(clado_B_df)
saveRDS(clado_B_df, "03_Results/candidatos_antismash/clado_B.rds")


#### escribir el csv con los genomas del clado que quiero que selccione ####
write.table(clado_B_df$assembly, sep = ",", "01_RawData/candidatos_genomas/cladoB/candidatos_cladoB.csv", 
row.names = FALSE, col.names = FALSE)




############################################
############################################
# Exploracion de los clados ################
############################################
############################################
colnames(clado_df_A)
clado_A_explorar <- select(clado_df_A, c("accession","especie","host", "isolation_source", "geo_loc_name"))
View(clado_A_explorar)


colnames(clado_B_df)
clado_B_explorar <- select(clado_B_df, c("accession","especie","host", "isolation_source", "geo_loc_name"))
View(clado_B_explorar)




#############################################
#############################################
### Especies finales ########################
#############################################

library(dplyr)
#####################
#### Clado A ########
#####################

candidatos_finales_cladoA <- c("WP_110044842.1", "WP_162357576.1", "WP_211020662.1",  "WP_068780233.1" ,"WP_233696727.1","WP_405452319.1", "WP_246061953.1", "WP_316561576.1", "WP_041047436.1", "WP_167056766.1", "WP_155701669.1",  "WP_230460689.1", "WP_144699461.1" )
candidatos_finales_cladoA <- as.list(candidatos_finales_cladoA)


c <- 1
clado_A_Final <- list()
clado_A_Final <- vector("list", length= length(candidatos_finales_cladoA))

for (i in 1:length(candidatos_finales_cladoA)){
  clado_final <- subset(clado_df_A, clado_df_A$accession == candidatos_finales_cladoA[c])
  clado_final <- as.data.frame(clado_final)
  print(i)

  clado_A_Final[[c]] <- clado_final

  c <- c + 1
}

clado_A_Final_df <- do.call(rbind, clado_A_Final)
View(clado_A_Final_df)


write.table(clado_A_Final_df$assembly, sep = ",", "01_RawData/candidatos_genomas/cladoA/candidatos_finales_cladoA.csv", 
row.names = FALSE, col.names = FALSE)

############################
####### Clado B ############
############################
candidatos_finales_cladoB <- c("MBP3040987.1", "WP_335424268.1", "WP_091587821.1", "WP_213104404.1", "WP_060206736.1", "WP_239432773.1", "WP_213106905.1", "WP_101569076.1", "WP_338753039.1", "WP_251525773.1", "MGM0852453.1", "WP_136351696.1")
candidatos_finales_cladoB <- as.list(candidatos_finales_cladoB)


d <- 1
clado_B_Final <- list()
clado_B_Final <- vector("list", length= length(candidatos_finales_cladoB))

for (i in 1:length(candidatos_finales_cladoB)){
  clado_finalB <- subset(clado_B_df, clado_B_df$accession == candidatos_finales_cladoB[d])
  clado_finalB <- as.data.frame(clado_finalB)
  print(i)

  clado_B_Final[[d]] <- clado_finalB

  d <- d + 1
}

clado_B_Final_df <- do.call(rbind, clado_B_Final)
View(clado_B_Final_df)

write.table(clado_B_Final_df$assembly, sep = ",", "01_RawData/candidatos_genomas/cladoB/candidatos_cladoB_finales.csv", 
row.names = FALSE, col.names = FALSE)
