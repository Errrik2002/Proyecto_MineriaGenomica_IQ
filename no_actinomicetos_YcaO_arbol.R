# El siguiente analisis sera para la enzima de YcaO de thiocillina, lactocillina y microccocina.
 
library(Biostrings)
library(msa)
library(seqinr)
library(ggmsa)
library(ape)
library(ggtree)
library(DECIPHER)
library(microseq)
library(heattree)


conjunto_noactinomicetos_Ycao <- readAAStringSet("01_RawData/conjunto_thiopeptidos_no_actinomicetos/no_actinomicetos_YcaO.txt")
conjunto_noactinomicetos_Ycao
#son 756 secuencias
summary(width(conjunto_noactinomicetos_Ycao))

#Eliminar secuencias repetidas en caso de que tenga
conjunto_noactinomicetos_Ycao_unicos <- conjunto_noactinomicetos_Ycao[!duplicated(as.character(conjunto_noactinomicetos_Ycao))]
conjunto_noactinomicetos_Ycao_unicos
#### Si habia repetidas, las secuencias finales son 542 secuencias
summary(width(conjunto_noactinomicetos_Ycao_unicos))
saveRDS(conjunto_noactinomicetos_Ycao_unicos,"03_Results/conjunto_noactinomicetos_YcaO/secuencias_unicas_YcaO.rds")


#########Hacer el MSA
conjunto_noactinomicetos_MSA <- msa(conjunto_noactinomicetos_Ycao_unicos, method = "Muscle")
conjunto_noactinomicetos_MSA_compatible <- unmasked(conjunto_noactinomicetos_MSA)
BrowseSeqs(conjunto_noactinomicetos_MSA_compatible)
writeXStringSet(conjunto_noactinomicetos_MSA_compatible, "03_Results/conjunto_noactinomicetos_YcaO/conjunto_MSA_original.fasta")

#filtrado con ClipKit
conjunto_noactinomicetos_MSA_ajustado <- readAAMultipleAlignment("03_Results/conjunto_noactinomicetos_YcaO/conjunto_MSA_ajustado_clipkit.fasta")
conjunto_noactinomicetos_MSA_ajustado_compatible <- unmasked(conjunto_noactinomicetos_MSA_ajustado)
BrowseSeqs(conjunto_noactinomicetos_MSA_ajustado_compatible)

#######Funcion para conocer porcentajes de Gaps
porcentaje_de_gaps <- function(x){
  print("Debes de introducir un objeto AAMultipleAlignment. Lee tu fasta alignment con readAAMultipleAlignment()")
  print("-------------------------------------------------------------------------------------------------------")
 resultado <- readline(prompt = "Si quieres ver el %GAPS por renglon, pon 1. Si quieres ver por renglon, pon 2. Si quieres el global del MSA, pon 3. Si quieres ver todos. Pon 4:")
 resultado <- as.numeric(resultado)
  
  #Renglones/Rows
x_gaps_compatible <- unmasked(x)
gaps_row <- letterFrequency(x_gaps_compatible, letters = "-")
##
porcen_gaps_row <- (gaps_row / width(x_gaps_compatible)) * 100
##

#Columnas/Col
mc_x <- consensusMatrix(x)
gaps_col_x <- mc_x["-",]  
seq_totales <- nrow(x)
#
col_gaps_porcentaje <- (gaps_col_x / seq_totales) * 100
# 
mcompleta_x <- as.matrix(x)
gaps_global <- (sum(mcompleta_x=="-") / length(mcompleta_x)) * 100
  
#Condicional
if (resultado == 1){
print(porcen_gaps_row)
} else if (resultado == 2) { 
print(col_gaps_porcentaje)
} else if (resultado == 3) {
  cat("El alineamiento tiene un", gaps_global,"% de GAPS")
} else if (resultado == 4) {
  print(porcen_gaps_row)
  print(col_gaps_porcentaje)
  cat("El alineamiento tiene un", gaps_global,"% global de GAPS")
} else {
  print("Solo introduce un numero, del 1 al 4, solo uno por favor")
}
}

#MSA original
porcentaje_de_gaps(conjunto_noactinomicetos_MSA) #25.93% de gaps
#MSA ajustado
porcentaje_de_gaps(conjunto_noactinomicetos_MSA_ajustado) #20.8974%


#Arbol por distancias 
arbol_conjunto_noactino <- msaConvert(conjunto_noactinomicetos_MSA_ajustado, type = "seqinr::alignment")

#matriz de distancia
dm_conjunto_noactino_ycao <- dist.alignment(arbol_conjunto_noactino, matrix = "identity")

###Arbol NJ

tree_conjunto_noactinomicetos_Ycao_NJ <- nj(dm_conjunto_noactino_ycao)

#Ver arbol de NJ

plot(tree_conjunto_noactinomicetos_Ycao_NJ)

heat_tree(tree = tree_conjunto_noactinomicetos_Ycao_NJ)


##### Arbol con IQTREE
metadata_no_actinos <- readRDS("03_Results/conjunto_noactinomicetos_YcaO/metadata_base_FASTA.rds")


iqtree_conjunto_noactinomicetos_Ycao_ML <- read.tree("03_Results/conjunto_noactinomicetos_YcaO/conjunto_MSA_ajustado_clipkit.fasta.treefile")
heat_tree(iqtree_conjunto_noactinomicetos_Ycao_ML, metadata = metadata_no_actinos)


