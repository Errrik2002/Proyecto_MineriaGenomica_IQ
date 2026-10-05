library(Biostrings)
library(msa)
library(seqinr)
library(ggmsa)
library(ape)
library(ggtree)
library(DECIPHER)
library(microseq)
library(heattree)


ycao_lactocillina_PI30_QC60 <- readAAStringSet("01_RawData/lactocillina_4_10_26/lactocillina_30PI_60QC_041026_YcaO.txt")
#
ycao_lactocillina_PI30_QC60 #son 183 secuencias
summary(width(ycao_lactocillina_PI30_QC60))


#Aqui no deberia de haber secuencias repetidas

lactocilina_msa <- msa(ycao_lactocillina_PI30_QC60, method="Muscle")
lactocilina_msa_compatible <- unmasked(lactocilina_msa)
BrowseSeqs(lactocilina_msa_compatible)
writeXStringSet(lactocilina_msa_compatible, "03_Results/lactocillina_MSA_YcaO/lactocillina_msa_original.fasta")

#Cargar con ajuste por ClipKit
lactocilina_msa_ajustado <- readAAMultipleAlignment("03_Results/lactocillina_MSA_YcaO/ajustado_lactocillina_MSA.fasta")
lactocilina_msa_ajustado_compatible <- unmasked(lactocilina_msa_ajustado)
BrowseSeqs(lactocilina_msa_ajustado_compatible)




######## Funcion
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


porcentaje_de_gaps(lactocilina_msa)
#la lactocillina msa_original tiene un 12.12% de GAPS
porcentaje_de_gaps(lactocilina_msa_ajustado)
#ya ajustado tiene un 4.75 de Gaps globales




#Arbol NJ 

lactocillina_arbol <- msaConvert(lactocilina_msa_ajustado, type="seqinr::alignment")

#matriz de distancias
lactocillina_matriz_arbol <- dist.alignment(lactocillina_arbol)


#####Arbol NJ

lactocillina_arbol_NJ <- nj(lactocillina_matriz_arbol)

heat_tree(tree=lactocillina_arbol_NJ)


###########Visualizar arborl ML IQTREE