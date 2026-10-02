library(Biostrings)
library(msa)
library(seqinr)
library(ggmsa)
library(ape)
library(ggtree)
library(DECIPHER)
library(microseq)
library(heattree)

#Leer la secuencia cargada en la carpeta de 
#RawData y dentro de la carpeta de macrociclacion

macrociclacion_junta <- readAAStringSet("01_RawData/macrociclacion_BLAST_02_Oct_2026/seq_PSIBLAST_conjunto_macrociclacion.txt")
#Son un total de 58 secuencias
summary(width(macrociclacion_junta))

#Hay secuencias repetidas, hay que eliminarlas
macrociclacion_junta_unicas <- macrociclacion_junta[!duplicated(as.character(macrociclacion_junta))]
macrociclacion_junta_unicas
#quedan un total de 42 seq
summary(width(macrociclacion_junta_unicas))




##Realiar el alineamiento multiple MSA

msa_macrociclacion_conjunta <- msa(macrociclacion_junta_unicas, method="Muscle")
#msa_macrociclacion_conjunta
msa_macrociclacion_conjunta_compatible <- unmasked(msa_macrociclacion_conjunta)
BrowseSeqs(msa_macrociclacion_conjunta_compatible)

#Guardar MSA
writeXStringSet(msa_macrociclacion_conjunta_compatible,"03_Results/macrociclacion_resultados_021026/alineamiento_macrociclacion_conjuntas_compatible.fasta")

#HACER EL TRIMMING EN LA TERMINAL CON CLIPKIT

msa_clipkit_macrociclacion <- readAAMultipleAlignment("03_Results/macrociclacion_resultados_021026/clipkit_triming_MSA_macrociclacion.fasta")
msa_clipkit_macrociclacion_compatible <- unmasked(msa_clipkit_macrociclacion)
BrowseSeqs(msa_clipkit_macrociclacion_compatible)

#APLICAR FUNCION PARA CONOCER PORCENTAJES DE GAPS

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
porcentaje_de_gaps(msa_macrociclacion_conjunta) #secuencia original con un 27.0033% de gaps
porcentaje_de_gaps(msa_clipkit_macrociclacion) #Con trimming mejoro al 6.1335% de GAPS


#arbol filogenetico basado en distancias
arbol_msa_clipkit_macrociclacion <- msaConvert(msa_clipkit_macrociclacion, type= "seqinr::alignment")

#matrix de distancia
dm_arbol_msa_clipkit_macrocicalcion <- dist.alignment(arbol_msa_clipkit_macrociclacion, matrix="identity")

##############
## Arbol NJ ##
##############
#arbol con DECIPHER
msa_decipher <- readAAStringSet("03_Results/macrociclacion_resultados_021026/clipkit_triming_MSA_macrociclacion.fasta")
class(msa_decipher)
DECIPHER_arbol_macrociclacion <- TreeLine(msa_decipher,myDistMatrix=dm_arbol_msa_clipkit_macrocicalcion , method= "NJ")


nj_arbol_macrociclacion <- nj(dm_arbol_msa_clipkit_macrocicalcion)

#vISUALIZAR ARBOL
#1
plot(nj_arbol_macrociclacion)
#2
ggtree(nj_arbol_macrociclacion, layout="circular", size=1)+
  geom_tiplab(size=2, aes(angle=angle))+
  geom_nodelab(geom='label')+
  hexpand(0.05)
#3
heat_tree(tree=nj_arbol_macrociclacion)

#DECIPHER
plot(DECIPHER_arbol_macrociclacion)

ggtree(DECIPHER_arbol_macrociclacion, layout="circular", size=1)+
  geom_tiplab(size=2, aes(angle=angle))+
  geom_nodelab(geom='label')+
  hexpand(0.05)
#3

################
## Arbol ML ####
################

iqtree_arbol_macrociclacion <- read.tree("03_Results/macrociclacion_resultados_021026/clipkit_triming_MSA_macrociclacion.fasta.treefile")

plot(iqtree_arbol_macrociclacion, main="Arbol enzima macrocilacion")

ggtree_arbol_macrociclacion <- ggtree(iqtree_arbol_macrociclacion, layout="circular", size=1)+
  geom_tiplab(size=2, aes(angle=angle))+
  geom_nodelab(geom='label')+
  hexpand(0.05)

ggtree_arbol_macrociclacion


heat_tree(tree= iqtree_arbol_macrociclacion, layout="circular")
heat_tree(tree= iqtree_arbol_macrociclacion)
