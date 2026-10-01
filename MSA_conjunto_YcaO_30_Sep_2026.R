#BiocManager::install("pacman")

#pacman::p_load(seqinr, msa, ggmsa, Biostrings, ape, ggtree, DECIPHER)


install.packages("seqinr") #manejar secuencias 
BiocManager::install("msa") #Hacer el alineamiento
BiocManager::install("ggmsa") #Visualizacion del alineamiento
BiocManager::install("Biostrings") #Manejar secuencias
install.packages("ape") #Realizar el arbol filogenetico
BiocManager::install("ggtree") #Visualizar el arbol filogenetico
BiocManager::install("DECIPHER") #Hacer alineamiento y curation
BiocManager::install("microseq")

?DECIPHER
library(Biostrings)
library(msa)
library(seqinr)
library(ggmsa)
library(ape)
library(ggtree)
library(DECIPHER)
library(microseq)

######################################
#### Leer las secuencias cargadas ####
######################################


secuencias_juntas_Ycao <- readAAStringSet("01_RawData/ycao_Thiopeptina_syomicina_GE37468_nosiheptide/YcaO_combinados_thiopeptina_syiomicina_GE_Nosiheptide.txt")
secuencias_juntas_Ycao
summary(width(secuencias_juntas_Ycao))


#################################################
#### Realizar los Multple Sequence Alignment ####
#################################################

?msa()

msa_ycao_juntas <- msa(secuencias_juntas_Ycao, "Muscle")
class(msa_ycao_juntas)
#Visalizar MSA
msa_ycao_juntas_compatible <- unmasked(msa_ycao_juntas)
class(msa_ycao_juntas_compatible)
BrowseSeqs(msa_ycao_juntas_compatible)

#Guardar MSA
writeXStringSet(msa_ycao_juntas_compatible,"03_Results/MSA_YcaO_juntas_original.fasta")

msa_ycaoJUNTAS_original <- readAAMultipleAlignment("03_Results/MSA_YcaO_juntas_original.fasta")

########
####TRIMMING 
#########
?msaTrim()


#######################################
#### Posterior al trimming clipkit ####
#######################################

msa_ycao_juntas_curadas <- readAAMultipleAlignment("03_Results/Msa_YcaO_juntas_arbol/MSA_Ycao_juntas_curadas_Clipkit.fasta")
msa_ycao_juntas_curadas_compatibles <- unmasked(msa_ycao_juntas_curadas)
BrowseSeqs(msa_ycao_juntas_curadas_compatibles)





###################################################################
##### Visualizar % de Gaps en cada MSA por secuencia FILA #########
###################################################################
#msa_ycao_juntas_curadas_compatibles <- unmasked(msa_ycao_juntas_curadas)
#msa_ycao_juntas_originales_compatibles <- unmasked(msa_ycaoJUNTAS_original)

#class(msa_ycao_juntas_curadas_compatibles)
#class(msa_ycao_juntas_originales_compatibles)

#gaps_original <- letterFrequency(msa_ycao_juntas_originales_compatibles, letters="-")
#porcen_gaps_original <- (gaps_original/width(msa_ycao_juntas_originales_compatibles))* 100
#print(porcen_gaps_original)


####################################################
#### pORCENTAJE DE GAPS POR COLUMNAS ###############
####################################################

#Hacer esto en una funcion

#matrix_consenso_Ycao_original <- consensusMatrix(msa_ycaoJUNTAS_original)
#gaps_columna <- matrix_consenso_Ycao_original["-",]
#gaps_columna

#seq_totales <- nrow(msa_ycaoJUNTAS_original)

#porcentajes_gaps_columns <- (gaps_columna / seq_totales) * 100

#Ver % de GAPS
#porcentajes_gaps_columns


###############################################################
#### pORCENTAJE DE GAPS GLOBAL del alineamiento ###############
###############################################################

#matrix_completa <- as.matrix(msa_ycaoJUNTAS_original)
#matrix_completa_curado <- as.matrix("03_Results/Msa_YcaO_juntas_arbol/MSA_Ycao_juntas_curadas_Clipkit.fasta")
#porcentaje global de gaps

#gaps_globales <- (sum(matrix_completa_curado=="-") / length(matrix_completa_curado)) * 100

#cat("El alineamiento total tiene un", gaps_globales,"% de gaps.")





################################################################
######################Funcion porcentaje de GAPS ###############
################################################################


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
  
#cat("El alineamiento tiene un", gaps_global,"% de GAPS")  

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


porcentaje_de_gaps(x1)

#gaps


##############################################################################
################################
#### Matrices de distancia #####
################################
class(msa_ycao_juntas_curadas)

tree_msa_ycao_juntas_curadas <- msaConvert(msa_ycao_juntas_curadas, type="seqinr::alignment")

class(tree_msa_ycao_juntas_curadas)

#####matriz de distancia
distance_matrix_msa_ycao_juntas <- dist.alignment(tree_msa_ycao_juntas_curadas, matrix="identity")

###########################
#### Arbol NJ #############
###########################

nj_tree_ycao_juntas <- nj(distance_matrix_msa_ycao_juntas)


ggtree(nj_tree_ycao_juntas, layout="circular", size=1)+
  geom_tiplab(size=2, aes(angle=angle))+
  geom_nodelab(geom='label')+
  hexpand(0.05)



#DECIPHER 
arbol_ycao_juntas_decipher <- Treeline(msa_ycao_juntas_curadas,
method="ML",
)
########################################
########### Arboles ML IQTREE ##########
########################################


ycao_juntas_arbol_ML <- read.tree("03_Results/Msa_YcaO_juntas_arbol/MSA_Ycao_juntas_curadas_Clipkit.fasta.treefile")

#Visualizar el arbol
plot(ycao_juntas_arbol_ML, main= "Arbol ML de YcaO juntas (Thiopeptina, syomicina, GE37468A, nosiheptide)")
