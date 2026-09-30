#BiocManager::install("pacman")

#pacman::p_load(seqinr, msa, ggmsa, Biostrings, ape, ggtree, DECIPHER)


install.packages("seqinr") #manejar secuencias 
BiocManager::install("msa") #Hacer el alineamiento
BiocManager::install("ggmsa") #Visualizacion del alineamiento
BiocManager::install("Biostrings") #Manejar secuencias
install.packages("ape") #Realizar el arbol filogenetico
BiocManager::install("ggtree") #Visualizar el arbol filogenetico
BiocManager::install("DECIPHER") #Hacer alineamiento y curation

?DECIPHER
library(Biostrings)
library(msa)
library(seqinr)
library(ggmsa)
library(ape)
library(ggtree)
library(DECIPHER)

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


#######################################
#### Posterior al trimming clipkit ####
#######################################

msa_ycao_juntas_curadas <- readAAMultipleAlignment("03_Results/MSA_Ycao_juntas_curadas_Clipkit.fasta")
msa_ycao_juntas_curadas_compatibles <- unmasked(msa_ycao_juntas_curadas)
BrowseSeqs(msa_ycao_juntas_curadas_compatibles)


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