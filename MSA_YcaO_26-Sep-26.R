#############################################
#ARBOLES FILOGENETICOS DE LA ENZIMA YCAO#
#############################################


#Guardar objetos
# saveRDS()
# readRDS()
#### Carga de paquetes ya instalados ####

library(pacman)
library(Biostrings)
library(msa)
library(seqinr)
library(ggmsa)
library(ape)
library(ggtree)
library(DECIPHER)
library(phangorn)



##############################
#### Cargar las secuencias ####
##############################

thiopeptina_YcaO <- readAAStringSet("01_RawData/thiopeptina_BGC0001474_blast_240926/ycao_blast_sequences/thiopeptina_PYC80225.1_YcaO.txt")
thiopeptina_YcaO #No debo pasarlo a fasta necesariamente. #403SECUENCIAS 
summary(width(thiopeptina_YcaO))

syomicina_YcaO <- readAAStringSet("01_RawData/syiomicina_BGC0000611_BLAST_26-09-26/siomicina_GAP53284.1_YCAO.txt")
syomicina_YcaO  #POCAS SECUENCIAS 19 #CLUSTALW
summary(width(syomicina_YcaO))

GE37468A_YcaO <- readAAStringSet("01_RawData/GE37468A_BGC0000605_BLAST_240926/GE37468A_AEM00626.1_YCAO.txt")
GE37468A_YcaO #273SECUENCIAS
summary(width(GE37468A_YcaO))

nosiheptide_YcaO <- readAAStringSet("01_RawData/nosiheptide_BGC0000610_BLAST_24-09-26/nosiheptide_ACR48336.1_YCAO.txt")
nosiheptide_YcaO #POCAS SECUENCIAS 37 #CLUSTALW
summary(width(nosiheptide_YcaO))


##################################################
#### Realizar los Multiple Sequence Alignment ####
##################################################


?msa()

###Thiopeptina #MUSCLE#
msa_thiopeptina_Ycao <- msa(thiopeptina_YcaO, method="Muscle")
#msa_thiopeptina_Ycao_Decipher <- AlignSeqs(thiopeptina_YcaO)
#DECIPHER no puedo modificar el algoritmo

#Visualizar MSA en Browser
#Pasar el msa original a XStringSet
msa_thiopeptina_Ycao_compatible <- unmasked(msa_thiopeptina_Ycao)
BrowseSeqs(msa_thiopeptina_Ycao_compatible, highlight = 0)

#### Guardar el MSA ####
writeXStringSet(msa_thiopeptina_Ycao_compatible, "03_Results/msa_thiopeptina_ycao.fasta")



###Syomicina #CLUSTALW#
#msa
msa_syiomicina_Ycao <- msa(syomicina_YcaO)
#msa_syiomicina_Ycao_decipher <- AlignSeqs(syomicina_YcaO)

#Visualizar en Browser, hacerlo compatible con BrowSeqs
msa_syiomicina_Ycao_compatible <- unmasked(msa_syiomicina_Ycao)
BrowseSeqs(msa_syiomicina_Ycao_compatible)

#Guardar MSA
writeXStringSet(msa_syiomicina_Ycao_compatible, "03_Results/msa_syomicina_ycao.fasta")




###GE37468A #MUSCLE#
msa_GE37468A_Ycao <- msa(GE37468A_YcaO,method = "Muscle")
#msa_GE37468A_Ycao_DECIPHER <- AlignSeqs(GE37468A_YcaO)

#Visualizar en Browser
msa_GE37468A_Ycao_compatible <- unmasked(msa_GE37468A_Ycao)
BrowseSeqs(msa_GE37468A_Ycao_compatible)

#Guardar MSA 
writeXStringSet(msa_GE37468A_Ycao_compatible, "03_Results/msa_GE37468A_Ycao.fasta")




###Nosiheptide #ClustalW#
msa_nosiheptide_Ycao <- msa(nosiheptide_YcaO)
msa_nosiheptide_Ycao_decipher <- AlignSeqs(nosiheptide_YcaO)

#Visualizar en Browser
msa_nosiheptide_Ycao_compatible <- unmasked(msa_nosiheptide_Ycao)
BrowseSeqs(msa_nosiheptide_Ycao_compatible)

#Guardar MSA
writeXStringSet(msa_nosiheptide_Ycao_compatible, "03_Results/msa_nosiheptide_ycao.fasta")


?TreeLine()
?ape
library(help = ape)
?nj()
?upgma()



#########################################

#### POSTERIOR AL AJUSTE POR CLIPKIT ####

#########################################


#Thiopeptina
msa_thiopeptina_ycao_AJclipkit <- readAAMultipleAlignment("03_Results/msa_thiopeptina_ycao_ajustado_clipkit.fasta")
msa_thiopeptina_ycao_AJclipkit <- unmasked(msa_thiopeptina_ycao_AJclipkit)

msa_thiopeptina_original <- readAAMultipleAlignment("03_Results/msa_thiopeptina.fasta")
msa_thiopeptina_original <- unmasked(msa_thiopeptina_original)

BrowseSeqs(msa_thiopeptina_original)
BrowseSeqs(msa_thiopeptina_ycao_AJclipkit)



#Syomicina#
msa_syiomicina_Ycao_Ajustado_clipkit <- readAAMultipleAlignment("03_Results/msa_syomicina_ycao_ajustado_Clipkit.fasta")
msa_syiomicina_Ycao_Ajustado_clipkit <- unmasked(msa_syiomicina_Ycao_Ajustado_clipkit)

msa_syomicina_original <- readAAMultipleAlignment("03_Results/msa_syomicina_ycao.fasta")
msa_syomicina_original <- unmasked(msa_syomicina_original)

BrowseSeqs(msa_syomicina_original)
BrowseSeqs(msa_syiomicina_Ycao_Ajustado_clipkit)


#GE37468A#
msa_GE37468A_Ycao_ajustado_Clipkit <- readAAMultipleAlignment("03_Results/msa_GE37468A_Ycao_ajustado_clipkit.fasta")
msa_GE37468A_Ycao_ajustado_Clipkit <- unmasked(msa_GE37468A_Ycao_ajustado_Clipkit)

msa_GE37468A_Ycao_original <- readAAMultipleAlignment("03_Results/msa_GE37468A_Ycao.fasta")
msa_GE37468A_Ycao_original <- unmasked(msa_GE37468A_Ycao_original)

BrowseSeqs(msa_GE37468A_Ycao_original)
BrowseSeqs(msa_GE37468A_Ycao_ajustado_Clipkit)



#nosiheptide
msa_nosiheptide_Ycao_ajustado_clipkit <- readAAMultipleAlignment("03_Results/nosiheptide_ycao_tree/msa_nosiheptide_ycao_ajustado_Clipkit.fasta")
msa_nosiheptide_Ycao_ajustado_clipkit <- unmasked(msa_nosiheptide_Ycao_ajustado_clipkit)

msa_nosiheptide_Ycao_original <- readAAMultipleAlignment("03_Results/msa_nosiheptide_ycao.fasta")
msa_nosiheptide_Ycao_original <- unmasked(msa_nosiheptide_Ycao_original)

BrowseSeqs(msa_nosiheptide_Ycao_original)
BrowseSeqs(msa_nosiheptide_Ycao_ajustado_clipkit)




###############################
#### Matrices de distancia ####
###############################


#Modificar clase y poder hacer la matriz de distancia

tree_msa_thiopeptina_ycao_AJclipkit <- readAAMultipleAlignment("03_Results/msa_thiopeptina_ycao_ajustado_clipkit.fasta")
tree_msa_syiomicina_Ycao_Ajustado_clipkit <- readAAMultipleAlignment("03_Results/msa_syomicina_ycao_ajustado_Clipkit.fasta")
tree_msa_GE37468A_Ycao_ajustado_Clipkit <- readAAMultipleAlignment("03_Results/msa_GE37468A_Ycao_ajustado_clipkit.fasta")
tree_msa_nosiheptide_Ycao_ajustado_clipkit <- readAAMultipleAlignment("03_Results/nosiheptide_ycao_tree/msa_nosiheptide_ycao_ajustado_Clipkit.fasta")

msa_tree_thiopeptina_ycao <- msaConvert(tree_msa_thiopeptina_ycao_AJclipkit, type="seqinr::alignment")
msa_tree_syomicina_ycao <- msaConvert(tree_msa_syiomicina_Ycao_Ajustado_clipkit, type="seqinr::alignment")
msa_tree_GE37468A_ycao <- msaConvert(tree_msa_GE37468A_Ycao_ajustado_Clipkit, type= "seqinr::alignment")
msa_tree_nosiheptide_ycao <- msaConvert(tree_msa_nosiheptide_Ycao_ajustado_clipkit, type= "seqinr::alignment")

#?dist.alignment()
dm_thiopeptina_ycao <- dist.alignment(msa_tree_thiopeptina_ycao, matrix="identity")
dm_syomicina_ycao <- dist.alignment(msa_tree_syomicina_ycao, "identity")
dm_GE37468A_ycao <- dist.alignment(msa_tree_GE37468A_ycao, "identity")
dm_nosiheptide_ycao <- dist.alignment(msa_tree_nosiheptide_ycao, "identity")


###############################
#### Arboles filogeneticos ####
###############################

####Neighbor joining
?saveRDS
class(nj_tree_thiopeptina_ycao)


nj_tree_thiopeptina_ycao <- nj(dm_thiopeptina_ycao)
saveRDS(nj_tree_thiopeptina_ycao, "03_Results/arbol_nj_thiopeptina_YCAO.rds")
#plot(nj_tree_thiopeptina_ycao)

nj_tree_thiopeptina_ycao_plot <- ggtree(nj_tree_thiopeptina_ycao, layout="circular", size=0.2)+
  geom_tiplab(size=0.5, aes(angle=angle))
print(nj_tree_thiopeptina_ycao_plot)



nj_tree_syomicina_ycao <- nj(dm_syomicina_ycao)
plot(nj_tree_syomicina_ycao)

nj_tree_GE37468A_ycao <- nj(dm_GE37468A_ycao)
#plot(nj_tree_GE37468A_ycao)
nj_tree_GE37468A_ycao_plot <- ggtree(nj_tree_GE37468A_ycao, layout="circular", size=0.2)+
  geom_tiplab(size=0.5, aes(angle=angle))
print(nj_tree_GE37468A_ycao_plot)


nj_tree_nosiheptide_ycao <- nj(dm_nosiheptide_ycao)
plot(nj_tree_nosiheptide_ycao)




#Maximum parsimony desde terminal IQtree

####syomicina

#Cargar mir arbol
syomicina_tree_ycao <- read.tree("03_Results/syomicina_ycao_tree/msa_syomicina_ycao_ajustado_Clipkit.fasta.treefile")

#Visaulizar el arbol
plot(syomicina_tree_ycao, main= "Arbol ML de Ycao de Syomicina")
nodelabels(syomicina_tree_ycao$node.label, cex=1.5, frame="none")


####Nosiheptide 

nosiheptide_tree_Ycao <- read.tree("03_Results/nosiheptide_ycao_tree/msa_nosiheptide_ycao_ajustado_Clipkit.fasta.treefile")

plot(nosiheptide_tree_Ycao, main="Arbol ML de YcaO de nosiheptide")
nodelabels(nosiheptide_tree_Ycao$node.label, cex=1.2, frame="none")
