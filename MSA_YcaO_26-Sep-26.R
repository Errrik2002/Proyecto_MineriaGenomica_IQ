#Guardar objetos
# saveRDS()
# readRDS()

install.packages("phangorn")
install.packages("igraph")

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


#Cargar las secuencias 
thiopeptina_YcaO <- readAAStringSet("01_RawData/thiopeptina_BGC0001474_blast_240926/ycao_blast_sequences/thiopeptina_PYC80225.1_YcaO.txt")
thiopeptina_YcaO #No debo pasarlo a fasta necesariamente. 

syomicina_YcaO <- readAAStringSet("01_RawData/syiomicina_BGC0000611_BLAST_26-09-26/siomicina_GAP53284.1_YCAO.txt")
syomicina_YcaO

GE37468A_YcaO <- readAAStringSet("01_RawData/GE37468A_BGC0000605_BLAST_240926/GE37468A_AEM00626.1_YCAO.txt")
GE37468A_YcaO

nosiheptide_YcaO <- readAAStringSet("01_RawData/nosiheptide_BGC0000610_BLAST_24-09-26/nosiheptide_ACR48336.1_YCAO.txt")
nosiheptide_YcaO



#### Realizar los Multiple Sequence Alignment ####
?msa()

###Thiopeptina
msa_thiopeptina_Ycao <- msa(thiopeptina_YcaO)
msa_thiopeptina_Ycao_Decipher <- AlignSeqs(thiopeptina_YcaO)
#DECIPHER no puedo modificar el algoritmo

#Visualizar MSA en Browser
#Pasar el msa original a XStringSet
msa_thiopeptina_Ycao_compatible <- unmasked(msa_thiopeptina_Ycao)
BrowseSeqs(msa_thiopeptina_Ycao_compatible, highlight = 0)

#### Guardar el MSA ####
writeXStringSet(msa_thiopeptina_Ycao_compatible, "03_Results/msa_thiopeptina.fasta")



###Syomicina
#msa
msa_syiomicina_Ycao <- msa(syomicina_YcaO)
msa_syiomicina_Ycao_decipher <- AlignSeqs(syomicina_YcaO)

#Visualizar en Browser, hacerlo compatible con BrowSeqs
msa_syiomicina_Ycao_compatible <- unmasked(msa_syiomicina_Ycao)
BrowseSeqs(msa_syiomicina_Ycao_compatible)

#Guardar MSA
writeXStringSet(msa_syiomicina_Ycao_compatible, "03_Results/msa_syomicina_ycao.fasta")




###GE37468A
msa_GE37468A_Ycao <- msa(GE37468A_YcaO)
msa_GE37468A_Ycao_DECIPHER <- AlignSeqs(GE37468A_YcaO)

#Visualizar en Browser
msa_GE37468A_Ycao_compatible <- unmasked(msa_GE37468A_Ycao)
BrowseSeqs(msa_GE37468A_Ycao_compatible)

#Guardar MSA 
writeXStringSet(msa_GE37468A_Ycao_compatible, "03_Results/msa_GE37468A_Ycao.fasta")


###Nosiheptide
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



#### POSTERIOR AL AJUSTE POR CLIPKIT ####

#Thiopeptina
msa_thiopeptina_ycao_AJclipkit <- readAAMultipleAlignment("03_Results/msa_thiopeptina_ycao_AJclipkit.fasta")
msa_thiopeptina_ycao_AJclipkit <- unmasked(msa_thiopeptina_ycao_AJclipkit)

msa_thiopeptina_original <- readAAMultipleAlignment("03_Results/msa_thiopeptina.fasta")
msa_thiopeptina_original <- unmasked(msa_thiopeptina_original)

BrowseSeqs(msa_thiopeptina)
BrowseSeqs(msa_thiopeptina_ycao_AJclipkit)



#Syomicina
msa_syiomicina_Ycao_Ajustado_clipkit <- readAAMultipleAlignment("03_Results/msa_syomicina_ycao_ajustado_Clipkit.fasta")
msa_syiomicina_Ycao_Ajustado_clipkit <- unmasked(msa_syiomicina_Ycao_Ajustado_clipkit)

msa_syomicina_original <- readAAMultipleAlignment("03_Results/msa_syomicina_ycao.fasta")
msa_syomicina_original <- unmasked(msa_syomicina_original)

BrowseSeqs(msa_syiomicina_Ycao_Ajustado_clipkit)
BrowseSeqs(msa_syomicina_original)



#GE37468A
msa_GE37468A_Ycao_ajustado_Clipkit <- readAAMultipleAlignment("03_Results/msa_GE37468A_Ycao_ajustado_Clipkit.fasta")
msa_GE37468A_Ycao_ajustado_Clipkit <- unmasked(msa_GE37468A_Ycao_ajustado_Clipkit)

msa_GE37468A_Ycao_original <- readAAMultipleAlignment("03_Results/msa_GE37468A_Ycao.fasta")
msa_GE37468A_Ycao_original <- unmasked(msa_GE37468A_Ycao_original)

BrowseSeqs(msa_GE37468A_Ycao_original)
BrowseSeqs(msa_GE37468A_Ycao_ajustado_Clipkit)



#nosiheptide
msa_nosiheptide_Ycao_ajustado_clipkit <- readAAMultipleAlignment("03_Results/msa_nosiheptide_ycao_ajustado_Clipkit.fasta")
msa_nosiheptide_Ycao_ajustado_clipkit <- unmasked(msa_nosiheptide_Ycao_ajustado_clipkit)

msa_nosiheptide_Ycao_original <- readAAMultipleAlignment("03_Results/msa_nosiheptide_ycao.fasta")
msa_nosiheptide_Ycao_original <- unmasked(msa_nosiheptide_Ycao_original)

BrowseSeqs(msa_nosiheptide_Ycao_original)
BrowseSeqs(msa_nosiheptide_Ycao_ajustado_clipkit)
