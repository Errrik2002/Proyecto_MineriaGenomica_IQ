library(pacman)
library(Biostrings)
library(msa)
library(seqinr)
library(ggmsa)
library(ape)
library(ggtree)
library(DECIPHER)


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


#Thiopeptina
msa_thiopeptina_Ycao <- msa(thiopeptina_YcaO)

msa_thiopeptina_Ycao_Decipher <- AlignSeqs(thiopeptina_YcaO)
#DECIPHER no puedo modificar el algoritmo


#Visualizar MSA en Browser

#Pasar el msa original a XStringSet
msa_thiopeptina_Ycao_compatible <- unmasked(msa_thiopeptina_Ycao)
BrowseSeqs(msa_thiopeptina_Ycao_compatible, highlight = 1)

#### Guardar el MSA ####

writeXStringSet(msa_thiopeptina_Ycao_compatible, "03_Results/msa_thiopeptina.fasta")

#Syomicina

msa_syiomicina_Ycao <- msa