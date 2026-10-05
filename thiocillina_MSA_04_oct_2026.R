library(Biostrings)
library(msa)
library(seqinr)
library(ggmsa)
library(ape)
library(ggtree)
library(DECIPHER)
library(microseq)
library(heattree)

thiocillina_ycao <- readAAStringSet("01_RawData/thiocillina_4_OCT_2026/thiocilina_YcaO_041026.txt")
thiocillina_ycao #son 500 seqs
summary(width(thiocillina_ycao))


#MSA

thiocillina_MSA_YcaO <- msa(thiocillina_ycao, method= "Muscle")
#No me deja hacer el alineamiento, MUSCLE finaliza sin razon alguna, 
#ando corriendo al mismo tiempo el IQTREE, sera que no hay suficiente CPU para el MSA
thiocillina_MSA_ycao_compatible <- unmasked(thiocillina_MSA_YcaO)
BrowseSeqs(thiocillina_MSA_ycao_compatible)
