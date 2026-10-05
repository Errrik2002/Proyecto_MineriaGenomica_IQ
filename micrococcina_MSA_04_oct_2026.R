library(Biostrings)
library(msa)
library(seqinr)
library(ggmsa)
library(ape)
library(ggtree)
library(DECIPHER)
library(microseq)
library(heattree)


microccina_ycao <- readAAStringSet("01_RawData/microccoina_4_10_26/microccocina_PI30_QC60_4_oct_2026.txt")
microccina_ycao #73 secuencias horizontales
summary(width(microccina_ycao))


microccocina_MSA_Ycao <- msa(microccina_ycao, method= "Muscle")
microccocina_MSA_Ycao_comptible <- unmasked(microccocina_MSA_Ycao)
BrowseSeqs(microccocina_MSA_Ycao_comptible)
writeXStringSet(microccocina_MSA_Ycao_comptible, "03_Results/micrococcina_ycao_041026/micrococina_MSA_ycao_original.fasta")

