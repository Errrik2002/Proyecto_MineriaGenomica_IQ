install.packages("baseq")
install.packages("GaMaBioMD")


library(tidyverse)
library(Biostrings)
library(msa)
library(seqinr)
library(ggmsa)
library(ape)
library(ggtree)
library(DECIPHER)
library(microseq)
library(heattree)
library(htmltools)
library(baseq)
library(phytools)

metadata_comparar <- readRDS("03_Results/conjunto_noactinomicetos_YcaO/metadata_homologos.rds")
View(metadata_comparar)
metadata_bacillus <- subset(metadata_comparar,metadata_comparar$query_homologos == "Bacillus cereus")
View(metadata_bacillus)
print(metadata_bacillus$accession)

secuencias_bacillus <- readRDS("03_Results/conjunto_noactinomicetos_YcaO/secuencias_unicas_YcaO.rds")
class(secuencias_bacillus)
View(secuencias_bacillus)

secuencias_bacillus <- as.data.frame(secuencias_bacillus)
nombres <- rownames(secuencias_bacillus)
rownames(secuencias_bacillus) <- NULL
secuencias_bacillus_data_2Cols <- cbind(nombres, secuencias_bacillus)
View(secuencias_bacillus_data_2Cols)

seq_accession <- cbind(secuencias_bacillus_data_2Cols, metadata_comparar$accession)
View(seq_accession)

baccillus_a_comparar <- seq_accession[,]
View(baccillus_a_comparar)

bacillus1 <- subset(baccillus_a_comparar, baccillus_a_comparar$`metadata_comparar$accession` == "AAP11954.1")
bacillus2 <- subset(baccillus_a_comparar, baccillus_a_comparar$`metadata_comparar$accession` == "WP_098594832.1")
bacillus3 <- subset(baccillus_a_comparar, baccillus_a_comparar$`metadata_comparar$accession` == "TKI10472.1")

bacillus_totales <- rbind(bacillus1, bacillus3, bacillus2)
View(bacillus_totales)


bacillus_seqs <- readAAStringSet("01_RawData/conjunto_thiopeptidos_no_actinomicetos/bacillus.txt")
bacillus_seqs

summary(width(bacillus_seqs))

msa_bacillus_seq <- msa(bacillus_seqs, method= "Muscle")
comp_msa_bacillus_seq <- unmasked(msa_bacillus_seq)
BrowseSeqs(comp_msa_bacillus_seq)

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

porcentaje_de_gaps(msa_bacillus_seq)
