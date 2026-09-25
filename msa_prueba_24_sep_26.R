#PRUEBAS PARA HACER UN MSA PATITO
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

install.packages("seqinr")
BiocManager::install("msa")
BiocManager::install("ggmsa")
BiocManager::install("Biostrings")
install.packages("ape")
BiocManager::install("ggtree")

library(Biostrings)
library(msa)
library(seqinr)
library(ggmsa)
library(ape)
library(ggtree)



secuencias_fasta <- readAAStringSet("01_RawData/txt-to-fasta_17-sequences.fasta")

alineamiento <- msa(secuencias_fasta)

class(alineamiento)

saveWidth <- getOption("width")
options(width=100)
sink("03_Results/myAlignment.txt")
print(alineamiento, show="complete", halfNrow=-1)
sink()
options(width=saveWidth)


class(alineamiento) <- "AAMultipleAlignment"

pdf("03_Results/alineamiento_prueba.pdf", width = 60, height = 10)
# 2. Create the plot
ggmsa(alineamiento, char_width = 0.5, seq_name = T)+
  geom_seqlogo()+
  geom_msaBar()
# 3. Close the file
dev.off()

dim(alineamiento)

?readAAMultipleAlignment
leer_msa <- readAAMultipleAlignment(alineamiento, "clustal")

class(alineamiento)

