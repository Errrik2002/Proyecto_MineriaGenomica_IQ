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



#Esto es para guardar el alineamiento, pero me lo da mal
saveWidth <- getOption("width")
options(width=100)
sink("03_Results/myAlignment.txt")
print(alineamiento, show="complete", halfNrow=-1)
sink()
options(width=saveWidth)

#correr este para hacer el ggmsa
class(alineamiento) <- "AAMultipleAlignment"

pdf("03_Results/alineamiento_prueba.pdf", width = 600, height = 75)
# 2. Create the plot
ggmsa(alineamiento, char_width = 0.5, seq_name = T)+
  geom_seqlogo()+
  geom_msaBar()
# 3. Close the file
dev.off()


observar_ggmsa <- ggmsa(alineamiento,color="Clustal",
font="DroidSansMono", char_width=0.5, seq_name = T)


#observar_ggmsa
class(alineamiento)

alineamiento <- msaConvert(alineamiento, type="seqinr::alignment")
#?msaConvert()

d <- dist.alignment(alineamiento, "identity")
as.matrix(d)

ycaotree <- nj(d)

plot(ycaotree)

png(file="03_Results/ycao_thiopeptin_.png",
width=595, height=600)
plot(ycaotree)
dev.off()

class(ycaotree)

arbolcon_ggtree <- ggtree(ycaotree)+
  geom_tiplab()+
  geom_nodelab(geom='label')+
  hexpand(0.05)

png(file="03_Results/ycao_thiopeptin_.png",
width=2000, height=350)
ggtree(ycaotree)+
  geom_tiplab()+
  geom_nodelab(geom='label')+
  hexpand(0.05)
dev.off()