if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes")

remotes::install_github("YuLab-SMU/iggtree")

BiocManager::install("ggimage")


library(Biostrings)
library(msa)
library(seqinr)
library(ggmsa)
library(ape)
library(ggtree)
library(DECIPHER)
library(microseq)
library(heattree)
library(iggtree)
library(ggplot2)
library(ggiraph)


set.seed(123)
tr <- rtree(20)
dt <- data.frame(id = c(36, 38), type=c("A", "B"))
p <- ggtree(
       tr, 
       mapping = aes(
           tooltip = round(branch.length, 2), 
           data_id = node
         )
     ) +
     geom_hilight(
        data = dt, 
        mapping = aes(
           node = id, 
           fill = type, 
           tooltip = paste0("clade of node ", id), 
           data_id = type
        ), 
        to.bottom = TRUE
     )

girafe(ggobj = p)
