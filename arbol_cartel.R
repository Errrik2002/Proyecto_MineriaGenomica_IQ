BiocManager::install(c("ggtree", "tidytree"))
install.packages("ggplot2")

library(tidyverse)
library(tidytree)
library(ggtree) # for plotting tree
library(phangorn) # for midpoint root of tree
library(ggnewscale) #for new color/fill scale
library(viridis)
library(ggtree) # for plotting tree
library(phangorn) # for midpoint root of tree
library(treeio) # for loading special tree 
library(ggplot2)


plottree <- read.tree("03_Results/conjunto_noactinomicetos_YcaO/conjunto_MSA_ajustado_clipkit.fasta.treefile")

#meta <- read.table('Tree_tutorial/Phylogeny_module_tutorial_meta_data.csv', sep=',', header=T)
metadata_claude <- read.table("03_Results/conjunto_noactinomicetos_YcaO/metadata_enriched.csv",sep=',', header=TRUE, fill=TRUE)

vshv.tree <- midpoint(plottree)
vshv.tree

vshv.tree <- read.iqtree('03_Results/conjunto_noactinomicetos_YcaO/conjunto_MSA_ajustado_clipkit.fasta.treefile')
vshv.tree@phylo <- midpoint(vshv.tree@phylo)

# Just tree
plot(vshv.tree@phylo)

t1 <- ggtree(vshv.tree) %<+% metadata_claude +
  geom_tippoint(aes(color=host))

plot(t1)
t1 <- ggtree(vshv.tree, layout="circular") %<+% metadata_claude +
  geom_tippoint(aes(color= host))

length(metadata_claude$host)


# See the first few tip labels of your tree
head(vshv.tree@phylo$tip.label)

# See the column names and first few rows of your metadata
head(metadata_claude)
