install.packages("devtools")
install.packages("geneviewer")


library(devtools)
library(geneviewer)


#Leer los resultados de antismash
seq_gbk <- read_gbk("03_Results/results_antismash/Gen_BankCP114205.1/CP114205.1.region003.gbk")
seq_gbk

seq_gbk[1]

m_seq_gbk <- as.matrix(seq_gbk)
m_seq_gbk <- m_seq_gbk[1]

#### Sacar los features es lo mas esencial #####
features_reg03 <- seq_gbk$CP114205.1.region003$FEATURES
features_reg03

####caracteristicas de cada una

db_feat_reg03 <- as.matrix(features_reg03)
db_feat_reg03 #marca los tipos de listas que hay dentro del resultado de antismash


#Aqui se puede ver el tipo de producto
reg_reg03 <- features_reg03$region
reg_reg03  <- as.data.frame(reg_reg03)
View(reg_reg03)


#cand_cluster 
candcluster <- as.data.frame(features_reg03$cand_cluster)
View(candcluster)

#protocluster
protocluster <- as.data.frame(features_reg03$protocluster)
View(protocluster)

#protocore 
protocore <- as.data.frame(features_reg03$proto_core)
View(protocore)

#gene
gene <- features_reg03$gene
gene
length(gene) #por cada 2 objetos de lista es un gen, si hay 44 objetos, hay 22 genes
#Aqui vienene las demas regiones que componen el BGC.
#Vienen en listas

