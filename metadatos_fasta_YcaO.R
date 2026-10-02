library(Biostrings)
library(dplyr)
library(tidyverse)


fasta_metadatos <- readAAStringSet("01_RawData/ycao_Thiopeptina_syomicina_GE37468_nosiheptide/YcaO_combinados_thiopeptina_syiomicina_GE_Nosiheptide.txt")

fasta_metadatos_M <- as.data.frame(fasta_metadatos)
View(fasta_metadatos_M)

meta <- rownames(fasta_metadatos_M)
class(meta)

meta <- as.data.frame(meta)
#meta <- as.matrix(meta)
metacopia <- meta


meta_separada <- separate_wider_delim(meta, meta, sep="/",into= c("access", "especies"))

str_meta <- str_extract(metacopia[,1], "(?<=\\[).*?(?=\\])")
especies <- str_meta
especies

acceso_number <- str_extract(metacopia[,1], "") 

dim(metacopia)

IA_separacion <-metacopia %>%
  separate_wider_regex(
    meta,
    patterns = c(
      prefijo = ".*_",       # Lo que está antes del guion bajo
      numero  = "\\d+",      # Extrae los dígitos principales
      version = "\\.\\d+",   # Extrae el punto y versión
      resto   = ".*"         # El resto de la cadena
    ),
    cols_remove = FALSE,      # Mantiene la columna original
    too_few="align_start"
  )

metacopia

#############
str_meta <- str_extract(metacopia[,1], "(?<=\\[).*?(?=\\])")
especies <- str_meta
df_especies <- as.data.frame(especies)
class(especies)
saveRDS(df_especies, "03_Results/Msa_YcaO_juntas_arbol/metadata.rds")
view(df_especies)
