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

#Extrae todo lo de dentro de corchetes
str_meta <- str_extract(metacopia[,1], "(?<=\\[).*?(?=\\])")
especies <- str_meta


#############
str_meta <- str_extract(metacopia[,1], "(?<=\\[).*?(?=\\])")
especies <- str_meta
df_especies <- as.data.frame(especies)

saveRDS(df_especies, "03_Results/Msa_YcaO_juntas_arbol/metadata.rds")




#####Separar base de datos por espacios####

max_words <- max(str_count(metacopia$meta, "\\s+")) +1
max_words

nombre_col <- paste0("palabra_", 1:max_words)

meta_separada <- metacopia |> 
  separate_wider_delim(
cols= meta,
delim= regex("\\s+"),
names= nombre_col,
too_few= "align_start"
  )
View(meta_separada)
dim(metadata_separada_acceso)
metadata_separada_acceso <- meta_separada[,1]


especies <- as.data.frame(especies)
meta_completa <- merge(especies, metadata_separada_acceso)
meta_completo_1 <- bind_cols(especies$especies, metadata_separada_acceso$palabra_1)
meta_completo_1


columnas <- c("especies", "numero_acceso_node_id")
colnames(meta_completo_1) <- columnas
head(meta_completo_1)
dim(meta_completo_1)

saveRDS(meta_completo_1, "03_Results/Msa_YcaO_juntas_arbol/metadatos_arbol_cao.rds")


head(meta_completo_1)

df_separado <- meta_completo_1 %>%
  separate_wider_regex(
    col = especies,
    patterns = c(
      primera_palabra = "^[^ ]+", # Captura todo hasta antes del primer espacio
      " +",                       # Coincide con el primer espacio (o espacios)
      resto_frase = ".*$"        # Captura todo lo que queda después
    ),
    too_few="align_start"
  )

df_separado <- df_separado[,-2]


metadata_genero_especie_acceso <- bind_cols(especies$especies, 
  df_separado$primera_palabra, 
  df_separado$numero_acceso_node_id)



metadata_genero_especie_acceso
columna_completa <- c("espcie","genero", "acces_number_node_id")

colnames(metadata_genero_especie_acceso) <- columna_completa

head(metadata_genero_especie_acceso)
saveRDS(metadata_genero_especie_acceso, "03_Results/Msa_YcaO_juntas_arbol/metadata_especie_genero_access_number.rds")
