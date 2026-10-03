library(Biostrings)
library(dplyr)
library(tidyverse)

macro_metadatos <- readRDS("03_Results/macrociclacion_resultados_021026/secuencias_fasta_unicas_macrociclacion.fasta")
macro_metadatos

macro_metadatos <- as.data.frame(macro_metadatos)
class(macro_metadatos)

macro_metadatos <- rownames(macro_metadatos)
macro_metadatos <- as.data.frame(macro_metadatos)

macro_metadatos_especies <- str_extract(macro_metadatos[,1], "(?<=\\[).*?(?=\\])")
##### Aqui ya estan las especies. #####
macro_metadatos_especies <- as.data.frame(macro_metadatos_especies)



##### obtener el valor del Access
max_words_macro <- max(str_count(macro_metadatos$macro_metadatos, "\\s+")) +1
max_words_macro

nombre_col_macro <- paste0("palabra_", 1:max_words_macro)

meta_macro_separada <- macro_metadatos |> 
  separate_wider_delim(
cols= macro_metadatos,
delim= regex("\\s+"),
names= nombre_col,
too_few= "align_start"
  )

meta_macro_separada
#### Numeros de acceso #####
numero_acceso <- meta_macro_separada[,1]




##### Separar el genero de la especie
macro_metadatos_generos <- macro_metadatos_especies |> 
  separate_wider_regex(
    col= macro_metadatos_especies,
    patterns = c(
      genero = "^[^ ]+", #todo hasta el primer espacio
    " +",
    resto = ".*$"
    ),
    too_few = "align_start"
  )

macro_metadatos_generos <- macro_metadatos_generos[,1]


###################### Crear la metadata completa #######################

macro_meta_completa <- bind_cols(macro_metadatos_generos$genero, 
  macro_metadatos_especies$macro_metadatos_especies,
numero_acceso$palabra_1)

macro_nombre_col <- c("genero", "especie", "number_access_node_id")

colnames(macro_meta_completa) <- macro_nombre_col

macro_meta_completa
saveRDS(macro_meta_completa, "03_Results/macrociclacion_resultados_021026/metadatos_completos_arbol.rds")
