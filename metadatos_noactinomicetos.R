#install.packages("tidyverse")

library(Biostrings)
library(dplyr)
library(tidyverse)


no_actino_metadatos <- readRDS("03_Results/conjunto_noactinomicetos_YcaO/secuencias_unicas_YcaO.rds")
no_actino_metadatos

no_actino_metadatos_df <- as.data.frame(no_actino_metadatos) 

especies_noactino <- rownames(no_actino_metadatos_df) |> 
   as.data.frame() 


especies_nombres <- str_extract(especies_noactino[,1], "(?<=\\[).*?(?=\\])")
especies_nombres <- as.data.frame(especies_nombres)
class(especies_nombres)
especies_nombres

#### LOS WP de acceso 

max_words_noactino <- max(str_count(especies_noactino$`rownames(no_actino_metadatos_df)`, "\\s+")) +1
max_words_noactino

nombre_col_noactino <- paste0("palabra_", 1:max_words_noactino)

meta_macro_noactino <- especies_noactino |> 
  separate_wider_delim(
cols=`rownames(no_actino_metadatos_df)` ,
delim= regex("\\s+"),
names= nombre_col_noactino,
too_few= "align_start"
  )

meta_macro_noactino
#### Numeros de acceso #####
numero_acceso_noactino <- meta_macro_noactino[,1]
numero_acceso_noactino

especies_nombres$
#########Separar genero de la especie
generos_noactino <- especies_nombres |> 
  separate_wider_regex(
    col= especies_nombres,
    patterns = c(
      genero = "^[^ ]+", #todo hasta el primer espacio
    " +",
    resto = ".*$"
    ),
    too_few = "align_start"
  )

generos_noactino <- generos_noactino[,1]
generos_noactino

#metadata completa

metadata_noactinos <- bind_cols(generos_noactino$genero,
   especies_nombres$especies_nombres,
 numero_acceso_noactino$palabra_1)

columnas_noactino <- c("genero", "especie","numero_access_node_id")

colnames(metadata_noactinos) <- columnas_noactino
View(metadata_noactinos)

saveRDS(metadata_noactinos, "03_Results/conjunto_noactinomicetos_YcaO/metadata_base_FASTA.rds")


write.csv(metadata_noactinos,"03_Results/conjunto_noactinomicetos_YcaO/metadata_noactinos.csv", row.names = FALSE)
metadata_hoy <- readRDS("03_Results/conjunto_noactinomicetos_YcaO/metadata_base_FASTA.rds")
View(metadata_hoy)


write.csv(metadata_hoy$numero_access_node_id, "03_Results/conjunto_noactinomicetos_YcaO/prueba.csv")


metadatos_completos_Claude <- read.csv("03_Results/conjunto_noactinomicetos_YcaO/metadata_enriched.csv", header= TRUE)
dim(metadatos_completos_Claude)
dim(metadata_hoy)


View(metadatos_completos_Claude)

metadata_conQuery <- add_column(metadatos_completos_Claude)
?add_column()

View(metadata_conQuery)

df <- df %>%
        mutate(status= if_else(.$points > 20, 'Good', 'Bad'))

#, "WP_172686359.1", "AAP11954.1"

metadata_conQuery <- metadata_conQuery |> 
  mutate(query = if_else(metadata_conQuery$accession %in% c("WP_259592873.1","WP_172686359.1", "AAP11954.1" ), 'query', 'no'))  

View(metadata_conQuery)
metadata_conQuery <- as.data.frame(metadata_conQuery)

saveRDS(metadata_conQuery, "03_Results/conjunto_noactinomicetos_YcaO/metadata_con_Query.rds")
