# ========================================
# WORDCLOUD2
# ========================================

library(officer)
library(dplyr)
library(stringr)
library(tibble)
library(tidytext)
library(stopwords)
library(wordcloud2)

# ---- DOCX ----
read_docx_text <- function(path) {
  doc <- read_docx(path)
  txt <- docx_summary(doc)$text
  paste(txt, collapse = " ")
}

texto <- read_docx_text("//Memoria/TFM_wordcloud.docx")

# ---- Stopwords ----
stop_es <- stopwords("es")
stop_en <- stopwords("en")

stop_custom <- c(stop_es,stop_en,"Fig","Figura", "Tabla", "número", "puede", "forma", 
  "compatible", "resultado", "porcentaje", "mientras", "estrategias", "trabajo", 
  "relacionados", "doi", "realizó", "procesos", "nature", "negativo","dcha",
  "human", "expression", "front", "chen", "aunque", "relacionados", "clúster", "tras", 
  "analysis", "mostró", "ello", "resultados", "dos", "immunol", "según", "tabla", "realizó",
  "mismo", "frente", "nivel", "junto", "partir", "además", "pageref", "cells", "así", 
  "cell", "mayor","distintos","gene","mayor","menos","menor","motivos","relevante",
  "métodos", "dentro", "pueden", "manera", "posteriormente", "apartado", "central",
  "infection", "lugar", "modo", "proc","ser", "tamaño","tsv","posterior", "usando",
  "diferentes", "global", "estudio", "grupo", "contexto","ausencia","servidor",
  "lado", "versión","obstante","sino","tablas","varios","vez","memory","eje",
  "basado", "caso","aociados","efecto","regulatory","up-regulados","Res","coherente",
  "principalmente", "resto","procesos","permite","presencia","especialmente",
  "adj", "Nat", "Zhang", "MERGEFORMAT","ARABIC","REF","SEQ","muestran","fila","baja",
  "obtenidos","asociados","tres","abajo","arriba","filas","Representación","correspondientes",
  "fracción","izda","basados","indica","método","posteriores","representativos","superior",
  "biológica","ambos","compartidos","detectados","hacia","izquierda","previamente","sugiere",
  "utilizó","reducción","total","base","función","relativa","Tem","fase", "diferencialmente",
  "representa","tipo","cada", "mediante")

# stopwords a minúscula, pero no para las palabras originales
stop_custom <- str_to_lower(stop_custom)

# ---- Frecuencias ----
freq_words <- tibble(text = texto) |> unnest_tokens(output = word, input = text,
    token = "regex", pattern = "[^\\p{L}\\p{N}\\-]+", to_lower = FALSE) |>
  mutate(word = case_when(
    str_to_lower(word) %in% c("vih", "hiv") ~ "VIH",
    str_to_lower(word) %in% c("vih-1", "vih1", "hiv-1", "hiv1") ~ "VIH-1",
    str_to_lower(word) == "scrna-seq" ~ "scRNA-seq",
    str_to_lower(word) == "rna-seq" ~ "RNA-seq",
    str_to_lower(word) == "cd4" ~ "CD4",
    TRUE ~ word)) |> filter(!str_to_lower(word) %in% stop_custom) |>
  filter(!str_detect(word, "^[0-9]+$")) |> filter(!str_detect(word, "^[0-9\\-\\+]+$")) |>
  filter(str_length(word) > 2) |> count(word, sort = TRUE)

# ---- Top 150 palabras ----
freq_top <- freq_words |> slice_max(n, n = 150)
freq_top$word

pal_10 <- c("#F8766D","#00BF7D","#00BFC4","#00B0F6","#9590FF","#E76BF3",
  "#D89000","#FF62BC","#A3A500","#39B600")

wordcloud(data = freq_top,size = 0.9,minSize = 8,
           color = rep(pal_10,length.out = nrow(freq_top)))

# ---- Nube completa ----
wordcloud2(data = freq_words, size = 0.7, minSize = 5)
