# ==============================================================================
# TFM Bioinformática - VIU - Mayo 2026
#
# Caracterización de la heterogeneidad celular, transcriptómica y reguladora
# de linfocitos T CD4+ en la reversión de latencia del VIH-1 mediante scRNA-seq
#
# Datos obtenidos de repositorio público: GSE210824 / PRJNA867681
#
# ==============================================================================

## ===============================
# Script: 03b_ann_clust_gsea.R
## ===============================

# Objetivo:
# Realizar anotación funcional de clusters mediante GSEA pre-rankeado a partir
# de genes marcadores globales, conservados y filtrados, generar tablas y figuras
# de enriquecimiento funcional, añadir anotaciones manuales de cluster al objeto
# Seurat integrado y exportar un archivo loom anotado para pySCENIC.

# Inputs en Proc/:
# - 01_data_prep_and_QC.rds
# - 03_02_int_ann.rds

# Inputs en Ann/:
# - BloodGen3Module.rda

# Inputs en Results/:
# - Tables/02_03_FindAllMarkers_clust_MAST.tsv
# - Tables/02_03_FindConservedMarkers_clust_MAST.tsv

# Outputs en Proc/:
# - 03b_04_int_ann.rds

# Outputs en Results/Figures y Results/Tables

# Outputs en SCENIC/:
# - 02_postprocessed_inputs/ann.loom
# - 01_preprocessed_inputs/metadata.tsv



#----------------------------
library(dplyr)

markers.all <- read.delim("Results/Tables/02_03_FindAllMarkers_clust_MAST.tsv")
markers.cons <- read.delim("Results/Tables/02_03_FindConservedMarkers_clust_MAST.tsv")

## ========================================================================
## BLOQUE 1 - LISTAS DE RANKS POR CADA ESTRATEGIA DE D.E.
## ========================================================================

## ---- 1.1 lista ranks GSEA de FindAllMarkers filtrando significativos y ranking con score propio ----

# Eliminar all markers con p-valor > 0.05 y calcular score de rankeado
markers.all.rank <- markers.all |> filter(p_val_adj <= 0.05)
markers.all.rank <- markers.all.rank |> mutate(
  score = case_when(avg_log2FC > 0 ~ avg_log2FC * pct.1 * (pct.1 / (pct.1 + pct.2)),
                    avg_log2FC < 0 ~ avg_log2FC * pct.2 * (pct.2 / (pct.1 + pct.2)),
                    TRUE ~ 0))

# separar por cluster
markers.all.clust <- split(markers.all.rank, markers.all.rank$cluster)

# rankear cada cluster por score (desempate por p-val)
markers.all.clust <- lapply(markers.all.clust, function(df) {
  df |> arrange(desc(score), p_val_adj)})

# convertir cada cluster en un vector rankeado
ranks.all <- lapply(markers.all.clust, function(df) {
  ranks <- df$score
  names(ranks) <- df$gene
  ranks})


# ---- 1.2 lista ranks GSEA de FindConservedMarkers, filtrado significativos y ranking con score propio ----

# Eliminar conserved markers con adj p-valor > 0.05 y calcular score de rankeado
markers.cons.rank <- markers.cons |> filter(Control_p_val_adj <= 0.05, HA15_p_val_adj <= 0.05)
markers.cons.rank <- markers.cons.rank |>
  mutate(
    score_control = case_when(
      Control_avg_log2FC > 0 ~ Control_avg_log2FC * Control_pct.1 * (Control_pct.1 / (Control_pct.1 + Control_pct.2)),
      Control_avg_log2FC < 0 ~ Control_avg_log2FC * Control_pct.2 * (Control_pct.2 / (Control_pct.1 + Control_pct.2)),
      TRUE ~ 0
    ),
    score_HA15 = case_when(
      HA15_avg_log2FC > 0 ~ HA15_avg_log2FC * HA15_pct.1 * (HA15_pct.1 / (HA15_pct.1 + HA15_pct.2)),
      HA15_avg_log2FC < 0 ~ HA15_avg_log2FC * HA15_pct.2 * (HA15_pct.2 / (HA15_pct.1 + HA15_pct.2)),
      TRUE ~ 0
    ),
    weight_control = Control_pct.1 / (Control_pct.1 + HA15_pct.1),
    weight_HA15 = HA15_pct.1 / (Control_pct.1 + HA15_pct.1),
    score = score_control * weight_control + score_HA15 * weight_HA15
  )

# separar por cluster
markers.cons.clust <- split(markers.cons.rank, markers.cons.rank$cluster)

# rankear cada cluster por score (desempate por p-val)
markers.cons.clust <- lapply(markers.cons.clust, function(df) {
  df |> arrange(desc(score), max_pval)})

# convertir cada cluster en un vector rankeado
ranks.cons <- lapply(markers.cons.clust, function(df) {
  ranks <- df$score
  names(ranks) <- df$gene
  ranks})
ranks.cons <- lapply(ranks.cons, function(x) x[!is.na(x)])


## ---- 1.3 lista ranks GSEA de FindAllMarkers filtrando no significativos y no conservados y ranking con score propio ----

# Eliminar conserved markers con p-valor > 0.05 en cualquier grupo
markers.cons.filt <- markers.cons |> filter(Control_p_val_adj <= 0.05, HA15_p_val_adj <= 0.05)

# Eliminar all markers con p-valor > 0.05
markers.all.filt <- markers.all |> filter(p_val_adj <= 0.05)

# comprobar intersección entre pares cluster-gen all y conserved significativos
# nº de pares únicos en cada tabla
markers.cons.filt  |> dplyr::distinct(cluster, gene) |> nrow()
markers.all.filt |> dplyr::distinct(cluster, gene) |> nrow()

# intersección y diferencias
pairs_all  <- markers.all.filt  |> dplyr::distinct(cluster, gene) |>
  dplyr::transmute(pair = paste(cluster, gene, sep = "_")) |> dplyr::pull(pair)
pairs_cons <- markers.cons.filt |> dplyr::distinct(cluster, gene) |>
  dplyr::transmute(pair = paste(cluster, gene, sep = "_")) |> dplyr::pull(pair)

length(intersect(pairs_all, pairs_cons))
length(setdiff(pairs_all, pairs_cons)) # pares solo en all
length(setdiff(pairs_cons, pairs_all)) # pares solo en conserved

# Eliminar pares gen-cluster de all markers no presentes en conserved markers
markers.all.filt <- markers.all.filt |> dplyr::mutate(pair = paste(cluster, gene, sep = "_")) |>
  dplyr::filter(pair %in% pairs_cons) |>
  dplyr::select(-pair)

# Calcular score de rankeado #1 según pct.1/2 y especificidad por cluster
markers.all.filt <- markers.all.filt |> mutate(
  score = case_when(avg_log2FC > 0 ~ avg_log2FC * pct.1 * (pct.1 / (pct.1 + pct.2)),
                    avg_log2FC < 0 ~ avg_log2FC * pct.2 * (pct.2 / (pct.1 + pct.2)),
                    TRUE ~ 0))

# separar por cluster
markers.filt.clust <- split(markers.all.filt, markers.all.filt$cluster)

# rankear cada cluster por score (desempate por p-val)
markers.filt.clust <- lapply(markers.filt.clust, function(df) {
  df |> arrange(desc(score), p_val_adj)})

# convertir cada cluster en un vector rankeado
ranks.filt <- lapply(markers.filt.clust, function(df) {
  ranks <- df$score
  names(ranks) <- df$gene
  ranks})

#length(ranks.filt)
#head(ranks.filt[[8]])
#tail(ranks.filt[[8]])


## =======================================
## BLOQUE 2 - PRE-RANKED GSEA
## =======================================

## ---- 2.1 preparar genesets de DB para GSEA ----
library(fgsea)
library(msigdbr)
library(BloodGen3Module)
load("Ann/BloodGen3Module.rda")

# genesets humanos MSigDB (Hallmark, C7/immune, GO:BP) y BloodGen3Module
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- split(x = hallmark$gene_symbol, f = hallmark$gs_name)

c7 <- msigdbr(species = "Homo sapiens", collection = "C7")
c7 <- split(x = c7$gene_symbol, f = c7$gs_name)

gobp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
gobp <- split(x = gobp$gene_symbol, f = gobp$gs_name)

bg3m <- Module_listGen3 |>
  mutate(pathway = paste(Module, gsub(" ", "_", Function), sep = "_")) |>
  group_by(pathway) |>
  summarise(genes = list(Gene), .groups = "drop") |> tibble::deframe()

## ---- 2.2 pre-ranked GSEA (para ranks.all) ----
# fgsea por cluster + guardar tabla

# hallmark
fgsea.hallmark.all <- lapply(names(ranks.all), function(cl) {
  res <- fgsea(pathways = hallmark, stats = ranks.all[[cl]])
  res$cluster <- cl
  res
})

fgsea.hallmark.all <- bind_rows(fgsea.hallmark.all)
fgsea.hallmark.all.out <- fgsea.hallmark.all |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.table(fgsea.hallmark.all.out,"Results/Tables/02_04_fgsea_hallmark_all.tsv", sep = "\t", row.names = F, quote = F)

# immune
fgsea.c7.all <- lapply(names(ranks.all), function(cl) {
  res <- fgsea(pathways = c7, stats = ranks.all[[cl]])
  res$cluster <- cl
  res
})

fgsea.c7.all <- bind_rows(fgsea.c7.all)
fgsea.c7.all.out <- fgsea.c7.all |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.table(fgsea.c7.all.out,"Results/Tables/02_04_fgsea_c7_all.tsv", sep = "\t", row.names = F, quote = F)

# BP
fgsea.gobp.all <- lapply(names(ranks.all), function(cl) {
  res <- fgsea(pathways = gobp, stats = ranks.all[[cl]])
  res$cluster <- cl
  res
})

fgsea.gobp.all <- bind_rows(fgsea.gobp.all)
fgsea.gobp.all.out <- fgsea.gobp.all |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.table(fgsea.gobp.all.out,"Results/Tables/02_04_fgsea_gobp_all.tsv", sep = "\t", row.names = F, quote = F)

# BloodGen3Module
fgsea.bg3m.all <- lapply(names(ranks.all), function(cl) {
  res <- fgsea(pathways = bg3m, stats = ranks.all[[cl]])
  res$cluster <- cl
  res
})

fgsea.bg3m.all <- bind_rows(fgsea.bg3m.all)
fgsea.bg3m.all.out <- fgsea.bg3m.all |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.table(fgsea.bg3m.all.out,"Results/Tables/02_04_fgsea_bg3m_all.tsv", sep = "\t", row.names = F, quote = F)

## ---- 2.3 pre-ranked GSEA (para ranks.filt) ----
# fgsea por cluster + guardar tabla

# hallmark
fgsea.hallmark.filt <- lapply(names(ranks.filt), function(cl) {
  res <- fgsea(pathways = hallmark, stats = ranks.filt[[cl]])
  res$cluster <- cl
  res
})

fgsea.hallmark.filt <- bind_rows(fgsea.hallmark.filt)
fgsea.hallmark.filt.out <- fgsea.hallmark.filt |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.table(fgsea.hallmark.filt.out,"Results/Tables/02_04_fgsea_hallmark_filt.tsv", sep = "\t", row.names = F, quote = F)

# immune
fgsea.c7.filt <- lapply(names(ranks.filt), function(cl) {
  res <- fgsea(pathways = c7, stats = ranks.filt[[cl]])
  res$cluster <- cl
  res
})

fgsea.c7.filt <- bind_rows(fgsea.c7.filt)
fgsea.c7.filt.out <- fgsea.c7.filt |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.table(fgsea.c7.filt.out,"Results/Tables/02_04_fgsea_c7_filt.tsv", sep = "\t", row.names = F, quote = F)

# BP
fgsea.gobp.filt <- lapply(names(ranks.filt), function(cl) {
  res <- fgsea(pathways = gobp, stats = ranks.filt[[cl]])
  res$cluster <- cl
  res
})

fgsea.gobp.filt <- bind_rows(fgsea.gobp.filt)
fgsea.gobp.filt.out <- fgsea.gobp.filt |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.table(fgsea.gobp.filt.out,"Results/Tables/02_04_fgsea_gobp_filt.tsv", sep = "\t", row.names = F, quote = F)

# BloodGen3Module
fgsea.bg3m.filt <- lapply(names(ranks.filt), function(cl) {
  res <- fgsea(pathways = bg3m, stats = ranks.filt[[cl]])
  res$cluster <- cl
  res
})

fgsea.bg3m.filt <- bind_rows(fgsea.bg3m.filt)
fgsea.bg3m.filt.out <- fgsea.bg3m.filt |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.table(fgsea.bg3m.filt.out,"Results/Tables/02_04_fgsea_bg3m_filt.tsv", sep = "\t", row.names = F, quote = F)

## ---- 2.4 pre-ranked GSEA (para ranks.cons) ----
# fgsea por cluster + guardar tabla

# hallmark
fgsea.hallmark.cons <- lapply(names(ranks.cons), function(cl) {
  res <- fgsea(pathways = hallmark, stats = ranks.cons[[cl]])
  res$cluster <- cl
  res
})

fgsea.hallmark.cons <- bind_rows(fgsea.hallmark.cons)
fgsea.hallmark.cons.out <- fgsea.hallmark.cons |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.table(fgsea.hallmark.cons.out,"Results/Tables/02_04_fgsea_hallmark_cons.tsv", sep = "\t", row.names = F, quote = F)

# immune
fgsea.c7.cons <- lapply(names(ranks.cons), function(cl) {
  res <- fgsea(pathways = c7, stats = ranks.cons[[cl]])
  res$cluster <- cl
  res
})

fgsea.c7.cons <- bind_rows(fgsea.c7.cons)
fgsea.c7.cons.out <- fgsea.c7.cons |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.table(fgsea.c7.cons.out,"Results/Tables/02_04_fgsea_c7_cons.tsv", sep = "\t", row.names = F, quote = F)

# BP
fgsea.gobp.cons <- lapply(names(ranks.cons), function(cl) {
  res <- fgsea(pathways = gobp, stats = ranks.cons[[cl]])
  res$cluster <- cl
  res
})

fgsea.gobp.cons <- bind_rows(fgsea.gobp.cons)
fgsea.gobp.cons.out <- fgsea.gobp.cons |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.table(fgsea.gobp.cons.out,"Results/Tables/02_04_fgsea_gobp_cons.tsv", sep = "\t", row.names = F, quote = F)

# BloodGen3Module
fgsea.bg3m.cons <- lapply(names(ranks.cons), function(cl) {
  res <- fgsea(pathways = bg3m, stats = ranks.cons[[cl]])
  res$cluster <- cl
  res
})

fgsea.bg3m.cons <- bind_rows(fgsea.bg3m.cons)
fgsea.bg3m.cons.out <- fgsea.bg3m.cons |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.table(fgsea.bg3m.cons.out,"Results/Tables/02_04_fgsea_bg3m_cons.tsv", sep = "\t", row.names = F, quote = F)


## =================================================
## BLOQUE 3 - RESULTADOS GSEA - TABLAS Y FIGURAS
## =================================================
library(dplyr)

## ---- 3.1 tablas resumen de resultados significativos y tendencias

fgsea_all_sig <- bind_rows(
  fgsea.hallmark.all  |> mutate(method = "all",  database = "Hallmark"),
  fgsea.c7.all        |> mutate(method = "all",  database = "C7"),
  fgsea.gobp.all      |> mutate(method = "all",  database = "GOBP"),
  fgsea.bg3m.all      |> mutate(method = "all",  database = "BG3M"),
  
  fgsea.hallmark.filt |> mutate(method = "filt", database = "Hallmark"),
  fgsea.c7.filt       |> mutate(method = "filt", database = "C7"),
  fgsea.gobp.filt     |> mutate(method = "filt", database = "GOBP"),
  fgsea.bg3m.filt     |> mutate(method = "filt", database = "BG3M"),
  
  fgsea.hallmark.cons |> mutate(method = "cons", database = "Hallmark"),
  fgsea.c7.cons       |> mutate(method = "cons", database = "C7"),
  fgsea.gobp.cons     |> mutate(method = "cons", database = "GOBP"),
  fgsea.bg3m.cons     |> mutate(method = "cons", database = "BG3M")
) |>
  filter(padj <= 0.05) |>
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) |>
  arrange(desc(NES)) |>
  select(method, database, cluster, pathway, pval, padj, NES, size, log2err, leadingEdge)

write.table(fgsea_all_sig,"Results/Tables/02_04_fgsea_all_methods_sig.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## tabla resumen de tendencias

fgsea_trends <- bind_rows(
  fgsea.hallmark.all  |> mutate(method = "all",  database = "Hallmark"),
  fgsea.c7.all        |> mutate(method = "all",  database = "C7"),
  fgsea.gobp.all      |> mutate(method = "all",  database = "GOBP"),
  fgsea.bg3m.all      |> mutate(method = "all",  database = "BG3M"),
  
  fgsea.hallmark.filt |> mutate(method = "filt", database = "Hallmark"),
  fgsea.c7.filt       |> mutate(method = "filt", database = "C7"),
  fgsea.gobp.filt     |> mutate(method = "filt", database = "GOBP"),
  fgsea.bg3m.filt     |> mutate(method = "filt", database = "BG3M"),
  
  fgsea.hallmark.cons |> mutate(method = "cons", database = "Hallmark"),
  fgsea.c7.cons       |> mutate(method = "cons", database = "C7"),
  fgsea.gobp.cons     |> mutate(method = "cons", database = "GOBP"),
  fgsea.bg3m.cons     |> mutate(method = "cons", database = "BG3M")
) |>
  filter(padj > 0.05, padj <= 0.2) |>
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) |>
  arrange(desc(NES)) |>
  select(method, database, cluster, pathway, pval, padj, NES, size, log2err, leadingEdge)

write.table(fgsea_trends, "Results/Tables/02_04_fgsea_all_methods_trends.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## ---- 3.2 figuras GSEA agrupadas de términos top (grupo significativo)

res.sig <- fgsea_all_sig |> filter(database %in% c("Hallmark", "GOBP")) |> split(~ cluster)
ranks_map <- list(all = ranks.all, filt = ranks.filt, cons = ranks.cons)

plot_fgsea_top <- function(res.sig, ranks_map) {
  for (cl in names(res.sig)) {
    df_cl <- res.sig[[cl]]
    methods <- unique(df_cl$method)
    for (m in methods) {
      fgseaRes <- df_cl |> filter(method == m)
      if (nrow(fgseaRes) == 0) next
      top_up <- fgseaRes |> filter(NES > 0) |> arrange(desc(NES)) |> slice_head(n = 10) |> pull(pathway)
      top_down <- fgseaRes |> filter(NES < 0) |> arrange(NES) |> slice_head(n = 10) |> pull(pathway)
      fgsea.top <- c(top_up, rev(top_down))
      if (length(fgsea.top) == 0) next
      p_fgsea.top <- plotGseaTable(pathways = c(hallmark, gobp)[fgsea.top], ranks_map[[m]][[cl]], fgseaRes, gseaParam = 0.5)
      ggsave(filename = paste0("Results/Figures/02_04_clust", cl, "_fgsea_", m, ".png"), p_fgsea.top, dpi = 300)
    }
  }
}

plot_fgsea_top(res.sig, ranks_map)

## ---- 3.3 figuras GSEA agrupadas de términos top (grupo de tendencias)
res.tend <- fgsea_trends |> filter(database %in% c("Hallmark", "GOBP")) |> split(~ cluster)

plot_fgsea_top_tend <- function(res.sig, ranks_map) {
  for (cl in names(res.sig)) {
    df_cl <- res.sig[[cl]]
    methods <- unique(df_cl$method)
    for (m in methods) {
      fgseaRes <- df_cl |> filter(method == m)
      if (nrow(fgseaRes) == 0) next
      top_up <- fgseaRes |> filter(NES > 0) |> arrange(desc(NES)) |> slice_head(n = 10) |> pull(pathway)
      top_down <- fgseaRes |> filter(NES < 0) |> arrange(NES) |> slice_head(n = 10) |> pull(pathway)
      fgsea.top <- c(top_up, rev(top_down))
      if (length(fgsea.top) == 0) next
      p_fgsea.top <- plotGseaTable(pathways = c(hallmark, gobp)[fgsea.top], ranks_map[[m]][[cl]], fgseaRes, gseaParam = 0.5)
      ggsave(filename = paste0("Results/Figures/02_04_clust", cl, "_fgsea_", m, "_trend.png"), p_fgsea.top, dpi = 300)
    }
  }
}

plot_fgsea_top_tend(res.tend, ranks_map)


## ---- 3.4 figuras plotenrichment representativas

plotEnrichment(hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], ranks.all[["1"]])



for (cl in as.character(1:11)) {
  p <- plotEnrichment(hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], ranks.all[[cl]]) +
    labs(title = paste0("Respuesta a IFN gamma (Cluster ", cl, ")"))
  ggsave(paste0("Results/Figures/02_04_clust", cl, "_IFNg.png"), p, dpi = 300)
}

## ---- 3.4 Dotplot de términos enriquecidos (a partir de fgsea) ----

library(dplyr)
library(ggplot2)
library(forcats)

# elegir subconjunto
dot.df <- fgsea_all_sig |> filter(method == "all", database %in% c("Hallmark", "GOBP")) |>
  mutate(cluster = factor(cluster, levels = sort(unique(as.numeric(as.character(cluster))))))

# selección de términos

terms.keep <- c(
  # inflamación / IFN / TNF
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_COMPLEMENT",
  "GOBP_INFLAMMATORY_RESPONSE",
  "GOBP_PYROPTOTIC_INFLAMMATORY_RESPONSE",
  "GOBP_INTERLEUKIN_10_PRODUCTION",
  "GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION",
  "GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
  
  # activación T / respuesta inmune efectora
  "HALLMARK_IL2_STAT5_SIGNALING",
  "GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
  "GOBP_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
  "GOBP_POSITIVE_REGULATION_OF_CELL_ACTIVATION",
  "GOBP_T_CELL_ACTIVATION",
  "GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_ADAPTIVE_IMMUNE_RESPONSE",
  "GOBP_IMMUNE_EFFECTOR_PROCESS",
  "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY",
  "GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY",
  
  # presentación antigénica
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_VIA_MHC_CLASS_IB",
  
  # adhesión / interacción / tráfico vesicular
  "GOBP_CELL_CELL_ADHESION",
  "GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION",
  "GOBP_POSITIVE_REGULATION_OF_EXOCYTOSIS",
  "GOBP_REGULATION_OF_ENDOCYTOSIS",
  "GOBP_POSITIVE_REGULATION_OF_PROTEOLYSIS",
  
  # diferenciación / memoria / crecimiento
  "GOBP_POSITIVE_REGULATION_OF_CELL_DIFFERENTIATION",
  "GOBP_POSITIVE_REGULATION_OF_MEMORY_T_CELL_DIFFERENTIATION",
  "GOBP_POSITIVE_REGULATION_OF_RAS_PROTEIN_SIGNAL_TRANSDUCTION",
  
  # traducción / ribosoma / biosíntesis
  "HALLMARK_HEME_METABOLISM",
  "HALLMARK_MYC_TARGETS_V1",
  "GOBP_CYTOPLASMIC_TRANSLATION",
  "GOBP_RIBOSOME_ASSEMBLY",
  
  # sincitio / fusión
  "GOBP_NEGATIVE_REGULATION_OF_SYNCYTIUM_FORMATION_BY_PLASMA_MEMBRANE_FUSION",
  
  # proliferación / ciclo celular general
  "GOBP_POSITIVE_REGULATION_OF_T_CELL_PROLIFERATION",
  "GOBP_CELL_CYCLE",
  "HALLMARK_E2F_TARGETS",
  "GOBP_DNA_METABOLIC_PROCESS",
  
  # división / mitosis
  "HALLMARK_G2M_CHECKPOINT",
  "GOBP_MITOTIC_CELL_CYCLE",
  "GOBP_MITOTIC_SISTER_CHROMATID_SEPARATION",
  
  # apoptosis / muerte / supervivencia
  "GOBP_NEGATIVE_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY",
  "GOBP_NEGATIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
  "GOBP_NEGATIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_BY_P53_CLASS_MEDIATOR",
  "GOBP_NEGATIVE_REGULATION_OF_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR"
)

dot.sel <- dot.df |> filter(pathway %in% terms.keep) |> distinct(cluster, pathway, .keep_all = TRUE)
path_order <- terms.keep[terms.keep %in% unique(dot.sel$pathway)]
dot.sel <- dot.sel |> mutate(pathway = factor(pathway, levels = path_order),
    pathway_label = pathway |> as.character() |> gsub("^HALLMARK_", "", x = _) |>
      gsub("^GOBP_", "", x = _) |> gsub("_", " ", x = _))

label_order <- dot.sel |> distinct(pathway, pathway_label) |> arrange(pathway) |>
  pull(pathway_label)

dot.sel <- dot.sel |> mutate(pathway_label = factor(pathway_label, levels = label_order))
dot.sel <- dot.sel |> mutate(NES_plot = pmax(pmin(NES, 2), -2), padj_plot = -log10(padj))

p_dot <- ggplot(dot.sel, aes(x = cluster, y = pathway_label)) +
  geom_point(aes(size = padj_plot, color = NES_plot)) +
  scale_size_continuous(name = expression(-log[10](FDR)), range = c(2, 8)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES") +
  labs(x = "Cluster", y = NULL) +
  theme_bw() + theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 10))

ggsave("Results/Figures/02_04_dotplot_fgsea_all.png", p_dot, width = 12, height = 16, dpi = 300)



## =================================================
## BLOQUE 4 - ANOTACIÓN MANUAL
## =================================================
library(Seurat)
library(dplyr)
library(ggplot2)

ann <- readRDS("Proc/03_02_int_ann.rds")

manual_ann <- c(
  "1"  = "01 - Activación II (polarización Th1)",
  "2"  = "02 - Activación I (reprogramación metabólica)",
  "3"  = "03 - Activación (replicación viral)",
  "4"  = "04 - Efector CTL",
  "5"  = "05 - Transición efector-memoria I (CTL)",
  "6"  = "06 - Transición efector-memoria II (polarización Treg/Th1)",
  "7"  = "07 - Treg memoria Th1-like",
  "8"  = "08 - Proliferación G2/M",
  "9"  = "09 - Proliferación S",
  "10" = "10 - Treg memoria Th1-like / Polarización Th2",
  "11" = "11 - Memoria efectora Th2"
)

# añadir anotación al objeto
Idents(ann) <- "seurat_clusters"
ann$functional_ann <- unname(manual_ann[as.character(ann$seurat_clusters)])

p_functional_ann <- DimPlot(ann, reduction = "umap", group.by = "functional_ann")
ggsave("Results/Figures/03b_04_plot_functional_ann.png", p_functional_ann, width = 10, height = 6, dpi = 300)


# versión 2
manual_ann2 <- c(
  "1"  = "Act (Th1)",
  "2"  = "Act (reprog)",
  "3"  = "Act (replic)",
  "4"  = "Tem/Teff CTL",
  "5"  = "Tem CTL",
  "6"  = "Tem Th1/Treg",
  "7"  = "mTreg Th1-like",
  "8"  = "Prolif G2/M",
  "9"  = "Prolif S",
  "10" = "mTreg Th1/Th2-like",
  "11" = "Tem Th2"
)

# añadir anotación al objeto
Idents(ann) <- "seurat_clusters"
ann$functional_ann2 <- unname(manual_ann2[as.character(ann$seurat_clusters)])

p_functional_ann2 <- DimPlot(ann, reduction = "umap", group.by = "functional_ann2")
ggsave("Results/Figures/03b_04_plot_functional_ann2.png", p_functional_ann2, width = 8, height = 6, dpi = 300)

saveRDS(ann, "Proc/03b_04_int_ann.rds")


## ========================================================
## BLOQUE 5 - EXPORTACIÓN DE DATOS PARA pySCENIC 
## ========================================================
library(Seurat)
library(SCopeLoomR)

obj <- readRDS("Proc/01_data_prep_and_QC.rds")
ann <- readRDS("Proc/03b_04_int_ann.rds")

# obtener archivo loom
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
counts <- counts[Matrix::rowSums(counts > 0) >= 3, ]

loom <- build_loom(
  file.name = "SCENIC/02_postprocessed_inputs/ann.loom",
  dgem = counts,
  title = "seurat input for notebook 2",
  genome = "Human",
  default.embedding = as.data.frame(ann@reductions$umap@cell.embeddings[colnames(counts), 1:2]),
  default.embedding.name = "Seurat_UMAP")

# añadir metadatos de anotación y clusters seurat
for (col in colnames(ann@meta.data)) {
  values <- ann@meta.data[colnames(counts), col]
  
  if (is.factor(values)) {values <- as.character(values)}
  
  n_unique <- length(unique(values[!is.na(values)]))
  use_as_annotation <- !is.numeric(values) && n_unique <= 245
  
  add_col_attr(loom = loom, key = col, value = values,
               as.annotation = use_as_annotation)
}


close_loom(loom)

