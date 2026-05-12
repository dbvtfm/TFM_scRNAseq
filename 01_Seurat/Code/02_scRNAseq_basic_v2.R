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
# Script: 02_scRNAseq_basic_v2.R
## ===============================

# Objetivo:
# Aplicar SCTransform v2 sobre el objeto Seurat filtrado, integración mediante RPCA,
# reducción de dimensionalidad, clustering Leiden, genes marcadores globales y 
# conservados por cluster mediante MAST, heatmaps y UpSet plots rankeados con
# score propio y caracterización preliminar de clusters

# Inputs en Proc/:
# - 01_data_prep_and_QC.rds

# Outputs en Proc/:
# - 02_01_sct.rds
# - 02_02_sct_int-rpca.rds
# - 02_02_sct_int-rpca_leiden_1.rds
# - 02_03_sct_prep_findmarkers.rds

# Outputs en Results/:
# - Figures/02_*.png
# - Tables/02_*.tsv


## ==========================
## BLOQUE 1 - SCTransform v2
## ==========================
obj <- readRDS("Proc/01_data_prep_and_QC.rds")
hvg <- obj

library(Seurat)
library(glmGamPoi)
library(patchwork)
library(ggplot2)

# split assay RNA en layers por muestra
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Sample)

# SCTransform
obj <- SCTransform(obj, vst.flavor = "v2", return.only.var.genes = F, 
                   vars.to.regress = "percent.mt", verbose = F)
saveRDS(obj, "Proc/02_01_sct.rds")

# gráfico ajuste varianza residuos SCTransform
library(sctransform)

obj_tmp <- JoinLayers(obj, assay = "RNA")
umi <- GetAssayData(obj_tmp, assay = "RNA", layer = "counts")
vst_out <- vst(umi, vst.flavor = "v2", return_cell_attr = T, return_gene_attr = T, verbosity = 0)

p_model_pars <- plot_model_pars(vst_out, show_theta = T, show_var = T)
p_model_gene <- plot_model(vst_out, umi, goi = c("MALAT1", "IL32", "STAT3"), show_fit = T, show_nr = T, plot_residual = T)

ggsave("Results/Figures/02_01_SCT_model_parameters.png", p_model_pars, width = 8, height = 6, dpi = 300)
ggsave("Results/Figures/02_01_SCT_model_genes.png", p_model_gene, width = 10, height = 10, dpi = 300)

# HVG
top10 <- head(VariableFeatures(obj), 10)
top10

# Figura HVG
DefaultAssay(hvg) <- "RNA"
hvg_list = SplitObject(hvg, split.by="Sample")
hvg_list <- lapply(X = hvg_list, FUN = SCTransform, vst.flavor = "v2", return.only.var.genes = F, vars.to.regress = "percent.mt")
hvg_list <- lapply(X = hvg_list, FUN = FindVariableFeatures, selection.method = "vst", nfeatures = 3000, assay="SCT")
top10_1 <- head(VariableFeatures(hvg_list[[5]]), 10)
top10_1
p_hvg <- VariableFeaturePlot(hvg_list[[5]])
p_hvg <- LabelPoints(plot = p_hvg, points = top10_1, repel = T)
ggsave("Results/Figures/02_01_SCT_HVG.png", p_hvg, width = 8, height = 4, dpi = 300) 


## ========================================================================
## BLOQUE 2 - Integración, reducción de dimensionalidad y clustering
## ========================================================================

## ---- 2.1 Integración (RPCA) ----
set.seed(123)
obj <- RunPCA(obj, npcs = 300, verbose = FALSE)
int <- IntegrateLayers(object = obj, method = RPCAIntegration, normalization.method = "SCT",
                            orig.reduction = "pca", new.reduction = "integrated.rpca")
saveRDS(int, "Proc/02_02_sct_int-rpca.rds")

## ---- 2.2 Elbow Plot ----
p_elbow <- ElbowPlot(int, ndims = 300, reduction = "pca")
ggsave("Results/Figures/02_02_elbow_rpca.png", p_elbow, width = 8, height = 5, dpi = 300)

## ---- 2.3 Reducción de dimensionalidad y clustering ----
set.seed(123)
int <- RunUMAP(int, reduction = "integrated.rpca", dims = 1:30, verbose = FALSE)
int <- FindNeighbors(int, reduction = "integrated.rpca", dims = 1:30, verbose = FALSE)
int <- FindClusters(int, resolution = 1, algorithm = 4, verbose = FALSE)

saveRDS(int, "Proc/02_02_sct_int-rpca_leiden_1.rds")

## ---- 2.4 UMAP ----
p_clust <- DimPlot(int, reduction = "umap", label = TRUE)
p_samp <- DimPlot(int, reduction = "umap", group.by = "Sample")
p_treat <- DimPlot(int, reduction = "umap", group.by = "Treatment")
p_rep <- DimPlot(int, reduction = "umap", group.by = "Replicate")
p_umap <- (p_clust | p_samp) / (p_treat | p_rep)

#ggsave("Results/Figures/02_02_UMAP_clusters.png", p_clust, width = 7, height = 5, dpi = 300)
#ggsave("Results/Figures/02_02_UMAP_sample.png", p_samp, width = 7, height = 5, dpi = 300)
#ggsave("Results/Figures/02_02_UMAP_treatment.png", p_treat, width = 7, height = 5, dpi = 300)
#ggsave("Results/Figures/02_02_UMAP_replicate.png", p_rep, width = 7, height = 5, dpi = 300)
ggsave("Results/Figures/02_02_UMAP_int-rpca_leiden-1.png", p_umap, width = 12, height = 10, dpi = 300)


## ---- 2.5 Figuras y tablas resumen del clustering ----
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

meta_clust <- int[[]] |> tibble::rownames_to_column("cell") |>
  select(cell, seurat_clusters, Sample, Treatment, Replicate) |>rename(cluster = seurat_clusters)

# nº células
tab_ncells_cluster <- meta_clust |> count(cluster, name = "n_cells") |> arrange(as.numeric(as.character(cluster))) 
p_ncells_cluster <- ggplot(tab_ncells_cluster, aes(x = factor(cluster), y = n_cells)) +
  geom_col(width = 0.8) + labs(x = "Cluster", y = "Número de células") + theme_bw()
#ggsave("Results/Figures/02_02b_ncells_by_cluster.png", p_ncells_cluster, width = 7, height = 5, dpi = 300)

# % tratamiento
tab_cluster_treatment <- meta_clust |> count(cluster, Treatment, name = "n_cells") |> 
  group_by(cluster) |> mutate(prop = n_cells / sum(n_cells)) |> ungroup()
p_cluster_treatment <- ggplot(tab_cluster_treatment, aes(x = factor(cluster), y = prop, fill = Treatment)) +
  geom_col(width = 0.8) + scale_y_continuous(labels = percent_format()) + labs(x = "Cluster") + theme_bw()
#ggsave("Results/Figures/02_02b_cluster_by_treat.png", p_cluster_treatment, width = 7, height = 5, dpi = 300)

# % donante (réplica)
tab_cluster_replicate <- meta_clust |> count(cluster, Replicate, name = "n_cells") |> group_by(cluster) |>
  mutate(prop = n_cells / sum(n_cells)) |> ungroup()
p_cluster_replicate <- ggplot(tab_cluster_replicate, aes(x = factor(cluster), y = prop, fill = Replicate)) +
  geom_col(width = 0.8) + scale_y_continuous(labels = percent_format()) + labs(x = "Cluster") + theme_bw()
#ggsave("Results/Figures/02_02b_cluster_by_rep.png", p_cluster_replicate, width = 7, height = 5, dpi = 300)

# % muestra
tab_cluster_sample <- meta_clust |> count(cluster, Sample, name = "n_cells") |> group_by(cluster) |>
  mutate(prop = n_cells / sum(n_cells)) |> ungroup()
p_cluster_sample <- ggplot(tab_cluster_sample, aes(x = factor(cluster), y = prop, fill = Sample)) +
  geom_col(width = 0.8) + scale_y_continuous(labels = percent_format()) + labs(x = "Cluster") + theme_bw()
#ggsave("Results/Figures/02_02b_cluster_by_samp.png", p_cluster_sample, width = 7, height = 5, dpi = 300)

# patchwork
p_summary <- (p_ncells_cluster | p_cluster_treatment) / (p_cluster_replicate | p_cluster_sample)
ggsave("Results/Figures/02_02b_clusters_summary.png", p_summary, width = 12, height = 10, dpi = 300)

# tablas
write.table(tab_ncells_cluster, "Results/Tables/02_02b_ncells_by_cluster.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tab_cluster_treatment, "Results/Tables/02_02b_cluster_by_treatment.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tab_cluster_replicate, "Results/Tables/02_02b_cluster_by_replicate.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tab_cluster_sample, "Results/Tables/02_02b_cluster_by_sample.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


## ======================================================
## BLOQUE 3 - IDENTIFICACIÓN DE GENES MARCADORES (MAST)
## ======================================================
library(Seurat)
library(dplyr)
library(ggplot2)

# Preparación del assay SCT para DE
DefaultAssay(int) <- "SCT"
prep <- PrepSCTFindMarkers(int, assay = "SCT")
Idents(prep) <- "seurat_clusters"
saveRDS(prep, "Proc/02_03_sct_prep_findmarkers.rds")
#prep <- readRDS("Proc/02_03_sct_prep_findmarkers.rds")

# Find all markers
markers.all <- FindAllMarkers(prep, assay = "SCT", test.use = "MAST")
write.table(markers.all,"Results/Tables/02_03_FindAllMarkers_clust_MAST.tsv", sep = "\t", row.names = F, quote = F)
#markers.all <- read.delim("Results/Tables/02_03_FindAllMarkers_clust_MAST.tsv")

# Find conserved markers
clusters <- levels(Idents(prep))
markers.cons <- lapply(clusters, function(cl) {
  df <- FindConservedMarkers(prep, ident.1 = cl, assay = "SCT", test.use = "MAST", grouping.var = "Treatment")
  df$gene <- rownames(df)
  df$cluster <- cl
  df
})
names(markers.cons) <- clusters
rm(clusters)
markers.cons <- bind_rows(markers.cons)
write.table(markers.cons,"Results/Tables/02_03_FindConservedMarkers_clust_MAST.tsv", sep = "\t", row.names = F, quote = F)
#markers.cons <- read.delim("Results/Tables/02_03_FindConservedMarkers_clust_MAST.tsv")

## ---- 3.1 HM de top10 de FindAllMarkers con filtrado no significativos y ranking con score propio ----

# Eliminar all markers con p-valor > 0.05 y calcular score de rankeado
markers.all.rank <- markers.all |> filter(p_val_adj <= 0.05)
markers.all.rank <- markers.all.rank |> mutate(
  score = case_when(avg_log2FC > 0 ~ avg_log2FC * pct.1 * (pct.1 / (pct.1 + pct.2)),
                    avg_log2FC < 0 ~ avg_log2FC * pct.2 * (pct.2 / (pct.1 + pct.2)),
                    TRUE ~ 0))

# top 10 markers por cluster (up y down) y 50 up
top10up.all <- markers.all.rank |> filter(score > 0) |> group_by(cluster) |>
  arrange(desc(score), p_val_adj, .by_group = T) |> slice_head(n = 10) |> ungroup()

top10down.all <- markers.all.rank |> filter(score < 0) |> group_by(cluster) |>
  arrange(score, p_val_adj, .by_group = T) |> slice_head(n = 10) |> ungroup()

top50up.all <- markers.all.rank |> filter(score > 0) |> group_by(cluster) |>
  arrange(desc(score), p_val_adj, .by_group = T) |> slice_head(n = 50) |> ungroup()

# heatmaps top 10 up y down
prep <- GetResidual(prep, features = top10up.all$gene)
hm_top10up.all <- DoHeatmap(prep, features = top10up.all$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  #+ NoLegend()
ggsave("Results/Figures/02_03_FindAllMarkers_top10up.png", hm_top10up.all, width = 9, height = 12)

prep <- GetResidual(prep, features = top10down.all$gene)
hm_top10down.all <- DoHeatmap(prep, features = top10down.all$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  #+ NoLegend()
ggsave("Results/Figures/02_03_FindAllMarkers_top10down.png", hm_top10down.all, width = 9, height = 12)

prep <- GetResidual(prep, features = top50up.all$gene)
hm_top50up.all <- DoHeatmap(prep, features = top50up.all$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  #+ NoLegend()
ggsave("Results/Figures/02_03_FindAllMarkers_top50up.png", hm_top50up.all, width = 9, height = 12)


## ---- 3.2 HM de top10 de FindConservedMarkers, filtrado no significativos y ranking con score propio ----

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

# top 10 markers conservados por cluster (up y down) y 50 up
top10up.cons <- markers.cons.rank |> filter(score > 0) |> group_by(cluster) |>
  arrange(desc(score), max_pval, .by_group = TRUE) |> slice_head(n = 10) |> ungroup()

top10down.cons <- markers.cons.rank |> filter(score < 0) |> group_by(cluster) |>
  arrange(score, max_pval, .by_group = TRUE) |> slice_head(n = 10) |> ungroup()

top50up.cons <- markers.cons.rank |> filter(score > 0) |> group_by(cluster) |>
  arrange(desc(score), max_pval, .by_group = TRUE) |> slice_head(n = 50) |> ungroup()

# heatmaps top 10 up y down
prep <- GetResidual(prep, features = top10up.cons$gene)
hm_top10up.cons <- DoHeatmap(prep, features = top10up.cons$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  #+ NoLegend()
ggsave("Results/Figures/02_03_FindConservedMarkers_top10up.png", hm_top10up.cons, width = 9, height = 12)

prep <- GetResidual(prep, features = top10down.cons$gene)
hm_top10down.cons <- DoHeatmap(prep, features = top10down.cons$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  #+ NoLegend()
ggsave("Results/Figures/02_03_FindConservedMarkers_top10down.png", hm_top10down.cons, width = 9, height = 12)

prep <- GetResidual(prep, features = top50up.cons$gene)
hm_top50up.cons <- DoHeatmap(prep, features = top50up.cons$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  #+ NoLegend()
ggsave("Results/Figures/02_03_FindConservedMarkers_top50up.png", hm_top50up.cons, width = 9, height = 12)


## ---- 3.3 HM de top10 de FindAllMarkers, filtrado no significativos y no conservados significativos y ranking con score propio ----

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

# top 10 markers filtrados por cluster (up y down) y 50 up
top10up.filt <- markers.all.filt |> filter(score > 0) |> group_by(cluster) |>
  arrange(desc(score), p_val_adj, .by_group = TRUE) |> slice_head(n = 10) |> ungroup()

top10down.filt <- markers.all.filt |> filter(score < 0) |> group_by(cluster) |>
  arrange(score, p_val_adj, .by_group = TRUE) |> slice_head(n = 10) |> ungroup()

top50up.filt <- markers.all.filt |> filter(score > 0) |> group_by(cluster) |>
  arrange(desc(score), p_val_adj, .by_group = TRUE) |> slice_head(n = 50) |> ungroup()

# heatmaps top 10 up y down
prep <- GetResidual(prep, features = top10up.filt$gene)
hm_top10up.filt <- DoHeatmap(prep, features = top10up.filt$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  #+ NoLegend()
ggsave("Results/Figures/02_03_FindAllMarkers_filt_top10up.png", hm_top10up.filt, width = 9, height = 12)

prep <- GetResidual(prep, features = top10down.filt$gene)
hm_top10down.filt <- DoHeatmap(prep, features = top10down.filt$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  #+ NoLegend()
ggsave("Results/Figures/02_03_FindAllMarkers_filt_top10down.png", hm_top10down.filt, width = 9, height = 12)

prep <- GetResidual(prep, features = top50up.filt$gene)
hm_top50up.filt <- DoHeatmap(prep, features = top50up.filt$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  #+ NoLegend()
ggsave("Results/Figures/02_03_FindAllMarkers_filt_top50up.png", hm_top50up.filt, width = 9, height = 12)


## ---- 3.4 Tabla resumen de top 10 genes DE por cluster según cada método ----

# tabla resumen de 10 genes up y down más diferenciales 
tab_top10 <- bind_rows(
  top10up.all   |> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop") |> mutate(method = "FindAllMarkers", direction = "up"),
  top10down.all |> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop") |> mutate(method = "FindAllMarkers", direction = "down"),
  top10up.filt  |> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop") |> mutate(method = "FindAllMarkers_filt", direction = "up"),
  top10down.filt|> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop") |> mutate(method = "FindAllMarkers_filt", direction = "down"),
  top10up.cons  |> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop") |> mutate(method = "FindConservedMarkers", direction = "up"),
  top10down.cons|> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop") |> mutate(method = "FindConservedMarkers", direction = "down")
) |> 
  select(method, direction, cluster, gene)

write.table(tab_top10, "Results/Tables/02_03_top10_genes.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
#top10 <- read.delim("Results/Tables/02_03_top10_genes.tsv")


## ---- 3.4 UpSet plot del total de genes DE conservados por cluster ----
library(dplyr)
library(tidyr)
library(ComplexUpset)
library(ggplot2)


# UpSet plot de genes conservados up-regulados compartidos entre clusters
sets_all_up.cons <- markers.cons |> filter(Control_p_val_adj <= 0.05, HA15_p_val_adj <= 0.05,
                                           Control_avg_log2FC > 0, HA15_avg_log2FC > 0) |>
  select(cluster, gene) |> distinct() |> mutate(set = paste0("cluster", cluster))

tab_upset_all_up.cons <- sets_all_up.cons |> mutate(value = 1) |> select(set, gene, value) |>
  distinct() |> pivot_wider(names_from = set, values_from = value, values_fill = 0)

set_order_all_up.cons <- paste0("cluster", 1:11)

p_upset_all_up.cons <- upset(tab_upset_all_up.cons, intersect = set_order_all_up.cons,
  name = "Genes", sort_sets = F, sort_intersections = "descending", min_size = 10,
  sort_intersections_by = c("cardinality", "degree"), min_degree = 2, set_sizes = upset_set_size(geom = geom_bar(fill = "red")), height_ratio = 3, width_ratio = 0.12,
  base_annotations = list("Intersect. size" = intersection_size(counts = TRUE,text = list(vjust = -0.3), fill ="red")))
ggsave("Results/Figures/02_03_UpSet_findconservedmarkers_up.png", p_upset_all_up.cons, width = 8, height = 4)


# UpSet plot de genes conservados down-regulados compartidos entre clusters
sets_all_down.cons <- markers.cons |> filter(Control_p_val_adj <= 0.05, HA15_p_val_adj <= 0.05,
                                             Control_avg_log2FC < 0, HA15_avg_log2FC < 0) |>
  select(cluster, gene) |> distinct() |> mutate(set = paste0("cluster", cluster))

tab_upset_all_down.cons <- sets_all_down.cons |> mutate(value = 1) |> select(set, gene, value) |>
  distinct() |> pivot_wider(names_from = set, values_from = value, values_fill = 0)

set_order_all_down.cons <- paste0("cluster", 1:11)

p_upset_all_down.cons <- upset(tab_upset_all_down.cons, intersect = set_order_all_down.cons,
  name = "Genes", sort_sets = F, sort_intersections = "descending", min_size = 8,
  sort_intersections_by = c("cardinality", "degree"), min_degree = 2, set_sizes = upset_set_size(geom = geom_bar(fill = "blue")), height_ratio = 3, width_ratio = 0.12,
  base_annotations = list("Intersect. \n size" = intersection_size(counts = TRUE,text = list(vjust = -0.3), fill ="blue")))
ggsave("Results/Figures/02_03_UpSet_findconservedmarkers_down.png", p_upset_all_down.cons, width = 8, height = 4)

# opuestos
sets_all_opp.cons <- markers.cons |>
  filter(Control_p_val_adj <= 0.05, HA15_p_val_adj <= 0.05) |>
  mutate(direction = case_when(
    Control_avg_log2FC > 0 & HA15_avg_log2FC > 0 ~ "up",
    Control_avg_log2FC < 0 & HA15_avg_log2FC < 0 ~ "down",
    TRUE ~ NA_character_
  )) |>
  filter(!is.na(direction)) |>
  select(cluster, gene, direction) |> distinct() |>
  mutate(set = paste0("cluster", cluster, "_", direction))

tab_upset_all_opp.cons <- sets_all_opp.cons |> mutate(value = 1) |> select(set, gene, value) |>
  distinct() |> pivot_wider(names_from = set, values_from = value, values_fill = 0) |>
  filter(rowSums(across(ends_with("_up"))) > 0 & rowSums(across(ends_with("_down"))) > 0)

set_order_all_opp.cons <- c(rbind(paste0("cluster", 1:11, "_up"), paste0("cluster", 1:11, "_down")))

p_upset_all_opp.cons <- upset(tab_upset_all_opp.cons, intersect = set_order_all_opp.cons,
                              name = "Genes", sort_sets = F, sort_intersections = "descending", min_size = 10,
                              sort_intersections_by = c("cardinality", "degree"), min_degree = 2, set_sizes = upset_set_size(),
                              height_ratio = 3, width_ratio = 0.12,
                              stripes = upset_stripes(colors = c("salmon","lightblue","salmon","lightblue","lightblue", "salmon", "lightblue", "salmon", "lightblue", "salmon", "salmon", "lightblue", "salmon", "lightblue")),
                              base_annotations = list("Intersect. size" = intersection_size(counts = TRUE, text = list(vjust = -0.3))))
ggsave("Results/Figures/02_03_UpSet_findconservedmarkers_opposite.png", p_upset_all_opp.cons, width = 8, height = 4)


## =============================================================================
## BLOQUE 4 - CARACTERIZACIÓN PRELIMINAR DE CLUSTERS Y POBLACIONES INMUNITARIAS
## =============================================================================
library(Seurat)
library(ggplot2)

int <- readRDS("Proc/02_02_sct_int-rpca_leiden_1.rds")
DefaultAssay(int) <- "SCT"

# selección negativa de CD8, CD14, CD16, CD19, CD20, CD36, CD56, CD66b, CD123, GlyA/CD235a, and TCRgd
MACS = c("CD8A", "CD8B","CD14", "FCGR3A", "FCGR3B", "CD19", "MS4A1", "CD36", "NCAM1",
         "CEACAM8", "IL3RA", "GYPA", "TRDC", "TRGC1", "TRGC2")
MACS_SCT <- intersect(MACS, rownames(int[["SCT"]]))

p_MACS = DotPlot(int, features = MACS_SCT, assay = "SCT", scale = F, scale.max = 100, scale.min = 0) + RotatedAxis()
ggsave("Results/Figures/02_04_dotplot_MACS.png", p_MACS, width = 6, height = 4)

# genes altamente específicos de poblaciones hematopoyéticas
panHemato = "PTPRC"
panT = c("CD3D", "CD3E", "CD3G")
Tab = c("TRAC", "TRBC1", "TRBC2", "CD4", "CD8A", "CD8B")
Tgd = c("TRDC", "TRGC1", "TRGC2")
B = c("MS4A1", "CD19", "CD79A", "CD79B", "JCHAIN", "MZB1", "SDC1", "IGKC")
Mielo_NK = c("LYZ", "FCN1", "S100A8", "S100A9", "CD14", "CEACAM8", "IL3RA",
             "TPSAB1", "TPSB2", "KIT", "CPA3", "FCGR3B","NCAM1", "FCGR3A")
eri_plaq = c("HBB", "HBA1", "HBA2", "ALAS2", "GYPA","PPBP", "PF4", "NRGN")

qc_CD4_SCT <- intersect(c(panHemato, panT, Tab, Tgd, B, Mielo_NK, eri_plaq), rownames(int[["SCT"]]))

p_qc_cd4 = DotPlot(int, features = qc_CD4_SCT, assay = "SCT", scale = F, scale.max = 100) + RotatedAxis()
ggsave("Results/Figures/02_04_dotplot_qc_cd4.png", p_qc_cd4, width = 10, height = 4)

qc_CD4_SCT_v2 <- intersect(c(panHemato, panT, Tab, Tgd), rownames(int[["SCT"]]))

p_qc_cd4_v2 = DotPlot(int, features = qc_CD4_SCT_v2, assay = "SCT", scale = F, scale.max = 100) + RotatedAxis()
ggsave("Results/Figures/02_04_dotplot_qc_cd4_v2.png", p_qc_cd4_v2, width = 10, height = 4)

# genes VIH
VIH = c("genomic-mRNA", "GagPol", "Rev", "Tat", "Vif", "Vpr", "Vpu", "EnvdGFPiresThy1.2insert")
VIH_SCT <- intersect(VIH, rownames(int[["SCT"]]))

p_VIH <- FeaturePlot(int, features = VIH_SCT, cols = c("grey","red"), ncol =3, order = T, reduction = "umap")
ggsave("Results/Figures/02_04_featplot_VIH.png", p_VIH, width = 12, height = 4)

p_VIH_2 <- FeaturePlot(int, features = VIH_SCT, cols = c("grey","red"), ncol =3, order = T, reduction = "umap", split.by = "Treatment")
ggsave("Results/Figures/02_04_featplot_VIH_2.png", p_VIH_2, width = 12, height = 8)

p_combined <- DotPlot(int, features = c(qc_CD4_SCT, VIH_SCT), assay = "SCT", cols = c("grey","red"), scale = T, scale.max = 100) + RotatedAxis()
ggsave("Results/Figures/02_04_dotplot_combined_VIH_scale.png", p_combined, width = 10, height = 5)

# particularidades sobre exclusión alélica del TCR-beta y linaje gamma-delta
p_tcr <- FeaturePlot(int, features = c("TRAC", "TRBC1", "TRBC2"), ncol = 3, order = T, reduction = "umap")
ggsave("Results/Figures/02_04_featplot_TCR.png", p_tcr, width = 12, height = 4)

# cuantificaciones

# % CD3
check <- WhichCells(int, expression = CD3D > 0)
length(check) / ncol(int)
check <- WhichCells(int, expression = CD3E > 0)
length(check) / ncol(int)
check <- WhichCells(int, expression = CD3G > 0)
length(check) / ncol(int)
check <- WhichCells(int, expression = (CD3D > 0 | CD3G > 0))
length(check) / ncol(int)
check <- WhichCells(int, expression = (CD3D > 0 & CD3G > 0))
length(check) / ncol(int)
check <- WhichCells(int, expression = (CD3D > 0 & CD3E > 0 & CD3G > 0))
length(check) / ncol(int)
check <- WhichCells(int, expression = (CD3D == 0 & CD3G == 0 & CD3E > 0))
length(check) / ncol(int)

# % alfa-beta
check <- WhichCells(int, expression = TRAC > 0)
length(check) / ncol(int)
check <- WhichCells(int, expression = TRBC1 > 0)
length(check) / ncol(int)
check <- WhichCells(int, expression = TRBC2 > 0)
length(check) / ncol(int)
check <- WhichCells(int, expression = (TRBC1 > 0 & TRBC2 > 0))
length(check) / ncol(int)
check <- WhichCells(int, expression = (TRAC > 0 | TRBC1 > 0 | TRBC2 > 0))
length(check) / ncol(int)
check <- WhichCells(int, expression = (TRAC > 0 & (TRBC1 > 0 | TRBC2 > 0)))
length(check) / ncol(int)
check <- WhichCells(int, expression = (TRAC > 0 & TRBC1 > 0))
length(check) / ncol(int)
check <- WhichCells(int, expression = (TRAC > 0 & TRBC2 > 0))
length(check) / ncol(int)
check <- WhichCells(int, expression = (TRAC > 0 & TRBC1 > 0 & TRBC2 == 0))
length(check) / ncol(int)
check <- WhichCells(int, expression = (TRAC > 0 & TRBC2 > 0 & TRBC1 == 0))
length(check) / ncol(int)
check <- WhichCells(int, expression = (TRAC > 0 & TRBC1 > 0 & TRBC2 > 0))
length(check) / ncol(int)
check <- WhichCells(int, expression = (TRAC > 0 & TRBC2 == 0 & TRBC1 == 0))
length(check) / ncol(int)

# % gamma-delta
check <- WhichCells(int, expression = TRDC > 0)
length(check) / ncol(int) 
check <- WhichCells(int, expression = TRGC1 > 0)
length(check) / ncol(int) 
check <- WhichCells(int, expression = TRGC2 > 0)
length(check) / ncol(int) 
check <- WhichCells(int, expression = (TRDC > 0 & (TRGC1 > 0 | TRGC2 > 0)))
length(check) / ncol(int) 
check <- WhichCells(int, expression = (TRAC > 0 & TRGC2 > 0))
length(check) / ncol(int)
check <- WhichCells(int, expression = (TRBC1 > 0 & TRGC2 > 0))
length(check) / ncol(int)
check <- WhichCells(int, expression = (TRBC2 > 0 & TRGC2 > 0))
length(check) / ncol(int)
check <- WhichCells(int, expression = (TRGC2 > 0 & (TRBC1 > 0 | TRBC2 > 0)))
length(check) / ncol(int)
check <- WhichCells(int, expression = (TRGC2 > 0 & (TRBC1 > 0 | TRBC2 > 0) & TRAC > 0))
length(check) / ncol(int)

check <- WhichCells(int, expression = (CD4 > 0 & TRAC > 0))
length(check) / ncol(int)

#conserved <- read.table("Results/Tables/03_02_FindConservedMarkers_clust.tsv", header = T, sep = "\t")
#cons.genelist <- read.delim("Results/Tables/03_02_FindConservedMarkers_clust_genelist_summary.tsv", sep = "\t", header = T, stringsAsFactors = F)

# featureplots de todos los marcadores up-regulados para exploración visual, seleccionando método

fp_summary <- function(clust_obj, tab, method, direction, genes_per_plot = 36) {
  lista <- unique(trimws(unlist(strsplit(
    paste(tab$gene[tab$method == method & tab$direction == direction], collapse = ","),","))))
  bloques <- split(lista, ceiling(seq_along(lista) / genes_per_plot))
  for (i in seq_along(bloques)) {
    p_all <- FeaturePlot(clust_obj, features = bloques[[i]], reduction = "umap", order = T, ncol = 4)
    ggsave(paste0("results/figures/02_04_markers_", method, "_", direction, "_part", i, ".png"),p_all, width = 14, height = 20, dpi = 300)
  }
}

# featureplots de todos los marcadores up-regulados únicos de entre todos los métodos

fp_summary2 <- function(clust_obj, tab, direction, genes_per_plot = 36) {
  lista <- unique(trimws(unlist(strsplit(
    paste(tab$gene[tab$direction == direction], collapse = ","), ","))))
  bloques <- split(lista, ceiling(seq_along(lista) / genes_per_plot))
  for (i in seq_along(bloques)) {
    p_all <- FeaturePlot(clust_obj, features = bloques[[i]], reduction = "umap", order = TRUE, ncol = 4)
    ggsave(paste0("results/figures/02_04_markers_allmethods_", direction, "_part", i, ".png"),
           p_all, width = 14, height = 20, dpi = 300)
  }
}

top10 <- read.delim("Results/Tables/02_03_top10_genes.tsv")
fp_summary2(int, top10, "up")


# dotplots para un determinado pct
get_list <- function(tab, direction_sel) {
  gene_string <- paste(tab$gene[tab$direction == direction_sel], collapse = ",")
  lista <- unique(trimws(unlist(strsplit(gene_string, ","))))
  return(lista)
}

dp_markers <- DotPlot(int, get_list(top10, "up")) + coord_flip() + 
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) 
ggsave("results/figures/02_04_dotplot_markers_relexp.png", dp_markers, width = 12, height = 16, dpi = 300)

