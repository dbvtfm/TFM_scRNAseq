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
# Script: 02_scRNAseq_basic_v1.R
## ===============================

# Objetivo:
# Aplicar SCTransform v2 sobre el objeto Seurat filtrado, integración mediante RPCA,
# reducción de dimensionalidad, clustering Leiden, genes marcadores globales y 
# conservados por cluster mediante Wilcoxon, heatmaps desglosados por porcentaje
# de células expresoras

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


## ---- 2.6 UMAPs exploratorios ----
DefaultAssay(int) <- "SCT"
# pan-T
plot_T <- FeaturePlot(int, features = c("CD3D", "CD3E", "TRAC"), ncol =3, cols = c("grey","red"), order= T, reduction = "umap")
ggsave("Results/Figures/02_02b_plot_T.png", plot_T, width = 12, height = 4)

# QC PBMCs
plot_qc <- FeaturePlot(int, features = c("LYZ", "FCN1", "MS4A1", "CD79A", "PPBP", "PF4", "HBB", "HBA1"), cols = c("grey","red"), reduction = "umap")
ggsave("Results/Figures/02_02b_plot_QC.png", plot_qc, width = 12, height = 12)

# VIH
plot_VIH <- FeaturePlot(int, features = c("genomic-mRNA", "GagPol", "Rev", "Tat", "Vif", "Vpr", "Vpu", "EnvdGFPiresThy1.2insert"),
                        cols = c("grey","red"), ncol =3, order = T, reduction = "umap")
ggsave("Results/Figures/02_02b_plot_VIH.png", plot_VIH, width = 12, height = 4)

plot_VIH_2 <- FeaturePlot(int, features = c("Rev", "Tat"),
                          cols = c("grey","red"), order = T, reduction = "umap")
ggsave("Results/Figures/02_02b_plot_VIH_2.png", plot_VIH_2, width = 9, height = 4)



## ================================================
## BLOQUE 3 - IDENTIFICACIÓN DE GENES MARCADORES
## ================================================
library(dplyr)
library(ggplot2)

# Preparación del assay SCT para DE
DefaultAssay(int) <- "SCT"
prep <- PrepSCTFindMarkers(int, assay = "SCT")
Idents(prep) <- "seurat_clusters"
saveRDS(prep, "Proc/02_03_sct_prep_findmarkers.rds")
#prep <- readRDS("Proc/02_03_sct_prep_findmarkers.rds")


## ---- 3.1 FINDALLMARKERS: MARCADORES GLOBALES POR CLUSTER (one-vs-rest) ----
markers.all <- FindAllMarkers(prep, assay = "SCT", logfc.threshold = 0.25, min.pct = 0.10, test.use = "wilcox")
write.table(markers.all,"Results/Tables/02_03_FindAllMarkers_clust.tsv", sep = "\t", row.names = F, quote = F)

markers.all_pos <- FindAllMarkers(prep, assay = "SCT", logfc.threshold = 0.25, min.pct = 0.25, test.use = "wilcox", only.pos = T)
write.table(markers.all_pos,"Results/Tables/02_03_FindAllMarkers_pos_clust.tsv", sep = "\t", row.names = F, quote = F)

# Top 50 (todos)
top50_markers.all <- markers.all |> filter(p_val_adj <= 0.05) |> group_by(cluster) |>
  arrange(desc(abs(avg_log2FC)), p_val_adj, .by_group = T) |> slice_head(n = 50) |> ungroup()
write.table(top50_markers.all,"Results/Tables/02_03_FindAllMarkers_clust_top50_all.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top50_markers.all$gene)
hm_top50_markers.all <- DoHeatmap(prep, features = top50_markers.all$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  + NoLegend()
ggsave("Results/Figures/02_03_FindAllMarkers_clust_top50_hm.png", hm_top50_markers.all, width = 9, height = 12)

# Top 10 (up/down irrespective) - pct.1 >= 0.9
top10_09_markers.all <- markers.all |> filter(p_val_adj <= 0.05, pct.1 >= 0.9) |> group_by(cluster) |>
  arrange(desc(abs(avg_log2FC)), p_val_adj, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10_09_markers.all,"Results/Tables/02_03_FindAllMarkers_clust_top10_09_all.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10_09_markers.all$gene)
hm_top10_09_markers.all <- DoHeatmap(prep, features = top10_09_markers.all$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  + NoLegend()
ggsave("Results/Figures/02_03_FindAllMarkers_clust_top10_09_hm.png", hm_top10_09_markers.all, width = 9, height = 12)

# Top 10 (up/down irrespective) - pct.1 >= 0.1
top10_01_markers.all <- markers.all |> filter(p_val_adj <= 0.05, pct.1 >= 0.1) |> group_by(cluster) |>
  arrange(desc(abs(avg_log2FC)), p_val_adj, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10_01_markers.all,"Results/Tables/02_03_FindAllMarkers_clust_top10_01_all.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10_01_markers.all$gene)
hm_top10_01_markers.all <- DoHeatmap(prep, features = top10_01_markers.all$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  + NoLegend()
ggsave("Results/Figures/02_03_FindAllMarkers_clust_top10_01_hm.png", hm_top10_01_markers.all, width = 9, height = 12)

# Top 10 (up) - pct.1 >= 0.9
top10up_09_markers.all <- markers.all |> filter(p_val_adj <= 0.05, avg_log2FC > 0, pct.1 >= 0.9) |> group_by(cluster) |>
  arrange(desc(avg_log2FC), p_val_adj, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10up_09_markers.all,"Results/Tables/02_03_FindAllMarkers_clust_top10_09_up.tsv", sep = "\t", row.names = F, quote = F)

tab_top10up_09_genes <- top10up_09_markers.all |> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop")
write.table(tab_top10up_09_genes, "Results/Tables/02_03_FindAllMarkers_clust_top10_up_09_genelist.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10up_09_markers.all$gene)
hm_top10up_09_markers.all <- DoHeatmap(prep, features = top10up_09_markers.all$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  + NoLegend()
ggsave("Results/Figures/02_03_FindAllMarkers_clust_top10up_09_hm.png", hm_top10up_09_markers.all, width = 9, height = 12)

# Top 10 (up) - pct.1 >= 0.1
top10up_01_markers.all <- markers.all |> filter(p_val_adj <= 0.05, avg_log2FC > 0, pct.1 >= 0.1) |> group_by(cluster) |>
  arrange(desc(avg_log2FC), p_val_adj, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10up_01_markers.all,"Results/Tables/02_03_FindAllMarkers_clust_top10_01_up.tsv", sep = "\t", row.names = F, quote = F)

tab_top10up_01_genes <- top10up_01_markers.all |> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop")
write.table(tab_top10up_01_genes, "Results/Tables/02_03_FindAllMarkers_clust_top10_up_01_genelist.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10up_01_markers.all$gene)
hm_top10up_01_markers.all <- DoHeatmap(prep, features = top10up_01_markers.all$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  + NoLegend()
ggsave("Results/Figures/02_03_FindAllMarkers_clust_top10up_01_hm.png", hm_top10up_01_markers.all, width = 9, height = 12)

# Top 10 (down) - pct.1 >= 0.9
top10down_09_markers.all <- markers.all |> filter(p_val_adj <= 0.05, avg_log2FC < 0, pct.1 >= 0.9) |> group_by(cluster) |>
  arrange(avg_log2FC, p_val_adj, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10down_09_markers.all,"Results/Tables/02_03_FindAllMarkers_clust_top10_09_down.tsv", sep = "\t", row.names = F, quote = F)

tab_top10down_09_genes <- top10down_09_markers.all |> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop")
write.table(tab_top10down_09_genes, "Results/Tables/02_03_FindAllMarkers_clust_top10_down_09_genelist.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10down_09_markers.all$gene)
hm_top10down_09_markers.all <- DoHeatmap(prep, features = top10down_09_markers.all$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  + NoLegend()
ggsave("Results/Figures/02_03_FindAllMarkers_clust_top10down_09_hm.png", hm_top10down_09_markers.all, width = 9, height = 12)

# Top 10 (down) - pct.1 >= 0.1
top10down_01_markers.all <- markers.all |> filter(p_val_adj <= 0.05, avg_log2FC < 0, pct.1 >= 0.1) |> group_by(cluster) |>
  arrange(avg_log2FC, p_val_adj, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10down_01_markers.all,"Results/Tables/02_03_FindAllMarkers_clust_top10_01_down.tsv", sep = "\t", row.names = F, quote = F)

tab_top10down_01_genes <- top10down_01_markers.all |> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop")
write.table(tab_top10down_01_genes, "Results/Tables/02_03_FindAllMarkers_clust_top10_down_01_genelist.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10down_01_markers.all$gene)
hm_top10down_01_markers.all <- DoHeatmap(prep, features = top10down_01_markers.all$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  + NoLegend()
ggsave("Results/Figures/02_03_FindAllMarkers_clust_top10down_01_hm.png", hm_top10down_01_markers.all, width = 9, height = 12)


# Top 10 (up/down irrespective) - pct.1 >= 0.5
top10_05_markers.all <- markers.all |> filter(p_val_adj <= 0.05, pct.1 >= 0.5) |> group_by(cluster) |>
  arrange(desc(abs(avg_log2FC)), p_val_adj, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10_05_markers.all,"Results/Tables/02_03_FindAllMarkers_clust_top10_05_all.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10_05_markers.all$gene)
hm_top10_05_markers.all <- DoHeatmap(prep, features = top10_05_markers.all$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  + NoLegend()
ggsave("Results/Figures/02_03_FindAllMarkers_clust_top10_05_hm.png", hm_top10_05_markers.all, width = 9, height = 12)

# Top 10 (up) - pct.1 >= 0.5
top10up_05_markers.all <- markers.all |> filter(p_val_adj <= 0.05, avg_log2FC > 0, pct.1 >= 0.5) |> group_by(cluster) |>
  arrange(desc(avg_log2FC), p_val_adj, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10up_05_markers.all,"Results/Tables/02_03_FindAllMarkers_clust_top10_05_up.tsv", sep = "\t", row.names = F, quote = F)

tab_top10up_05_genes <- top10up_05_markers.all |> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop")
write.table(tab_top10up_05_genes, "Results/Tables/02_03_FindAllMarkers_clust_top10_up_05_genelist.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10up_05_markers.all$gene)
hm_top10up_05_markers.all <- DoHeatmap(prep, features = top10up_05_markers.all$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  + NoLegend()
ggsave("Results/Figures/02_03_FindAllMarkers_clust_top10up_05_hm.png", hm_top10up_05_markers.all, width = 9, height = 12)

# Top 10 (down) - pct.1 >= 0.5
top10down_05_markers.all <- markers.all |> filter(p_val_adj <= 0.05, avg_log2FC < 0, pct.1 >= 0.5) |> group_by(cluster) |>
  arrange(avg_log2FC, p_val_adj, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10down_05_markers.all,"Results/Tables/02_03_FindAllMarkers_clust_top10_05_down.tsv", sep = "\t", row.names = F, quote = F)

tab_top10down_05_genes <- top10down_05_markers.all |> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop")
write.table(tab_top10down_05_genes, "Results/Tables/02_03_FindAllMarkers_clust_top10_down_05_genelist.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10down_05_markers.all$gene)
hm_top10down_05_markers.all <- DoHeatmap(prep, features = top10down_05_markers.all$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12))  + NoLegend()
ggsave("Results/Figures/02_03_FindAllMarkers_clust_top10down_05_hm.png", hm_top10down_05_markers.all, width = 9, height = 12)

# Tabla resumen de los gene lists
tab_summary_genelists <- full_join(tab_top10up_09_genes, tab_top10up_05_genes, by = "cluster", suffix = c("_up_09", "_up_05")) |>
  full_join(tab_top10up_01_genes, by = "cluster") |> rename(gene_up_01 = gene) |>
  full_join(tab_top10down_09_genes, by = "cluster") |> rename(gene_down_09 = gene) |>
  full_join(tab_top10down_05_genes, by = "cluster") |> rename(gene_down_05 = gene) |>
  full_join(tab_top10down_01_genes, by = "cluster") |> rename(gene_down_01 = gene) |>
  arrange(as.numeric(as.character(cluster)))
write.table(tab_summary_genelists,"Results/Tables/02_03_FindAllMarkers_clust_genelist_summary.tsv", sep = "\t", row.names = F, quote = F)


# UpSet plot de marcadores up-regulados compartidos entre clusters
library(dplyr)
library(tidyr)
library(ComplexUpset)
library(ggplot2)

# colapsar genes upregulated de todos los umbrales por cluster
sets_all_up_collapsed <- bind_rows(tab_top10up_09_genes |> mutate(source = "09"), tab_top10up_05_genes |> mutate(source = "05"),
                                   tab_top10up_01_genes |> mutate(source = "01")) |> select(cluster, gene) |> separate_rows(gene, sep = ",\\s*") |> distinct() |>
  mutate(set = paste0("cluster", cluster))

# tabla binaria genes x clusters
tab_upset_all_up_collapsed <- sets_all_up_collapsed |> mutate(value = 1) |> select(set, gene, value) |> distinct() |>
  pivot_wider(names_from = set, values_from = value, values_fill = 0)
set_order_all_up_collapsed <- paste0("cluster", 1:11)

p_upset_all_up_collapsed <- upset(tab_upset_all_up_collapsed, intersect = set_order_all_up_collapsed,
  name = "Genes", sort_sets = F, sort_intersections = "descending", sort_intersections_by = c("cardinality", "degree"),
  min_degree = 2, set_sizes = F, height_ratio = 3,
  base_annotations = list("Intersection size" = intersection_size(counts = T,text = list(vjust = -0.3))))
ggsave("Results/Figures/02_03_UpSet_findallmarkers_up_collapsed.png", p_upset_all_up_collapsed, width = 10, height = 5)


# UpSet plot de marcadores down-regulados compartidos entre clusters
# colapsar genes downregulated de todos los umbrales por cluster
sets_all_down_collapsed <- bind_rows(tab_top10down_09_genes |> mutate(source = "09"), tab_top10down_05_genes |> mutate(source = "05"),
                                     tab_top10down_01_genes |> mutate(source = "01")) |> select(cluster, gene) |> separate_rows(gene, sep = ",\\s*") |> distinct() |>
  mutate(set = paste0("cluster", cluster))

# tabla binaria genes x clusters
tab_upset_all_down_collapsed <- sets_all_down_collapsed |> mutate(value = 1) |> select(set, gene, value) |> distinct() |>
  pivot_wider(names_from = set, values_from = value, values_fill = 0)
set_order_all_down_collapsed <- paste0("cluster", 1:11)

p_upset_all_down_collapsed <- upset(tab_upset_all_down_collapsed, intersect = set_order_all_down_collapsed,
  name = "Genes", sort_sets = F, sort_intersections = "descending", sort_intersections_by = c("cardinality", "degree"),
  min_degree = 2, set_sizes = F, height_ratio = 3,
  base_annotations = list("Intersection size" = intersection_size(counts = T, text = list(vjust = -0.3))))
ggsave("Results/Figures/02_03_UpSet_findallmarkers_down_collapsed.png", p_upset_all_down_collapsed, width = 10, height = 5)



## ---- 2.2 FindConservedMarkers: Marcadores conservados por cluster independientemente del tratamiento ----
clusters <- levels(Idents(prep))
conserved.all <- lapply(clusters, function(cl) {
  df <- FindConservedMarkers(prep, ident.1 = cl, assay = "SCT", grouping.var = "Treatment")
  df$gene <- rownames(df)
  df$cluster <- cl
  df
})
names(conserved.all) <- clusters
conserved.all <- bind_rows(conserved.all)
write.table(conserved.all,"Results/Tables/02_03_FindConservedMarkers_clust.tsv", sep = "\t", row.names = F, quote = F)

clusters_pos <- levels(Idents(prep))
conserved.all_pos <- lapply(clusters_pos, function(cl) {
  df <- FindConservedMarkers(prep, ident.1 = cl, assay = "SCT", grouping.var = "Treatment", only.pos = T)
  df$gene <- rownames(df)
  df$cluster <- cl
  df
})
names(conserved.all_pos) <- clusters_pos
conserved.all_pos <- bind_rows(conserved.all_pos)
write.table(conserved.all_pos,"Results/Tables/02_03_FindConservedMarkers_pos_clust.tsv", sep = "\t", row.names = F, quote = F)

# Top 50 (todos)
top50_markers.cons <- conserved.all |> filter(max_pval <= 0.05) |> group_by(cluster) |>
  arrange(desc(abs((Control_avg_log2FC + HA15_avg_log2FC) / 2)), max_pval, .by_group = T) |> slice_head(n = 50) |> ungroup()
write.table(top50_markers.cons,"Results/Tables/02_03_FindConservedMarkers_clust_top50_all.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top50_markers.cons$gene)
hm_top50_markers.cons <- DoHeatmap(prep, features = top50_markers.cons$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12)) + NoLegend()
ggsave("Results/Figures/02_03_FindConservedMarkers_clust_top50_hm.png", hm_top50_markers.cons, width = 9, height = 12)

# Top 10 (up/down irrespective) - pct.1 >= 0.9 en ambas condiciones
top10_09_markers.cons <- conserved.all |> filter(max_pval <= 0.05, Control_pct.1 >= 0.9, HA15_pct.1 >= 0.9) |> group_by(cluster) |>
  arrange(desc(abs((Control_avg_log2FC + HA15_avg_log2FC) / 2)), max_pval, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10_09_markers.cons,"Results/Tables/02_03_FindConservedMarkers_clust_top10_09_all.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10_09_markers.cons$gene)
hm_top10_09_markers.cons <- DoHeatmap(prep, features = top10_09_markers.cons$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12)) + NoLegend()
ggsave("Results/Figures/02_03_FindConservedMarkers_clust_top10_09_hm.png", hm_top10_09_markers.cons, width = 9, height = 12)

# Top 10 (up/down irrespective) - pct.1 >= 0.1 en ambas condiciones
top10_01_markers.cons <- conserved.all |> filter(max_pval <= 0.05, Control_pct.1 >= 0.1, HA15_pct.1 >= 0.1) |> group_by(cluster) |>
  arrange(desc(abs((Control_avg_log2FC + HA15_avg_log2FC) / 2)), max_pval, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10_01_markers.cons,"Results/Tables/02_03_FindConservedMarkers_clust_top10_01_all.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10_01_markers.cons$gene)
hm_top10_01_markers.cons <- DoHeatmap(prep, features = top10_01_markers.cons$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12)) + NoLegend()
ggsave("Results/Figures/02_03_FindConservedMarkers_clust_top10_01_hm.png", hm_top10_01_markers.cons, width = 9, height = 12)

# Top 10 (up) - pct.1 >= 0.9 en ambas condiciones
top10up_09_markers.cons <- conserved.all |> filter(max_pval <= 0.05, Control_avg_log2FC > 0, HA15_avg_log2FC > 0, Control_pct.1 >= 0.9, HA15_pct.1 >= 0.9) |> group_by(cluster) |>
  arrange(desc((Control_avg_log2FC + HA15_avg_log2FC) / 2), max_pval, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10up_09_markers.cons,"Results/Tables/02_03_FindConservedMarkers_clust_top10_09_up.tsv", sep = "\t", row.names = F, quote = F)

tab_top10up_09_genes.cons <- top10up_09_markers.cons |> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop")
write.table(tab_top10up_09_genes.cons, "Results/Tables/02_03_FindConservedMarkers_clust_top10_up_09_genelist.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10up_09_markers.cons$gene)
hm_top10up_09_markers.cons <- DoHeatmap(prep, features = top10up_09_markers.cons$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12)) + NoLegend()
ggsave("Results/Figures/02_03_FindConservedMarkers_clust_top10up_09_hm.png", hm_top10up_09_markers.cons, width = 9, height = 12)

# Top 10 (up) - pct.1 >= 0.1 en ambas condiciones
top10up_01_markers.cons <- conserved.all |> filter(max_pval <= 0.05, Control_avg_log2FC > 0, HA15_avg_log2FC > 0, Control_pct.1 >= 0.1, HA15_pct.1 >= 0.1) |> group_by(cluster) |>
  arrange(desc((Control_avg_log2FC + HA15_avg_log2FC) / 2), max_pval, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10up_01_markers.cons,"Results/Tables/02_03_FindConservedMarkers_clust_top10_01_up.tsv", sep = "\t", row.names = F, quote = F)

tab_top10up_01_genes.cons <- top10up_01_markers.cons |> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop")
write.table(tab_top10up_01_genes.cons, "Results/Tables/02_03_FindConservedMarkers_clust_top10_up_01_genelist.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10up_01_markers.cons$gene)
hm_top10up_01_markers.cons <- DoHeatmap(prep, features = top10up_01_markers.cons$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12)) + NoLegend()
ggsave("Results/Figures/02_03_FindConservedMarkers_clust_top10up_01_hm.png", hm_top10up_01_markers.cons, width = 9, height = 12)

# Top 10 (down) - pct.1 >= 0.9 en ambas condiciones
top10down_09_markers.cons <- conserved.all |> filter(max_pval <= 0.05, Control_avg_log2FC < 0, HA15_avg_log2FC < 0, Control_pct.1 >= 0.9, HA15_pct.1 >= 0.9) |> group_by(cluster) |>
  arrange((Control_avg_log2FC + HA15_avg_log2FC) / 2, max_pval, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10down_09_markers.cons,"Results/Tables/02_03_FindConservedMarkers_clust_top10_09_down.tsv", sep = "\t", row.names = F, quote = F)

tab_top10down_09_genes.cons <- top10down_09_markers.cons |> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop")
write.table(tab_top10down_09_genes.cons, "Results/Tables/02_03_FindConservedMarkers_clust_top10_down_09_genelist.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10down_09_markers.cons$gene)
hm_top10down_09_markers.cons <- DoHeatmap(prep, features = top10down_09_markers.cons$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12)) + NoLegend()
ggsave("Results/Figures/02_03_FindConservedMarkers_clust_top10down_09_hm.png", hm_top10down_09_markers.cons, width = 9, height = 12)

# Top 10 (down) - pct.1 >= 0.1 en ambas condiciones
top10down_01_markers.cons <- conserved.all |> filter(max_pval <= 0.05, Control_avg_log2FC < 0, HA15_avg_log2FC < 0, Control_pct.1 >= 0.1, HA15_pct.1 >= 0.1) |> group_by(cluster) |>
  arrange((Control_avg_log2FC + HA15_avg_log2FC) / 2, max_pval, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10down_01_markers.cons,"Results/Tables/02_03_FindConservedMarkers_clust_top10_01_down.tsv", sep = "\t", row.names = F, quote = F)

tab_top10down_01_genes.cons <- top10down_01_markers.cons |> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop")
write.table(tab_top10down_01_genes.cons, "Results/Tables/02_03_FindConservedMarkers_clust_top10_down_01_genelist.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10down_01_markers.cons$gene)
hm_top10down_01_markers.cons <- DoHeatmap(prep, features = top10down_01_markers.cons$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12)) + NoLegend()
ggsave("Results/Figures/02_03_FindConservedMarkers_clust_top10down_01_hm.png", hm_top10down_01_markers.cons, width = 9, height = 12)

# Top 10 (up/down irrespective) - pct.1 >= 0.5 en ambas condiciones
top10_05_markers.cons <- conserved.all |> filter(max_pval <= 0.05, Control_pct.1 >= 0.5, HA15_pct.1 >= 0.5) |> group_by(cluster) |>
  arrange(desc(abs((Control_avg_log2FC + HA15_avg_log2FC) / 2)), max_pval, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10_05_markers.cons,"Results/Tables/02_03_FindConservedMarkers_clust_top10_05_all.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10_05_markers.cons$gene)
hm_top10_05_markers.cons <- DoHeatmap(prep, features = top10_05_markers.cons$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12)) + NoLegend()
ggsave("Results/Figures/02_03_FindConservedMarkers_clust_top10_05_hm.png", hm_top10_05_markers.cons, width = 9, height = 12)

# Top 10 (up) - pct.1 >= 0.5 en ambas condiciones
top10up_05_markers.cons <- conserved.all |> filter(max_pval <= 0.05, Control_avg_log2FC > 0, HA15_avg_log2FC > 0, Control_pct.1 >= 0.5, HA15_pct.1 >= 0.5) |> group_by(cluster) |>
  arrange(desc((Control_avg_log2FC + HA15_avg_log2FC) / 2), max_pval, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10up_05_markers.cons,"Results/Tables/02_03_FindConservedMarkers_clust_top10_05_up.tsv", sep = "\t", row.names = F, quote = F)

tab_top10up_05_genes.cons <- top10up_05_markers.cons |> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop")
write.table(tab_top10up_05_genes.cons, "Results/Tables/02_03_FindConservedMarkers_clust_top10_up_05_genelist.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10up_05_markers.cons$gene)
hm_top10up_05_markers.cons <- DoHeatmap(prep, features = top10up_05_markers.cons$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12)) + NoLegend()
ggsave("Results/Figures/02_03_FindConservedMarkers_clust_top10up_05_hm.png", hm_top10up_05_markers.cons, width = 9, height = 12)

# Top 10 (down) - pct.1 >= 0.5 en ambas condiciones
top10down_05_markers.cons <- conserved.all |> filter(max_pval <= 0.05, Control_avg_log2FC < 0, HA15_avg_log2FC < 0, Control_pct.1 >= 0.5, HA15_pct.1 >= 0.5) |> group_by(cluster) |>
  arrange((Control_avg_log2FC + HA15_avg_log2FC) / 2, max_pval, .by_group = T) |> slice_head(n = 10) |> ungroup()
write.table(top10down_05_markers.cons,"Results/Tables/02_03_FindConservedMarkers_clust_top10_05_down.tsv", sep = "\t", row.names = F, quote = F)

tab_top10down_05_genes.cons <- top10down_05_markers.cons |> group_by(cluster) |> summarise(gene = paste(gene, collapse = ", "), .groups = "drop")
write.table(tab_top10down_05_genes.cons, "Results/Tables/02_03_FindConservedMarkers_clust_top10_down_05_genelist.tsv", sep = "\t", row.names = F, quote = F)

prep <- GetResidual(prep, features = top10down_05_markers.cons$gene)
hm_top10down_05_markers.cons <- DoHeatmap(prep, features = top10down_05_markers.cons$gene, assay = "SCT", size = 3) +
  theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12)) + NoLegend()
ggsave("Results/Figures/02_03_FindConservedMarkers_clust_top10down_05_hm.png", hm_top10down_05_markers.cons, width = 9, height = 12)

# Tabla resumen de los gene lists
tab_summary_genelists.cons <- full_join(tab_top10up_09_genes.cons, tab_top10up_05_genes.cons, by = "cluster", suffix = c("_up_09", "_up_05")) |>
  full_join(tab_top10up_01_genes.cons, by = "cluster") |> rename(gene_up_01 = gene) |>
  full_join(tab_top10down_09_genes.cons, by = "cluster") |> rename(gene_down_09 = gene) |>
  full_join(tab_top10down_05_genes.cons, by = "cluster") |> rename(gene_down_05 = gene) |>
  full_join(tab_top10down_01_genes.cons, by = "cluster") |> rename(gene_down_01 = gene) |>
  arrange(as.numeric(as.character(cluster)))
write.table(tab_summary_genelists.cons,"Results/Tables/02_03_FindConservedMarkers_clust_genelist_summary.tsv", sep = "\t", row.names = F, quote = F)


# UpSet plot de marcadores up-regulados conservados compartidos entre clusters
library(dplyr)
library(tidyr)
library(ComplexUpset)
library(ggplot2)

sets_all_up_collapsed.cons <- bind_rows(tab_top10up_09_genes.cons |> mutate(source = "09"), tab_top10up_05_genes.cons |> mutate(source = "05"), tab_top10up_01_genes.cons |> mutate(source = "01")) |> select(cluster, gene) |> separate_rows(gene, sep = ",\\s*") |> distinct() |>
  mutate(set = paste0("cluster", cluster))

tab_upset_all_up_collapsed.cons <- sets_all_up_collapsed.cons |> mutate(value = 1) |> select(set, gene, value) |> distinct() |>
  pivot_wider(names_from = set, values_from = value, values_fill = 0)
set_order_all_up_collapsed.cons <- paste0("cluster", 1:11)

p_upset_all_up_collapsed.cons <- upset(tab_upset_all_up_collapsed.cons, intersect = set_order_all_up_collapsed.cons,
  name = "Genes", sort_sets = F, sort_intersections = "descending", sort_intersections_by = c("cardinality", "degree"),
  min_degree = 2, set_sizes = F, height_ratio = 3,
  base_annotations = list("Intersection size" = intersection_size(counts = T, text = list(vjust = -0.3))))
ggsave("Results/Figures/02_03_UpSet_findconservedmarkers_up_collapsed.png", p_upset_all_up_collapsed.cons, width = 10, height = 5)

# UpSet plot de marcadores down-regulados conservados compartidos entre clusters
sets_all_down_collapsed.cons <- bind_rows(tab_top10down_09_genes.cons |> mutate(source = "09"), tab_top10down_05_genes.cons |> mutate(source = "05"),
                                          tab_top10down_01_genes.cons |> mutate(source = "01")) |> select(cluster, gene) |> separate_rows(gene, sep = ",\\s*") |> distinct() |>
  mutate(set = paste0("cluster", cluster))

tab_upset_all_down_collapsed.cons <- sets_all_down_collapsed.cons |> mutate(value = 1) |> select(set, gene, value) |> distinct() |>
  pivot_wider(names_from = set, values_from = value, values_fill = 0)
set_order_all_down_collapsed.cons <- paste0("cluster", 1:11)

p_upset_all_down_collapsed.cons <- upset(tab_upset_all_down_collapsed.cons, intersect = set_order_all_down_collapsed.cons,
  name = "Genes", sort_sets = F, sort_intersections = "descending", sort_intersections_by = c("cardinality", "degree"),
  min_degree = 2, set_sizes = F, height_ratio = 3,
  base_annotations = list("Intersection size" = intersection_size(counts = T, text = list(vjust = -0.3))))
ggsave("Results/Figures/02_03_UpSet_findconservedmarkers_down_collapsed.png", p_upset_all_down_collapsed.cons, width = 10, height = 5)


