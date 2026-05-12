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
# Script: 00_benchmark_int_clust.R
## ===============================

# Objetivo:
# Comparar estrategias de normalización, integración y clustering sobre el objeto
# Seurat filtrado, incluyendo SCTransform v2 sin integración, integración Seurat v4
# y Seurat v5, split por muestra o réplica, CCA/RPCA y resoluciones Leiden 0.25,
# 0.5 y 1, generando objetos, UMAPs, elbow plots y tabla resumen de benchmarking.

# Inputs en Proc/:
# - 01_data_prep_and_QC.rds

# Outputs en Proc/:
# - 02_01_*.rds
# - 02_02_*.rds

# Outputs en Results/:
# - Figures/02_02_*.png
# - Tables/02_02_benchmark_summary.csv


## ==========================
## BLOQUE 1 - SCTransform v2
## ==========================
library(Seurat)
library(glmGamPoi)
library(patchwork)
library(ggplot2)

obj <- readRDS("Proc/01_data_prep_and_QC.rds")
obj_byrep <- obj

DefaultAssay(obj) <- "RNA"
DefaultAssay(obj_byrep) <- "RNA"

#split por muestra en lista de objetos (para tests #1,2), Seurat v4
list <- SplitObject(obj, split.by = "Sample")

# SCTransform por muestra (para tests #1,2), Seurat v4
list <- lapply(X = list, FUN = SCTransform, vst.flavor = "v2", return.only.var.genes = FALSE,
                   vars.to.regress = "percent.mt", verbose = FALSE)
saveRDS(list, "Proc/02_01_sct_bysample_v4.rds")

# paso adicional de FindVariableFeatures tras SCTransform (para test #1, paper original), Seurat v4
list2 <- lapply(X = list, FUN = FindVariableFeatures, selection.method = "vst", nfeatures = 3000, assay = "SCT")
saveRDS(list2, "Proc/02_01_sct_bysample_v4_paper.rds")

#split por réplica en lista de objetos (para tests #6), Seurat v4
list_byrep <- SplitObject(obj, split.by = "Replicate")

# SCTransform por réplica (para test #6), Seurat v4
list_byrep <- lapply(X = list_byrep, FUN = SCTransform, vst.flavor = "v2", return.only.var.genes = FALSE,
               vars.to.regress = "percent.mt", verbose = FALSE)
saveRDS(list_byrep, "Proc/02_01_sct_byrep_v4.rds")

# control de SCTransform sin split ni integrar (para test #3)
obj_ctrl <- SCTransform(obj, vst.flavor = "v2", return.only.var.genes = FALSE, 
                        vars.to.regress = "percent.mt", verbose = FALSE)
saveRDS(obj_ctrl, "Proc/02_01_sct_nonsplit.rds")

# split assay RNA en layers por muestra para SCTransform (para tests #4,7), Seurat v5
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Sample)

# SCTransform por muestra (para tests #4,7), Seurat v5
obj <- SCTransform(obj, vst.flavor = "v2", return.only.var.genes = FALSE, 
                   vars.to.regress = "percent.mt", verbose = FALSE)
saveRDS(obj, "Proc/02_01_sct_bysample.rds")

# split assay RNA en layers por réplica para SCTransform (para test #5), Seurat v5
obj_byrep[["RNA"]] <- split(obj_byrep[["RNA"]], f = obj_byrep$Replicate)

# SCTransform por réplica (para test #5), Seurat v5
obj_byrep <- SCTransform(obj_byrep, vst.flavor = "v2", return.only.var.genes = FALSE, 
                         vars.to.regress = "percent.mt", verbose = FALSE)
saveRDS(obj_byrep, "Proc/02_01_sct_byrep.rds")


#top10 <- head(VariableFeatures(obj), 10)
#top10
#VariableFeaturePlot(obj)



## ========================================================================
## BLOQUE 2 - Integración, reducción de dimensionalidad y clustering
## ========================================================================

## ---- 2.1 Integración ----

# Integración (Seurat v4) por muestra (para test #2)
all.features_list <- SelectIntegrationFeatures(object.list = list, nfeatures = 5000)
list <- PrepSCTIntegration(object.list = list, anchor.features = all.features_list)
all.anchors_list <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT", anchor.features = all.features_list)
int_list <- IntegrateData(anchorset = all.anchors_list, normalization.method = "SCT")

saveRDS(int_list, "Proc/02_02_sct_integrated_v4_bysample.rds")

# Integración (Seurat v4) por muestra (para test #1)
all.features_list2 <- SelectIntegrationFeatures(object.list = list2, nfeatures = 5000)
list2 <- PrepSCTIntegration(object.list = list2, anchor.features = all.features_list2)
all.anchors_list2 <- FindIntegrationAnchors(object.list = list2, normalization.method = "SCT", anchor.features = all.features_list2)
int_list2 <- IntegrateData(anchorset = all.anchors_list2, normalization.method = "SCT")

saveRDS(int_list2, "Proc/02_02_sct_integrated_v4_bysample_paper.rds")

# Integración (Seurat v4) por réplica (para test #6)
all.features_byrep <- SelectIntegrationFeatures(object.list = list_byrep, nfeatures = 5000)
list_byrep <- PrepSCTIntegration(object.list = list_byrep, anchor.features = all.features_byrep)
all.anchors_byrep <- FindIntegrationAnchors(object.list = list_byrep, normalization.method = "SCT", anchor.features = all.features_byrep)
int_list_byrep <- IntegrateData(anchorset = all.anchors_byrep, normalization.method = "SCT")

saveRDS(int_list_byrep, "Proc/02_02_sct_integrated_v4_byrep.rds")

# Integración (Seurat v5) por muestra (para test #4)
set.seed(123)
obj <- RunPCA(obj, npcs = 300, verbose = FALSE)
int <- IntegrateLayers(object = obj, method = CCAIntegration, normalization.method = "SCT",
                       orig.reduction = "pca", new.reduction = "integrated.dr")
saveRDS(int, "Proc/02_02_sct_integrated.rds")

# Integración (Seurat v5) por muestra, RPCA (para test #7)
int_rpca <- IntegrateLayers(object = obj, method = RPCAIntegration, normalization.method = "SCT",
                            orig.reduction = "pca", new.reduction = "integrated.rpca")
saveRDS(int_rpca, "Proc/02_02_sct_integrated_rpca.rds")

# Integración (Seurat v5) por réplica, (para test #5)
set.seed(123)
obj_byrep <- RunPCA(obj_byrep, npcs = 300, verbose = FALSE)
int_byrep <- IntegrateLayers(object = obj_byrep, method = CCAIntegration, normalization.method = "SCT",
                             orig.reduction = "pca", new.reduction = "integrated.dr")
saveRDS(int_byrep, "Proc/02_02_sct_integrated_byrep.rds")



## ---- 2.2 Reducción de dimensionalidad y clustering ----

# (Seurat v4) por muestra (para test #2)
set.seed(123)
int_list <- RunPCA(int_list, verbose = FALSE, npcs = 300)
p_elbow_int_list <- ElbowPlot(int_list, ndims = 300, reduction = "pca")
ggsave("Results/Figures/02_02_ElbowPlot_PCA_v4.png", p_elbow_int_list, width = 8, height = 5, dpi = 300)

set.seed(123)
int_list <- RunUMAP(int_list, dims = 1:30, verbose = FALSE)
int_list <- FindNeighbors(int_list, dims = 1:30, verbose = FALSE)
set.seed(123)
int_list <- FindClusters(int_list, verbose = FALSE, resolution = 0.25, algorithm = 4)
saveRDS(int_list, "Proc/02_02_sct_integrated_v4_clust0.25.rds")

# (Seurat v4) por muestra (para test #1)
set.seed(123)
int_list2 <- RunPCA(int_list2, verbose = FALSE, npcs = 300)
p_elbow_int_list2 <- ElbowPlot(int_list2, ndims = 300, reduction = "pca")
ggsave("Results/Figures/02_02_ElbowPlot_PCA_v4_paper.png", p_elbow_int_list2, width = 8, height = 5, dpi = 300)

set.seed(123)
int_list2 <- RunUMAP(int_list2, dims = 1:30, verbose = FALSE)
int_list2 <- FindNeighbors(int_list2, dims = 1:30, verbose = FALSE)
set.seed(123)
int_list2 <- FindClusters(int_list2, verbose = FALSE, resolution = 0.25, algorithm = 4)
saveRDS(int_list2, "Proc/02_02_sct_integrated_v4_paper_clust0.25.rds")

# (Seurat v4) por réplica (para test #6)
set.seed(123)
int_list_byrep <- RunPCA(int_list_byrep, verbose = FALSE, npcs = 300)
p_elbow_int_list_byrep <- ElbowPlot(int_list_byrep, ndims = 300, reduction = "pca")
ggsave("Results/Figures/02_02_ElbowPlot_PCA_v4_byrep.png", p_elbow_int_list_byrep, width = 8, height = 5, dpi = 300)

set.seed(123)
int_list_byrep <- RunUMAP(int_list_byrep, dims = 1:30, verbose = FALSE)
int_list_byrep <- FindNeighbors(int_list_byrep, dims = 1:30, verbose = FALSE)
set.seed(123)
int_list_byrep <- FindClusters(int_list_byrep, verbose = FALSE, resolution = 0.25, algorithm = 4)
saveRDS(int_list_byrep, "Proc/02_02_sct_integrated_v4_byrep_clust0.25.rds")

# split por muestra + integración
set.seed(123)
int <- RunUMAP(int, reduction = "integrated.dr", dims = 1:30, verbose = FALSE)
int <- FindNeighbors(int, reduction = "integrated.dr", dims = 1:30, verbose = FALSE)
set.seed(123)
int <- FindClusters(int, resolution = 0.25, algorithm = 4, verbose = FALSE)

saveRDS(int, "Proc/02_02_sct_integrated_clust0.25.rds")

# split por muestra + integración, RPCA
set.seed(123)
int_rpca <- RunUMAP(int_rpca, reduction = "integrated.rpca", dims = 1:30, verbose = FALSE)
int_rpca <- FindNeighbors(int_rpca, reduction = "integrated.rpca", dims = 1:30, verbose = FALSE)
set.seed(123)
int_rpca <- FindClusters(int_rpca, resolution = 0.25, algorithm = 4, verbose = FALSE)

saveRDS(int_rpca, "Proc/02_02_sct_integrated_rpca_clust0.25.rds")

# split por réplica + integración
set.seed(123)
int_byrep <- RunUMAP(int_byrep, reduction = "integrated.dr", dims = 1:30, verbose = FALSE)
int_byrep <- FindNeighbors(int_byrep, reduction = "integrated.dr", dims = 1:30, verbose = FALSE)
set.seed(123)
int_byrep <- FindClusters(int_byrep, resolution = 0.25, algorithm = 4, verbose = FALSE)

saveRDS(int_byrep, "Proc/02_02_sct_integrated_byrep_clust0.25.rds")

# control sin split ni integración
set.seed(123)
obj_ctrl <- RunPCA(obj_ctrl, verbose = FALSE, npcs = 300)
obj_ctrl <- RunUMAP(obj_ctrl, dims = 1:30, verbose = FALSE)
obj_ctrl <- FindNeighbors(obj_ctrl, dims = 1:30, verbose = FALSE)
set.seed(123)
obj_ctrl <- FindClusters(obj_ctrl, verbose = FALSE, resolution = 0.25, algorithm = 4)

saveRDS(obj_ctrl, "Proc/02_02_sct_nonintegrated_clust0.25.rds")



## ---- Resoluciones adicionales ----

# Seurat v4 por muestra (test #2)
set.seed(123)
int_list_res0.5 <- FindClusters(int_list, resolution = 0.5, algorithm = 4, verbose = FALSE)
set.seed(123)
int_list_res1 <- FindClusters(int_list, resolution = 1, algorithm = 4, verbose = FALSE)

saveRDS(int_list_res0.5, "Proc/02_02_sct_integrated_v4_clust0.5.rds")
saveRDS(int_list_res1, "Proc/02_02_sct_integrated_v4_clust1.rds")

# Seurat v4 por muestra con paso adicional de FindVariableFeatures post-SCT (test #1)
set.seed(123)
int_list2_res0.5 <- FindClusters(int_list2, resolution = 0.5, algorithm = 4, verbose = FALSE)
set.seed(123)
int_list2_res1 <- FindClusters(int_list2, resolution = 1, algorithm = 4, verbose = FALSE)

saveRDS(int_list2_res0.5, "Proc/02_02_sct_integrated_v4_paper_clust0.5.rds")
saveRDS(int_list2_res1, "Proc/02_02_sct_integrated_v4_paper_clust1.rds")

# Seurat v4 por réplica (test #6)
set.seed(123)
int_list_byrep_res0.5 <- FindClusters(int_list_byrep, resolution = 0.5, algorithm = 4, verbose = FALSE)
set.seed(123)
int_list_byrep_res1 <- FindClusters(int_list_byrep, resolution = 1, algorithm = 4, verbose = FALSE)

saveRDS(int_list_byrep_res0.5, "Proc/02_02_sct_integrated_v4_byrep_clust0.5.rds")
saveRDS(int_list_byrep_res1, "Proc/02_02_sct_integrated_v4_byrep_clust1.rds")

# split por muestra + integración
set.seed(123)
int_res0.5 <- FindClusters(int, resolution = 0.5, algorithm = 4, verbose = FALSE)
set.seed(123)
int_res1 <- FindClusters(int, resolution = 1, algorithm = 4, verbose = FALSE)

saveRDS(int_res0.5, "Proc/02_02_sct_integrated_clust0.5.rds")
saveRDS(int_res1, "Proc/02_02_sct_integrated_clust1.rds")

# split por muestra + integración, RPCA
set.seed(123)
int_rpca_res0.5 <- FindClusters(int_rpca, resolution = 0.5, algorithm = 4, verbose = FALSE)
set.seed(123)
int_rpca_res1 <- FindClusters(int_rpca, resolution = 1, algorithm = 4, verbose = FALSE)

saveRDS(int_rpca_res0.5, "Proc/02_02_sct_integrated_rpca_clust0.5.rds")
saveRDS(int_rpca_res1, "Proc/02_02_sct_integrated_rpca_clust1.rds")

# split por réplica + integración
set.seed(123)
int_byrep_res0.5 <- FindClusters(int_byrep, resolution = 0.5, algorithm = 4, verbose = FALSE)
set.seed(123)
int_byrep_res1 <- FindClusters(int_byrep, resolution = 1, algorithm = 4, verbose = FALSE)

saveRDS(int_byrep_res0.5, "Proc/02_02_sct_integrated_byrep_clust0.5.rds")
saveRDS(int_byrep_res1, "Proc/02_02_sct_integrated_byrep_clust1.rds")

# control sin split ni integración
set.seed(123)
obj_ctrl_res0.5 <- FindClusters(obj_ctrl, resolution = 0.5, algorithm = 4, verbose = FALSE)
set.seed(123)
obj_ctrl_res1 <- FindClusters(obj_ctrl, resolution = 1, algorithm = 4, verbose = FALSE)

saveRDS(obj_ctrl_res0.5, "Proc/02_02_sct_nonintegrated_clust0.5.rds")
saveRDS(obj_ctrl_res1, "Proc/02_02_sct_nonintegrated_clust1.rds")


## ---- 2.3 Elbow Plots ----

# split por muestra + integración
p_elbow <- ElbowPlot(int, ndims = 300, reduction = "pca")
ggsave("Results/Figures/02_02_elbow_PCA.png", p_elbow, width = 8, height = 5, dpi = 300)

# split por muestra + integración, RPCA
p_elbow_rpca <- ElbowPlot(int_rpca, ndims = 300, reduction = "pca")
ggsave("Results/Figures/02_02_elbow_rpca_PCA.png", p_elbow_rpca, width = 8, height = 5, dpi = 300)

# split por réplica + integración
p_elbow_byrep <- ElbowPlot(int_byrep, ndims = 300, reduction = "pca")
ggsave("Results/Figures/02_02_elbow_PCA_byrep.png", p_elbow_byrep, width = 8, height = 5, dpi = 300)

# control sin split ni integración
p_elbow_ctrl <- ElbowPlot(obj_ctrl, ndims = 300, reduction = "pca") + ggtitle("SCTransform sin split ni integración")
ggsave("Results/Figures/02_02_elbow_PCA_ctrl.png", p_elbow_ctrl, width = 8, height = 5, dpi = 300)


## ---- 2.4 Figuras UMAP ----

# split por muestra + integración, Seurat v4 (test #2)
p_clust_v4 <- DimPlot(int_list, reduction = "umap", label = TRUE)
p_samp_v4 <- DimPlot(int_list, reduction = "umap", group.by = "Sample")
p_treat_v4 <- DimPlot(int_list, reduction = "umap", group.by = "Treatment")
p_rep_v4 <- DimPlot(int_list, reduction = "umap", group.by = "Replicate")
p_umap_int_v4 <- (p_clust_v4 | p_samp_v4) / (p_treat_v4 | p_rep_v4)

ggsave("Results/Figures/02_02_UMAP_combined_integrated_v4.png", p_umap_int_v4, width = 12, height = 10, dpi = 300)

# split por muestra + integración con paso adicional de FindVariableFeatures post-SCT, Seurat v4 (test #1)
p_clust_v4_paper <- DimPlot(int_list2, reduction = "umap", label = TRUE)
p_samp_v4_paper <- DimPlot(int_list2, reduction = "umap", group.by = "Sample")
p_treat_v4_paper <- DimPlot(int_list2, reduction = "umap", group.by = "Treatment")
p_rep_v4_paper <- DimPlot(int_list2, reduction = "umap", group.by = "Replicate")
p_umap_int_v4_paper <- (p_clust_v4_paper | p_samp_v4_paper) / (p_treat_v4_paper | p_rep_v4_paper)

ggsave("Results/Figures/02_02_UMAP_combined_integrated_v4_paper.png", p_umap_int_v4_paper, width = 12, height = 10, dpi = 300)

# split por réplica + integración, Seurat v4 (test #6)
p_clust_v4_byrep <- DimPlot(int_list_byrep, reduction = "umap", label = TRUE)
p_samp_v4_byrep <- DimPlot(int_list_byrep, reduction = "umap", group.by = "Sample")
p_treat_v4_byrep <- DimPlot(int_list_byrep, reduction = "umap", group.by = "Treatment")
p_rep_v4_byrep <- DimPlot(int_list_byrep, reduction = "umap", group.by = "Replicate")
p_umap_int_v4_byrep <- (p_clust_v4_byrep | p_samp_v4_byrep) / (p_treat_v4_byrep | p_rep_v4_byrep)

ggsave("Results/Figures/02_02_UMAP_combined_integrated_v4_byrep.png", p_umap_int_v4_byrep, width = 12, height = 10, dpi = 300)

# split por muestra + integración (test #4)
p_clust <- DimPlot(int, reduction = "umap", label = TRUE)
p_samp <- DimPlot(int, reduction = "umap", group.by = "Sample")
p_treat <- DimPlot(int, reduction = "umap", group.by = "Treatment")
p_rep <- DimPlot(int, reduction = "umap", group.by = "Replicate")
p_umap_int <- (p_clust | p_samp) / (p_treat | p_rep)

#ggsave("Results/Figures/02_02_UMAP_clusters.png", p_clust, width = 7, height = 5, dpi = 300)
#ggsave("Results/Figures/02_02_UMAP_sample.png", p_samp, width = 7, height = 5, dpi = 300)
#ggsave("Results/Figures/02_02_UMAP_treatment.png", p_treat, width = 7, height = 5, dpi = 300)
#ggsave("Results/Figures/02_02_UMAP_replicate.png", p_rep, width = 7, height = 5, dpi = 300)
ggsave("Results/Figures/02_02_UMAP_combined_integrated.png", p_umap_int, width = 12, height = 10, dpi = 300)

# split por muestra + integración, RPCA (test #7)
p_clust_rpca <- DimPlot(int_rpca, reduction = "umap", label = TRUE)
p_samp_rpca <- DimPlot(int_rpca, reduction = "umap", group.by = "Sample")
p_treat_rpca <- DimPlot(int_rpca, reduction = "umap", group.by = "Treatment")
p_rep_rpca <- DimPlot(int_rpca, reduction = "umap", group.by = "Replicate")
p_umap_int_rpca <- (p_clust_rpca | p_samp_rpca) / (p_treat_rpca | p_rep_rpca)

#ggsave("Results/Figures/02_02_UMAP_clusters_rpca.png", p_clust_rpca, width = 7, height = 5, dpi = 300)
#ggsave("Results/Figures/02_02_UMAP_sample_rpca.png", p_samp_rpca, width = 7, height = 5, dpi = 300)
#ggsave("Results/Figures/02_02_UMAP_treatment_rpca.png", p_treat_rpca, width = 7, height = 5, dpi = 300)
#ggsave("Results/Figures/02_02_UMAP_replicate_rpca.png", p_rep_rpca, width = 7, height = 5, dpi = 300)
ggsave("Results/Figures/02_02_UMAP_combined_integrated_rpca.png", p_umap_int_rpca, width = 12, height = 10, dpi = 300)

# split por réplica + integración (test #5)
p_clust_byrep <- DimPlot(int_byrep, reduction = "umap", label = TRUE)
p_samp_byrep <- DimPlot(int_byrep, reduction = "umap", group.by = "Sample")
p_treat_byrep <- DimPlot(int_byrep, reduction = "umap", group.by = "Treatment")
p_rep_byrep <- DimPlot(int_byrep, reduction = "umap", group.by = "Replicate")
p_umap_int_byrep <- (p_clust_byrep | p_samp_byrep) / (p_treat_byrep | p_rep_byrep)

#ggsave("Results/Figures/02_02_UMAP_clusters_byrep.png", p_clust_byrep, width = 7, height = 5, dpi = 300)
#ggsave("Results/Figures/02_02_UMAP_sample_byrep.png", p_samp_byrep, width = 7, height = 5, dpi = 300)
#ggsave("Results/Figures/02_02_UMAP_treatment_byrep.png", p_treat_byrep, width = 7, height = 5, dpi = 300)
#ggsave("Results/Figures/02_02_UMAP_replicate_byrep.png", p_rep_byrep, width = 7, height = 5, dpi = 300)
ggsave("Results/Figures/02_02_UMAP_combined_integrated_byrep.png", p_umap_int_byrep, width = 12, height = 10, dpi = 300)

# control sin split ni integración (test #3)
p_clust_ctrl <- DimPlot(obj_ctrl, reduction = "umap", label = TRUE)
p_samp_ctrl <- DimPlot(obj_ctrl, reduction = "umap", group.by = "Sample")
p_treat_ctrl <- DimPlot(obj_ctrl, reduction = "umap", group.by = "Treatment")
p_rep_ctrl <- DimPlot(obj_ctrl, reduction = "umap", group.by = "Replicate")
p_umap_ctrl <- (p_clust_ctrl | p_samp_ctrl) / (p_treat_ctrl | p_rep_ctrl)

#ggsave("Results/Figures/02_02_UMAP_clusters_ctrl.png", p_clust_ctrl, width = 7, height = 5, dpi = 300)
#ggsave("Results/Figures/02_02_UMAP_sample_ctrl.png", p_samp_ctrl, width = 7, height = 5, dpi = 300)
#ggsave("Results/Figures/02_02_UMAP_treatment_ctrl.png", p_treat_ctrl, width = 7, height = 5, dpi = 300)
#ggsave("Results/Figures/02_02_UMAP_replicate_ctrl.png", p_rep_ctrl, width = 7, height = 5, dpi = 300)
ggsave("Results/Figures/02_02_UMAP_combined_ctrl.png", p_umap_ctrl, width = 12, height = 10, dpi = 300)

## ---- Patchworks multirresolución (3 columnas x 2 filas) ----

make_multires_patchwork <- function(obj_res0.25, obj_res0.5, obj_res1, reduction_name = "umap") {
  # fila 1: clustering a distintas resoluciones
  p_clust_0.25 <- DimPlot(obj_res0.25, reduction = reduction_name, label = TRUE) +
    ggtitle("Clusters (res = 0.25)")
  
  p_clust_0.5 <- DimPlot(obj_res0.5, reduction = reduction_name, label = TRUE) +
    ggtitle("Clusters (res = 0.5)")
  
  p_clust_1 <- DimPlot(obj_res1, reduction = reduction_name, label = TRUE) +
    ggtitle("Clusters (res = 1)")
  
  # fila 2: metadatos sobre el mismo embedding
  p_samp <- DimPlot(obj_res0.25, reduction = reduction_name, group.by = "Sample") +
    ggtitle("Sample")
  
  p_treat <- DimPlot(obj_res0.25, reduction = reduction_name, group.by = "Treatment") +
    ggtitle("Treatment")
  
  p_rep <- DimPlot(obj_res0.25, reduction = reduction_name, group.by = "Replicate") +
    ggtitle("Replicate")
  
  (p_clust_0.25 | p_clust_0.5 | p_clust_1) /
    (p_samp | p_treat | p_rep)
}

# patchworks finales v4
p_umap_int_v4_multires <- make_multires_patchwork(int_list, int_list_res0.5, int_list_res1, reduction_name = "umap")
p_umap_int_v4_paper_multires <- make_multires_patchwork(int_list2, int_list2_res0.5, int_list2_res1, reduction_name = "umap")
p_umap_int_v4_byrep_multires <- make_multires_patchwork(int_list_byrep, int_list_byrep_res0.5, int_list_byrep_res1, reduction_name = "umap")

ggsave("Results/Figures/02_02_UMAP_multires_integrated_v4.png", p_umap_int_v4_multires, width = 16, height = 10, dpi = 300)
ggsave("Results/Figures/02_02_UMAP_multires_integrated_v4_paper.png", p_umap_int_v4_paper_multires, width = 16, height = 10, dpi = 300)
ggsave("Results/Figures/02_02_UMAP_multires_integrated_v4_byrep.png", p_umap_int_v4_byrep_multires, width = 16, height = 10, dpi = 300)

# patchworks finales v5
p_umap_int_multires <- make_multires_patchwork(int, int_res0.5, int_res1, reduction_name = "umap")
p_umap_int_rpca_multires <- make_multires_patchwork(int_rpca, int_rpca_res0.5, int_rpca_res1, reduction_name = "umap")
p_umap_int_byrep_multires <- make_multires_patchwork(int_byrep, int_byrep_res0.5, int_byrep_res1, reduction_name = "umap")
p_umap_ctrl_multires <- make_multires_patchwork(obj_ctrl, obj_ctrl_res0.5, obj_ctrl_res1, reduction_name = "umap")

ggsave("Results/Figures/02_02_UMAP_multires_integrated.png", p_umap_int_multires, width = 16, height = 10, dpi = 300)
ggsave("Results/Figures/02_02_UMAP_multires_integrated_rpca.png", p_umap_int_rpca_multires, width = 16, height = 10, dpi = 300)
ggsave("Results/Figures/02_02_UMAP_multires_integrated_byrep.png", p_umap_int_byrep_multires, width = 16, height = 10, dpi = 300)
ggsave("Results/Figures/02_02_UMAP_multires_ctrl.png", p_umap_ctrl_multires, width = 16, height = 10, dpi = 300)



## ---- 2.5 Tabla resumen del benchmarking ----

summarise_clusters <- function(obj, strategy, resolution) {
  cl_sizes <- as.integer(table(Idents(obj)))
  data.frame(
    Strategy = strategy,
    Resolution = resolution,
    N_cells = ncol(obj),
    N_clusters = length(cl_sizes),
    Min_cells_per_cluster = min(cl_sizes),
    Median_cells_per_cluster = median(cl_sizes),
    Max_cells_per_cluster = max(cl_sizes)
  )
}

benchmark_summary <- rbind(
  summarise_clusters(int_list, "integrated_v4", 0.25),
  summarise_clusters(int_list_res0.5, "integrated_v4", 0.5),
  summarise_clusters(int_list_res1, "integrated_v4", 1),
  
  summarise_clusters(int_list2, "integrated_v4_paper", 0.25),
  summarise_clusters(int_list2_res0.5, "integrated_v4_paper", 0.5),
  summarise_clusters(int_list2_res1, "integrated_v4_paper", 1),
  
  summarise_clusters(int_list_byrep, "integrated_v4_byrep", 0.25),
  summarise_clusters(int_list_byrep_res0.5, "integrated_v4_byrep", 0.5),
  summarise_clusters(int_list_byrep_res1, "integrated_v4_byrep", 1),
  
  summarise_clusters(int, "integrated", 0.25),
  summarise_clusters(int_res0.5, "integrated", 0.5),
  summarise_clusters(int_res1, "integrated", 1),
  
  summarise_clusters(int_rpca, "integrated_rpca", 0.25),
  summarise_clusters(int_rpca_res0.5, "integrated_rpca", 0.5),
  summarise_clusters(int_rpca_res1, "integrated_rpca", 1),
  
  summarise_clusters(int_byrep, "integrated_byrep", 0.25),
  summarise_clusters(int_byrep_res0.5, "integrated_byrep", 0.5),
  summarise_clusters(int_byrep_res1, "integrated_byrep", 1),
  
  summarise_clusters(obj_ctrl, "ctrl_nonsplit", 0.25),
  summarise_clusters(obj_ctrl_res0.5, "ctrl_nonsplit", 0.5),
  summarise_clusters(obj_ctrl_res1, "ctrl_nonsplit", 1)
)

write.csv(benchmark_summary, "Results/Tables/02_02_benchmark_summary.csv", row.names = FALSE)

