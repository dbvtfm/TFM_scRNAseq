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
# Script: 00_benchmark_ann_cell.R
## ===============================

# Objetivo:
# Comparar de forma ampliada métodos y referencias de anotación celular automática,
# incluyendo SingleR, ScType, EasyCellType y mapeo a atlas CD4 basado en Symphony,
# usando distintas bases de referencia y excluyendo genes de ciclo celular

# Inputs en Proc/:
# - 02_02_sct_int-rpca_leiden_1.rds
# - 01_data_prep_and_QC.rds

# Inputs en Ann/:
# - pbmc_multimodal_2023.rds (https://zenodo.org/records/7779017/files/pbmc_multimodal_2023.rds)
# - ScTypeDB_full.xlsx
# - Cell_marker_Human.xlsx
# - Cell_marker_Seq.xlsx
# - screfmapping/ (https://github.com/yyoshiaki/screfmapping)
# - screfmapping/data/ref_Reference_Mapping_20220525.RData (https://doi.org/10.6084/m9.figshare.25052648)

# Inputs en Results/:
# - Tables/03_02_FindConservedMarkers_clust_top10_05_up.tsv
# - Tables/02_03_FindAllMarkers_pos_clust.tsv

# Outputs en Results/

#----------------------------
library(Seurat)
library(SingleCellExperiment)

# inputs from script 02b_scRNAseq_basic.R
ann <- readRDS("Proc/02_02_sct_int-rpca_leiden_1.rds")
DefaultAssay(ann) <- "SCT"

tab.genes <- read.delim("Results/Tables/03_02_FindConservedMarkers_clust_top10_05_up.tsv",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gene.markers <- tab.genes$gene
gene.markers <- unique(na.omit(gene.markers))

# genes de ciclo celular a excluir
cycle.markers <- c(
  "BIRC5","DLGAP5","UBE2C","KIF14","CDCA2","ESCO2","MKI67","NEK2","PLK1",
  "KIF15","SPC24","HJURP","HMMR","ANLN","TROAP","SPC25","KIF23","CKAP2L",
  "GTSE1","CENPF","RRM2","ASPM","CDCA8","CDCA3","MYBL2","TYMS","KIF4A",
  "KIFC1","STMN1","TOP2A","DTL","MCM10","CDC20","KIF20A","CCNB2","SKA1",
  "CCNE2","DEPDC1","CLSPN","DIAPH3","DSCC1","CEP55","ZWINT","BUB1B","CIT",
  "ASF1B","CENPA","TK1","SKA3","PBK","RAD51AP1","CDT1","SHCBP1","STIL",
  "PCLAF","KIF2C","TRIP13","E2F2","CENPW","CDC6","CDK1","E2F8","CDCA5",
  "CDC25C","PIMREG","KIF18B","AURKB","PSRC1","FBXO43","GINS2","E2F1","MCM2",
  "XRCC2","ORC1","CDC25A","MCM4","UNG","UHRF1","AUNIP","CDC45","RMI2",
  "GINS1","CHEK1","CENPU","MCM7","GINS3","POLE2","ZNF367","MCM5",
  "H1-5","H2AC14","H2AC12","H2BC3","H2AC4","H3C2","H3C15","H3C8","H3C3","H3C7")

nocycle.markers <- gene.markers[!gene.markers %in% cycle.markers]


## ==========================================
## BLOQUE 1 - ANOTACIÓN AUTOMÁTICA: SingleR
## ==========================================
library(SingleR)
library(celldex)
library(scrapper)

## ---- 1.1 SingleR basado en bases de datos integradas ----

# bases de datos
ref_dice <- celldex::DatabaseImmuneCellExpressionData()
ref_mi <- celldex::MonacoImmuneData()
ref_nhd <- celldex::NovershternHematopoieticData()
ref_bp <- celldex::BlueprintEncodeData()

# dataframe de resultados de las etiquetas celulares mapeadas
# el objeto seurat se tiene que pasar como SingleCellExperiment
# observar df ref y las labels con colData(ref)
res_dice <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_dice, labels = ref_dice$label.fine,
                    restrict = nocycle.markers)
res_dice_main <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_dice, labels = ref_dice$label.main,
                         restrict = nocycle.markers)

res_mi <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_mi, labels = ref_mi$label.fine, 
                  restrict = nocycle.markers)
res_mi_main <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_mi, labels = ref_mi$label.main,
                       restrict = nocycle.markers)

res_nhd <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_nhd, labels = ref_nhd$label.fine,
                   restrict = nocycle.markers)
res_nhd_main <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_nhd, labels = ref_nhd$label.main,
                        restrict = nocycle.markers)

res_dice_clust <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_dice, labels = ref_dice$label.fine,
                          clusters = ann$seurat_clusters, restrict = nocycle.markers)
res_dice_main_clust <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_dice,
                               labels = ref_dice$label.main, clusters = ann$seurat_clusters, restrict = nocycle.markers)

res_bp <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_bp, labels = ref_bp$label.fine,
                  restrict = nocycle.markers)
res_bp_main <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_bp, labels = ref_bp$label.main,
                       restrict = nocycle.markers)


# transferencia de etiquetas de res (singleR) a metadatos de objeto seurat (columna nueva)
ann$singleR_dice <- res_dice$labels
ann$singleR_dice_main <- res_dice_main$labels

ann$singleR_mi <- res_mi$labels
ann$singleR_mi_main <- res_mi_main$labels

ann$singleR_nhd <- res_nhd$labels
ann$singleR_nhd_main <- res_nhd_main$labels

lab_map <- setNames(res_dice_clust$labels, rownames(res_dice_clust))
ann$singleR_dice_clust <- unname(lab_map[as.character(ann$seurat_clusters)])
rm(lab_map)

lab_map <- setNames(res_dice_main_clust$labels, rownames(res_dice_main_clust))
ann$singleR_dice_main_clust <- unname(lab_map[as.character(ann$seurat_clusters)])
rm(lab_map)

ann$singleR_bp <- res_bp$labels
ann$singleR_bp_main <- res_bp_main$labels


# dimplots con los nuevos metadatos de las anotaciones
p_dice <- DimPlot(ann, reduction = "umap", group.by = "singleR_dice", order = T)
ggsave("Results/Figures/03_01_ann_dice.png", p_dice, width = 8, height = 6, dpi = 300)

p_dice_main <- DimPlot(ann, reduction = "umap", group.by = "singleR_dice_main", order = T)
ggsave("Results/Figures/03_01_ann_dice_main.png", p_dice_main, width = 8, height = 6, dpi = 300)

p_mi <- DimPlot(ann, reduction = "umap", group.by = "singleR_mi", order = T)
ggsave("Results/Figures/03_01_ann_mi.png", p_mi, width = 8, height = 6, dpi = 300)

p_mi_main <- DimPlot(ann, reduction = "umap", group.by = "singleR_mi_main", order = T)
ggsave("Results/Figures/03_01_ann_mi_main.png", p_mi_main, width = 8, height = 6, dpi = 300)

p_nhd <- DimPlot(ann, reduction = "umap", group.by = "singleR_nhd", order = T)
ggsave("Results/Figures/03_01_ann_nhd.png", p_nhd, width = 8, height = 6, dpi = 300)

p_nhd_main <- DimPlot(ann, reduction = "umap", group.by = "singleR_nhd_main", order = T)
ggsave("Results/Figures/03_01_ann_nhd_main.png", p_nhd_main, width = 8, height = 6, dpi = 300)

p_dice_clust <- DimPlot(ann, reduction = "umap", group.by = "singleR_dice_clust", order = T)
ggsave("Results/Figures/03_01_ann_dice_clust.png", p_dice_clust, width = 8, height = 6, dpi = 300)

p_dice_main_clust <- DimPlot(ann, reduction = "umap", group.by = "singleR_dice_main_clust", order = T)
ggsave("Results/Figures/03_01_ann_dice_main_clust.png", p_dice_main_clust, width = 8, height = 6, dpi = 300)

p_bp <- DimPlot(ann, reduction = "umap", group.by = "singleR_bp", order = T)
ggsave("Results/Figures/03_01_ann_bp.png", p_bp, width = 8, height = 6, dpi = 300)

p_bp_main <- DimPlot(ann, reduction = "umap", group.by = "singleR_bp_main", order = T)
ggsave("Results/Figures/03_01_ann_bp_main.png", p_bp_main, width = 8, height = 6, dpi = 300)


## ---- 1.2 SingleR basado en referencias single-cell externas ----

library(scRNAseq)
library(scuttle)

all.ds <- listDatasets()
all.ds$Part

#"peripheral blood mononuclear cells" x4 --> 1 es de ratón
#"T cells"
#"T helper cells" --> es de ratón

sub.ds <- all.ds[
  grepl("9606", all.ds$Taxonomy) & (grepl("T cells", all.ds$Part) | grepl("T helper cells", all.ds$Part) |
        grepl("peripheral blood mononuclear cells", all.ds$Part)),]
sub.ds


# A. singleR con referencias externas de scRNAseq (Biconductor): T cells
Tcells <- BacherTCellData()
Tcells <- logNormCounts(Tcells)

sce <- as.SingleCellExperiment(ann, assay = "SCT")
assay(sce, "logcounts") <- NULL
sce <- logNormCounts(sce)

res_tcells <- SingleR(test = sce, ref = Tcells,
                      labels = Tcells$new_cluster_names, de.method = "wilcox",
                      assay.type.test = "logcounts", assay.type.ref = "logcounts",
                      restrict = nocycle.markers)
ann$singleR_tcells <- res_tcells$labels

res_tcells_clust <- SingleR(test = sce, ref = Tcells,
                      labels = Tcells$new_cluster_names, de.method = "wilcox",
                      assay.type.test = "logcounts", assay.type.ref = "logcounts",
                      restrict = nocycle.markers, clusters = ann$seurat_clusters)
lab_map <- setNames(res_tcells_clust$labels, rownames(res_tcells_clust))
ann$singleR_tcells_clust <- unname(lab_map[as.character(ann$seurat_clusters)])
rm(lab_map)

p_tcells <- DimPlot(ann, reduction = "umap", group.by = "singleR_tcells", order = TRUE)
ggsave("Results/Figures/03_01_ann_tcells.png", p_tcells, width = 8, height = 6, dpi = 300)

p_tcells_clust <- DimPlot(ann, reduction = "umap", group.by = "singleR_tcells_clust", order = T)
ggsave("Results/Figures/03_01_ann_tcells_clust.png", p_tcells_clust, width = 8, height = 6, dpi = 300)


# lista de anotaciones SingleR
singler_res_list <- list(dice = res_dice, dice_main = res_dice_main, mi = res_mi,
                         mi_main = res_mi_main, nhd = res_nhd, nhd_main = res_nhd_main, dice_clust = res_dice_clust,
                         dice_main_clust = res_dice_main_clust, bp = res_bp, bp_main = res_bp_main, tcells = res_tcells,
                         tcells_clust = res_tcells_clust)

# tablas abundancia etiqueta por cluster
#table(ann$seurat_clusters, unname(setNames(res$labels, rownames(res))[as.character(ann$seurat_clusters)]))

for (nm in names(singler_res_list)) {
  res <- singler_res_list[[nm]]
  if (length(res$labels) == length(ann$seurat_clusters)) {lab <- res$labels}
  else {lab <- unname(setNames(res$labels, rownames(res))[as.character(ann$seurat_clusters)])}
  tab <- table(cluster = ann$seurat_clusters, label = lab)
  write.table(as.data.frame.matrix(tab),paste0("Results/Tables/03_01_ann_", nm, "_cluster_vs_label.tsv"),
              sep = "\t", quote = FALSE, col.names = NA)
}

## ---- ann diagnostics ----
res_dice$scores[1:10,]

# heatmap del score
for (nm in names(singler_res_list)) {
  p_qc <- plotScoreHeatmap(singler_res_list[[nm]])
  ggsave(paste0("Results/Figures/03_01_ann_", nm, "_qc.png"), p_qc,
         width = 10, height = 6, dpi = 300)
}

# distribución delta
for (nm in names(singler_res_list)) {
  p_delta <- plotDeltaDistribution(singler_res_list[[nm]])
  ggsave(paste0("Results/Figures/03_01_ann_", nm, "_qc-delta.png"), p_delta,
         width = 10, height = 6, dpi = 300)
}


# comparing to unsupervised clustering
library(bluster)
for (nm in names(singler_res_list)) {
  res <- singler_res_list[[nm]]
  if (length(res$labels) != length(ann$seurat_clusters)) next
  ari <- pairwiseRand(ann$seurat_clusters, res$labels, mode = "index")
  tab <- table(cluster = ann$seurat_clusters, label = res$labels)
  hm <- pheatmap::pheatmap(log10(tab + 10))
  ggsave(paste0("Results/Figures/03_01_ann_", nm, "_qc-comp.png"),hm,
         width = 7, height = 10, dpi = 300)
}
rm(list = ls(pattern = "^p_"), hm, res, ari, nm, tab)


# marcadores por etiqueta
#p_dice_markers <- plotMarkerHeatmap(res_dice, test = ann_sce)
#ggsave("Results/Figures/03_01_ann_dice_qc-markers.png", p_dice_markers, width = 10, height = 6, dpi = 300)


#library(bluster)
#pairwiseRand(ann$seurat_clusters, res_dice$labels, mode="index")
#tab <- table(cluster=ann$seurat_clusters, label=res_dice$labels) 
#hm_comp <- pheatmap::pheatmap(log10(tab+10)) # using a larger pseudo-count for smoothing. 
#ggsave("Results/Figures/03_01_ann_dice_qc-comp.png", hm_comp, width = 7, height = 10, dpi = 300)


# adjusted rand index (ARI) debe ser > 0.5
#[1] 0.01708637   


#to.remove <- pruneScores(res_dice, min.diff.med=0.2)
#table(Label=res_dice$labels, Removed=to.remove)
#                                  Removed
# Label                              TRUE
# NK cells                           75
# T cells, CD4+, memory TREG       3815
# T cells, CD4+, naive, stimulated   22
# T cells, CD4+, TFH                 47
# T cells, CD4+, Th1                701
# T cells, CD4+, Th1_17               3
# T cells, CD4+, Th17               525
# T cells, CD4+, Th2                117
# T cells, CD8+, naive, stimulated    7

#to.remove_clust <- pruneScores(res_dice_clust, min.diff.med=0.2)
#table(Label=res_dice_clust$labels, Removed=to.remove_clust)
#                             Removed
#Label                        TRUE
# T cells, CD4+, memory TREG    4
# T cells, CD4+, Th1            4
# T cells, CD4+, Th17           3


#tab <- table(cluster=sceG$cluster, label=pred.grun$labels) 
#pheatmap::pheatmap(log10(tab+10)) # using a larger pseudo-count for smoothing. 


# score heatmap
#p_tcells_qc <- plotScoreHeatmap(res_tcells)
#ggsave("Results/Figures/03_01_ann_tcells_qc.png", p_tcells_qc, width = 10, height = 6, dpi = 300)

# delta distribution
#p_tcells_delta <- plotDeltaDistribution(res_tcells)
#ggsave("Results/Figures/03_01_ann_tcells_qc-delta.png", p_tcells_delta, width = 10, height = 6, dpi = 300)

# pruning de asignaciones ambiguas
#to.remove_tcells <- pruneScores(res_tcells, min.diff.med = 0.2)
#table(Label = res_tcells$labels, Removed = to.remove_tcells)

# marker heatmap global
#p_tcells_markers <- plotMarkerHeatmap(res_tcells_t, test = sce, label = res_tcells_t)
#ggsave("Results/Figures/03_01_ann_tcells_qc-markers.png", p_tcells_markers, width = 10, height = 6, dpi = 300)

# comparación con clustering no supervisado
#library(bluster)
#pairwiseRand(ann$seurat_clusters, res_tcells$labels, mode = "index")
#tab_tcells <- table(cluster = ann$seurat_clusters, label = res_tcells$labels)
#hm_tcells <- pheatmap::pheatmap(log10(tab_tcells + 10))
#ggsave("Results/Figures/03_01_ann_tcells_qc-comp.png", hm_tcells, width = 7, height = 10, dpi = 300)



# -----> tras la larga espera, lo asigna todo como cycling (pero MKI67, TOP2A y MCM5 sólo están en clusters 8-9)
#FeaturePlot(ann, c("MKI67", "TOP2A", "MCM5"))  # genes de ciclo localizados en 2 clusters. no todo el dataset



# B. singleR con referencias externas de scRNAseq (Biconductor): pbmc1
#pbmc1 <- MairPBMCData()
#pbmc1 <- logNormCounts(pbmc1)

#res_pbmc1 <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = pbmc1, labels = pbmc1$Sample_Tag)
#ann$singleR_pbmc1 <- res_pmbc1$labels

#p_pbmc1 <- DimPlot(ann, reduction = "umap", group.by = "singleR_pbmc1", order = T)
#ggsave("Results/Figures/03_01_ann_pbmc1.png", p_pbmc1, width = 8, height = 6, dpi = 300)

# no encuentro anotaciones

#> head(rownames(pbmc1))
#[1] "ADAM15_RNA"     "ADAM17_RNA"     "ADAM28_PolyA_1" "ADA_RNA"        "ADGRE1_RNA"     "ADGRG3_RNA"    
#> head(rownames(Tcells))
#[1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3"  "AL627309.2" 
#> head(rownames(ref_bp))
#[1] "TSPAN6"   "TNMD"     "DPM1"     "SCYL3"    "C1orf112" "FGR"     
#> head(rownames(ref_dice))
#[1] "5S_rRNA"  "7SK"      "A1BG"     "A1BG-AS1" "A1CF"     "A2M"     
#> head(rownames(ref_mi))
#[1] "A1BG"     "A1BG-AS1" "A1CF"     "A2M"      "A2M-AS1"  "A2ML1"   
#> head(rownames(ref_nhd))
#[1] "13CDNA73" "15E1.2"   "2'-PDE"   "3.8-1"    "76P"      "A2BP1"   
#> head(rownames(ann))
#[1] "7SK"      "A1BG"     "A1BG-AS1" "AAAS"     "AACS"     "AAGAB"   
#> length(intersect(rownames(ann), rownames(ref_dice))) / length(rownames(ann))
#[1] 0.7796013
#> length(intersect(rownames(ann), rownames(Tcells))) / length(rownames(ann))
#[1] 0.8620598
#> length(intersect(rownames(ann), rownames(pbmc1))) / length(rownames(ann))
#[1] 0
#> length(intersect(rownames(ann), rownames(ref_mi))) / length(rownames(ann))
#[1] 0.9533555
#> length(intersect(rownames(ann), rownames(ref_nhd))) / length(rownames(ann))
#[1] 0.4230565


# B. singleR con referencias externas de scRNAseq (Biconductor): pbmc2
#pbmc2 <- KotliarovPBMCData()
#pbmc2 <- logNormCounts(pbmc2)

#res_pbmc2 <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = pbmc2, labels = pbmc2$)
#ann$singleR_pbmc2 <- res_pmbc2$labels

#p_pbmc2 <- DimPlot(ann, reduction = "umap", group.by = "singleR_pbmc2", order = T)
#ggsave("Results/Figures/03_01_ann_pbmc2.png", p_pbmc2, width = 8, height = 6, dpi = 300)
# no encuentro anotaciones

# B. singleR con referencias externas de scRNAseq (Biconductor): pbmc3
#pbmc3 <- StoeckiusHashingData()
#pbmc3 <- logNormCounts(pbmc3)

#res_pbmc3 <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = pbmc3, labels = pbmc3$sizeF)
# no encuentro anotaciones



# Users can also supply their own custom marker lists to SingleR(), facilitating
#incorporation of prior biological knowledge into the annotation process. 



# ---- 1.3 SingleR con atlas PBMCs multimodal ----

pbmc <- readRDS("Ann/pbmc_multimodal_2023.rds")
pbmc_ref <- as.SingleCellExperiment(pbmc, assay = "SCT")

query <- as.SingleCellExperiment(ann, assay = "SCT")

pbmc1 <- SingleR(test = query, ref = pbmc_ref, labels = pbmc_ref$celltype.l1,
                  assay.type.test = "logcounts", assay.type.ref = "logcounts",
                  de.method = "wilcox")

pbmc2 <- SingleR(test = query, ref = pbmc_ref, labels = pbmc_ref$celltype.l2,
                  assay.type.test = "logcounts", assay.type.ref = "logcounts",
                  de.method = "wilcox")

pbmc3 <- SingleR(test = query, ref = pbmc_ref, labels = pbmc_ref$celltype.l3,
                  assay.type.test = "logcounts", assay.type.ref = "logcounts",
                  de.method = "wilcox")

## !! ha tardado horas en local
save.image(file = "Proc/03a_singleR_multimodal_PBMC.RData")

ann$pbmc1 <- pbmc1$labels
ann$pbmc2 <- pbmc2$labels
ann$pbmc3 <- pbmc3$labels

p_pbmc1 <- DimPlot(ann, reduction = "umap", group.by = "pbmc1", order = TRUE)
ggsave("Results/Figures/03_01_ann_pbmc1.png", p_pbmc1, width = 8, height = 6, dpi = 300)

p_pbmc2 <- DimPlot(ann, reduction = "umap", group.by = "pbmc2", order = TRUE)
ggsave("Results/Figures/03_01_ann_pbmc2.png", p_pbmc2, width = 8, height = 6, dpi = 300)

p_pbmc3 <- DimPlot(ann, reduction = "umap", group.by = "pbmc3", order = TRUE)
ggsave("Results/Figures/03_01_ann_pbmc3.png", p_pbmc3, width = 10, height = 6, dpi = 300)


## !! ha tardado horas en local

# filtrando fuera anotaciones no CD4






# restringiendo a marcadores identificados con FindConservedMarkers
# restrict = nocycle.markers





# restringiendo a marcadores + filtrando fuera anotaciones no CD4





# restringiendo a marcadores + filtrando fuera genes de ciclo celular





# restringiendo a marcadores + filtrando etiquetas no CD4 + filtrando genes cico celular


saveRDS(ann, "Proc/03_01_ann_pbmc.rds")
#ann <- readRDS("Proc/03_01_ann_pbmc.rds")


## =================================================
## BLOQUE Z - MAPEO A ATLAS CD4 (BASADO EN SYMPHONY)
## =================================================
## ----library and source-----
library(Seurat)
library(SeuratData)
#library(Azimuth)
library(patchwork)
library(tidyverse)
library(sctransform)
library(Matrix)


source('Ann/screfmapping/ref_mapping_seuratobj.R')
source('Ann/screfmapping/utils_seurat.R')

## ----Load　reference----
# Symphony
load("Ann/screfmapping/data/ref_Reference_Mapping_20220525.RData")
file.copy(from = 'Ann/screfmapping/data/cache_symphony_sct.uwot', 
          to = 'cache_symphony_sct.uwot',
          overwrite = FALSE)

## ----parameter setting (change here)----
#project.name <- "TFM"
#prefix <- paste0("./output/", project.name, "/", project.name)
#dir.create(paste0("./output/", project.name), recursive = T)

## ----reading data (raw seurat) ----
q <- readRDS("Proc/01_data_prep_and_QC.rds")

## ----run Symphony----
reference_mapping_seuratobj(ref, q, "prueba/")

## anotación
pred <- read.csv("prueba/_Reference_Mapping.csv", row.names = 1)

ann$clusterL1_ref      <- pred[colnames(ann), "clusterL1"]
ann$clusterL1_prob_ref <- pred[colnames(ann), "clusterL1_prob"]
ann$clusterL2_ref      <- pred[colnames(ann), "clusterL2"]
ann$clusterL2_prob_ref <- pred[colnames(ann), "clusterL2_prob"]

p_L1 <- DimPlot(ann, reduction = "umap", group.by = "clusterL1_ref", order = T)
ggsave("Results/Figures/03_02_ann_cd4_L1.png", p_L1, width = 8, height = 6, dpi = 300)

p_L2 <- DimPlot(ann, reduction = "umap", group.by = "clusterL2_ref", order = T)
ggsave("Results/Figures/03_02_ann_cd4_L2.png", p_L2, width = 8, height = 6, dpi = 300)


saveRDS(ann, "Proc/03_01_ann_cd4.rds")


## ======================================================
## BLOQUE 2 - ANOTACIÓN AUTOMÁTICA: ScType con ScTypeDB
## ======================================================
library(Seurat)
library(HGNChelper)
library(openxlsx)
library(dplyr)
library(ggplot2)

## ---- 2.1 Cargar funciones, bases de datos, tejido y genes  ----
cellmarkerDB_hu <- "Ann/Cell_marker_Human.xlsx"
cellmarkerDB_sc <- "Ann/Cell_marker_Seq.xlsx"
sctypeDB <- "Ann/ScTypeDB_full.xlsx"

# scType is not on CRAN or Bioconductor — functions directly from GitHub
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
#   gene_sets_prepare() - loads and processes ScTypeDB for your tissue
#   sctype_score() - scores your Seurat object against the DB

# función auxiliar
sctype_results <- function(obj_seurat, es_max) {
  clusters <- unique(obj_seurat$seurat_clusters)
  summary <- lapply(clusters, function(cl) {
    # Get cell barcodes for this cluster
    cell_ids <- WhichCells(obj_seurat, idents = cl)
    # Sum scores across cells in cluster
    if (length(cell_ids) == 1) {
      cluster_scores <- es_max[, cell_ids, drop = FALSE]
      scores_sum <- rowSums(cluster_scores)
    } else {
      cluster_scores <- es_max[, cell_ids]
      scores_sum <- rowSums(cluster_scores)
    }
    # Sort and take top 3
    top3 <- sort(scores_sum, decreasing = TRUE)[1:3]
    # Confidence = gap between 1st and 2nd normalised by total
    confidence <- ifelse(
      length(top3) >= 2 && top3[2] > 0,
      round((top3[1] - top3[2]) / abs(top3[1]) * 100, 1),
      100
    )
    data.frame(
      cluster = cl,
      rank1_type = names(top3)[1],
      rank1_score = round(top3[1], 2),
      rank2_type = names(top3)[2],
      rank2_score = round(top3[2], 2),
      rank3_type = names(top3)[3],
      rank3_score = round(top3[3], 2),
      confidence_gap = confidence,
      stringsAsFactors = FALSE
    )
  }) %>%
    dplyr::bind_rows()
  return(summary)
}


# ScTypeDB is an Excel file hosted on GitHub
#"https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
# Columns: tissueType | CellName | geneSymbolmore1 (positive) | geneSymbolmore2 (negative)

# available tissues in the DB
db_raw <- read.xlsx(sctypeDB, sheet = 1)
print(unique(db_raw$tissueType))

# Preparar genes para un tejido 
tissue_type <- "Immune system"
gs_list <- gene_sets_prepare(sctypeDB, tissue_type)

# components:
print(names(gs_list$gs_positive))
print(gs_list$gs_positive[["Naive CD4+ T cells"]])
print(gs_list$gs_negative[["Naive CD4+ T cells"]])

## ---- 2.2 Score clusters with scType  ----

# sctype_score() needs:
#     scRNAseqData : scaled expression matrix (genes x cells)
#     scaled       : TRUE if already scaled
#     gs           : positive marker gene sets
#     gs2          : negative marker gene sets

# Extract scaled expression matrix from Seurat
scaled_matrix <- GetAssayData(ann, assay = "SCT", layer = "scale.data")
dim(scaled_matrix)
scaled_matrix[1:3,1:3]

# Run scType scoring
es_max <- sctype_score(scaled_matrix, scaled = T, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
# es_max is a matrix: cell types x cells
# Each value is the scType score for that cell type in that cell
dim(es_max)      # rows = cell types, cols = cells
head(es_max[,1:5])


## ---- 2.3 Assign cell types per cluster  ----

#  For each Seurat cluster:
#    - Sum scType scores across all cells in the cluster
#    - Cell type with highest sum = cluster annotation
#    - Also compute a confidence: top score / (top score + 2nd score)

sctype.res <- sctype_results(ann, es_max)

## ---- 2.4 Add scType labels to Seurat object ----
ann$scType1 <- sctype.res$rank1_type
ann$scType2 <- sctype.res$rank2_type
ann$scType3 <- sctype.res$rank3_type

## ---- 2.5 Figuras ----
p_scType1 <- DimPlot(ann, reduction = "umap", group.by = "scType1", order = T)
ggsave("Results/Figures/03_02_ann_scType1.png", p_scType1, width = 8, height = 6, dpi = 300)

p_scType2 <- DimPlot(ann, reduction = "umap", group.by = "scType2", order = T)
ggsave("Results/Figures/03_02_ann_scType2.png", p_scType2, width = 8, height = 6, dpi = 300)

p_scType3 <- DimPlot(ann, reduction = "umap", group.by = "scType3", order = T)
ggsave("Results/Figures/03_02_ann_scType3.png", p_scType3, width = 8, height = 6, dpi = 300)


## ================================================
## BLOQUE 4 - ANOTACIÓN AUTOMÁTICA: EasyCellType
## ================================================

library(EasyCellType)
library(stats)
library(ggplot2)

markers.pos <- read.delim("Results/Tables/02_03_FindAllMarkers_pos_clust.tsv", sep = "\t")

library("org.Hs.eg.db")
library(AnnotationDbi)
markers.pos$entrezid <- mapIds(org.Hs.eg.db,
                           keys=markers.pos$gene, #Column containing Ensembl gene ids
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
markers.pos <- na.omit(markers.pos)

library(dplyr)
markers_sort <- data.frame(gene=markers.pos$entrezid, cluster=markers.pos$cluster, 
                           score=markers.pos$avg_log2FC) %>% 
  group_by(cluster) %>% 
  mutate(rank = rank(score),  ties.method = "random") %>% 
  arrange(desc(rank)) 
input.d <- as.data.frame(markers_sort[, 1:3])
write.csv(input.d,"Results/Tables/03_03_easycelltype_input.csv", row.names = F, quote = F)

annot.GSEA <- easyct(input.d, db="cellmarker", species="Human", 
                     tissue=c("Blood", "Peripheral blood"), p_cut=0.3,
                     test="GSEA", scoretype = "pos")

annot.GSEA_panglao <- easyct(input.d, db="panglao", species="Human", 
                     tissue=c("Blood", "Immune system"), p_cut=0.3,
                     test="GSEA", scoretype = "pos")

annot.GSEA_mole <- easyct(input.d, db="clustermole", species="Human", 
                             tissue=c("Blood", "Lymph node", "Peripheral blood", "Immune system"), p_cut=0.3,
                             test="GSEA", scoretype = "pos")


GSEA_dot <- plot_dot(test="GSEA", annot.GSEA)
ggsave("Results/Figures/03_03_ann_easycelltype_dot.png",GSEA_dot, width = 8, height = 6, dpi = 300)

pdf("Results/Figures/03_03_ann_easycelltype_bar.pdf", width = 8, height = 6)
plot_bar(test = "GSEA", annot.GSEA)
dev.off()


GSEA_dot_panglao <- plot_dot(test = "GSEA", annot.GSEA_panglao)
ggsave("Results/Figures/03_03_ann_easycelltype_dot_panglao.png", GSEA_dot_panglao, width = 8, height = 6, dpi = 300)

pdf("Results/Figures/03_03_ann_easycelltype_bar_panglao.pdf", width = 8, height = 6)
plot_bar(test = "GSEA", annot.GSEA_panglao)
dev.off()


GSEA_dot_mole <- plot_dot(test = "GSEA", annot.GSEA_mole)
ggsave("Results/Figures/03_03_ann_easycelltype_dot_mole.png", GSEA_dot_mole, width = 8, height = 6, dpi = 300)

pdf("Results/Figures/03_03_ann_easycelltype_bar_mole.pdf", width = 8, height = 6)
plot_bar(test = "GSEA", annot.GSEA_mole)
dev.off()



coremarkers("GSEA", annot.GSEA, "Human")

data(markers.pos)

process_results("GSEA", annot.GSEA)

summarycelltype("GSEA", annot.GSEA, 1)


saveRDS(ann, "Proc/03_sct_int_ann2.rds")



cat("Porcentaje respecto al objeto padre:", pct_sub, "%\n")




# funcion para generar todos los featureplot
#fp_summary <- function(clust_obj, tab, gene_col) {
#  for (i in seq_len(nrow(tab))) {
#    cluster_id <- tab$cluster[i]
#    pct_label <- sub(".*_(\\d+)$", "\\1", gene_col)
#    gene_string <- tab[[gene_col]][i]
#    lista <- trimws(unlist(strsplit(gene_string, ",")))
#    p_clust <- FeaturePlot(clust_obj, features = lista, reduction = "umap", cols = c("grey", "blue"), order = T, ncol = 4)
#    ggsave(paste0("results/figures/03_01_markers_cluster", cluster_id, "_up_", pct_label, ".png"), p_clust, width = 14, height = 9, dpi = 300)
#  }
#}

#fp_summary(clust, cons.genelist, "gene_up_05")
#fp_summary(clust, cons.genelist, "gene_up_09")
#fp_summary(clust, cons.genelist, "gene_up_01")


#################################




####################################

ann <- readRDS("Proc/02_02_sct_int-rpca_leiden_1.rds")

library(Seurat)
library(ggplot2)
library(scGate)

#Manually define a simple scGate gating model to purify
#eg. Natural Killer (NK) cells, using a positive marker KLRD1 and negative marker CD3D
my_scGate_model <- gating_model(name = "CCR", signature = c("CD3D","CD3E","CD3G"))  

#scGate it!
ann <- scGate(data = ann, model = my_scGate_model)

#Use Seurat to visualize "Pure" and "Impure" cells
DimPlot(ann, group.by = "is.pure")

#Use Seurat to subset pure cells
seurat_object_purified <- subset(seurat_object, subset = `is.pure` == "Pure" )


#Get scGate database of pre-defined gating models
scGate_models_DB <- get_scGateDB()

#For example, filter abT cells using one of scGate pre-defined gating models
ann <- scGate(ann, model = scGate_models_DB$human$generic$Tcell.alphabeta)

DimPlot(ann)

ann <- scGate(ann, model = scGate_models_DB$human$generic$CD4T)

DimPlot(ann)



scGate::plot_tree(scGate_models_DB$human$PBMC$)


# scGate as a multi-class classifier
models.list <- scGate_models_DB$human$generic[c("Bcell","MoMacDC","CD8T","CD4T","Erythrocyte")]
ann <- scGate(ann, model = models.list)
DimPlot(ann)


models.TME <- scGate_models_DB$human$TME_HiRes
obj <- scGate(obj, model=models.TME)
table(ann$scGate_multi)

