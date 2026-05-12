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
# Script: 03a_ann_cell_auto.R
## ===============================

# Objetivo:
# Realizar anotación celular automática mediante SingleR, ScType, CellTypist y
# mapeo a un atlas CD4 basado en Symphony, incorporando las etiquetas al objeto
# Seurat integrado y generando figuras y tablas de evaluación por cluster.

# Inputs en Proc/:
# - 02_02_sct_int-rpca_leiden_1.rds
# - 01_data_prep_and_QC.rds (para mapeo con Symphony)

# Inputs en Ann/:
# - pbmc_multimodal_2023.rds (https://zenodo.org/records/7779017/files/pbmc_multimodal_2023.rds)
# - ScTypeDB_full.xlsx
# - PanglaoDB_markers_all.tsv
# - PanglaoDB_markers_Tcells.tsv
# - Cell_marker_Human.xlsx
# - Cell_marker_Seq.xlsx
# - celltypist_metadata.csv
# - screfmapping/ (https://github.com/yyoshiaki/screfmapping)
# - screfmapping/data/ref_Reference_Mapping_20220525.RData (https://doi.org/10.6084/m9.figshare.25052648)

# Outputs en Proc/:
# - 03_01_ann_singleR.rds
# - 03_02_int_ann.rds

# Outputs en Results/:
# - Figures/03_*.png
# - Tables/03_*.tsv




library(Seurat)
library(SingleCellExperiment)

# inputs from script 02b_scRNAseq_basic.R
ann <- readRDS("Proc/02_02_sct_int-rpca_leiden_1.rds")
DefaultAssay(ann) <- "SCT"

## ==========================================
## BLOQUE 1 - ANOTACIÓN AUTOMÁTICA: SingleR
## ==========================================
library(SingleR)
library(celldex)
library(scrapper)

## ---- 1.1 SingleR modo clásico: bases de datos bulk purificado por tipo celular ----

# bases de datos built-in
ref_dice <- celldex::DatabaseImmuneCellExpressionData()
ref_mi <- celldex::MonacoImmuneData()
ref_nhd <- celldex::NovershternHematopoieticData()
ref_bp <- celldex::BlueprintEncodeData()

# dataframe de resultados de las etiquetas celulares mapeadas
# el objeto seurat se tiene que pasar como SingleCellExperiment
# observar df ref y las labels con colData(ref)
res_dice <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_dice, labels = ref_dice$label.fine)
res_dice_main <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_dice, labels = ref_dice$label.main)

res_mi <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_mi, labels = ref_mi$label.fine)
res_mi_main <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_mi, labels = ref_mi$label.main)

res_nhd <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_nhd, labels = ref_nhd$label.fine)
res_nhd_main <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_nhd, labels = ref_nhd$label.main)

res_bp <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_bp, labels = ref_bp$label.fine)
res_bp_main <- SingleR(test = as.SingleCellExperiment(ann, assay = "SCT"), ref = ref_bp, labels = ref_bp$label.main)


# transferencia de etiquetas de res (singleR) a metadatos de objeto seurat (columna nueva)
ann$SR_dice <- res_dice$labels
ann$SR_dice_main <- res_dice_main$labels

ann$SR_mi <- res_mi$labels
ann$SR_mi_main <- res_mi_main$labels

ann$SR_nhd <- res_nhd$labels
ann$singleR_nhd_main <- res_nhd_main$labels

ann$SR_bp <- res_bp$labels
ann$SR_bp_main <- res_bp_main$labels


# dimplots con los nuevos metadatos de las anotaciones
p_dice <- DimPlot(ann, reduction = "umap", group.by = "SR_dice", order = T)
ggsave("Results/Figures/03_01_ann_dice.png", p_dice, width = 8, height = 6, dpi = 300)

p_dice_main <- DimPlot(ann, reduction = "umap", group.by = "SR_dice_main", order = T)
ggsave("Results/Figures/03_01_ann_dice_main.png", p_dice_main, width = 8, height = 6, dpi = 300)

p_mi <- DimPlot(ann, reduction = "umap", group.by = "SR_mi", order = T)
ggsave("Results/Figures/03_01_ann_mi.png", p_mi, width = 8, height = 6, dpi = 300)

p_mi_main <- DimPlot(ann, reduction = "umap", group.by = "SR_mi_main", order = T)
ggsave("Results/Figures/03_01_ann_mi_main.png", p_mi_main, width = 8, height = 6, dpi = 300)

p_nhd <- DimPlot(ann, reduction = "umap", group.by = "SR_nhd", order = T)
ggsave("Results/Figures/03_01_ann_nhd.png", p_nhd, width = 8, height = 6, dpi = 300)

p_nhd_main <- DimPlot(ann, reduction = "umap", group.by = "SR_nhd_main", order = T)
ggsave("Results/Figures/03_01_ann_nhd_main.png", p_nhd_main, width = 8, height = 6, dpi = 300)

p_bp <- DimPlot(ann, reduction = "umap", group.by = "SR_bp", order = T)
ggsave("Results/Figures/03_01_ann_bp.png", p_bp, width = 8, height = 6, dpi = 300)

p_bp_main <- DimPlot(ann, reduction = "umap", group.by = "SR_bp_main", order = T)
ggsave("Results/Figures/03_01_ann_bp_main.png", p_bp_main, width = 8, height = 6, dpi = 300)


## ---- 1.2 SingleR basado en referencias single-cell externas: BatcherTCelldata (scRNAseq) ----
library(scRNAseq)
library(scuttle)

# inspección de bases de datos del paquete scRNAseq de Bioconductor
all.ds <- listDatasets()
all.ds$Part

sub.ds <- all.ds[
  grepl("9606", all.ds$Taxonomy) & (grepl("T cells", all.ds$Part) | grepl("T helper cells", all.ds$Part) |
        grepl("peripheral blood mononuclear cells", all.ds$Part)),]
sub.ds

#"peripheral blood mononuclear cells" x4 --> 1 es de ratón, de los otros no encuentro anotaciones
#"T helper cells" --> es de ratón
#"T cells" -> BatcherTCelldata

# BatcherTCelldata
Tcells <- BacherTCellData()
Tcells <- logNormCounts(Tcells)

sce <- as.SingleCellExperiment(ann, assay = "SCT")
assay(sce, "logcounts") <- NULL
sce <- logNormCounts(sce)

res_tcells <- SingleR(test = sce, ref = Tcells,
                      labels = Tcells$new_cluster_names, de.method = "wilcox",
                      assay.type.test = "logcounts", assay.type.ref = "logcounts")
ann$SR_tcells <- res_tcells$labels

p_tcells <- DimPlot(ann, reduction = "umap", group.by = "SR_tcells", order = TRUE)
ggsave("Results/Figures/03_01_ann_tcells.png", p_tcells, width = 8, height = 6, dpi = 300)


## ---- 1.3 SingleR basado en referencias single-cell externas: atlas PBMC multimodal ----
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

#save.image(file = "Proc/03a_singleR_multimodal_PBMC.RData")

ann$SR_pbmc1 <- pbmc1$labels
ann$SR_pbmc2 <- pbmc2$labels
ann$SR_pbmc3 <- pbmc3$labels

p_pbmc1 <- DimPlot(ann, reduction = "umap", group.by = "SR_pbmc1", order = TRUE)
ggsave("Results/Figures/03_01_ann_pbmc1.png", p_pbmc1, width = 8, height = 6, dpi = 300)

p_pbmc2 <- DimPlot(ann, reduction = "umap", group.by = "SR_pbmc2", order = TRUE)
ggsave("Results/Figures/03_01_ann_pbmc2.png", p_pbmc2, width = 8, height = 6, dpi = 300)

p_pbmc3 <- DimPlot(ann, reduction = "umap", group.by = "SR_pbmc3", order = TRUE)
ggsave("Results/Figures/03_01_ann_pbmc3.png", p_pbmc3, width = 10, height = 6, dpi = 300)


saveRDS(ann, "Proc/03_01_ann_singleR.rds")


## ---- 1.4 SingleR: diagnostics ----

# lista de anotaciones SingleR
singler_res_list <- list(dice = res_dice, dice_main = res_dice_main, mi = res_mi,
                         mi_main = res_mi_main, nhd = res_nhd, nhd_main = res_nhd_main,
                         bp = res_bp, bp_main = res_bp_main, tcells = res_tcells,
                         pbmc1 = pbmc1, pbmc2 = pbmc2, pbmc3 = pbmc3)

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

# marcadores
for (nm in names(singler_res_list)) {
  res <- singler_res_list[[nm]]
  p_markers <- plotMarkerHeatmap(res,test = sce,label = res$labels)
  ggsave(paste0("Results/Figures/03_01_ann_", nm, "_qc-markers.png"),
    p_markers, width = 10, height = 6, dpi = 300)
}

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



## ================================================================================
## BLOQUE 2 - ANOTACIÓN AUTOMÁTICA: ScType con ScTypeDB, PanglaoDB y CellMarker
## ================================================================================
library(Seurat)
library(HGNChelper)
library(openxlsx)
library(dplyr)
library(ggplot2)

## ---- 2.1 Cargar funciones ----

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
    cell_ids <- colnames(obj_seurat)[obj_seurat$seurat_clusters == cl]
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
    confidence <- ifelse(length(top3) >= 2 && top3[2] > 0,
      round((top3[1] - top3[2]) / abs(top3[1]) * 100, 1), 100)
    data.frame(cluster = cl, rank1_type = names(top3)[1],
      rank1_score = round(top3[1], 2), rank2_type = names(top3)[2],
      rank2_score = round(top3[2], 2), rank3_type = names(top3)[3],
      rank3_score = round(top3[3], 2), confidence_gap = confidence,
      stringsAsFactors = FALSE)
  }) %>%
    dplyr::bind_rows()
  return(summary)
}

## ---- 2.2 ScType con base de datos ScTypeDB ----

# base de datos (ScTypeDB is an Excel file hosted on GitHub)
#"https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
# Columns: tissueType | CellName | geneSymbolmore1 (positive) | geneSymbolmore2 (negative)
sctypeDB <- "Ann/ScTypeDB_full.xlsx"

# available tissues in the DB
print(unique(read.xlsx(sctypeDB, sheet = 1)$tissueType))

# Preparar genes para un tejido 
gs_list <- gene_sets_prepare(sctypeDB, "Immune system")

# componentes para CD4: naive, memory, effector, NKT-like. + ISG expressing cells
print(names(gs_list$gs_positive))
print(gs_list$gs_positive[["Naive CD4+ T cells"]])
print(gs_list$gs_negative[["Naive CD4+ T cells"]])

## Score clusters with scType

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


## Assign cell types per cluster

#  For each Seurat cluster:
#    - Sum scType scores across all cells in the cluster
#    - Cell type with highest sum = cluster annotation
#    - Also compute a confidence: top score / (top score + 2nd score)

sctype.res <- sctype_results(ann, es_max)

## Add scType labels to Seurat object 
Idents(ann) <- "seurat_clusters"

# mapas cluster -> etiqueta
map_scType1 <- setNames(sctype.res$rank1_type, sctype.res$cluster)
map_scType2 <- setNames(sctype.res$rank2_type, sctype.res$cluster)
map_scType3 <- setNames(sctype.res$rank3_type, sctype.res$cluster)
map_conf    <- setNames(sctype.res$confidence_gap, sctype.res$cluster)

# asignar a cada célula la etiqueta de su cluster
ann$scType1 <- unname(map_scType1[as.character(ann$seurat_clusters)])
ann$scType2 <- unname(map_scType2[as.character(ann$seurat_clusters)])
ann$scType3 <- unname(map_scType3[as.character(ann$seurat_clusters)])
ann$scType_conf <- unname(map_conf[as.character(ann$seurat_clusters)])

##   Figuras 
p_scType1 <- DimPlot(ann, reduction = "umap", group.by = "scType1", order = T)
ggsave("Results/Figures/03_02_ann_scType1.png", p_scType1, width = 8, height = 6, dpi = 300)

p_scType2 <- DimPlot(ann, reduction = "umap", group.by = "scType2", order = T)
ggsave("Results/Figures/03_02_ann_scType2.png", p_scType2, width = 10, height = 6, dpi = 300)

p_scType3 <- DimPlot(ann, reduction = "umap", group.by = "scType3", order = T)
ggsave("Results/Figures/03_02_ann_scType3.png", p_scType3, width = 8, height = 6, dpi = 300)


## ---- 2.3 ScType con base de datos ScTypeDB restringida a CD4 ----

cd4_labels <- c("Naive CD4+ T cells", "Memory CD4+ T cells", "Effector CD4+ T cells",
                "CD4+ NKT-like cells", "ISG expressing immune cells")

gs_list_cd4 <- list(gs_positive = gs_list$gs_positive[cd4_labels],
                    gs_negative = gs_list$gs_negative[cd4_labels])

es_max_cd4 <- sctype_score(scaled_matrix, scaled = T, gs = gs_list_cd4$gs_positive,
                           gs2 = gs_list_cd4$gs_negative)
sctype.res.cd4 <- sctype_results(ann, es_max_cd4)


Idents(ann) <- "seurat_clusters"

map_scType1_cd4 <- setNames(sctype.res.cd4$rank1_type, sctype.res.cd4$cluster)
map_scType2_cd4 <- setNames(sctype.res.cd4$rank2_type, sctype.res.cd4$cluster)
map_scType3_cd4 <- setNames(sctype.res.cd4$rank3_type, sctype.res.cd4$cluster)
map_conf_cd4    <- setNames(sctype.res.cd4$confidence_gap, sctype.res.cd4$cluster)

# asignar a cada célula la etiqueta de su cluster
ann$scType1_cd4 <- unname(map_scType1_cd4[as.character(ann$seurat_clusters)])
ann$scType2_cd4 <- unname(map_scType2_cd4[as.character(ann$seurat_clusters)])
ann$scType3_cd4 <- unname(map_scType3_cd4[as.character(ann$seurat_clusters)])
ann$scType_conf_cd4 <- unname(map_conf_cd4[as.character(ann$seurat_clusters)])

##   Figuras 
p_scType1_cd4 <- DimPlot(ann, reduction = "umap", group.by = "scType1_cd4", order = T)
ggsave("Results/Figures/03_02_ann_scType1_cd4.png", p_scType1_cd4, width = 8, height = 6, dpi = 300)

p_scType2_cd4 <- DimPlot(ann, reduction = "umap", group.by = "scType2_cd4", order = T)
ggsave("Results/Figures/03_02_ann_scType2_cd4.png", p_scType2_cd4, width = 10, height = 6, dpi = 300)

p_scType3_cd4 <- DimPlot(ann, reduction = "umap", group.by = "scType3_cd4", order = T)
ggsave("Results/Figures/03_02_ann_scType3_cd4.png", p_scType3_cd4, width = 8, height = 6, dpi = 300)


## ---- 2.4 ScType con base de datos PanglaoDB ----

##  Format needed:
##    gs_positive: named list of positive marker vectors
##    gs_negative: named list of negative marker vectors (can be empty)

# cargar PanglaoDB
panglaoDB <- read.delim("Ann/PanglaoDB_markers_all.tsv", stringsAsFactors = F) %>%
  filter(grepl("Hs", species), organ %in% c("Blood", "Immune system")) %>%
  dplyr::select(gene = official.gene.symbol, cell_type = cell.type, spec = specificity_human)

panglaoDB_t <- read.delim("Ann/PanglaoDB_markers_Tcells.tsv", stringsAsFactors = F) %>%
  filter(grepl("Hs", species), organ %in% c("Blood", "Immune system")) %>%
  dplyr::select(gene = official.gene.symbol, cell_type = cell.type, spec = specificity_human)

# Build positive marker list from PanglaoDB (no tiene negativos)
panglao_positive <- panglaoDB %>% group_by(cell_type) %>% summarise(genes = list(unique(gene)), .groups = "drop")
gs_panglao_positive <- setNames(panglao_positive$genes,panglao_positive$cell_type)

panglao_positive_t <- panglaoDB_t %>% group_by(cell_type) %>% summarise(genes = list(unique(gene)), .groups = "drop")
gs_panglao_positive_t <- setNames(panglao_positive_t$genes,panglao_positive_t$cell_type)

# Create empty negative list of same length
gs_panglao_negative <- lapply(gs_panglao_positive, function(x) c())
gs_panglao_negative_t <- lapply(gs_panglao_positive_t, function(x) c())

# Score with PanglaoDB markers using scType engine
es_panglao <- sctype_score(scaled_matrix, scaled = T, gs = gs_panglao_positive, gs2 = gs_panglao_negative)
sctype.res.panglao <- sctype_results(ann, es_panglao)
sctype.res.panglao

# Tranferir anotaciones
Idents(ann) <- "seurat_clusters"

map_scType1_panglao <- setNames(sctype.res.panglao$rank1_type, sctype.res.panglao$cluster)
ann$scType1_panglao <- unname(map_scType1_panglao[as.character(ann$seurat_clusters)])
p_scType1_panglao <- DimPlot(ann, reduction = "umap", group.by = "scType1_panglao", order = T)
ggsave("Results/Figures/03_02_ann_scType1_panglao.png", p_scType1_panglao, width = 8, height = 6, dpi = 300)

map_scType2_panglao <- setNames(sctype.res.panglao$rank2_type, sctype.res.panglao$cluster)
ann$scType2_panglao <- unname(map_scType2_panglao[as.character(ann$seurat_clusters)])
p_scType2_panglao <- DimPlot(ann, reduction = "umap", group.by = "scType2_panglao", order = T)
ggsave("Results/Figures/03_02_ann_scType2_panglao.png", p_scType2_panglao, width = 8, height = 6, dpi = 300)

map_scType3_panglao <- setNames(sctype.res.panglao$rank3_type, sctype.res.panglao$cluster)
ann$scType3_panglao <- unname(map_scType3_panglao[as.character(ann$seurat_clusters)])
p_scType3_panglao <- DimPlot(ann, reduction = "umap", group.by = "scType3_panglao", order = T)
ggsave("Results/Figures/03_02_ann_scType3_panglao.png", p_scType3_panglao, width = 8, height = 6, dpi = 300)


## ---- 2.5 ScType con base de datos PanglaoDB restringida a linfocitos T ----
es_panglao_t <- sctype_score(scaled_matrix, scaled = T, gs = gs_panglao_positive_t, gs2 = gs_panglao_negative_t)
sctype.res.panglao_t <- sctype_results(ann, es_panglao_t)

Idents(ann) <- "seurat_clusters"

map_scType1_panglao_t <- setNames(sctype.res.panglao_t$rank1_type, sctype.res.panglao_t$cluster)
ann$scType1_panglao_t <- unname(map_scType1_panglao_t[as.character(ann$seurat_clusters)])
p_scType1_panglao_t <- DimPlot(ann, reduction = "umap", group.by = "scType1_panglao_t", order = T)
ggsave("Results/Figures/03_02_ann_scType1_panglao_t.png", p_scType1_panglao_t, width = 8, height = 6, dpi = 300)

map_scType2_panglao_t <- setNames(sctype.res.panglao_t$rank2_type, sctype.res.panglao_t$cluster)
ann$scType2_panglao_t <- unname(map_scType2_panglao_t[as.character(ann$seurat_clusters)])
p_scType2_panglao_t <- DimPlot(ann, reduction = "umap", group.by = "scType2_panglao_t", order = T)
ggsave("Results/Figures/03_02_ann_scType2_panglao_t.png", p_scType2_panglao_t, width = 8, height = 6, dpi = 300)

map_scType3_panglao_t <- setNames(sctype.res.panglao_t$rank3_type, sctype.res.panglao_t$cluster)
ann$scType3_panglao_t <- unname(map_scType3_panglao_t[as.character(ann$seurat_clusters)])
p_scType3_panglao_t <- DimPlot(ann, reduction = "umap", group.by = "scType3_panglao_t", order = T)
ggsave("Results/Figures/03_02_ann_scType3_panglao_t.png", p_scType3_panglao_t, width = 8, height = 6, dpi = 300)


## ---- 2.6 ScType con base de datos Cell Marker ----
cellmarkerDB_hu <- "Ann/Cell_marker_Human.xlsx"
cellmarkerDB_sc <- "Ann/Cell_marker_Seq.xlsx"

cellmarker_hu <- read.xlsx("Ann/Cell_marker_Human.xlsx", sep = ",") %>%
  dplyr::select(tissue = tissue_type, cell_type = cell_name, gene = Symbol) %>%
  dplyr::filter(grepl("blood", tissue, ignore.case = TRUE)) %>%
  dplyr::filter(!is.na(gene), gene != "", gene != "NA") %>%
  dplyr::mutate(gene = strsplit(gene, "[,;]+")) %>% tidyr::unnest(gene) %>%
  dplyr::mutate(gene = trimws(gene)) %>% dplyr::filter(gene != "")

cellmarker_sc <- read.xlsx("Ann/Cell_marker_Seq.xlsx", sep = ",") %>%
  dplyr::select(tissue = tissue_type, cell_type = cell_name, gene = Symbol) %>%
  dplyr::filter(grepl("blood", tissue, ignore.case = TRUE)) %>%
  dplyr::filter(!is.na(gene), gene != "", gene != "NA") %>%
  dplyr::mutate(gene = strsplit(gene, "[,;]+")) %>% tidyr::unnest(gene) %>%
  dplyr::mutate(gene = trimws(gene)) %>% dplyr::filter(gene != "")

gs_cm_positive_hu <- cellmarker_hu %>% group_by(cell_type) %>%
  summarise(genes = list(unique(gene)), .groups = "drop") %>%
  { setNames(.$genes, .$cell_type) }

gs_cm_negative_hu <- lapply(gs_cm_positive_hu, function(x) c())

gs_cm_positive_sc <- cellmarker_sc %>% group_by(cell_type) %>%
  summarise(genes = list(unique(gene)), .groups = "drop") %>%
  { setNames(.$genes, .$cell_type) }

gs_cm_negative_sc <- lapply(gs_cm_positive_sc, function(x) c())

es_cellmarker_hu <- sctype_score(scaled_matrix, scaled = T, gs = gs_cm_positive_hu, gs2 = gs_cm_negative_hu)
es_cellmarker_sc <- sctype_score(scaled_matrix, scaled = T, gs = gs_cm_positive_sc, gs2 = gs_cm_negative_sc)

sctype.res.cellmarker_hu <- sctype_results(ann, es_cellmarker_hu)
sctype.res.cellmarker_sc <- sctype_results(ann, es_cellmarker_sc)

Idents(ann) <- "seurat_clusters"

map_scType1_cm_hu <- setNames(sctype.res.cellmarker_hu$rank1_type, sctype.res.cellmarker_hu$cluster)
ann$scType1_cm_hu <- unname(map_scType1_cm_hu[as.character(ann$seurat_clusters)])
p_scType1_cm_hu <- DimPlot(ann, reduction = "umap", group.by = "scType1_cm_hu", order = T)
ggsave("Results/Figures/03_02_ann_scType1_cm_hu.png", p_scType1_cm_hu, width = 8, height = 6, dpi = 300)

map_scType2_cm_hu <- setNames(sctype.res.cellmarker_hu$rank2_type, sctype.res.cellmarker_hu$cluster)
ann$scType2_cm_hu <- unname(map_scType2_cm_hu[as.character(ann$seurat_clusters)])
p_scType2_cm_hu <- DimPlot(ann, reduction = "umap", group.by = "scType2_cm_hu", order = T)
ggsave("Results/Figures/03_02_ann_scType2_cm_hu.png", p_scType2_cm_hu, width = 8, height = 6, dpi = 300)

map_scType3_cm_hu <- setNames(sctype.res.cellmarker_hu$rank3_type, sctype.res.cellmarker_hu$cluster)
ann$scType3_cm_hu <- unname(map_scType3_cm_hu[as.character(ann$seurat_clusters)])
p_scType3_cm_hu <- DimPlot(ann, reduction = "umap", group.by = "scType3_cm_hu", order = T)
ggsave("Results/Figures/03_02_ann_scType3_cm_hu.png", p_scType3_cm_hu, width = 8, height = 6, dpi = 300)

map_scType1_cm_sc <- setNames(sctype.res.cellmarker_sc$rank1_type, sctype.res.cellmarker_sc$cluster)
ann$scType1_cm_sc <- unname(map_scType1_cm_sc[as.character(ann$seurat_clusters)])
p_scType1_cm_sc <- DimPlot(ann, reduction = "umap", group.by = "scType1_cm_sc", order = T)
ggsave("Results/Figures/03_02_ann_scType1_cm_sc.png", p_scType1_cm_sc, width = 8, height = 6, dpi = 300)

map_scType2_cm_sc <- setNames(sctype.res.cellmarker_sc$rank2_type, sctype.res.cellmarker_sc$cluster)
ann$scType2_cm_sc <- unname(map_scType2_cm_sc[as.character(ann$seurat_clusters)])
p_scType2_cm_sc <- DimPlot(ann, reduction = "umap", group.by = "scType2_cm_sc", order = T)
ggsave("Results/Figures/03_02_ann_scType2_cm_sc.png", p_scType2_cm_sc, width = 8, height = 6, dpi = 300)

map_scType3_cm_sc <- setNames(sctype.res.cellmarker_sc$rank3_type, sctype.res.cellmarker_sc$cluster)
ann$scType3_cm_sc <- unname(map_scType3_cm_sc[as.character(ann$seurat_clusters)])
p_scType3_cm_sc <- DimPlot(ann, reduction = "umap", group.by = "scType3_cm_sc", order = T)
ggsave("Results/Figures/03_02_ann_scType3_cm_sc.png", p_scType3_cm_sc, width = 8, height = 6, dpi = 300)


## =================================================
## BLOQUE 3 - CELLTYPIST (REALIZADO EN PYTHON)
## =================================================

meta_ct <- read.csv("Ann/celltypist_metadata.csv", row.names = 1)

# comprobar solapamiento
all(rownames(meta_ct) == colnames(ann))
#[1] TRUE

colnames(meta_ct)[colnames(meta_ct) == "predicted_labels"] <- "celltypist_predicted"
colnames(meta_ct)[colnames(meta_ct) == "majority_voting"] <- "celltypist_majority"
colnames(meta_ct)[colnames(meta_ct) == "conf_score"] <- "celltypist_confscore"

# añadir metadatos
ann <- AddMetaData(ann, metadata = meta_ct)

p_maj <- DimPlot(ann, reduction = "umap", group.by = "majority_voting")
ggsave("results/figures/03_01_celltypist_maj.png", p_maj, width = 10, height = 8, dpi = 300)

p_pred <- DimPlot(ann, reduction = "umap", group.by = "predicted_labels")
ggsave("results/figures/03_01_celltypist_pred.png", p_pred, width = 15, height = 8, dpi = 300)


## =================================================
## BLOQUE 4 - MAPEO A ATLAS CD4 (BASADO SYMPHONY)
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
reference_mapping_seuratobj(ref, q, "Ann/")

## anotación
pred <- read.csv("Ann/_Reference_Mapping.csv", row.names = 1)

ann$clusterL1_ref      <- pred[colnames(ann), "clusterL1"]
ann$clusterL1_prob_ref <- pred[colnames(ann), "clusterL1_prob"]
ann$clusterL2_ref      <- pred[colnames(ann), "clusterL2"]
ann$clusterL2_prob_ref <- pred[colnames(ann), "clusterL2_prob"]

p_L1 <- DimPlot(ann, reduction = "umap", group.by = "clusterL1_ref", order = T)
ggsave("Results/Figures/03_02_ann_cd4_L1.png", p_L1, width = 8, height = 6, dpi = 300)

p_L2 <- DimPlot(ann, reduction = "umap", group.by = "clusterL2_ref", order = T)
ggsave("Results/Figures/03_02_ann_cd4_L2.png", p_L2, width = 8, height = 6, dpi = 300)


saveRDS(ann, "Proc/03_02_int_ann.rds")



## =================================================
## BLOQUE 5 - ANÁLISIS DE RESULTADOS
## =================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

ann <- readRDS("Proc/03_02_int_ann.rds")

## ---- Evaluación de anotaciones automáticas por cluster ----
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

## columnas de anotación automática a evaluar
ann_cols <- c("pbmc1", "pbmc2", "pbmc3","cd4_L1", "cd4_L2","SR50_dice", "SR50_dice_main",
  "SR50_bp", "SR50_bp_main", "SR50_nhd", "SR50_nhd_main","SR50_tcells","SR10_dice", 
  "SR10_nhd", "SR10_nhd_main", "SR10_bp", "SR10_bp_main","SR10_tcells","SR50_mi",
  "SR_dice", "SR_dice_main", "SR_mi", "SR_mi_main","SR_nhd", "SR_nhd_main", "SR_bp", 
  "celltypist_predicted", "celltypist_majority","scType1", "scType2", "scType3",
  "scType1_cd4", "scType2_cd4", "scType3_cd4","scType1_panglao", "scType2_panglao",
  "scType1_panglao_t", "scType2_panglao_t", "scType3_panglao_t","SR50_mi_main",
  "scType1_cm_hu", "scType2_cm_hu", "scType3_cm_hu","scType3_panglao","SR_bp_main", "SR_tcells",
  "scType1_cm_sc", "scType2_cm_sc", "scType3_cm_sc","SR10_dice_main", "SR10_mi", "SR10_mi_main"
)

meta_ann <- ann[[]] |> tibble::rownames_to_column("cell") |>
  dplyr::select(cell, seurat_clusters, dplyr::all_of(ann_cols)) |>
  dplyr::rename(cluster = seurat_clusters)

## función para stacked plot de una anotación concreta
plot_ann_by_cluster <- function(meta, ann_col, out_prefix = "03_02_ann_eval") {
  tab <- meta |> count(cluster, label = .data[[ann_col]], name = "n_cells") |>
    group_by(cluster) |> mutate(prop = n_cells / sum(n_cells)) |> ungroup()
  
  p <- ggplot(tab, aes(x = cluster, y = prop, fill = label)) + geom_col(width = 0.8) +
    scale_y_continuous(labels = percent_format()) + labs(x = "Cluster", y = "% células", fill = ann_col) +
    theme_bw() + theme(legend.position = "right", axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  ggsave(paste0("Results/Figures/", "03a_05_", out_prefix, "_", ann_col, "_stacked.png"),
    p, width = 10, height = 6, dpi = 300)
  
  write.table(tab,paste0("Results/Tables/", "03a_05_", out_prefix, "_", ann_col, "_by_cluster.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
  return(p)
}

plots_ann <- lapply(ann_cols, function(x) {plot_ann_by_cluster(meta_ann, x)})

