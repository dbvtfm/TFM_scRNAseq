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
# Script: 01_data_prep_and_QC.R
## ===============================

# Objetivo:
# Importar matrices GEX y HTO procedentes de alevin-fry, construir el objeto
# Seurat inicial, realizar demultiplexado HTO, anotar metadatos experimentales,
# aplicar control de calidad celular y exportar los inputs para pySCENIC.

# Inputs en Raw/:
# - GEX_quants_mat.mtx
# - GEX_quants_mat_rows.txt
# - GEX_quants_mat_cols.txt
# - HTO_quants_mat.mtx
# - HTO_quants_mat_rows.txt
# - HTO_quants_mat_cols.txt
# - gencode-v36_withExtraMito_HIV_p379delta6_tx2gene.tsv
# - translation_3M-february-2018.txt
#(https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/translation_3M-february-2018.txt.gz)

# Outputs:
# - Proc/01_data_prep_and_QC.rds
# - SCENIC/01_preprocessed_inputs/metadata.tsv
# - SCENIC/01_preprocessed_inputs/counts.loom
# - Results/Figures/01_*.png
# - Results/Tables/01_*.tsv
  

## ===============================
## BLOQUE 1 - IMPORTACIÓN DE DATOS Y CREACIÓN DEL OBJETO SEURAT
## ===============================


## ---- 1.1 CONSTRUCCIÓN DE MATRICES GEX Y HTO A PARTIR DE OUTPUTS DE ALEVIN ---- 

# función para leer outputs de alevin y prepararlos para Seurat: 
# lee matriz dispersa de UMIs y tabla de cell barcodes (rows) y features (cols)
# asigna nombres a columnas y filas y traspone la matriz 
# (en Seurat: cells = cols y features = rows, en alevin al revés)

library(Matrix)

read_alevin <- function(mtx, cells, feats) {
  M <- readMM(mtx)
  c <- readLines(cells)
  f <- readLines(feats)
  stopifnot(nrow(M) == length(c), ncol(M) == length(f))
  rownames(M) <- c
  colnames(M) <- f
  t(M)
}

GEX <- read_alevin("Raw/GEX_quants_mat.mtx",
                   cells = "Raw/GEX_quants_mat_rows.txt",
                   feats = "Raw/GEX_quants_mat_cols.txt")

HTO <- read_alevin("Raw/HTO_quants_mat.mtx",
                   cells = "Raw/HTO_quants_mat_rows.txt",
                   feats = "Raw/HTO_quants_mat_cols.txt")

#dim(GEX)
#[1] 182010  10577
#dim(HTO)
#[1]    10 10893


## ---- 1.2 PREPARACIÓN DE LAS MATRICES GEX y HTO PARA SEURAT ---- 

# función para colapsar formato USA a recuentos S+A, sumando y quitando unspliced
# las features spliced no llevan sufijo, las ambiguas -A, y las unspliced -U
# sólo aplica para la matriz GEX, la HTO no se cuantificó en modo USA
# también mapea ENS a gene symbol para evitar filas duplicadas

count_SA <- function(read_alevin_output) {
  # filtrar -U fuera, conservar feats sin sufijo (S) y con -A, eliminar sufijo -A
  sa <- !grepl("-U$", rownames(read_alevin_output))
  M <- read_alevin_output[sa,]
  f <- sub("-A$", "", rownames(M))
  #matriz de pertenencia de feats con mismo nombre para preparar multiplicación matricial
  f_SA <- t(sparse.model.matrix(~ 0 + f))
  M_SA <- f_SA %*% M    # multiplicación matricial (suma de feats S + A)
  #eliminar el prefijo "f" que añade la matriz de pertenencia
  rownames(M_SA) <- sub("^f", "", rownames(M_SA))
  M_SA
}

GEX <- count_SA(GEX)

#dim(GEX)
#[1] 60670 10577  ----> 1/3 de las features originales (formato USA), correcto


# tabla de equivalencias ENSG - gene symbol del genoma quimérico humano-VIH
ens2gene <- read.delim("Raw/gencode-v36_withExtraMito_HIV_p379delta6_tx2gene.tsv",
                       header = FALSE, col.names = c("ens", "gene"))

# función para mapear ENSG a nombres de genes y colapsar los duplicados
# necesita otra multiplicación matricial para sumar filas con mismo nombre de gen

gene_name_dedup <- function(count_SA_output, ens2gene) {
  # asignar gene symbol a cada ENSG
  g <- ens2gene$gene[match(rownames(count_SA_output), ens2gene$ens)]
  # matriz de pertenencia de feats con el mismo gene symbol
  g_mat <- t(sparse.model.matrix(~ 0 + g)) # multiplicación matricial
  out <- g_mat %*% count_SA_output
  rownames(out) <- sub("^g", "", rownames(out)) # eliminar prefijo "g"
  out
}

GEX <- gene_name_dedup(GEX, ens2gene)
rm(ens2gene)

#dim(GEX)
#[1] 59437 10577
#> anyDuplicated(rownames(GEX))
#[1] 0
#sum(is.na(rownames(GEX)))
#[1] 0

# es necesaria una tabla interna de CellRanger para mapear los barcodes de HTO con los barcodes de GEX
#https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3fb.html

map <- read.table("Raw/translation_3M-february-2018.txt")
rownames(map) <- map$V1
colnames(GEX) <- map[colnames(GEX), 2]
rm(map)

#length(intersect(colnames(GEX), colnames(HTO)))
# [1] 9951

# los autores sólo usaron 9 HTO, si no se elimina el resto, habrá dobletes
# nota: de los 9 HTO, sólo 6 corresponden a las muestras de este estudio
HTO <- HTO[1:9,]

# intersección de cell barcodes entre GEX y HTO (células en común)
common_cells <- intersect(colnames(GEX), colnames(HTO))

GEX <- GEX[, common_cells]
HTO <- HTO[, common_cells]
rm(common_cells)

#dim(GEX)
#[1] 59437  9951
#dim(HTO)
#[1]    9 9951


## ---- 1.4 CREACIÓN DEL OBJETO SEURAT ---- 
library(Seurat)

obj <- CreateSeuratObject(counts = GEX, assay = "RNA")
obj[["HTO"]] <- CreateAssayObject(counts = HTO)



## ========================================================
## BLOQUE 2 - DEMULTIPLEXADO HTO Y ANOTACIÓN DE METADATOS
## ========================================================
library(Seurat)
library(dplyr)
library(ggplot2)

## ---- 2.1 Normalización y demultiplexado ----

obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR")
obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.99)

## ---- 2.2 Visualización del demultiplexado ----

# Ridge Plot
Idents(obj) <- "HTO_maxID"
HTO_ridge <- RidgePlot(obj, assay = "HTO", features = rownames(obj[["HTO"]]), ncol = 3)
ggsave("Results/Figures/01_02_HTO_ridgePlot.png", HTO_ridge, width = 12, height = 8, dpi = 300)


# UMAP de maxID y hashID
obj <- ScaleData(obj, assay = "HTO")
VariableFeatures(obj, assay = "HTO") <- rownames(obj[["HTO"]])
obj <- RunPCA(obj,assay = "HTO",reduction.name = "hto.pca",
              reduction.key = "HPC_",approx = FALSE)
obj <- RunUMAP(obj,reduction = "hto.pca",dims = 1:9,reduction.name = "hto.umap",
               reduction.key = "HUMAP_")

HTO_umap_max <- DimPlot(obj, reduction = "hto.umap", group.by = "HTO_maxID", label = TRUE)
ggsave("Results/Figures/01_02_HTO__UMAP_maxID.png", HTO_umap_max, width = 7, height = 5, dpi = 300)

HTO_umap_ID <- DimPlot(obj, reduction = "hto.umap", group.by = "hash.ID", label = TRUE)
ggsave("Results/Figures/01_02_HTO__UMAP_hashID.png", HTO_umap_ID, width = 7, height = 5, dpi = 300)


## ---- 2.3 Tablas resumen demultiplexado ----

# resumen global del demultiplexado
tab_demux <- as.data.frame(table(obj$HTO_classification.global)) |>
  rename(HTO_classification.global = Var1, n_cells = Freq) |>
  mutate(prop_total = n_cells / sum(n_cells))

write.table(tab_demux, "Results/Tables/01_02_HTO_demux_global.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

## 2) Singlets por hash.ID
tab_singlets <- as.data.frame(table(obj$hash.ID[obj$HTO_classification.global == "Singlet"])) |>
  rename(hash.ID = Var1, n_singlets = Freq) |>
  mutate(prop_singlets = n_singlets / sum(n_singlets)) |>
  arrange(dplyr::desc(n_singlets))

write.table(tab_singlets, "Results/Tables/01_02_HTO_singlets_by_hashID.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

rm(tab_demux, tab_singlets)


## ---- 2.4 Filtrado de dobletes y hashtags y anotación de metadatos ----
library(stringr)

obj <- subset(obj, subset = hash.ID != "Doublet")
obj <- subset(obj, subset = hash.ID %in% c("Hashtag-1","Hashtag-2","Hashtag-3","Hashtag-4","Hashtag-5","Hashtag-6"))

annot = as.data.frame(cbind(
  "HTO" = sort(unique(as.character(obj$hash.ID))),
  "Sample" = c("Control_Rep1","Control_Rep2","Control_Rep3","HA15_Rep1","HA15_Rep2","HA15_Rep3")))
rownames(annot) = annot$HTO
annot$Treatment = cbind(str_replace_all(annot$Sample,"(.+)_.+","\\1"))
annot$Replicate = cbind(str_replace_all(annot$Sample,".+_(.+)","\\1"))

obj <- AddMetaData(obj,metadata=as.character(annot[as.character(obj$hash.ID),2]),col.name="Sample")
obj <- AddMetaData(obj,metadata=as.character(annot[as.character(obj$hash.ID),3]),col.name="Treatment")
obj <- AddMetaData(obj,metadata=as.character(annot[as.character(obj$hash.ID),4]),col.name="Replicate")
rm(annot)


## ========================================================
## BLOQUE 3 - FILTRADO DE DATOS GEX Y CONTROL DE CALIDAD 
## ========================================================
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)


library(gt)

## ---- 3.1 QC prefiltrado ----

DefaultAssay(obj) <- "RNA"
obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

obj_preqc <- obj

# función para histograma + densidad 
hist_qc <- function(df, x, title, bins = 100) {
  ggplot(df, aes(x = .data[[x]])) +
    geom_histogram(aes(y = after_stat(density)), bins = bins) +
    geom_density() + labs(title = title, x = x, y = "Frequency") + theme_bw()}

# histogramas
df_preqc <- data.frame(n_genes = obj_preqc$nFeature_RNA, n_counts = obj_preqc$nCount_RNA,
                       Percent_mito = obj_preqc$percent.mt / 100)

p1 <- hist_qc(df_preqc, "n_genes", "Number of genes expressed per cell", bins = 100)
p2 <- hist_qc(df_preqc, "n_counts", "Counts per cell", bins = 100)
p3 <- hist_qc(df_preqc, "Percent_mito", "Mitochondrial read fraction per cell", bins = 100)
p_preqc <- p1 + p2 + p3 + plot_layout(ncol = 3)

ggsave("Results/Figures/01_03_QC_histo_prefilter.png", p_preqc, width = 12, height = 4, dpi = 300)
rm(p1,p2,p3,p_preqc,df_preqc)

# violin y scatter plots by sample
vln_pre <- VlnPlot(obj_preqc, features = c("nCount_RNA", "nFeature_RNA", "percent.mt",
                                   "percent.ribo"), group.by = "Sample", pt.size = 0.1, ncol = 2)
sc1_pre <- FeatureScatter(obj_preqc, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample")
sc2_pre <- FeatureScatter(obj_preqc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample")

ggsave("Results/Figures/01_03_QC_viol_bysample_prefilter.png", vln_pre, width = 12, height = 8, dpi = 300)
ggsave("Results/Figures/01_03_QC_sc_nCount-mt_bysample_prefilter.png", sc1_pre, width = 7, height = 5, dpi = 300)
ggsave("Results/Figures/01_03_QC_sc_nCount-nFeat_bysample_prefilter.png", sc2_pre, width = 7, height = 5, dpi = 300)
rm(vln_pre,sc1_pre,sc2_pre)

# violin y scatter plots by treatment
vln_pre <- VlnPlot(obj_preqc, features = c("nCount_RNA", "nFeature_RNA", "percent.mt",
                                           "percent.ribo"), group.by = "Treatment", pt.size = 0.1, ncol = 2)
sc1_pre <- FeatureScatter(obj_preqc, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Treatment")
sc2_pre <- FeatureScatter(obj_preqc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Treatment")

ggsave("Results/Figures/01_03_QC_viol_bytreat_prefilter.png", vln_pre, width = 12, height = 8, dpi = 300)
ggsave("Results/Figures/01_03_QC_sc_nCount-mt_bytreat_prefilter.png", sc1_pre, width = 7, height = 5, dpi = 300)
ggsave("Results/Figures/01_03_QC_sc_nCount-nFeat_bytreat_prefilter.png", sc2_pre, width = 7, height = 5, dpi = 300)
rm(vln_pre,sc1_pre,sc2_pre)

#tabla resumen

tb_preqc <- data.frame(n_cells = ncol(obj_preqc),
  median_nCount = median(obj_preqc$nCount_RNA),
  median_nFeature = median(obj_preqc$nFeature_RNA),
  median_percent_mt = median(obj_preqc$percent.mt),
  median_percent_ribo = median(obj_preqc$percent.ribo))

write.table(tb_preqc,"Results/Tables/01_03_qc_prefilter.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
rm(tb_preqc)

## ---- 3.2 Filtrado ----
library(miQC)
library(SeuratWrappers)
library(flexmix)

# miQC
obj <- RunMiQC(object = obj, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA",
  posterior.cutoff = 0.75, model.slot = "flexmix_model")

miqc_meta <- obj@meta.data

# filtrado miQC + umbrales RNA mínimos
obj <- subset(obj, subset = miQC.keep == "keep")
obj <- subset(obj, subset = nCount_RNA > 2000 & nFeature_RNA > 1000)

#dim(obj)
#[1] 59437  5312

## ---- 3.3 QC Postfiltrado ----

# figura miQC
p_miqc <- ggplot(miqc_meta, aes(x = nFeature_RNA, y = percent.mt, color = miQC.keep)) +
  geom_point(size = 0.4, alpha = 0.7) + theme_bw() + labs(title = "miQC classification",
    x = "Detected genes per cell", y = "Percent mitochondrial reads")

ggsave("Results/Figures/01_03_miQC_class.png", p_miqc, width = 7, height = 5, dpi = 300)
rm(p_miqc)

# tabla miQC
tb_miqc <- as.data.frame(table(miqc_meta$miQC.keep))
colnames(tb_miqc) <- c("miQC_keep", "n_cells")

write.table(tb_miqc,"Results/Tables/01_03_miQC_keep_counts.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
rm(tb_miqc)

# tabla cell counts pre vs post
qc_compare <- data.frame(stage = c("post_HTO_singlets", "post_miQC_and_thresholds"),
  n_cells = c(ncol(obj_preqc), ncol(obj)))

write.table(qc_compare,
  "Results/Tables/01_03_QC_cell_counts_pre_vs_post.tsv", sep = "\t",
  row.names = FALSE, quote = FALSE)
rm(qc_compare)

# histogramas post
df_postqc <- data.frame(n_genes = obj$nFeature_RNA, n_counts = obj$nCount_RNA,
  Percent_mito = obj$percent.mt / 100)

p1_post <- hist_qc(df_postqc, "n_genes", "Number of genes expressed per cell", bins = 100)
p2_post <- hist_qc(df_postqc, "n_counts", "Counts per cell", bins = 100)
p3_post <- hist_qc(df_postqc, "Percent_mito", "Mitochondrial read fraction per cell", bins = 100)

p_postqc <- p1_post + p2_post + p3_post + plot_layout(ncol = 3)
ggsave("Results/Figures/01_03_QC_histo_postfilter.png", p_postqc, width = 12, height = 4, dpi = 300)
rm(df_postqc, p1_post, p2_post, p3_post, p_postqc)

# violin y scatter plots by sample
vln_post <- VlnPlot(obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt",
                                      "percent.ribo"), group.by = "Sample", pt.size = 0.1, ncol = 2)
sc1_post <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample")
sc2_post <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample")

ggsave("Results/Figures/01_03_QC_viol_bysample_postfilter.png", vln_post, width = 12, height = 8, dpi = 300)
ggsave("Results/Figures/01_03_QC_sc_nCount-mt_bysample_postfilter.png", sc1_post, width = 7, height = 5, dpi = 300)
ggsave("Results/Figures/01_03_QC_sc_nCount-nFeat_bysample_postfilter.png", sc2_post, width = 7, height = 5, dpi = 300)
rm(vln_post,sc1_post,sc2_post)

# violin y scatter plots by treatment
vln_post <- VlnPlot(obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt",
                                      "percent.ribo"), group.by = "Treatment", pt.size = 0.1, ncol = 2)
sc1_post <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Treatment")
sc2_post <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Treatment")

ggsave("Results/Figures/01_03_QC_viol_bytreat_postfilter.png", vln_post, width = 12, height = 8, dpi = 300)
ggsave("Results/Figures/01_03_QC_sc_nCount-mt_bytreat_postfilter.png", sc1_post, width = 7, height = 5, dpi = 300)
ggsave("Results/Figures/01_03_QC_sc_nCount-nFeat_bytreat_postfilter.png", sc2_post, width = 7, height = 5, dpi = 300)
rm(vln_post,sc1_post,sc2_post)

# tabla resumen QC post
tb_postqc <- data.frame(n_cells = ncol(obj),
  median_nCount = median(obj$nCount_RNA),
  median_nFeature = median(obj$nFeature_RNA),
  median_percent_mt = median(obj$percent.mt),
  median_percent_ribo = median(obj$percent.ribo))

write.table(tb_postqc,"Results/Tables/01_03_qc_postfilter.tsv",sep = "\t",
            row.names = FALSE, quote = FALSE)
rm(tb_postqc)



saveRDS(obj, "Proc/01_data_prep_and_QC.rds")


## ========================================================
## BLOQUE 4 - EXPORTACIÓN DE DATOS PARA pySCENIC 
## ========================================================
library(Seurat)
library(SCopeLoomR)

# guardar metadatos de Seurat (para obtener RSS de una de las columnas)
write.table(obj@meta.data, "SCENIC/01_preprocessed_inputs/metadata.tsv",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)


# obtener archivo loom
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
counts <- counts[Matrix::rowSums(counts > 0) >= 3, ]

build_loom(
  file.name = "SCENIC/01_preprocessed_inputs/counts.loom",
  dgem = counts,
  title = "input loom file from seurat",
  genome = "Human"
)



