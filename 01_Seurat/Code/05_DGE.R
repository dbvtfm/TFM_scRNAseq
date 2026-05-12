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
# Script: 05_DGE.R
## ===============================

# Objetivo:
# Realizar análisis de expresión diferencial mediante pseudobulk con edgeR,
# comparando HA15 frente a control a nivel global y por anotación funcional,
# tanto en el conjunto total como en células VIH+, y comparando células VIH+
# frente a VIH- dentro de HA15 y de forma global.

# Inputs en Proc/:
# - 03b_04_int_ann.rds

# Outputs en Results/:
# - Figures/05_*.png
# - Tables/04_02_*.tsv


#----------------------------
library(Seurat)
library(Matrix)
library(edgeR)
library(dplyr)

ann <- readRDS("Proc/03b_04_int_ann.rds")

hivPos = colnames(ann) %in% WhichCells(ann, slot = 'counts', expression = Rev > 0 | Tat > 0 )
hivPos[hivPos==TRUE] <- "HIV+"
hivPos[hivPos==FALSE] <- "HIV-"
ann = AddMetaData(ann, metadata=hivPos,col.name="HIV_status")

## ============================================================
## BLOQUE 1 - HA15 vs control, total y por población celular
## ============================================================

run_edger_pb <- function(pb_counts, sample_info) {
  sample_info <- sample_info[colnames(pb_counts), , drop = FALSE]
  
  y <- DGEList(counts = pb_counts)
  design <- model.matrix(~ Replicate + Treatment, data = sample_info)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = "TreatmentHA15")
  
  res <- topTags(qlf, n = Inf)$table
  res$gene <- rownames(res)
  res <- res[, c("gene", setdiff(colnames(res), "gene"))]
  isde <- decideTests(qlf, adjust.method = "BH", p.value = 0.05, lfc = 0.5)
  
  list(res = res, y = y, isde = isde)
}

## ---- pseudobulk total por muestra ----
DefaultAssay(ann) <- "RNA"
ann <- JoinLayers(ann, assay = "RNA")
counts <- GetAssayData(ann, assay = "RNA", layer = "counts")

sample_info <- unique(ann[[]][, c("Sample", "Treatment", "Replicate")])
sample_info <- sample_info[order(sample_info$Sample), ]
rownames(sample_info) <- sample_info$Sample
sample_info$Treatment <- factor(sample_info$Treatment, levels = c("Control", "HA15"))
sample_info$Replicate <- factor(sample_info$Replicate)

sample <- factor(ann$Sample, levels = sample_info$Sample)
mm <- sparse.model.matrix(~ 0 + sample)
colnames(mm) <- levels(sample)

pb_counts <- counts %*% mm
pb_counts <- as.matrix(pb_counts)

out_pb_total <- run_edger_pb(pb_counts, sample_info)

res_pb_total <- out_pb_total$res
y_pb_total <- out_pb_total$y
isde_pb_total <- out_pb_total$isde
summary(isde_pb_total)

write.table(res_pb_total, "Results/Tables/04_02_edgeR_pseudobulk_total_HA15_vs_Control.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE)

res_pb_total_sig <- res_pb_total %>% dplyr::filter(FDR <= 0.05, abs(logFC) >= 0.5)
write.table(res_pb_total_sig, "Results/Tables/04_02_edgeR_pseudobulk_total_HA15_vs_Control_sig.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)


## ---- pseudobulk por anotación funcional ----

res_pb_ann <- list()
y_pb_ann <- list()
isde_pb_ann <- list()
sample_info_pb_ann <- list()

for (lb in sort(unique(ann$functional_ann2))) {
  cells_use <- colnames(ann)[ann$functional_ann2 == lb]
  sample <- factor(ann$Sample[cells_use], levels = sample_info$Sample)
  
  mm <- sparse.model.matrix(~ 0 + sample)
  colnames(mm) <- levels(sample)
  pb_counts <- counts[, cells_use] %*% mm
  pb_counts <- as.matrix(pb_counts)
  
  if (any(colSums(pb_counts) == 0)) next
  out <- run_edger_pb(pb_counts, sample_info)
  res <- out$res
  res$Label <- lb
  
  res_pb_ann[[lb]] <- res
  y_pb_ann[[lb]] <- out$y
  isde_pb_ann[[lb]] <- out$isde
  sample_info_pb_ann[[lb]] <- sample_info[colnames(pb_counts), , drop = FALSE]
}

res_pb_ann_tbl <- bind_rows(res_pb_ann)
write.table(res_pb_ann_tbl, "Results/Tables/04_02_edgeR_pseudobulk_by_functional_ann_HA15_vs_Control.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE)

res_pb_ann_sig <- res_pb_ann_tbl %>% dplyr::filter(FDR <= 0.05, abs(logFC) >= 0.5)
write.table(res_pb_ann_sig,"Results/Tables/04_02_edgeR_pseudobulk_by_functional_ann_HA15_vs_Control_sig.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE)

lapply(isde_pb_ann, summary)


## ===========================================================================
## BLOQUE 2 - HA15 vs control (dentro de VIH+), total y por población celular
## ===========================================================================

ann_hiv <- subset(ann, subset = HIV_status == "HIV+")

counts_hiv <- counts[, colnames(ann_hiv)]

sample_info_hiv <- unique(ann_hiv[[]][, c("Sample", "Treatment", "Replicate")])
sample_info_hiv <- sample_info_hiv[order(sample_info_hiv$Sample), ]
rownames(sample_info_hiv) <- sample_info_hiv$Sample
sample_info_hiv$Treatment <- factor(sample_info_hiv$Treatment, levels = c("Control", "HA15"))
sample_info_hiv$Replicate <- factor(sample_info_hiv$Replicate)

## ---- pseudobulk total HIV+ por muestra ----

sample <- factor(ann_hiv$Sample, levels = sample_info_hiv$Sample)
mm <- Matrix::sparse.model.matrix(~ 0 + sample)
colnames(mm) <- levels(sample)

pb_counts <- counts_hiv %*% mm
pb_counts <- as.matrix(pb_counts)

out_pb_total_hiv <- run_edger_pb(pb_counts, sample_info_hiv)
res_pb_total_hiv <- out_pb_total_hiv$res
y_pb_total_hiv <- out_pb_total_hiv$y
isde_pb_total_hiv <- out_pb_total_hiv$isde
summary(isde_pb_total_hiv)

write.table(res_pb_total_hiv, "Results/Tables/04_02_edgeR_pseudobulk_total_HIVpos_HA15_vs_Control.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE)

res_pb_total_hiv_sig <- res_pb_total_hiv %>% dplyr::filter(FDR <= 0.05, abs(logFC) >= 0.5)
write.table(res_pb_total_hiv_sig, "Results/Tables/04_02_edgeR_pseudobulk_total_HIVpos_HA15_vs_Control_sig.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)


## ---- pseudobulk HIV+ por anotación funcional ----

res_pb_ann_hiv <- list()
y_pb_ann_hiv <- list()
isde_pb_ann <- list()
sample_info_pb_ann_hiv <- list()

for (lb in sort(unique(ann_hiv$functional_ann2))) {
  
  cells_use <- colnames(ann_hiv)[ann_hiv$functional_ann2 == lb]
  sample <- factor(ann_hiv$Sample[cells_use], levels = sample_info_hiv$Sample)
  
  mm <- Matrix::sparse.model.matrix(~ 0 + sample)
  colnames(mm) <- levels(sample)
  pb_counts <- counts_hiv[, cells_use] %*% mm
  pb_counts <- as.matrix(pb_counts)
  
  if (any(colSums(pb_counts) == 0)) next
  
  out <- run_edger_pb(pb_counts, sample_info_hiv)
  res <- out$res
  res$Label <- lb
  res_pb_ann_hiv[[lb]] <- res
  y_pb_ann_hiv[[lb]] <- out$y
  isde_pb_ann[[lb]] <- out$isde
  sample_info_pb_ann_hiv[[lb]] <- sample_info_hiv[colnames(pb_counts), , drop = FALSE]
}

res_pb_ann_hiv_tbl <- dplyr::bind_rows(res_pb_ann_hiv)

write.table(res_pb_ann_hiv_tbl,"Results/Tables/04_02_edgeR_pseudobulk_by_functional_ann_HIVpos_HA15_vs_Control.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE)

res_pb_ann_hiv_sig <- res_pb_ann_hiv_tbl %>% dplyr::filter(FDR <= 0.05, abs(logFC) >= 0.5)
write.table(res_pb_ann_hiv_sig,"Results/Tables/04_02_edgeR_pseudobulk_by_functional_ann_HIVpos_HA15_vs_Control_sig.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

lapply(isde_pb_ann, summary)

## ====================================================
## BLOQUE 3 - VIH+ vs VIH- (dentro de HA15), total
## ====================================================

## ---- 1. HIV+ tratadas vs HIV- tratadas ----

ann_ha15 <- subset(ann, subset = Treatment == "HA15")
counts_ha15 <- counts[, colnames(ann_ha15)]

sample_info_HA15_HIV <- unique(ann_ha15[[]][, c("Sample", "HIV_status", "Replicate")])
sample_info_HA15_HIV <- sample_info_HA15_HIV[order(sample_info_HA15_HIV$Sample, sample_info_HA15_HIV$HIV_status), ]
sample_info_HA15_HIV$group <- paste(sample_info_HA15_HIV$Sample, sample_info_HA15_HIV$HIV_status, sep = "_")
rownames(sample_info_HA15_HIV) <- sample_info_HA15_HIV$group
sample_info_HA15_HIV$HIV_status <- factor(sample_info_HA15_HIV$HIV_status, levels = c("HIV-", "HIV+"))
sample_info_HA15_HIV$Replicate <- factor(sample_info_HA15_HIV$Replicate)

sample <- factor(paste(ann_ha15$Sample, ann_ha15$HIV_status, sep = "_"),
  levels = sample_info_HA15_HIV$group)

mm <- Matrix::sparse.model.matrix(~ 0 + sample)
colnames(mm) <- levels(sample)

pb_counts <- counts_ha15 %*% mm
pb_counts <- as.matrix(pb_counts)

sample_info_HA15_HIV <- sample_info_HA15_HIV[colnames(pb_counts), , drop = FALSE]

y <- DGEList(counts = pb_counts)
design <- model.matrix(~ Replicate + HIV_status, data = sample_info_HA15_HIV)

keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]

y <- calcNormFactors(y)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "HIV_statusHIV+")

y_pb_HA15_HIVpos_vs_HIVneg <- y

res_pb_HA15_HIVpos_vs_HIVneg <- topTags(qlf, n = Inf)$table
res_pb_HA15_HIVpos_vs_HIVneg$gene <- rownames(res_pb_HA15_HIVpos_vs_HIVneg)
res_pb_HA15_HIVpos_vs_HIVneg <- res_pb_HA15_HIVpos_vs_HIVneg[, c("gene", setdiff(colnames(res_pb_HA15_HIVpos_vs_HIVneg), "gene"))]

write.table(res_pb_HA15_HIVpos_vs_HIVneg, "Results/Tables/04_02_edgeR_pseudobulk_total_HA15_HIVpos_vs_HIVneg.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE)

isde_3 <- decideTests(qlf, adjust.method = "BH", p.value = 0.05, lfc = 0.5)

res_pb_HA15_HIVpos_vs_HIVneg_sig <- res_pb_HA15_HIVpos_vs_HIVneg %>% dplyr::filter(FDR <= 0.05, abs(logFC) >= 0.5)
write.table(res_pb_HA15_HIVpos_vs_HIVneg_sig, "Results/Tables/4_02_edgeR_pseudobulk_total_HA15_HIVpos_vs_HIVneg_sig.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

dim(res_pb_HA15_HIVpos_vs_HIVneg_sig)
summary(isde_3)

## ====================================================
## ## ---- HIV+ vs HIV- global ----
## ====================================================

sample_info_global_HIV <- unique(ann[[]][, c("Sample", "HIV_status", "Treatment", "Replicate")])
sample_info_global_HIV <- sample_info_global_HIV[order(sample_info_global_HIV$Sample, sample_info_global_HIV$HIV_status), ]
sample_info_global_HIV$group <- paste(sample_info_global_HIV$Sample, sample_info_global_HIV$HIV_status, sep = "_")
rownames(sample_info_global_HIV) <- sample_info_global_HIV$group
sample_info_global_HIV$HIV_status <- factor(sample_info_global_HIV$HIV_status, levels = c("HIV-", "HIV+"))
sample_info_global_HIV$Treatment <- factor(sample_info_global_HIV$Treatment, levels = c("Control", "HA15"))
sample_info_global_HIV$Replicate <- factor(sample_info_global_HIV$Replicate)

sample <- factor(paste(ann$Sample, ann$HIV_status, sep = "_"), levels = sample_info_global_HIV$group)

mm <- Matrix::sparse.model.matrix(~ 0 + sample)
colnames(mm) <- levels(sample)
pb_counts <- counts %*% mm
pb_counts <- as.matrix(pb_counts)
sample_info_global_HIV <- sample_info_global_HIV[colnames(pb_counts), , drop = FALSE]

y <- DGEList(counts = pb_counts)
design <- model.matrix(~ Replicate + Treatment + HIV_status, data = sample_info_global_HIV)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "HIV_statusHIV+")
y_pb_global_HIVpos_vs_HIVneg <- y

res_pb_global_HIVpos_vs_HIVneg <- topTags(qlf, n = Inf)$table
res_pb_global_HIVpos_vs_HIVneg$gene <- rownames(res_pb_global_HIVpos_vs_HIVneg)
res_pb_global_HIVpos_vs_HIVneg <- res_pb_global_HIVpos_vs_HIVneg[, c("gene", setdiff(colnames(res_pb_global_HIVpos_vs_HIVneg), "gene"))]

write.table(res_pb_global_HIVpos_vs_HIVneg, "Results/Tables/04_02_edgeR_pseudobulk_total_global_HIVpos_vs_HIVneg.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE)


isde_4 <- decideTests(qlf, adjust.method = "BH", p.value = 0.05, lfc = 0.5)

res_pb_global_HIVpos_vs_HIVneg_sig <- res_pb_global_HIVpos_vs_HIVneg %>% dplyr::filter(FDR <= 0.05, abs(logFC) >= 0.5)
write.table(res_pb_global_HIVpos_vs_HIVneg_sig, "Results/Tables/04_02_edgeR_pseudobulk_total_global_HIVpos_vs_HIVneg_sig.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

dim(res_pb_global_HIVpos_vs_HIVneg_sig)
summary(isde_4)

## ==========================================
## BLOQUE 4 - FIGURAS
## ==========================================

library(edgeR)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(ggrepel)

## ---- 4.1. HIV+ tratadas vs HIV- tratadas ----

plot_volcano_edger <- function(res, title = NULL) {
  res <- as.data.frame(res)
  
  if (!"gene" %in% colnames(res)) {res$gene <- rownames(res)}
  
  res$Significance <- "NS"
  res$Significance[res$FDR <= 0.05 & res$logFC >= 0.5] <- "Up"
  res$Significance[res$FDR <= 0.05 & res$logFC <= -0.5] <- "Down"
  
  top_lab <- res %>% 
    filter(FDR <= 0.05, abs(logFC) >= 0.5) %>%
    arrange(desc(abs(logFC))) %>% 
    slice_head(n = 20)
  
  ggplot(res, aes(x = logFC, y = -log10(FDR), color = Significance)) +
    geom_point(size = 1.2, alpha = 0.8) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, colour = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
    geom_text_repel(data = top_lab, aes(label = gene), size = 3, max.overlaps = 20, show.legend = F) +
    scale_color_manual(values = c("Down" = "royalblue3", "NS" = "grey70", "Up" = "firebrick3")) +
    theme_bw() + labs(title = title, x = "log2FC", y = "-log10(FDR)")
}

plot_heatmap_edger <- function(y, res, sample_info, title = NULL, n_up = 10, n_down = 10) {
  res <- as.data.frame(res)
  
  if (!"gene" %in% colnames(res)) {res$gene <- rownames(res)}

  genes_up <- res %>% filter(FDR <= 0.05, logFC >= 0.5) %>% arrange(FDR, desc(logFC)) %>%
    slice_head(n = n_up) %>% pull(gene)
  genes_down <- res %>% filter(FDR <= 0.05, logFC <= -0.5) %>% arrange(FDR, logFC) %>%
    slice_head(n = n_down) %>% pull(gene)
  genes_use <- unique(c(genes_up, genes_down))
  
  if (length(genes_use) == 0) {
    message("No hay genes significativos para el heatmap")
    return(NULL)
  }
  
  mat <- cpm(y, log = TRUE)[genes_use, , drop = FALSE]
  mat <- t(scale(t(mat)))
  
  id_col <- NULL
  if ("Sample" %in% colnames(sample_info) && all(colnames(mat) %in% sample_info$Sample)) {
    id_col <- "Sample"
  }
  if ("group" %in% colnames(sample_info) && all(colnames(mat) %in% sample_info$group)) {
    id_col <- "group"
  }
  
  sample_info <- sample_info[match(colnames(mat), sample_info[[id_col]]), , drop = FALSE]
  rownames(sample_info) <- sample_info[[id_col]]
  ann_col <- sample_info[, setdiff(colnames(sample_info), c(id_col, "Sample")), drop = FALSE]
  
  ann_colors <- list(Treatment = c("Control" = "grey70", "HA15" = "darkseagreen3"),
    HIV_status = c("HIV-" = "grey70", "HIV+" = "indianred3"),
    Replicate = c("Rep1" = "grey20", "Rep2" = "grey4", "Rep3" = "grey40"))
  ann_colors <- ann_colors[names(ann_colors) %in% colnames(ann_col)]
  
  pheatmap(mat, annotation_col = ann_col, annotation_colors = ann_colors,
    cellwidth = 16, cluster_rows = TRUE, cluster_cols = TRUE,
    fontsize_row = 8, fontsize_col = 8, main = title, silent = FALSE)
}

## ---- 4.2 Volcano + heatmap total HA15 vs Control ----

sample_info_total <- sample_info[match(colnames(y_pb_total$counts), sample_info$Sample), , drop = FALSE]
rownames(sample_info_total) <- NULL

p_vol_total <- plot_volcano_edger(res_pb_total, title = "HA15 vs Control - total")
ggsave("Results/Figures/05_01_volcano_total_HA15_vs_Control.png", p_vol_total, width = 6, height = 5, dpi = 300)

png("Results/Figures/05_01_heatmap_total_HA15_vs_Control.png",width = 1800, height = 1600, res = 300)

plot_heatmap_edger(y = y_pb_total,res = res_pb_total,sample_info = sample_info_total,
  title = "HA15 vs Control - total",n_up = 10, n_down = 10)
dev.off()

## ---- 4.3 Volcano + heatmap por anotación funcional HA15 vs Control ----

for (nm in names(res_pb_ann)) {
  
  p_vol <- plot_volcano_edger(res_pb_ann[[nm]],title = paste0(nm, " - HA15 vs Control"))
  
  ggsave(paste0("Results/Figures/05_02_volcano_", make.names(nm), "_HA15_vs_Control.png"),
    p_vol, width = 6, height = 5, dpi = 300)
  
  png(paste0("Results/Figures/05_02_heatmap_", make.names(nm), "_HA15_vs_Control.png"),
    width = 1800, height = 1600, res = 300)
  
  plot_heatmap_edger(y = y_pb_ann[[nm]],res = res_pb_ann[[nm]],sample_info = sample_info_pb_ann[[nm]],
    title = paste0(nm, " - HA15 vs Control"),n_up = 10, n_down = 10)
  
  dev.off()
}

## ---- 4.4 Volcano + heatmap total HIV+ HA15 vs HIV+ Control ----

sample_info_total_hiv <- sample_info_hiv[match(colnames(y_pb_total_hiv$counts), sample_info_hiv$Sample), , drop = FALSE]
rownames(sample_info_total_hiv) <- NULL

p_vol_total_hiv <- plot_volcano_edger(res_pb_total_hiv,title = "HIV+ HA15 vs HIV+ Control - total")
ggsave("Results/Figures/05_03_volcano_total_HIVpos_HA15_vs_Control.png",p_vol_total_hiv,
       width = 6, height = 5, dpi = 300)

png("Results/Figures/05_03_heatmap_total_HIVpos_HA15_vs_Control.png",width = 1800, height = 1600, res = 300)

plot_heatmap_edger(y = y_pb_total_hiv,res = res_pb_total_hiv,sample_info = sample_info_total_hiv,
  title = "HIV+ HA15 vs HIV+ Control - total",n_up = 10, n_down = 10)
dev.off()

## ---- 4.5 Volcano + heatmap por anotación funcional HIV+ HA15 vs HIV+ Control ----

for (nm in names(res_pb_ann_hiv)) {
  
  p_vol <- plot_volcano_edger(res_pb_ann_hiv[[nm]],title = paste0(nm, " - HIV+ HA15 vs HIV+ Control"))
  
  ggsave(paste0("Results/Figures/05_04_volcano_", make.names(nm), "_HIVpos_HA15_vs_Control.png"),
    p_vol, width = 6, height = 5, dpi = 300)
  
  png(paste0("Results/Figures/05_04_heatmap_", make.names(nm), "_HIVpos_HA15_vs_Control.png"),
    width = 1800, height = 1600, res = 300)
  
  plot_heatmap_edger(y = y_pb_ann_hiv[[nm]],res = res_pb_ann_hiv[[nm]],sample_info = sample_info_pb_ann_hiv[[nm]],
    title = paste0(nm, " - HIV+ HA15 vs HIV+ Control"),n_up = 10, n_down = 10)
  
  dev.off()
}

## ---- 4.6 Volcano + heatmap total HIV+ vs HIV- en HA15 ----

p_vol_HA15_HIV <- plot_volcano_edger(res_pb_HA15_HIVpos_vs_HIVneg,title = "HIV+ vs HIV- en HA15 - total")
ggsave("Results/Figures/05_05_volcano_total_HA15_HIVpos_vs_HIVneg.png",p_vol_HA15_HIV,
       width = 6, height = 5, dpi = 300)

png("Results/Figures/05_05_heatmap_total_HA15_HIVpos_vs_HIVneg.png",width = 1800, height = 1600, res = 300)

plot_heatmap_edger(y = y_pb_HA15_HIVpos_vs_HIVneg,res = res_pb_HA15_HIVpos_vs_HIVneg,
  sample_info = sample_info_HA15_HIV,title = "HIV+ vs HIV- en HA15 - total",n_up = 10, n_down = 10)
dev.off()

## ---- 4.7 Volcano + heatmap total HIV+ vs HIV- global ----

p_vol_global_HIV <- plot_volcano_edger(res_pb_global_HIVpos_vs_HIVneg,title = "HIV+ vs HIV- global - total")
ggsave("Results/Figures/05_06_volcano_total_global_HIVpos_vs_HIVneg.png",p_vol_global_HIV,
       width = 6, height = 5, dpi = 300)

png("Results/Figures/05_06_heatmap_total_global_HIVpos_vs_HIVneg.png",width = 1800, height = 1600, res = 300)

plot_heatmap_edger(y = y_pb_global_HIVpos_vs_HIVneg,res = res_pb_global_HIVpos_vs_HIVneg,
  sample_info = sample_info_global_HIV,title = "HIV+ vs HIV- global - total",n_up = 10, n_down = 10)
dev.off()

