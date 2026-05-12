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
# Script: 04_cell_abundance.R
## ===============================

# Objetivo:
# Analizar abundancia celular diferencial por tratamiento, composición de células
# VIH+ y diferencias de composición entre células VIH+ y VIH- en el grupo HA15,
# usando anotaciones funcionales de cluster y modelos por muestra/réplica.

# Inputs en Proc/:
# - 03b_04_int_ann.rds

# Outputs en Results/:
# - Figures/04_*.png
# - Tables/04_*.tsv


ann <- readRDS("Proc/03b_04_int_ann.rds")

## ==========================================
## BLOQUE 1 - ABUNDANCIA CELULAR
## ==========================================
library(speckle)
library(SingleCellExperiment)
library(CellBench)
library(limma)
library(ggplot2)
library(scater)
library(patchwork)
library(edgeR)
library(statmod)


## ---- 1.1 Efecto del tratamiento (HA15) en la composición celular ----

# DA de anotación
props <- getTransformedProps(clusters = ann$functional_ann2, sample = ann$Sample, transform = "logit")

sample_info <- unique(ann[[]][, c("Sample", "Treatment", "Replicate")])
sample_info <- sample_info[match(colnames(props$TransformedProps), sample_info$Sample), ]
rownames(sample_info) <- NULL
pair <- as.numeric(factor(sample_info$Replicate))

design <- model.matrix(~ 0 + Treatment + pair, data = sample_info)
contrasts <- makeContrasts(HA15_vs_Control = TreatmentHA15 - TreatmentControl, levels = design)

res_prop <- propeller.ttest(prop.list = props, design = design, contrasts = contrasts,
  robust = TRUE, trend = FALSE, sort = TRUE)

write.table(res_prop,"Results/Tables/04_01_res_prop_ann.tsv",sep = "\t", row.names = T,
  col.names = NA, quote = FALSE)

# figura
proptable <- as.data.frame(table(ann$functional_ann2, ann$Sample)) %>%
  tibble::as_tibble() %>% dplyr::rename(Sample = Var2) %>% dplyr::left_join(
    tibble::as_tibble(as.data.frame(table(ann$Sample))) %>% dplyr::rename(Sample = Var1),
    by = "Sample") %>% dplyr::mutate(Proportion = Freq.x / Freq.y) %>% dplyr::select(-Freq.x, -Freq.y) %>%
  tidyr::pivot_wider(names_from = Sample, values_from = Proportion) %>% dplyr::rename(Cluster = Var1)

prop_long <- proptable %>%
  tidyr::pivot_longer(cols = !Cluster, names_to = "Sample", values_to = "Proportion") %>%
  tidyr::extract(Sample, c("Treatment", "Replicate"), "(^.+)_(.+)") %>%
  dplyr::mutate(Treatment = factor(Treatment, levels = c("Control", "HA15")))
sig <- data.frame(Cluster = rownames(res_prop),FDR = res_prop$FDR)
sig$label <- ifelse(sig$FDR < 0.001, "***", ifelse(sig$FDR < 0.01, "**", ifelse(sig$FDR < 0.05, "*", "")))
sig <- sig[sig$label != "", ]
ymax <- tapply(prop_long$Proportion, prop_long$Cluster, max)
sig$y <- ymax[match(sig$Cluster, names(ymax))] + 0.02

p_prop <- ggplot(prop_long, aes(x = Cluster, y = Proportion, fill = Treatment, color = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.75)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), size = 1.5) +
  geom_text(data = sig, aes(x = Cluster, y = y, label = label), inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = c("Control" = "grey70", "HA15" = "darkseagreen3")) +
  scale_color_manual(values = c("Control" = "grey40", "HA15" = "darkgreen")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + labs(x = "Annotation", y = "Proportion")
ggsave("Results/Figures/04_functional_ann_by_treat.png", p_prop, width = 9, height = 4, dpi = 300)


## ---- 1.2 Efecto del tratamiento (HA15) en la composición celular VIH+ ----
## (efecto en la reactivación viral)

# asignación de estatus VIH

hivPos = colnames(ann) %in% WhichCells(ann, slot = 'counts', expression = Rev > 0 | Tat > 0 )
hivPos[hivPos==TRUE] <- "HIV+"
hivPos[hivPos==FALSE] <- "HIV-"
ann = AddMetaData(ann, metadata=hivPos,col.name="HIV_status")

# DA de anotación en células HIV+, composición celular
ann_hiv <- subset(ann, subset = HIV_status == "HIV+")

props <- getTransformedProps(clusters = ann_hiv$functional_ann2, sample = ann_hiv$Sample, transform = "logit")

sample_info <- unique(ann_hiv[[]][, c("Sample", "Treatment", "Replicate")])
sample_info <- sample_info[match(colnames(props$TransformedProps), sample_info$Sample), ]
rownames(sample_info) <- NULL

pair <- as.numeric(factor(sample_info$Replicate))
design <- model.matrix(~ 0 + Treatment + pair, data = sample_info)
contrasts <- makeContrasts(HA15_vs_Control = TreatmentHA15 - TreatmentControl, levels = design)

res_prop <- propeller.ttest(prop.list = props, design = design, contrasts = contrasts,
                            robust = TRUE, trend = FALSE, sort = TRUE)

write.table(res_prop,"Results/Tables/04_01_res_prop_ann_HIVpos.tsv",sep = "\t", row.names = T,
            col.names = NA, quote = FALSE)

# figura
proptable <- as.data.frame(table(ann_hiv$functional_ann2, ann_hiv$Sample)) %>%
  tibble::as_tibble() %>% dplyr::rename(Sample = Var2) %>% dplyr::left_join(
    tibble::as_tibble(as.data.frame(table(ann_hiv$Sample))) %>% dplyr::rename(Sample = Var1),
    by = "Sample") %>% dplyr::mutate(Proportion = Freq.x / Freq.y) %>% dplyr::select(-Freq.x, -Freq.y) %>%
  tidyr::pivot_wider(names_from = Sample, values_from = Proportion) %>% dplyr::rename(Cluster = Var1)

prop_long <- proptable %>%
  tidyr::pivot_longer(cols = !Cluster, names_to = "Sample", values_to = "Proportion") %>%
  tidyr::extract(Sample, c("Treatment", "Replicate"), "(^.+)_(.+)") %>%
  dplyr::mutate(Treatment = factor(Treatment, levels = c("Control", "HA15")))
sig <- data.frame(Cluster = rownames(res_prop),FDR = res_prop$FDR)
sig$label <- ifelse(sig$FDR < 0.001, "***", ifelse(sig$FDR < 0.01, "**", ifelse(sig$FDR < 0.05, "*", "")))
sig <- sig[sig$label != "", ]
ymax <- tapply(prop_long$Proportion, prop_long$Cluster, max)
sig$y <- ymax[match(sig$Cluster, names(ymax))] + 0.02

p_prop <- ggplot(prop_long, aes(x = Cluster, y = Proportion, fill = Treatment, color = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.75)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), size = 1.5) +
  geom_text(data = sig, aes(x = Cluster, y = y, label = label), inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = c("Control" = "grey70", "HA15" = "darkseagreen3")) +
  scale_color_manual(values = c("Control" = "grey40", "HA15" = "darkgreen")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + labs(x = "Annotation", y = "Proportion")
ggsave("Results/Figures/04_functional_ann_HIVpos_by_treat.png", p_prop, width = 9, height = 4, dpi = 300)


# DA de porcentaje total HIV+ / HIV-
props <- getTransformedProps(clusters = ann$HIV_status, sample = ann$Sample, transform = "logit")

sample_info <- unique(ann[[]][, c("Sample", "Treatment", "Replicate")])
sample_info <- sample_info[match(colnames(props$TransformedProps), sample_info$Sample), ]
rownames(sample_info) <- NULL

pair <- as.numeric(factor(sample_info$Replicate))
design <- model.matrix(~ 0 + Treatment + pair, data = sample_info)
contrasts <- makeContrasts(HA15_vs_Control = TreatmentHA15 - TreatmentControl, levels = design)

res_prop <- propeller.ttest(prop.list = props, design = design, contrasts = contrasts,
                              robust = TRUE, trend = FALSE, sort = TRUE)

write.table(res_prop, "Results/Tables/04_01_res_prop_HIV_status.tsv", sep = "\t",
            row.names = TRUE, col.names = NA, quote = FALSE)

# figura: solo proporción HIV+
proptable <- as.data.frame(table(ann$HIV_status, ann$Sample)) %>%
  tibble::as_tibble() %>%
  dplyr::rename(Sample = Var2) %>%
  dplyr::left_join(
    tibble::as_tibble(as.data.frame(table(ann$Sample))) %>% dplyr::rename(Sample = Var1),
    by = "Sample"
  ) %>%
  dplyr::mutate(Proportion = Freq.x / Freq.y) %>%
  dplyr::select(-Freq.x, -Freq.y) %>%
  tidyr::pivot_wider(names_from = Sample, values_from = Proportion) %>%
  dplyr::rename(HIV_status = Var1)

prop_long <- proptable %>%
  dplyr::filter(HIV_status == "HIV+") %>%
  tidyr::pivot_longer(cols = !HIV_status, names_to = "Sample", values_to = "Proportion") %>%
  tidyr::extract(Sample, c("Treatment", "Replicate"), "(^.+)_(.+)") %>%
  dplyr::mutate(Treatment = factor(Treatment, levels = c("Control", "HA15")))

p_prop <- ggplot(prop_long, aes(x = Treatment, y = Proportion, fill = Treatment, color = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_point(position = position_jitter(width = 0.08), size = 1.8) +
  scale_fill_manual(values = c("Control" = "grey70", "HA15" = "darkseagreen3")) +
  scale_color_manual(values = c("Control" = "grey40", "HA15" = "darkgreen")) +
  theme_bw() +
  labs(x = NULL, y = "Proportion HIV+")

ggsave("Results/Figures/04_HIVpos_total_by_treat.png", p_prop, width = 4, height = 4, dpi = 300)


## ---- 1.3 Efecto de la reactivación viral (VIH+) del grupo tratado en la composición celular  ----

# DA de anotación en células tratadas: HIV+ vs HIV-
ann_hiv <- subset(ann, subset = Treatment == "HA15")

props <- getTransformedProps(clusters = ann_hiv$functional_ann2,
  sample = paste(ann_hiv$Sample, ann_hiv$HIV_status, sep = "_"), transform = "logit")

sample_info <- unique(data.frame(Sample = paste(ann_hiv$Sample, ann_hiv$HIV_status, sep = "_"),
  HIV_status = ann_hiv$HIV_status, Replicate = ann_hiv$Replicate))

sample_info <- sample_info[match(colnames(props$TransformedProps), sample_info$Sample), ]
rownames(sample_info) <- NULL
sample_info$HIV_status <- factor(sample_info$HIV_status, levels = c("HIV-", "HIV+"), labels = c("HIVneg", "HIVpos"))

pair <- as.numeric(factor(sample_info$Replicate))
design <- model.matrix(~ 0 + HIV_status + pair, data = sample_info)
contrasts <- makeContrasts(HIVpos_vs_HIVneg = HIV_statusHIVpos - HIV_statusHIVneg, levels = design)

res_prop <- propeller.ttest(prop.list = props, design = design, contrasts = contrasts,
                            robust = TRUE, trend = FALSE, sort = TRUE)

write.table(res_prop, "Results/Tables/04_01_res_prop_ann_HA15_HIVpos_vs_HIVneg.tsv",
            sep = "\t", row.names = T, col.names = NA, quote = FALSE)

# figura
proptable <- as.data.frame(table(ann_hiv$functional_ann2, paste(ann_hiv$Sample, ann_hiv$HIV_status, sep = "_"))) %>%
  tibble::as_tibble() %>% dplyr::rename(Sample = Var2) %>% dplyr::left_join(
    tibble::as_tibble(as.data.frame(table(paste(ann_hiv$Sample, ann_hiv$HIV_status, sep = "_")))) %>%
      dplyr::rename(Sample = Var1),
    by = "Sample") %>% dplyr::mutate(Proportion = Freq.x / Freq.y) %>% dplyr::select(-Freq.x, -Freq.y) %>%
  tidyr::pivot_wider(names_from = Sample, values_from = Proportion) %>% dplyr::rename(Cluster = Var1)

prop_long <- proptable %>%
  tidyr::pivot_longer(cols = !Cluster, names_to = "Sample", values_to = "Proportion") %>%
  tidyr::extract(Sample, c("Sample", "HIV_status"), "(.+)_(HIV[+-])$") %>%
  dplyr::mutate(HIV_status = factor(HIV_status, levels = c("HIV-", "HIV+")))

sig <- data.frame(Cluster = rownames(res_prop), FDR = res_prop$FDR)
sig$label <- ifelse(sig$FDR < 0.001, "***", ifelse(sig$FDR < 0.01, "**", ifelse(sig$FDR < 0.05, "*", "")))
sig <- sig[sig$label != "", ]
ymax <- tapply(prop_long$Proportion, prop_long$Cluster, max)
sig$y <- ymax[match(sig$Cluster, names(ymax))] + 0.02

p_prop <- ggplot(prop_long, aes(x = Cluster, y = Proportion, fill = HIV_status, color = HIV_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.75)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), size = 1.5) +
  geom_text(data = sig, aes(x = Cluster, y = y, label = label), inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = c("HIV-" = "grey70", "HIV+" = "indianred3")) +
  scale_color_manual(values = c("HIV-" = "grey40", "HIV+" = "darkred")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Annotation", y = "Proportion")

ggsave("Results/Figures/04_functional_ann_HA15_HIVpos_vs_HIVneg.png",
       p_prop, width = 9, height = 4, dpi = 300)
