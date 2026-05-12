## ---- 3.4 figuras plotenrichment representativas

plotEnrichment(hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], ranks.all[["1"]])

for (cl in as.character(1:11)) {
  p <- plotEnrichment(hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], ranks.all[[cl]]) +
    labs(title = paste0("Respuesta a IFN gamma (Cluster ", cl, ")"))
  ggsave(paste0("Results/Figures/02_04_clust", cl, "_IFNg.png"), p, dpi = 300)
}


# Cluster 1

# HALLMARK

p <- plotEnrichment(hallmark[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]], ranks.all[["1"]]) + labs(title="Respuesta a IFN alfa (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_IFNa_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]], ranks.cons[["1"]]) + labs(title="Respuesta a IFN alfa (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_IFNa_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]], ranks.filt[["1"]]) + labs(title="Respuesta a IFN alfa (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_IFNa_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], ranks.all[["1"]]) + labs(title="Respuesta a IFN gamma (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_IFNg_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], ranks.cons[["1"]]) + labs(title="Respuesta a IFN gamma (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_IFNg_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], ranks.filt[["1"]]) + labs(title="Respuesta a IFN gamma (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_IFNg_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]], ranks.all[["1"]]) + labs(title="Señalización TNFa vía NF-kB (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_TNFa_NFkB_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]], ranks.cons[["1"]]) + labs(title="Señalización TNFa vía NF-kB (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_TNFa_NFkB_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]], ranks.filt[["1"]]) + labs(title="Señalización TNFa vía NF-kB (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_TNFa_NFkB_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_APICAL_JUNCTION"]], ranks.all[["1"]]) + labs(title="Apical junction (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_ApicalJunction_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_APICAL_JUNCTION"]], ranks.cons[["1"]]) + labs(title="Apical junction (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_ApicalJunction_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_APICAL_JUNCTION"]], ranks.filt[["1"]]) + labs(title="Apical junction (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_ApicalJunction_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_HEME_METABOLISM"]], ranks.all[["1"]]) + labs(title="Metabolismo del hemo (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_HemeMetabolism_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_HEME_METABOLISM"]], ranks.cons[["1"]]) + labs(title="Metabolismo del hemo (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_HemeMetabolism_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_HEME_METABOLISM"]], ranks.filt[["1"]]) + labs(title="Metabolismo del hemo (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_HemeMetabolism_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_MYC_TARGETS_V1"]], ranks.all[["1"]]) + labs(title="Dianas de MYC v1 (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_MYCtargetsV1_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_MYC_TARGETS_V1"]], ranks.cons[["1"]]) + labs(title="Dianas de MYC v1 (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_MYCtargetsV1_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_MYC_TARGETS_V1"]], ranks.filt[["1"]]) + labs(title="Dianas de MYC v1 (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_MYCtargetsV1_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.all[["1"]]) + labs(title="Checkpoint G2/M (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_G2M_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.cons[["1"]]) + labs(title="Checkpoint G2/M (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_G2M_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.filt[["1"]]) + labs(title="Checkpoint G2/M (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_G2M_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.all[["1"]]) + labs(title="Dianas de E2F (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_E2F_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.cons[["1"]]) + labs(title="Dianas de E2F (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_E2F_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.filt[["1"]]) + labs(title="Dianas de E2F (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_E2F_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_MITOTIC_SPINDLE"]], ranks.all[["1"]]) + labs(title="Huso mitótico (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_MitoticSpindle_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_MITOTIC_SPINDLE"]], ranks.cons[["1"]]) + labs(title="Huso mitótico (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_MitoticSpindle_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_MITOTIC_SPINDLE"]], ranks.filt[["1"]]) + labs(title="Huso mitótico (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_MitoticSpindle_filt.png", p, dpi = 300)

# GOBP

p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"]], ranks.all[["1"]]) + labs(title="Procesamiento y presentación de antígeno (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_AgProcessingPresentation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"]], ranks.cons[["1"]]) + labs(title="Procesamiento y presentación de antígeno (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_AgProcessingPresentation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"]], ranks.filt[["1"]]) + labs(title="Procesamiento y presentación de antígeno (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_AgProcessingPresentation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I"]], ranks.all[["1"]]) + labs(title="Presentación Ag vía MHC-I (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_MHCI_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I"]], ranks.cons[["1"]]) + labs(title="Presentación Ag vía MHC-I (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_MHCI_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I"]], ranks.filt[["1"]]) + labs(title="Presentación Ag vía MHC-I (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_MHCI_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_T_CELL_MEDIATED_CYTOTOXICITY"]], ranks.all[["1"]]) + labs(title="Citotoxicidad mediada por T (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_Tcytotoxicity_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_T_CELL_MEDIATED_CYTOTOXICITY"]], ranks.cons[["1"]]) + labs(title="Citotoxicidad mediada por T (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_Tcytotoxicity_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_T_CELL_MEDIATED_CYTOTOXICITY"]], ranks.filt[["1"]]) + labs(title="Citotoxicidad mediada por T (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_Tcytotoxicity_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_EXOCYTOSIS"]], ranks.all[["1"]]) + labs(title="Exocitosis (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_Exocytosis_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_EXOCYTOSIS"]], ranks.cons[["1"]]) + labs(title="Exocitosis (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_Exocytosis_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_EXOCYTOSIS"]], ranks.filt[["1"]]) + labs(title="Exocitosis (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_Exocytosis_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_CELL_CELL_ADHESION"]], ranks.all[["1"]]) + labs(title="Adhesión célula-célula (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_CellCellAdhesion_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_CELL_ADHESION"]], ranks.cons[["1"]]) + labs(title="Adhesión célula-célula (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_CellCellAdhesion_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_CELL_ADHESION"]], ranks.filt[["1"]]) + labs(title="Adhesión célula-célula (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_CellCellAdhesion_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]], ranks.all[["1"]]) + labs(title="Respuesta inmune adaptativa (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_AdaptiveImmuneResponse_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]], ranks.cons[["1"]]) + labs(title="Respuesta inmune adaptativa (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_AdaptiveImmuneResponse_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]], ranks.filt[["1"]]) + labs(title="Respuesta inmune adaptativa (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_AdaptiveImmuneResponse_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR"]], ranks.all[["1"]]) + labs(title="Respuesta a TNF (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_TNFresponse_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR"]], ranks.cons[["1"]]) + labs(title="Respuesta a TNF (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_TNFresponse_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR"]], ranks.filt[["1"]]) + labs(title="Respuesta a TNF (Cluster 1)")
ggsave("Results/Figures/02_04_fgsea_clust1_TNFresponse_filt.png", p, dpi = 300)


# Cluster 2

# HALLMARK

p <- plotEnrichment(hallmark[["HALLMARK_MYC_TARGETS_V1"]], ranks.all[["2"]]) + labs(title="MYC targets v1 (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_MYCtargetsV1_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_MYC_TARGETS_V1"]], ranks.cons[["2"]]) + labs(title="MYC targets v1 (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_MYCtargetsV1_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_MYC_TARGETS_V1"]], ranks.filt[["2"]]) + labs(title="MYC targets v1 (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_MYCtargetsV1_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_ALLOGRAFT_REJECTION"]], ranks.all[["2"]]) + labs(title="Allograft rejection (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Allograft_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_ALLOGRAFT_REJECTION"]], ranks.cons[["2"]]) + labs(title="Allograft rejection (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Allograft_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_ALLOGRAFT_REJECTION"]], ranks.filt[["2"]]) + labs(title="Allograft rejection (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Allograft_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.all[["2"]]) + labs(title="E2F targets (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_E2F_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.cons[["2"]]) + labs(title="E2F targets (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_E2F_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.filt[["2"]]) + labs(title="E2F targets (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_E2F_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.all[["2"]]) + labs(title="G2M checkpoint (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_G2M_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.cons[["2"]]) + labs(title="G2M checkpoint (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_G2M_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.filt[["2"]]) + labs(title="G2M checkpoint (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_G2M_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_APOPTOSIS"]], ranks.all[["2"]]) + labs(title="Apoptosis (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Apoptosis_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_APOPTOSIS"]], ranks.cons[["2"]]) + labs(title="Apoptosis (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Apoptosis_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_APOPTOSIS"]], ranks.filt[["2"]]) + labs(title="Apoptosis (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Apoptosis_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], ranks.all[["2"]]) + labs(title="IFN gamma (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_IFNg_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], ranks.cons[["2"]]) + labs(title="IFN gamma (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_IFNg_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], ranks.filt[["2"]]) + labs(title="IFN gamma (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_IFNg_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_IL2_STAT5_SIGNALING"]], ranks.all[["2"]]) + labs(title="IL2-STAT5 signaling (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_IL2STAT5_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_IL2_STAT5_SIGNALING"]], ranks.cons[["2"]]) + labs(title="IL2-STAT5 signaling (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_IL2STAT5_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_IL2_STAT5_SIGNALING"]], ranks.filt[["2"]]) + labs(title="IL2-STAT5 signaling (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_IL2STAT5_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_COMPLEMENT"]], ranks.all[["2"]]) + labs(title="Complement (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Complement_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_COMPLEMENT"]], ranks.cons[["2"]]) + labs(title="Complement (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Complement_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_COMPLEMENT"]], ranks.filt[["2"]]) + labs(title="Complement (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Complement_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_MITOTIC_SPINDLE"]], ranks.all[["2"]]) + labs(title="Mitotic spindle (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_MitoticSpindle_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_MITOTIC_SPINDLE"]], ranks.cons[["2"]]) + labs(title="Mitotic spindle (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_MitoticSpindle_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_MITOTIC_SPINDLE"]], ranks.filt[["2"]]) + labs(title="Mitotic spindle (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_MitoticSpindle_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_GLYCOLYSIS"]], ranks.all[["2"]]) + labs(title="Glycolysis (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Glycolysis_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_GLYCOLYSIS"]], ranks.cons[["2"]]) + labs(title="Glycolysis (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Glycolysis_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_GLYCOLYSIS"]], ranks.filt[["2"]]) + labs(title="Glycolysis (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Glycolysis_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_ESTROGEN_RESPONSE_LATE"]], ranks.all[["2"]]) + labs(title="Estrogen response late (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Estrogen_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_ESTROGEN_RESPONSE_LATE"]], ranks.cons[["2"]]) + labs(title="Estrogen response late (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Estrogen_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_ESTROGEN_RESPONSE_LATE"]], ranks.filt[["2"]]) + labs(title="Estrogen response late (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Estrogen_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_HYPOXIA"]], ranks.all[["2"]]) + labs(title="Hypoxia (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Hypoxia_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_HYPOXIA"]], ranks.cons[["2"]]) + labs(title="Hypoxia (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Hypoxia_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_HYPOXIA"]], ranks.filt[["2"]]) + labs(title="Hypoxia (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Hypoxia_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_MTORC1_SIGNALING"]], ranks.all[["2"]]) + labs(title="mTORC1 signaling (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_mTORC1_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_MTORC1_SIGNALING"]], ranks.cons[["2"]]) + labs(title="mTORC1 signaling (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_mTORC1_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_MTORC1_SIGNALING"]], ranks.filt[["2"]]) + labs(title="mTORC1 signaling (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_mTORC1_filt.png", p, dpi = 300)


# Hallmark

p <- plotEnrichment(gobp[["GOBP_CELL_CYCLE"]], ranks.all[["2"]]) + labs(title="Cell cycle (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_CellCycle_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_CYCLE"]], ranks.cons[["2"]]) + labs(title="Cell cycle (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_CellCycle_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_CYCLE"]], ranks.filt[["2"]]) + labs(title="Cell cycle (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_CellCycle_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_INFLAMMATORY_RESPONSE"]], ranks.all[["2"]]) + labs(title="Inflammatory response (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_InflammatoryResponse_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_INFLAMMATORY_RESPONSE"]], ranks.cons[["2"]]) + labs(title="Inflammatory response (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_InflammatoryResponse_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_INFLAMMATORY_RESPONSE"]], ranks.filt[["2"]]) + labs(title="Inflammatory response (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_InflammatoryResponse_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_TRANSLATION"]], ranks.all[["2"]]) + labs(title="Translation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Translation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_TRANSLATION"]], ranks.cons[["2"]]) + labs(title="Translation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Translation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_TRANSLATION"]], ranks.filt[["2"]]) + labs(title="Translation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_Translation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_TRANSLATIONAL_ELONGATION"]], ranks.all[["2"]]) + labs(title="Translational elongation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_TranslationalElongation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_TRANSLATIONAL_ELONGATION"]], ranks.cons[["2"]]) + labs(title="Translational elongation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_TranslationalElongation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_TRANSLATIONAL_ELONGATION"]], ranks.filt[["2"]]) + labs(title="Translational elongation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_TranslationalElongation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_RIBOSOME_ASSEMBLY"]], ranks.all[["2"]]) + labs(title="Ribosome assembly (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RibosomeAssembly_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_RIBOSOME_ASSEMBLY"]], ranks.cons[["2"]]) + labs(title="Ribosome assembly (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RibosomeAssembly_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_RIBOSOME_ASSEMBLY"]], ranks.filt[["2"]]) + labs(title="Ribosome assembly (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RibosomeAssembly_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_RIBOSOME_BIOGENESIS"]], ranks.all[["2"]]) + labs(title="Ribosome biogenesis (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RibosomeBiogenesis_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_RIBOSOME_BIOGENESIS"]], ranks.cons[["2"]]) + labs(title="Ribosome biogenesis (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RibosomeBiogenesis_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_RIBOSOME_BIOGENESIS"]], ranks.filt[["2"]]) + labs(title="Ribosome biogenesis (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RibosomeBiogenesis_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS"]], ranks.all[["2"]]) + labs(title="Ribonucleoprotein complex biogenesis (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RNPBiogenesis_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS"]], ranks.cons[["2"]]) + labs(title="Ribonucleoprotein complex biogenesis (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RNPBiogenesis_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS"]], ranks.filt[["2"]]) + labs(title="Ribonucleoprotein complex biogenesis (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RNPBiogenesis_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS"]], ranks.all[["2"]]) + labs(title="Ribosomal small subunit biogenesis (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_SmallSubunitBiogenesis_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS"]], ranks.cons[["2"]]) + labs(title="Ribosomal small subunit biogenesis (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_SmallSubunitBiogenesis_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS"]], ranks.filt[["2"]]) + labs(title="Ribosomal small subunit biogenesis (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_SmallSubunitBiogenesis_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_RRNA_METABOLIC_PROCESS"]], ranks.all[["2"]]) + labs(title="rRNA metabolic process (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_rRNAMetabolism_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_RRNA_METABOLIC_PROCESS"]], ranks.cons[["2"]]) + labs(title="rRNA metabolic process (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_rRNAMetabolism_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_RRNA_METABOLIC_PROCESS"]], ranks.filt[["2"]]) + labs(title="rRNA metabolic process (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_rRNAMetabolism_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_SYNCYTIUM_FORMATION"]], ranks.all[["2"]]) + labs(title="Syncytium formation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_SyncytiumFormation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_SYNCYTIUM_FORMATION"]], ranks.cons[["2"]]) + labs(title="Syncytium formation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_SyncytiumFormation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_SYNCYTIUM_FORMATION"]], ranks.filt[["2"]]) + labs(title="Syncytium formation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_SyncytiumFormation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_SYNCYTIUM_FORMATION_BY_PLASMA_MEMBRANE_FUSION"]], ranks.all[["2"]]) + labs(title="Negative regulation of syncytium formation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_NegSyncytium_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_SYNCYTIUM_FORMATION_BY_PLASMA_MEMBRANE_FUSION"]], ranks.cons[["2"]]) + labs(title="Negative regulation of syncytium formation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_NegSyncytium_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_SYNCYTIUM_FORMATION_BY_PLASMA_MEMBRANE_FUSION"]], ranks.filt[["2"]]) + labs(title="Negative regulation of syncytium formation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_NegSyncytium_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_MUSCLE_CELL_DIFFERENTIATION"]], ranks.all[["2"]]) + labs(title="Muscle cell differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_MuscleDifferentiation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_MUSCLE_CELL_DIFFERENTIATION"]], ranks.cons[["2"]]) + labs(title="Muscle cell differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_MuscleDifferentiation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_MUSCLE_CELL_DIFFERENTIATION"]], ranks.filt[["2"]]) + labs(title="Muscle cell differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_MuscleDifferentiation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_STRIATED_MUSCLE_CELL_DIFFERENTIATION"]], ranks.all[["2"]]) + labs(title="Striated muscle cell differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_StriatedMuscleDifferentiation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_STRIATED_MUSCLE_CELL_DIFFERENTIATION"]], ranks.cons[["2"]]) + labs(title="Striated muscle cell differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_StriatedMuscleDifferentiation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_STRIATED_MUSCLE_CELL_DIFFERENTIATION"]], ranks.filt[["2"]]) + labs(title="Striated muscle cell differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_StriatedMuscleDifferentiation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_MITOTIC_CELL_CYCLE"]], ranks.all[["2"]]) + labs(title="Mitotic cell cycle (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_MitoticCellCycle_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_MITOTIC_CELL_CYCLE"]], ranks.cons[["2"]]) + labs(title="Mitotic cell cycle (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_MitoticCellCycle_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_MITOTIC_CELL_CYCLE"]], ranks.filt[["2"]]) + labs(title="Mitotic cell cycle (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_MitoticCellCycle_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_CELL_DIVISION"]], ranks.all[["2"]]) + labs(title="Cell division (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_CellDivision_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_DIVISION"]], ranks.cons[["2"]]) + labs(title="Cell division (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_CellDivision_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_DIVISION"]], ranks.filt[["2"]]) + labs(title="Cell division (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_CellDivision_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_DIVISION"]], ranks.all[["2"]]) + labs(title="Positive regulation of cell division (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCellDivision_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_DIVISION"]], ranks.cons[["2"]]) + labs(title="Positive regulation of cell division (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCellDivision_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_DIVISION"]], ranks.filt[["2"]]) + labs(title="Positive regulation of cell division (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCellDivision_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION"]], ranks.all[["2"]]) + labs(title="Positive regulation of leukocyte proliferation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosLeukocyteProlif_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION"]], ranks.cons[["2"]]) + labs(title="Positive regulation of leukocyte proliferation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosLeukocyteProlif_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION"]], ranks.filt[["2"]]) + labs(title="Positive regulation of leukocyte proliferation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosLeukocyteProlif_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_POPULATION_PROLIFERATION"]], ranks.all[["2"]]) + labs(title="Positive regulation of cell population proliferation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosPopulationProlif_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_POPULATION_PROLIFERATION"]], ranks.cons[["2"]]) + labs(title="Positive regulation of cell population proliferation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosPopulationProlif_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_POPULATION_PROLIFERATION"]], ranks.filt[["2"]]) + labs(title="Positive regulation of cell population proliferation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosPopulationProlif_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_T_CELL_PROLIFERATION"]], ranks.all[["2"]]) + labs(title="Positive regulation of T cell proliferation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosTcellProlif_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_T_CELL_PROLIFERATION"]], ranks.cons[["2"]]) + labs(title="Positive regulation of T cell proliferation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosTcellProlif_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_T_CELL_PROLIFERATION"]], ranks.filt[["2"]]) + labs(title="Positive regulation of T cell proliferation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosTcellProlif_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE"]], ranks.all[["2"]]) + labs(title="Cell activation involved in immune response (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_CellActivationImmune_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE"]], ranks.cons[["2"]]) + labs(title="Cell activation involved in immune response (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_CellActivationImmune_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE"]], ranks.filt[["2"]]) + labs(title="Cell activation involved in immune response (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_CellActivationImmune_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_ACTIVATION"]], ranks.all[["2"]]) + labs(title="Positive regulation of cell activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCellActivation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_ACTIVATION"]], ranks.cons[["2"]]) + labs(title="Positive regulation of cell activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCellActivation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_ACTIVATION"]], ranks.filt[["2"]]) + labs(title="Positive regulation of cell activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCellActivation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_LYMPHOCYTE_ACTIVATION"]], ranks.all[["2"]]) + labs(title="Lymphocyte activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_LymphocyteActivation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_LYMPHOCYTE_ACTIVATION"]], ranks.cons[["2"]]) + labs(title="Lymphocyte activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_LymphocyteActivation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_LYMPHOCYTE_ACTIVATION"]], ranks.filt[["2"]]) + labs(title="Lymphocyte activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_LymphocyteActivation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_REGULATION_OF_LYMPHOCYTE_ACTIVATION"]], ranks.all[["2"]]) + labs(title="Regulation of lymphocyte activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RegLymphocyteActivation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_REGULATION_OF_LYMPHOCYTE_ACTIVATION"]], ranks.cons[["2"]]) + labs(title="Regulation of lymphocyte activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RegLymphocyteActivation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_REGULATION_OF_LYMPHOCYTE_ACTIVATION"]], ranks.filt[["2"]]) + labs(title="Regulation of lymphocyte activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RegLymphocyteActivation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_T_CELL_ACTIVATION"]], ranks.all[["2"]]) + labs(title="T cell activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_TcellActivation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_T_CELL_ACTIVATION"]], ranks.cons[["2"]]) + labs(title="T cell activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_TcellActivation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_T_CELL_ACTIVATION"]], ranks.filt[["2"]]) + labs(title="T cell activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_TcellActivation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION"]], ranks.all[["2"]]) + labs(title="CD4 alpha-beta T cell activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_CD4TcellActivation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION"]], ranks.cons[["2"]]) + labs(title="CD4 alpha-beta T cell activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_CD4TcellActivation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION"]], ranks.filt[["2"]]) + labs(title="CD4 alpha-beta T cell activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_CD4TcellActivation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_ALPHA_BETA_T_CELL_ACTIVATION"]], ranks.all[["2"]]) + labs(title="Alpha-beta T cell activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_AlphaBetaTcellActivation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ALPHA_BETA_T_CELL_ACTIVATION"]], ranks.cons[["2"]]) + labs(title="Alpha-beta T cell activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_AlphaBetaTcellActivation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ALPHA_BETA_T_CELL_ACTIVATION"]], ranks.filt[["2"]]) + labs(title="Alpha-beta T cell activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_AlphaBetaTcellActivation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_REGULATION_OF_ALPHA_BETA_T_CELL_ACTIVATION"]], ranks.all[["2"]]) + labs(title="Regulation of alpha-beta T cell activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RegAlphaBetaTcellActivation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_REGULATION_OF_ALPHA_BETA_T_CELL_ACTIVATION"]], ranks.cons[["2"]]) + labs(title="Regulation of alpha-beta T cell activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RegAlphaBetaTcellActivation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_REGULATION_OF_ALPHA_BETA_T_CELL_ACTIVATION"]], ranks.filt[["2"]]) + labs(title="Regulation of alpha-beta T cell activation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RegAlphaBetaTcellActivation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY"]], ranks.all[["2"]]) + labs(title="T cell receptor signaling pathway (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_TCRsignaling_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY"]], ranks.cons[["2"]]) + labs(title="T cell receptor signaling pathway (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_TCRsignaling_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY"]], ranks.filt[["2"]]) + labs(title="T cell receptor signaling pathway (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_TCRsignaling_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]], ranks.all[["2"]]) + labs(title="Adaptive immune response (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_AdaptiveImmune_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]], ranks.cons[["2"]]) + labs(title="Adaptive immune response (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_AdaptiveImmune_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]], ranks.filt[["2"]]) + labs(title="Adaptive immune response (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_AdaptiveImmune_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_IMMUNE_EFFECTOR_PROCESS"]], ranks.all[["2"]]) + labs(title="Immune effector process (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_ImmuneEffector_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_IMMUNE_EFFECTOR_PROCESS"]], ranks.cons[["2"]]) + labs(title="Immune effector process (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_ImmuneEffector_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_IMMUNE_EFFECTOR_PROCESS"]], ranks.filt[["2"]]) + labs(title="Immune effector process (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_ImmuneEffector_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_LEUKOCYTE_MEDIATED_IMMUNITY"]], ranks.all[["2"]]) + labs(title="Leukocyte-mediated immunity (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_LeukocyteImmunity_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_LEUKOCYTE_MEDIATED_IMMUNITY"]], ranks.cons[["2"]]) + labs(title="Leukocyte-mediated immunity (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_LeukocyteImmunity_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_LEUKOCYTE_MEDIATED_IMMUNITY"]], ranks.filt[["2"]]) + labs(title="Leukocyte-mediated immunity (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_LeukocyteImmunity_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY"]], ranks.all[["2"]]) + labs(title="Lymphocyte-mediated immunity (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_LymphocyteImmunity_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY"]], ranks.cons[["2"]]) + labs(title="Lymphocyte-mediated immunity (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_LymphocyteImmunity_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY"]], ranks.filt[["2"]]) + labs(title="Lymphocyte-mediated immunity (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_LymphocyteImmunity_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE"]], ranks.all[["2"]]) + labs(title="Positive regulation of immune response (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosImmuneResponse_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE"]], ranks.cons[["2"]]) + labs(title="Positive regulation of immune response (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosImmuneResponse_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE"]], ranks.filt[["2"]]) + labs(title="Positive regulation of immune response (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosImmuneResponse_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_ADHESION"]], ranks.all[["2"]]) + labs(title="Positive regulation of cell adhesion (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCellAdhesion_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_ADHESION"]], ranks.cons[["2"]]) + labs(title="Positive regulation of cell adhesion (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCellAdhesion_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_ADHESION"]], ranks.filt[["2"]]) + labs(title="Positive regulation of cell adhesion (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCellAdhesion_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_CELL_ADHESION"]], ranks.all[["2"]]) + labs(title="Positive regulation of cell-cell adhesion (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCellCellAdhesion_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_CELL_ADHESION"]], ranks.cons[["2"]]) + labs(title="Positive regulation of cell-cell adhesion (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCellCellAdhesion_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_CELL_ADHESION"]], ranks.filt[["2"]]) + labs(title="Positive regulation of cell-cell adhesion (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCellCellAdhesion_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION"]], ranks.all[["2"]]) + labs(title="Positive regulation of leukocyte cell-cell adhesion (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosLeukocyteCellCellAdhesion_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION"]], ranks.cons[["2"]]) + labs(title="Positive regulation of leukocyte cell-cell adhesion (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosLeukocyteCellCellAdhesion_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION"]], ranks.filt[["2"]]) + labs(title="Positive regulation of leukocyte cell-cell adhesion (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosLeukocyteCellCellAdhesion_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_LEUKOCYTE_DIFFERENTIATION"]], ranks.all[["2"]]) + labs(title="Leukocyte differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_LeukocyteDifferentiation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_LEUKOCYTE_DIFFERENTIATION"]], ranks.cons[["2"]]) + labs(title="Leukocyte differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_LeukocyteDifferentiation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_LEUKOCYTE_DIFFERENTIATION"]], ranks.filt[["2"]]) + labs(title="Leukocyte differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_LeukocyteDifferentiation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_DIFFERENTIATION"]], ranks.all[["2"]]) + labs(title="Positive regulation of cell differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCellDifferentiation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_DIFFERENTIATION"]], ranks.cons[["2"]]) + labs(title="Positive regulation of cell differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCellDifferentiation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_DIFFERENTIATION"]], ranks.filt[["2"]]) + labs(title="Positive regulation of cell differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCellDifferentiation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_REGULATION_OF_ALPHA_BETA_T_CELL_DIFFERENTIATION"]], ranks.all[["2"]]) + labs(title="Regulation of alpha-beta T cell differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RegAlphaBetaTcellDifferentiation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_REGULATION_OF_ALPHA_BETA_T_CELL_DIFFERENTIATION"]], ranks.cons[["2"]]) + labs(title="Regulation of alpha-beta T cell differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RegAlphaBetaTcellDifferentiation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_REGULATION_OF_ALPHA_BETA_T_CELL_DIFFERENTIATION"]], ranks.filt[["2"]]) + labs(title="Regulation of alpha-beta T cell differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RegAlphaBetaTcellDifferentiation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_MEMORY_T_CELL_DIFFERENTIATION"]], ranks.all[["2"]]) + labs(title="Positive regulation of memory T cell differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosMemoryTcellDifferentiation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_MEMORY_T_CELL_DIFFERENTIATION"]], ranks.cons[["2"]]) + labs(title="Positive regulation of memory T cell differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosMemoryTcellDifferentiation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_MEMORY_T_CELL_DIFFERENTIATION"]], ranks.filt[["2"]]) + labs(title="Positive regulation of memory T cell differentiation (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosMemoryTcellDifferentiation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION"]], ranks.all[["2"]]) + labs(title="Positive regulation of cytokine production (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCytokineProduction_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION"]], ranks.cons[["2"]]) + labs(title="Positive regulation of cytokine production (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCytokineProduction_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION"]], ranks.filt[["2"]]) + labs(title="Positive regulation of cytokine production (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_PosCytokineProduction_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_REGULATION_OF_PROGRAMMED_CELL_DEATH"]], ranks.all[["2"]]) + labs(title="Regulation of programmed cell death (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RegProgrammedCellDeath_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_REGULATION_OF_PROGRAMMED_CELL_DEATH"]], ranks.cons[["2"]]) + labs(title="Regulation of programmed cell death (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RegProgrammedCellDeath_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_REGULATION_OF_PROGRAMMED_CELL_DEATH"]], ranks.filt[["2"]]) + labs(title="Regulation of programmed cell death (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_RegProgrammedCellDeath_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_GRANZYME_MEDIATED_PROGRAMMED_CELL_DEATH_SIGNALING_PATHWAY"]], ranks.all[["2"]]) + labs(title="Granzyme-mediated programmed cell death signaling pathway (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_GranzymeCellDeath_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_GRANZYME_MEDIATED_PROGRAMMED_CELL_DEATH_SIGNALING_PATHWAY"]], ranks.cons[["2"]]) + labs(title="Granzyme-mediated programmed cell death signaling pathway (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_GranzymeCellDeath_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_GRANZYME_MEDIATED_PROGRAMMED_CELL_DEATH_SIGNALING_PATHWAY"]], ranks.filt[["2"]]) + labs(title="Granzyme-mediated programmed cell death signaling pathway (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_GranzymeCellDeath_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II"]], ranks.all[["2"]]) + labs(title="Antigen processing and presentation via MHC-II (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_MHCII_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II"]], ranks.cons[["2"]]) + labs(title="Antigen processing and presentation via MHC-II (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_MHCII_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II"]], ranks.filt[["2"]]) + labs(title="Antigen processing and presentation via MHC-II (Cluster 2)")
ggsave("Results/Figures/02_04_fgsea_clust2_MHCII_filt.png", p, dpi = 300)

# cluster 3

p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_SYNCYTIUM_FORMATION_BY_PLASMA_MEMBRANE_FUSION"]], ranks.all[["3"]]) + labs(title="Negative regulation of syncytium formation (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_NegSyncytium_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_SYNCYTIUM_FORMATION_BY_PLASMA_MEMBRANE_FUSION"]], ranks.cons[["3"]]) + labs(title="Negative regulation of syncytium formation (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_NegSyncytium_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_SYNCYTIUM_FORMATION_BY_PLASMA_MEMBRANE_FUSION"]], ranks.filt[["3"]]) + labs(title="Negative regulation of syncytium formation (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_NegSyncytium_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_MYOBLAST_FUSION"]], ranks.all[["3"]]) + labs(title="Myoblast fusion (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_MyoblastFusion_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_MYOBLAST_FUSION"]], ranks.cons[["3"]]) + labs(title="Myoblast fusion (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_MyoblastFusion_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_MYOBLAST_FUSION"]], ranks.filt[["3"]]) + labs(title="Myoblast fusion (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_MyoblastFusion_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS"]], ranks.all[["3"]]) + labs(title="Ribosomal small subunit biogenesis (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_SmallSubunitBiogenesis_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS"]], ranks.cons[["3"]]) + labs(title="Ribosomal small subunit biogenesis (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_SmallSubunitBiogenesis_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS"]], ranks.filt[["3"]]) + labs(title="Ribosomal small subunit biogenesis (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_SmallSubunitBiogenesis_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_VIA_MHC_CLASS_IB"]], ranks.all[["3"]]) + labs(title="Antigen processing and presentation via MHC class Ib (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_MHCIb_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_VIA_MHC_CLASS_IB"]], ranks.cons[["3"]]) + labs(title="Antigen processing and presentation via MHC class Ib (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_MHCIb_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_VIA_MHC_CLASS_IB"]], ranks.filt[["3"]]) + labs(title="Antigen processing and presentation via MHC class Ib (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_MHCIb_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_TRANSLATION"]], ranks.all[["3"]]) + labs(title="Translation (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_Translation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_TRANSLATION"]], ranks.cons[["3"]]) + labs(title="Translation (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_Translation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_TRANSLATION"]], ranks.filt[["3"]]) + labs(title="Translation (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_Translation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_STEM_CELL_DIFFERENTIATION"]], ranks.all[["3"]]) + labs(title="Negative regulation of stem cell differentiation (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_NegStemDiff_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_STEM_CELL_DIFFERENTIATION"]], ranks.cons[["3"]]) + labs(title="Negative regulation of stem cell differentiation (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_NegStemDiff_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_STEM_CELL_DIFFERENTIATION"]], ranks.filt[["3"]]) + labs(title="Negative regulation of stem cell differentiation (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_NegStemDiff_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_MITOTIC_CELL_CYCLE"]], ranks.all[["3"]]) + labs(title="Mitotic cell cycle (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_MitoticCellCycle_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_MITOTIC_CELL_CYCLE"]], ranks.cons[["3"]]) + labs(title="Mitotic cell cycle (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_MitoticCellCycle_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_MITOTIC_CELL_CYCLE"]], ranks.filt[["3"]]) + labs(title="Mitotic cell cycle (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_MitoticCellCycle_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_IL2_STAT5_SIGNALING"]], ranks.all[["3"]]) + labs(title="IL2-STAT5 signaling (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_IL2STAT5_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_IL2_STAT5_SIGNALING"]], ranks.cons[["3"]]) + labs(title="IL2-STAT5 signaling (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_IL2STAT5_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_IL2_STAT5_SIGNALING"]], ranks.filt[["3"]]) + labs(title="IL2-STAT5 signaling (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_IL2STAT5_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_CELL_DIVISION"]], ranks.all[["3"]]) + labs(title="Cell division (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_CellDivision_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_DIVISION"]], ranks.cons[["3"]]) + labs(title="Cell division (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_CellDivision_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_DIVISION"]], ranks.filt[["3"]]) + labs(title="Cell division (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_CellDivision_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.all[["3"]]) + labs(title="E2F targets (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_E2F_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.cons[["3"]]) + labs(title="E2F targets (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_E2F_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.filt[["3"]]) + labs(title="E2F targets (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_E2F_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_CELL_CYCLE"]], ranks.all[["3"]]) + labs(title="Cell cycle (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_CellCycle_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_CYCLE"]], ranks.cons[["3"]]) + labs(title="Cell cycle (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_CellCycle_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_CYCLE"]], ranks.filt[["3"]]) + labs(title="Cell cycle (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_CellCycle_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_CHROMOSOME_SEGREGATION"]], ranks.all[["3"]]) + labs(title="Chromosome segregation (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_ChromosomeSegregation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CHROMOSOME_SEGREGATION"]], ranks.cons[["3"]]) + labs(title="Chromosome segregation (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_ChromosomeSegregation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CHROMOSOME_SEGREGATION"]], ranks.filt[["3"]]) + labs(title="Chromosome segregation (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_ChromosomeSegregation_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.all[["3"]]) + labs(title="G2M checkpoint (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_G2M_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.cons[["3"]]) + labs(title="G2M checkpoint (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_G2M_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.filt[["3"]]) + labs(title="G2M checkpoint (Cluster 3)")
ggsave("Results/Figures/02_04_fgsea_clust3_G2M_filt.png", p, dpi = 300)


# cluster 4
p <- plotEnrichment(hallmark[["HALLMARK_COMPLEMENT"]], ranks.all[["4"]]) + labs(title="Complement (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_Complement_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_COMPLEMENT"]], ranks.cons[["4"]]) + labs(title="Complement (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_Complement_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_COMPLEMENT"]], ranks.filt[["4"]]) + labs(title="Complement (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_Complement_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_IL2_STAT5_SIGNALING"]], ranks.all[["4"]]) + labs(title="IL2-STAT5 signaling (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_IL2STAT5_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_IL2_STAT5_SIGNALING"]], ranks.cons[["4"]]) + labs(title="IL2-STAT5 signaling (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_IL2STAT5_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_IL2_STAT5_SIGNALING"]], ranks.filt[["4"]]) + labs(title="IL2-STAT5 signaling (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_IL2STAT5_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_ALLOGRAFT_REJECTION"]], ranks.all[["4"]]) + labs(title="Allograft rejection (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_Allograft_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_ALLOGRAFT_REJECTION"]], ranks.cons[["4"]]) + labs(title="Allograft rejection (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_Allograft_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_ALLOGRAFT_REJECTION"]], ranks.filt[["4"]]) + labs(title="Allograft rejection (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_Allograft_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], ranks.all[["4"]]) + labs(title="IFN gamma (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_IFNg_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], ranks.cons[["4"]]) + labs(title="IFN gamma (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_IFNg_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], ranks.filt[["4"]]) + labs(title="IFN gamma (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_IFNg_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_UNFOLDED_PROTEIN_RESPONSE"]], ranks.all[["4"]]) + labs(title="Unfolded protein response (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_UPR_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_UNFOLDED_PROTEIN_RESPONSE"]], ranks.cons[["4"]]) + labs(title="Unfolded protein response (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_UPR_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_UNFOLDED_PROTEIN_RESPONSE"]], ranks.filt[["4"]]) + labs(title="Unfolded protein response (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_UPR_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.all[["4"]]) + labs(title="E2F targets (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_E2F_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.cons[["4"]]) + labs(title="E2F targets (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_E2F_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.filt[["4"]]) + labs(title="E2F targets (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_E2F_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_MYC_TARGETS_V1"]], ranks.all[["4"]]) + labs(title="MYC targets v1 (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_MYCtargetsV1_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_MYC_TARGETS_V1"]], ranks.cons[["4"]]) + labs(title="MYC targets v1 (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_MYCtargetsV1_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_MYC_TARGETS_V1"]], ranks.filt[["4"]]) + labs(title="MYC targets v1 (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_MYCtargetsV1_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], ranks.all[["4"]]) + labs(title="Oxidative phosphorylation (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_OXPHOS_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], ranks.cons[["4"]]) + labs(title="Oxidative phosphorylation (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_OXPHOS_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], ranks.filt[["4"]]) + labs(title="Oxidative phosphorylation (Cluster 4)")
ggsave("Results/Figures/02_04_fgsea_clust4_OXPHOS_filt.png", p, dpi = 300)

#BP
terms <- c(
  "GOBP_POSITIVE_REGULATION_OF_CELL_ACTIVATION",
  "GOBP_INFLAMMATORY_RESPONSE",
  "GOBP_LEUKOCYTE_CELL_CELL_ADHESION",
  "GOBP_REGULATION_OF_T_CELL_ACTIVATION",
  "GOBP_T_CELL_ACTIVATION",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN",
  "GOBP_REGULATION_OF_LEUKOCYTE_DIFFERENTIATION",
  "GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION",
  "GOBP_NEGATIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
  "GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY",
  "GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_PROTEIN_COMPLEX",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "GOBP_ADAPTIVE_IMMUNE_RESPONSE",
  "GOBP_LYMPHOCYTE_DIFFERENTIATION",
  "GOBP_POSITIVE_REGULATION_OF_MAPK_CASCADE",
  "GOBP_POSITIVE_REGULATION_OF_PROTEOLYSIS",
  "GOBP_REGULATION_OF_ALPHA_BETA_T_CELL_ACTIVATION",
  "GOBP_ADAPTIVE_IMMUNE_RESPONSE_BASED_ON_SOMATIC_RECOMBINATION_OF_IMMUNE_RECEPTORS_BUILT_FROM_IMMUNOGLOBULIN_SUPERFAMILY_DOMAINS",
  "GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_CLASS_II_PROTEIN_COMPLEX",
  "GOBP_POSITIVE_REGULATION_OF_CELL_POPULATION_PROLIFERATION",
  "GOBP_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
  "GOBP_T_CELL_DIFFERENTIATION",
  "GOBP_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION",
  "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_IMMUNE_EFFECTOR_PROCESS",
  "GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION",
  "GOBP_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION",
  "GOBP_RESPONSE_TO_GROWTH_FACTOR",
  "GOBP_REGULATION_OF_CELL_DIFFERENTIATION",
  "GOBP_REGULATION_OF_IMMUNE_RESPONSE",
  "GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
  "GOBP_POSITIVE_REGULATION_OF_RAS_PROTEIN_SIGNAL_TRANSDUCTION",
  "GOBP_NEGATIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_BY_P53_CLASS_MEDIATOR",
  "GOBP_NEGATIVE_REGULATION_OF_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR",
  "GOBP_RIBOSOME_ASSEMBLY",
  "GOBP_OXIDATIVE_PHOSPHORYLATION",
  "GOBP_RRNA_METABOLIC_PROCESS",
  "GOBP_SYNCYTIUM_FORMATION",
  "GOBP_RIBOSOME_BIOGENESIS",
  "GOBP_REGULATION_OF_SYNCYTIUM_FORMATION_BY_PLASMA_MEMBRANE_FUSION",
  "GOBP_TRANSLATION",
  "GOBP_NEGATIVE_REGULATION_OF_SYNCYTIUM_FORMATION_BY_PLASMA_MEMBRANE_FUSION",
  "GOBP_CYTOPLASMIC_TRANSLATION"
)

for (term in terms) {
  
  p <- plotEnrichment(gobp[[term]], ranks.all[["4"]]) + 
    labs(title = paste0(gsub("GOBP_", "", term) |> gsub("_", " ", x = _), " (Cluster 4)"))
  ggsave(paste0("Results/Figures/02_04_fgsea_clust4_", gsub("GOBP_", "", term), "_all.png"), p, dpi = 300)
  
  p <- plotEnrichment(gobp[[term]], ranks.cons[["4"]]) + 
    labs(title = paste0(gsub("GOBP_", "", term) |> gsub("_", " ", x = _), " (Cluster 4)"))
  ggsave(paste0("Results/Figures/02_04_fgsea_clust4_", gsub("GOBP_", "", term), "_cons.png"), p, dpi = 300)
  
  p <- plotEnrichment(gobp[[term]], ranks.filt[["4"]]) + 
    labs(title = paste0(gsub("GOBP_", "", term) |> gsub("_", " ", x = _), " (Cluster 4)"))
  ggsave(paste0("Results/Figures/02_04_fgsea_clust4_", gsub("GOBP_", "", term), "_filt.png"), p, dpi = 300)
  
}


# cluster 5

p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.all[["5"]]) + labs(title="G2M checkpoint (Cluster 5)")
ggsave("Results/Figures/02_04_fgsea_clust5_G2M_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.cons[["5"]]) + labs(title="G2M checkpoint (Cluster 5)")
ggsave("Results/Figures/02_04_fgsea_clust5_G2M_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.filt[["5"]]) + labs(title="G2M checkpoint (Cluster 5)")
ggsave("Results/Figures/02_04_fgsea_clust5_G2M_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_CELL_DIVISION"]], ranks.all[["5"]]) + labs(title="Cell division (Cluster 5)")
ggsave("Results/Figures/02_04_fgsea_clust5_CellDivision_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_DIVISION"]], ranks.cons[["5"]]) + labs(title="Cell division (Cluster 5)")
ggsave("Results/Figures/02_04_fgsea_clust5_CellDivision_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_DIVISION"]], ranks.filt[["5"]]) + labs(title="Cell division (Cluster 5)")
ggsave("Results/Figures/02_04_fgsea_clust5_CellDivision_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_CELL_CYCLE"]], ranks.all[["5"]]) + labs(title="Cell cycle (Cluster 5)")
ggsave("Results/Figures/02_04_fgsea_clust5_CellCycle_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_CYCLE"]], ranks.cons[["5"]]) + labs(title="Cell cycle (Cluster 5)")
ggsave("Results/Figures/02_04_fgsea_clust5_CellCycle_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_CYCLE"]], ranks.filt[["5"]]) + labs(title="Cell cycle (Cluster 5)")
ggsave("Results/Figures/02_04_fgsea_clust5_CellCycle_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.all[["5"]]) + labs(title="E2F targets (Cluster 5)")
ggsave("Results/Figures/02_04_fgsea_clust5_E2F_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.cons[["5"]]) + labs(title="E2F targets (Cluster 5)")
ggsave("Results/Figures/02_04_fgsea_clust5_E2F_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.filt[["5"]]) + labs(title="E2F targets (Cluster 5)")
ggsave("Results/Figures/02_04_fgsea_clust5_E2F_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_CHROMOSOME_SEGREGATION"]], ranks.all[["5"]]) + labs(title="Chromosome segregation (Cluster 5)")
ggsave("Results/Figures/02_04_fgsea_clust5_ChromosomeSegregation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CHROMOSOME_SEGREGATION"]], ranks.cons[["5"]]) + labs(title="Chromosome segregation (Cluster 5)")
ggsave("Results/Figures/02_04_fgsea_clust5_ChromosomeSegregation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CHROMOSOME_SEGREGATION"]], ranks.filt[["5"]]) + labs(title="Chromosome segregation (Cluster 5)")
ggsave("Results/Figures/02_04_fgsea_clust5_ChromosomeSegregation_filt.png", p, dpi = 300)


# cluster 6

p <- plotEnrichment(hallmark[["HALLMARK_ANDROGEN_RESPONSE"]], ranks.all[["6"]]) + labs(title="Androgen response (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_Androgen_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_ANDROGEN_RESPONSE"]], ranks.cons[["6"]]) + labs(title="Androgen response (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_Androgen_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_ANDROGEN_RESPONSE"]], ranks.filt[["6"]]) + labs(title="Androgen response (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_Androgen_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_MITOTIC_SPINDLE"]], ranks.all[["6"]]) + labs(title="Mitotic spindle (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_MitoticSpindle_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_MITOTIC_SPINDLE"]], ranks.cons[["6"]]) + labs(title="Mitotic spindle (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_MitoticSpindle_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_MITOTIC_SPINDLE"]], ranks.filt[["6"]]) + labs(title="Mitotic spindle (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_MitoticSpindle_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.all[["6"]]) + labs(title="E2F targets (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_E2F_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.cons[["6"]]) + labs(title="E2F targets (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_E2F_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.filt[["6"]]) + labs(title="E2F targets (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_E2F_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.all[["6"]]) + labs(title="G2M checkpoint (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_G2M_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.cons[["6"]]) + labs(title="G2M checkpoint (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_G2M_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.filt[["6"]]) + labs(title="G2M checkpoint (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_G2M_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POST_TRANSLATIONAL_PROTEIN_TARGETING_TO_ENDOPLASMIC_RETICULUM_MEMBRANE"]], ranks.all[["6"]]) + labs(title="Post-translational protein targeting to ER membrane (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PostTranslERmem_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POST_TRANSLATIONAL_PROTEIN_TARGETING_TO_ENDOPLASMIC_RETICULUM_MEMBRANE"]], ranks.cons[["6"]]) + labs(title="Post-translational protein targeting to ER membrane (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PostTranslERmem_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POST_TRANSLATIONAL_PROTEIN_TARGETING_TO_ENDOPLASMIC_RETICULUM_MEMBRANE"]], ranks.filt[["6"]]) + labs(title="Post-translational protein targeting to ER membrane (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PostTranslERmem_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POST_TRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE_TRANSLOCATION"]], ranks.all[["6"]]) + labs(title="Post-translational protein targeting to membrane translocation (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PostTranslMemTranslocation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POST_TRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE_TRANSLOCATION"]], ranks.cons[["6"]]) + labs(title="Post-translational protein targeting to membrane translocation (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PostTranslMemTranslocation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POST_TRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE_TRANSLOCATION"]], ranks.filt[["6"]]) + labs(title="Post-translational protein targeting to membrane translocation (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PostTranslMemTranslocation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_CYTOPLASMIC_TRANSLATION"]], ranks.all[["6"]]) + labs(title="Cytoplasmic translation (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_CytoplasmicTranslation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CYTOPLASMIC_TRANSLATION"]], ranks.cons[["6"]]) + labs(title="Cytoplasmic translation (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_CytoplasmicTranslation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CYTOPLASMIC_TRANSLATION"]], ranks.filt[["6"]]) + labs(title="Cytoplasmic translation (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_CytoplasmicTranslation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY"]], ranks.all[["6"]]) + labs(title="Negative regulation of intrinsic apoptotic signaling pathway (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_NegIntrinsicApoptosis_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY"]], ranks.cons[["6"]]) + labs(title="Negative regulation of intrinsic apoptotic signaling pathway (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_NegIntrinsicApoptosis_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY"]], ranks.filt[["6"]]) + labs(title="Negative regulation of intrinsic apoptotic signaling pathway (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_NegIntrinsicApoptosis_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]], ranks.all[["6"]]) + labs(title="Adaptive immune response (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_AdaptiveImmune_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]], ranks.cons[["6"]]) + labs(title="Adaptive immune response (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_AdaptiveImmune_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]], ranks.filt[["6"]]) + labs(title="Adaptive immune response (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_AdaptiveImmune_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION"]], ranks.all[["6"]]) + labs(title="Positive regulation of leukocyte cell-cell adhesion (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PosLeukocyteCellCellAdhesion_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION"]], ranks.cons[["6"]]) + labs(title="Positive regulation of leukocyte cell-cell adhesion (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PosLeukocyteCellCellAdhesion_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION"]], ranks.filt[["6"]]) + labs(title="Positive regulation of leukocyte cell-cell adhesion (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PosLeukocyteCellCellAdhesion_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_ACTIVATION"]], ranks.all[["6"]]) + labs(title="Positive regulation of cell activation (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PosCellActivation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_ACTIVATION"]], ranks.cons[["6"]]) + labs(title="Positive regulation of cell activation (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PosCellActivation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_ACTIVATION"]], ranks.filt[["6"]]) + labs(title="Positive regulation of cell activation (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PosCellActivation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_ADHESION"]], ranks.all[["6"]]) + labs(title="Positive regulation of cell adhesion (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PosCellAdhesion_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_ADHESION"]], ranks.cons[["6"]]) + labs(title="Positive regulation of cell adhesion (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PosCellAdhesion_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CELL_ADHESION"]], ranks.filt[["6"]]) + labs(title="Positive regulation of cell adhesion (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PosCellAdhesion_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_CLASS_II_PROTEIN_COMPLEX"]], ranks.all[["6"]]) + labs(title="Peptide antigen assembly with MHC-II complex (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_MHCIIassembly_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_CLASS_II_PROTEIN_COMPLEX"]], ranks.cons[["6"]]) + labs(title="Peptide antigen assembly with MHC-II complex (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_MHCIIassembly_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_CLASS_II_PROTEIN_COMPLEX"]], ranks.filt[["6"]]) + labs(title="Peptide antigen assembly with MHC-II complex (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_MHCIIassembly_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE_BASED_ON_SOMATIC_RECOMBINATION_OF_IMMUNE_RECEPTORS_BUILT_FROM_IMMUNOGLOBULIN_SUPERFAMILY_DOMAINS"]], ranks.all[["6"]]) + labs(title="Adaptive immune response based on somatic recombination (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_AdaptiveImmuneSomaticRecomb_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE_BASED_ON_SOMATIC_RECOMBINATION_OF_IMMUNE_RECEPTORS_BUILT_FROM_IMMUNOGLOBULIN_SUPERFAMILY_DOMAINS"]], ranks.cons[["6"]]) + labs(title="Adaptive immune response based on somatic recombination (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_AdaptiveImmuneSomaticRecomb_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE_BASED_ON_SOMATIC_RECOMBINATION_OF_IMMUNE_RECEPTORS_BUILT_FROM_IMMUNOGLOBULIN_SUPERFAMILY_DOMAINS"]], ranks.filt[["6"]]) + labs(title="Adaptive immune response based on somatic recombination (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_AdaptiveImmuneSomaticRecomb_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_PROTEIN_COMPLEX"]], ranks.all[["6"]]) + labs(title="Peptide antigen assembly with MHC complex (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_MHCassembly_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_PROTEIN_COMPLEX"]], ranks.cons[["6"]]) + labs(title="Peptide antigen assembly with MHC complex (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_MHCassembly_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_PROTEIN_COMPLEX"]], ranks.filt[["6"]]) + labs(title="Peptide antigen assembly with MHC complex (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_MHCassembly_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"]], ranks.all[["6"]]) + labs(title="Antigen processing and presentation (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_AgProcessingPresentation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"]], ranks.cons[["6"]]) + labs(title="Antigen processing and presentation (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_AgProcessingPresentation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"]], ranks.filt[["6"]]) + labs(title="Antigen processing and presentation (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_AgProcessingPresentation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN"]], ranks.all[["6"]]) + labs(title="Antigen processing and presentation of peptide antigen (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PeptideAgPresentation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN"]], ranks.cons[["6"]]) + labs(title="Antigen processing and presentation of peptide antigen (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PeptideAgPresentation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN"]], ranks.filt[["6"]]) + labs(title="Antigen processing and presentation of peptide antigen (Cluster 6)")
ggsave("Results/Figures/02_04_fgsea_clust6_PeptideAgPresentation_filt.png", p, dpi = 300)





# cluster 7
p <- plotEnrichment(gobp[["GOBP_POST_TRANSLATIONAL_PROTEIN_TARGETING_TO_ENDOPLASMIC_RETICULUM_MEMBRANE"]], ranks.all[["7"]]) + labs(title="Post-translational protein targeting to ER membrane (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_PostTranslERmem_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POST_TRANSLATIONAL_PROTEIN_TARGETING_TO_ENDOPLASMIC_RETICULUM_MEMBRANE"]], ranks.cons[["7"]]) + labs(title="Post-translational protein targeting to ER membrane (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_PostTranslERmem_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POST_TRANSLATIONAL_PROTEIN_TARGETING_TO_ENDOPLASMIC_RETICULUM_MEMBRANE"]], ranks.filt[["7"]]) + labs(title="Post-translational protein targeting to ER membrane (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_PostTranslERmem_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POST_TRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE_TRANSLOCATION"]], ranks.all[["7"]]) + labs(title="Post-translational protein targeting to membrane translocation (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_PostTranslMemTranslocation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POST_TRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE_TRANSLOCATION"]], ranks.cons[["7"]]) + labs(title="Post-translational protein targeting to membrane translocation (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_PostTranslMemTranslocation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POST_TRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE_TRANSLOCATION"]], ranks.filt[["7"]]) + labs(title="Post-translational protein targeting to membrane translocation (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_PostTranslMemTranslocation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_CYTOPLASMIC_TRANSLATION"]], ranks.all[["7"]]) + labs(title="Cytoplasmic translation (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_CytoplasmicTranslation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CYTOPLASMIC_TRANSLATION"]], ranks.cons[["7"]]) + labs(title="Cytoplasmic translation (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_CytoplasmicTranslation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CYTOPLASMIC_TRANSLATION"]], ranks.filt[["7"]]) + labs(title="Cytoplasmic translation (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_CytoplasmicTranslation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_STEM_CELL_DIVISION"]], ranks.all[["7"]]) + labs(title="Stem cell division (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_StemCellDivision_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_STEM_CELL_DIVISION"]], ranks.cons[["7"]]) + labs(title="Stem cell division (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_StemCellDivision_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_STEM_CELL_DIVISION"]], ranks.filt[["7"]]) + labs(title="Stem cell division (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_StemCellDivision_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY"]], ranks.all[["7"]]) + labs(title="Negative regulation of apoptotic signaling pathway (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_NegApoptoticSignaling_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY"]], ranks.cons[["7"]]) + labs(title="Negative regulation of apoptotic signaling pathway (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_NegApoptoticSignaling_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY"]], ranks.filt[["7"]]) + labs(title="Negative regulation of apoptotic signaling pathway (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_NegApoptoticSignaling_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_REGULATION_OF_ENDOCYTOSIS"]], ranks.all[["7"]]) + labs(title="Regulation of endocytosis (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_Endocytosis_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_REGULATION_OF_ENDOCYTOSIS"]], ranks.cons[["7"]]) + labs(title="Regulation of endocytosis (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_Endocytosis_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_REGULATION_OF_ENDOCYTOSIS"]], ranks.filt[["7"]]) + labs(title="Regulation of endocytosis (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_Endocytosis_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION"]], ranks.all[["7"]]) + labs(title="Positive regulation of cytokine production (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_PosCytokineProduction_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION"]], ranks.cons[["7"]]) + labs(title="Positive regulation of cytokine production (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_PosCytokineProduction_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION"]], ranks.filt[["7"]]) + labs(title="Positive regulation of cytokine production (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_PosCytokineProduction_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE"]], ranks.all[["7"]]) + labs(title="T cell activation involved in immune response (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_TcellActivationImmune_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE"]], ranks.cons[["7"]]) + labs(title="T cell activation involved in immune response (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_TcellActivationImmune_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE"]], ranks.filt[["7"]]) + labs(title="T cell activation involved in immune response (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_TcellActivationImmune_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_T_CELL_ACTIVATION"]], ranks.all[["7"]]) + labs(title="T cell activation (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_TcellActivation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_T_CELL_ACTIVATION"]], ranks.cons[["7"]]) + labs(title="T cell activation (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_TcellActivation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_T_CELL_ACTIVATION"]], ranks.filt[["7"]]) + labs(title="T cell activation (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_TcellActivation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE"]], ranks.all[["7"]]) + labs(title="Cell activation involved in immune response (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_CellActivationImmune_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE"]], ranks.cons[["7"]]) + labs(title="Cell activation involved in immune response (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_CellActivationImmune_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE"]], ranks.filt[["7"]]) + labs(title="Cell activation involved in immune response (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_CellActivationImmune_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY"]], ranks.all[["7"]]) + labs(title="Lymphocyte-mediated immunity (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_LymphocyteImmunity_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY"]], ranks.cons[["7"]]) + labs(title="Lymphocyte-mediated immunity (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_LymphocyteImmunity_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY"]], ranks.filt[["7"]]) + labs(title="Lymphocyte-mediated immunity (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_LymphocyteImmunity_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION"]], ranks.all[["7"]]) + labs(title="Positive regulation of leukocyte cell-cell adhesion (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_PosLeukocyteCellCellAdhesion_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION"]], ranks.cons[["7"]]) + labs(title="Positive regulation of leukocyte cell-cell adhesion (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_PosLeukocyteCellCellAdhesion_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION"]], ranks.filt[["7"]]) + labs(title="Positive regulation of leukocyte cell-cell adhesion (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_PosLeukocyteCellCellAdhesion_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_IMMUNE_EFFECTOR_PROCESS"]], ranks.all[["7"]]) + labs(title="Immune effector process (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_ImmuneEffector_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_IMMUNE_EFFECTOR_PROCESS"]], ranks.cons[["7"]]) + labs(title="Immune effector process (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_ImmuneEffector_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_IMMUNE_EFFECTOR_PROCESS"]], ranks.filt[["7"]]) + labs(title="Immune effector process (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_ImmuneEffector_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]], ranks.all[["7"]]) + labs(title="Adaptive immune response (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_AdaptiveImmune_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]], ranks.cons[["7"]]) + labs(title="Adaptive immune response (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_AdaptiveImmune_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]], ranks.filt[["7"]]) + labs(title="Adaptive immune response (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_AdaptiveImmune_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE"]], ranks.all[["7"]]) + labs(title="Lymphocyte activation involved in immune response (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_LymphocyteActivationImmune_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE"]], ranks.cons[["7"]]) + labs(title="Lymphocyte activation involved in immune response (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_LymphocyteActivationImmune_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE"]], ranks.filt[["7"]]) + labs(title="Lymphocyte activation involved in immune response (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_LymphocyteActivationImmune_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE_BASED_ON_SOMATIC_RECOMBINATION_OF_IMMUNE_RECEPTORS_BUILT_FROM_IMMUNOGLOBULIN_SUPERFAMILY_DOMAINS"]], ranks.all[["7"]]) + labs(title="Adaptive immune response based on somatic recombination (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_AdaptiveImmuneSomaticRecomb_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE_BASED_ON_SOMATIC_RECOMBINATION_OF_IMMUNE_RECEPTORS_BUILT_FROM_IMMUNOGLOBULIN_SUPERFAMILY_DOMAINS"]], ranks.cons[["7"]]) + labs(title="Adaptive immune response based on somatic recombination (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_AdaptiveImmuneSomaticRecomb_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ADAPTIVE_IMMUNE_RESPONSE_BASED_ON_SOMATIC_RECOMBINATION_OF_IMMUNE_RECEPTORS_BUILT_FROM_IMMUNOGLOBULIN_SUPERFAMILY_DOMAINS"]], ranks.filt[["7"]]) + labs(title="Adaptive immune response based on somatic recombination (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_AdaptiveImmuneSomaticRecomb_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_CYTOKINE_PRODUCTION"]], ranks.all[["7"]]) + labs(title="Cytokine production (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_CytokineProduction_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CYTOKINE_PRODUCTION"]], ranks.cons[["7"]]) + labs(title="Cytokine production (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_CytokineProduction_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CYTOKINE_PRODUCTION"]], ranks.filt[["7"]]) + labs(title="Cytokine production (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_CytokineProduction_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE"]], ranks.all[["7"]]) + labs(title="Positive regulation of immune response (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_PosImmuneResponse_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE"]], ranks.cons[["7"]]) + labs(title="Positive regulation of immune response (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_PosImmuneResponse_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE"]], ranks.filt[["7"]]) + labs(title="Positive regulation of immune response (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_PosImmuneResponse_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II"]], ranks.all[["7"]]) + labs(title="Antigen processing and presentation via MHC-II (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_MHCII_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II"]], ranks.cons[["7"]]) + labs(title="Antigen processing and presentation via MHC-II (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_MHCII_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II"]], ranks.filt[["7"]]) + labs(title="Antigen processing and presentation via MHC-II (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_MHCII_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"]], ranks.all[["7"]]) + labs(title="Antigen processing and presentation (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_AgProcessingPresentation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"]], ranks.cons[["7"]]) + labs(title="Antigen processing and presentation (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_AgProcessingPresentation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"]], ranks.filt[["7"]]) + labs(title="Antigen processing and presentation (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_AgProcessingPresentation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_CELL_CYCLE"]], ranks.all[["7"]]) + labs(title="Cell cycle (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_CellCycle_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_CYCLE"]], ranks.cons[["7"]]) + labs(title="Cell cycle (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_CellCycle_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_CYCLE"]], ranks.filt[["7"]]) + labs(title="Cell cycle (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_CellCycle_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_MITOTIC_SPINDLE"]], ranks.all[["7"]]) + labs(title="Mitotic spindle (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_MitoticSpindle_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_MITOTIC_SPINDLE"]], ranks.cons[["7"]]) + labs(title="Mitotic spindle (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_MitoticSpindle_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_MITOTIC_SPINDLE"]], ranks.filt[["7"]]) + labs(title="Mitotic spindle (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_MitoticSpindle_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.all[["7"]]) + labs(title="E2F targets (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_E2F_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.cons[["7"]]) + labs(title="E2F targets (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_E2F_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_E2F_TARGETS"]], ranks.filt[["7"]]) + labs(title="E2F targets (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_E2F_filt.png", p, dpi = 300)

p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.all[["7"]]) + labs(title="G2M checkpoint (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_G2M_all.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.cons[["7"]]) + labs(title="G2M checkpoint (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_G2M_cons.png", p, dpi = 300)
p <- plotEnrichment(hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks.filt[["7"]]) + labs(title="G2M checkpoint (Cluster 7)")
ggsave("Results/Figures/02_04_fgsea_clust7_G2M_filt.png", p, dpi = 300)



# cluster 10

p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_GENE_EXPRESSION"]], ranks.all[["10"]]) + labs(title="Positive regulation of gene expression (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_PosGeneExpression_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_GENE_EXPRESSION"]], ranks.cons[["10"]]) + labs(title="Positive regulation of gene expression (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_PosGeneExpression_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_POSITIVE_REGULATION_OF_GENE_EXPRESSION"]], ranks.filt[["10"]]) + labs(title="Positive regulation of gene expression (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_PosGeneExpression_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY"]], ranks.all[["10"]]) + labs(title="Myeloid leukocyte-mediated immunity (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_MyeloidImmunity_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY"]], ranks.cons[["10"]]) + labs(title="Myeloid leukocyte-mediated immunity (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_MyeloidImmunity_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY"]], ranks.filt[["10"]]) + labs(title="Myeloid leukocyte-mediated immunity (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_MyeloidImmunity_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_LEUKOCYTE_DEGRANULATION"]], ranks.all[["10"]]) + labs(title="Leukocyte degranulation (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_Degranulation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_LEUKOCYTE_DEGRANULATION"]], ranks.cons[["10"]]) + labs(title="Leukocyte degranulation (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_Degranulation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_LEUKOCYTE_DEGRANULATION"]], ranks.filt[["10"]]) + labs(title="Leukocyte degranulation (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_Degranulation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY"]], ranks.all[["10"]]) + labs(title="Negative regulation of intrinsic apoptotic signaling pathway (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_NegIntrinsicApoptosis_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY"]], ranks.cons[["10"]]) + labs(title="Negative regulation of intrinsic apoptotic signaling pathway (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_NegIntrinsicApoptosis_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY"]], ranks.filt[["10"]]) + labs(title="Negative regulation of intrinsic apoptotic signaling pathway (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_NegIntrinsicApoptosis_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_TRANSLATION"]], ranks.all[["10"]]) + labs(title="Translation (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_Translation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_TRANSLATION"]], ranks.cons[["10"]]) + labs(title="Translation (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_Translation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_TRANSLATION"]], ranks.filt[["10"]]) + labs(title="Translation (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_Translation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II"]], ranks.all[["10"]]) + labs(title="Antigen processing and presentation via MHC-II (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_MHCII_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II"]], ranks.cons[["10"]]) + labs(title="Antigen processing and presentation via MHC-II (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_MHCII_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II"]], ranks.filt[["10"]]) + labs(title="Antigen processing and presentation via MHC-II (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_MHCII_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY"]], ranks.all[["10"]]) + labs(title="Negative regulation of apoptotic signaling pathway (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_NegApoptoticSignaling_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY"]], ranks.cons[["10"]]) + labs(title="Negative regulation of apoptotic signaling pathway (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_NegApoptoticSignaling_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY"]], ranks.filt[["10"]]) + labs(title="Negative regulation of apoptotic signaling pathway (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_NegApoptoticSignaling_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"]], ranks.all[["10"]]) + labs(title="Antigen processing and presentation (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_AgProcessingPresentation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"]], ranks.cons[["10"]]) + labs(title="Antigen processing and presentation (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_AgProcessingPresentation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"]], ranks.filt[["10"]]) + labs(title="Antigen processing and presentation (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_AgProcessingPresentation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_CYTOPLASMIC_TRANSLATION"]], ranks.all[["10"]]) + labs(title="Cytoplasmic translation (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_CytoplasmicTranslation_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CYTOPLASMIC_TRANSLATION"]], ranks.cons[["10"]]) + labs(title="Cytoplasmic translation (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_CytoplasmicTranslation_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CYTOPLASMIC_TRANSLATION"]], ranks.filt[["10"]]) + labs(title="Cytoplasmic translation (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_CytoplasmicTranslation_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_SYNCYTIUM_FORMATION_BY_PLASMA_MEMBRANE_FUSION"]], ranks.all[["10"]]) + labs(title="Negative regulation of syncytium formation (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_NegSyncytium_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_SYNCYTIUM_FORMATION_BY_PLASMA_MEMBRANE_FUSION"]], ranks.cons[["10"]]) + labs(title="Negative regulation of syncytium formation (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_NegSyncytium_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_NEGATIVE_REGULATION_OF_SYNCYTIUM_FORMATION_BY_PLASMA_MEMBRANE_FUSION"]], ranks.filt[["10"]]) + labs(title="Negative regulation of syncytium formation (Cluster 10)")
ggsave("Results/Figures/02_04_fgsea_clust10_NegSyncytium_filt.png", p, dpi = 300)


# cluster 11

p <- plotEnrichment(gobp[["GOBP_PYROPTOTIC_INFLAMMATORY_RESPONSE"]], ranks.all[["11"]]) + labs(title="Pyroptotic inflammatory response (Cluster 11)")
ggsave("Results/Figures/02_04_fgsea_clust11_PyroptoticInflammatoryResponse_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_PYROPTOTIC_INFLAMMATORY_RESPONSE"]], ranks.cons[["11"]]) + labs(title="Pyroptotic inflammatory response (Cluster 11)")
ggsave("Results/Figures/02_04_fgsea_clust11_PyroptoticInflammatoryResponse_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_PYROPTOTIC_INFLAMMATORY_RESPONSE"]], ranks.filt[["11"]]) + labs(title="Pyroptotic inflammatory response (Cluster 11)")
ggsave("Results/Figures/02_04_fgsea_clust11_PyroptoticInflammatoryResponse_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_INTERLEUKIN_10_PRODUCTION"]], ranks.all[["11"]]) + labs(title="Interleukin-10 production (Cluster 11)")
ggsave("Results/Figures/02_04_fgsea_clust11_IL10Production_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_INTERLEUKIN_10_PRODUCTION"]], ranks.cons[["11"]]) + labs(title="Interleukin-10 production (Cluster 11)")
ggsave("Results/Figures/02_04_fgsea_clust11_IL10Production_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_INTERLEUKIN_10_PRODUCTION"]], ranks.filt[["11"]]) + labs(title="Interleukin-10 production (Cluster 11)")
ggsave("Results/Figures/02_04_fgsea_clust11_IL10Production_filt.png", p, dpi = 300)

p <- plotEnrichment(gobp[["GOBP_CELL_KILLING"]], ranks.all[["11"]]) + labs(title="Cell killing (Cluster 11)")
ggsave("Results/Figures/02_04_fgsea_clust11_CellKilling_all.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_KILLING"]], ranks.cons[["11"]]) + labs(title="Cell killing (Cluster 11)")
ggsave("Results/Figures/02_04_fgsea_clust11_CellKilling_cons.png", p, dpi = 300)
p <- plotEnrichment(gobp[["GOBP_CELL_KILLING"]], ranks.filt[["11"]]) + labs(title="Cell killing (Cluster 11)")
ggsave("Results/Figures/02_04_fgsea_clust11_CellKilling_filt.png", p, dpi = 300)






