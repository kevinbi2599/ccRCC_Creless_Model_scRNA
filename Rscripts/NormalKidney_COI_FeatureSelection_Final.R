library(Seurat)

setwd("/broad/VanAllenLab/xchip/cga_home/kevinbi/StranskyLinehan_Mouse")

seur_ni_n = readRDS("Manuscript_Objects/NonImmune_NormalKidney_AmbFilt_SeuratObject_Final.RDS")
seur_ni_n = subset(seur_ni_n, idents = c("Immune_Contaminant","Unidentified","Unidentified_PT"), invert = T)
table(Idents(seur_ni_n))

markers = FindAllMarkers(seur_ni_n, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA",test.use = "LR",latent.vars = "Sex")

saveRDS(markers,file = "Manuscript_Objects/NormalKidney_COI_FeatureSelection_Output_Final.RDS")