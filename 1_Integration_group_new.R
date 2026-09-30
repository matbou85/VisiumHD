
#----[ load packages ]------#
{
  library("tibble")
  library("cowplot")
  library("Seurat")
  library("ggplot2")
  library("patchwork")
  library("dplyr")
  # library("SeuratData")
  library("tidyverse")
  library("SeuratWrappers")
  library("Banksy")
  options(future.globals.maxSize = 40000 * 1024^2)
}


#----[ load samples ]------#
{
  TNF <- readRDS(file= "/work_beegfs/sukmb687/FT_DPP_UC_VisiumHD_run2_2758/Analysis_TNF/UC_biopsy_TNF_merge_untransformed_Region.RDS")
  # TNF <- readRDS(file= "UC_biopsy_TNF_merge_untransformed_Region.RDS")
  TNF@meta.data$Condition <- TNF@meta.data$orig.ident
  TNF@meta.data$Condition <- gsub(".*TNF[1-6]_Infx_", "TNF_", TNF@meta.data$Condition)
  TNF@meta.data$Condition <- gsub(".*TNF[1-6]_Ada_", "TNF_", TNF@meta.data$Condition)
  names(TNF@images) <- c("TNF1", "TNF2", "TNF4", "TNF5", "TNF6")

  VEDO <- readRDS(file= "/work_beegfs/sukmb687/FT_DPP_UC_VisiumHD_run2_2758/Analysis_VEDO/UC_biopsy_VEDO_merge_untransformed_Region.RDS")
  # VEDO <- readRDS(file= "UC_biopsy_VEDO_merge_untransformed_Region.RDS")
  VEDO@meta.data$Condition <- VEDO@meta.data$orig.ident
  VEDO@meta.data$Condition <- gsub(".*VED[1-6]_", "VEDO_", VEDO@meta.data$Condition)
  names(VEDO@images) <- c("VED1", "VED2", "VED3", "VED4", "VED5", "VED6")

  Jaki <- readRDS(file= "/work_beegfs/sukmb687/FT_DPP_UC_VisiumHD_run2_2758/Analysis_Jaki/UC_biopsy_JAKi_merge_untransformed_Region_Cluster.RDS")
  # Jaki <- readRDS(file= "UC_biopsy_JAKi_merge_untransformed_Region.RDS")
  Jaki@meta.data$Condition <- Jaki@meta.data$orig.ident
  # Jaki@meta.data$Condition <- paste0("Jaki_", Jaki@meta.data$Condition)
  Jaki@meta.data$Condition <- gsub(".*JAK[1-6]_", "JAK_", Jaki@meta.data$Condition, ignore.case = T)
  Jaki@meta.data$Condition <- gsub(".*JAK_Upa_", "JAK_", Jaki@meta.data$Condition, ignore.case = T)
  Jaki@meta.data$Condition <- gsub(".*JAK_Tofa_", "JAK_", Jaki@meta.data$Condition, ignore.case = T)
  Jaki@meta.data$Condition <- gsub(".*JAK_Filgo_", "JAK_", Jaki@meta.data$Condition, ignore.case = T)
  names(Jaki@images) <- c("Jak1", "Jak2", "Jak3", "Jak4", "Jak5")

}


#----[ Merge samples ]------#
{
  obj <- merge(x = TNF, y = list(VEDO, Jaki))
  saveRDS(file= "UC_biopsy_merge_untransformed_Region_new.RDS", obj)

  obj <- readRDS(file= "UC_biopsy_merge_untransformed_Region_new.RDS")

}


#----[ Merge samples and re-define  ]------#
{
  DefaultAssay(obj) <- "RNA" ## Define layers to work with
  obj@assays[["RNA"]]

  obj[["RNA"]] <- JoinLayers(obj[["RNA"]]) ## Join all RNA layers
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Condition) ## Split object by Condition Jaki_R_W0, TNF_NR_W14...

  VlnPlot(obj, features = "nFeature_RNA", ncol = 1, split.by = "Condition") ## PLot QC
  ggsave("Rplot_violin_QC.tiff", units="in", width=20, height=21, dpi=300)
}


#----[ Standard workflow  ]------#
{
  obj <- SCTransform(obj, vst.flavor = "v2") ## Normalize
  DefaultAssay(obj) <- "SCT" ## Switch to Normalized counts layer

  VariableFeatures(obj[["SCT"]]) <- rownames(obj[["SCT"]]@scale.data) ## Define variable features for further analyses

  obj <- RunPCA(obj, verbose = TRUE)
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj, verbose = TRUE)
  obj <- RunUMAP(obj, dims = 1:30)
  DimPlot(obj, reduction = "umap", group.by = c("ident", "orig.ident"))
  ggsave("Rplot_obj_UMAP_new.tiff", units="in", width=10, height=21, dpi=300)

  saveRDS(file= "/work_beegfs/sukmb687/FT_DPP_UC_VisiumHD_run2_2758/Analysis_group/UC_biopsy_merge_clean_new.RDS", obj)
}


#----[ Remove batch effect using harmony reduction  ]------#
{
  obj <- readRDS(file= "UC_biopsy_merge_clean_new.RDS")

  obj@active.assay = "SCT"

  obj <- IntegrateLayers(object = obj,
                         method = HarmonyIntegration, assay = "SCT",
                         orig.reduction = "pca", new.reduction = "harmony",
                         batch = obj$cell.ids,
                         verbose = TRUE
  )

  saveRDS(file= "/work_beegfs/sukmb687/FT_DPP_UC_VisiumHD_run2_2758/Analysis_group/UC_biopsy_merge_clean_integrated_new.RDS", obj)

}


#----[ Find clusters for several resolutions  ]------#
{
  obj <- readRDS(file="UC_biopsy_merge_clean_integrated_new.RDS")

  obj <- RunPCA(obj, verbose = TRUE)
  obj <- FindNeighbors(obj, dims = 1:30)

  res = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

  for(i in 1:11){
    obj  <- FindClusters(obj, verbose = TRUE, resolution = res[i])
  }

  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
  DimPlot(obj, reduction = "umap.harmony", group.by = c("ident", "Condition")) + ggtitle("harmony_integrated") + theme_classic()
  ggsave("Rplot_obj_UMAP_harmony_new.tiff", units="in", width=20, height=21, dpi=300)

  saveRDS(file= "/work_beegfs/sukmb687/FT_DPP_UC_VisiumHD_run2_2758/Analysis_group/UC_biopsy_merge_clean_integrated_harmony_new.RDS", obj)

}










