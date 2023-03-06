
  
  # Integrating Bone marrow and Spleen Seurats

  
Seurats_Analysis_SS <- readRDS("~/Documents/PhD University of Glasgow/First Rotation/Data Analysis tutorial/Seurat analysis/Seurats_Analysis_SS.rds")
Seurats_Analysis_BM <- readRDS("~/Documents/PhD University of Glasgow/First Rotation/Data Analysis tutorial/Seurat analysis/Seurats_Analysis_BM.rds")

Seurats_Analysis_SS[["Tissue"]] <- "Spleen"
Seurats_Analysis_BM[["Tissue"]] <- "Bone_Marrow"

#Integrating the datasets

#Mix.list.b <- list(Seurats_Analysis_SS, Seurats_Analysis_BM)
# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list(Seurats_Analysis_SS, Seurats_Analysis_BM))
# Perform integration
Mix.anchors.b <- FindIntegrationAnchors(object.list = list(Seurats_Analysis_SS, Seurats_Analysis_BM), anchor.features = features)
# this command creates an 'integrated' data assay  of all common features
Mix_anchors.combined.b <- IntegrateData(anchorset = Mix.anchors.b)

DefaultAssay(Mix_anchors.combined.b) <- "integrated"

# Run the standard workflow for visualization and clustering
Mix_anchors.combined.b <- ScaleData(Mix_anchors.combined.b, verbose = FALSE)
Mix_anchors.combined.b <- RunPCA(Mix_anchors.combined.b, npcs = 50 , verbose = FALSE)
Mix_anchors.combined.b <- RunUMAP(Mix_anchors.combined.b, reduction = "pca", dims = 1:32)
Mix_anchors.combined.b <- FindNeighbors(Mix_anchors.combined.b, reduction = "pca", dims = 1:32)
Mix_anchors.combined.b <- FindClusters(Mix_anchors.combined.b, resolution = 0.15)

ElbowPlot(Mix_anchors.combined.b, ndims = 50)


# Visualization
Mix_an = subset(x = Mix_anchors.combined.b, 
                idents = c(0:8))

DimPlot(Mix_an, reduction = "umap", pt.size = 0.1, group.by = "Tissue")
DimPlot(Mix_an, reduction = "umap", pt.size = 0.1, split.by ="Tissue")
DimPlot(Mix_an, reduction = "umap", label = TRUE, pt.size = 0.1, repel = TRUE)
DimPlot(Mix_an, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()


Mix_an <- RenameIdents(Mix_an, 
                       `0` = "Pb_genes", `1` = "Erynthroblasts",`2` = "Monocytes", `3` = "B-Cells", 
                       `4` = "Macrophages", `5` = "T cells", `6` = "Proerynthroblasts", `7` = "Dcs",
                       `8` = "Late Pro-Bcells")

DimPlot(Mix_an, reduction = "umap", pt.size = 0.1, group.by = "Tissue")
DimPlot(Mix_an, reduction = "umap", label = TRUE, pt.size = 0.1, group.by = "Tissue")
DimPlot(Mix_an, reduction = "umap", pt.size = 0.1, split.by ="Tissue")
DimPlot(Mix_an, reduction = "umap", label = TRUE, pt.size = 0.1, repel = TRUE)
DimPlot(Mix_an, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()

#subset(x=Mix_an, subset = Pb_Status == "Infected")

DimPlot(subset(x=Mix_an, subset = Pb_Status == "Infected"), reduction = "umap", pt.size = 0.1, group.by = "Tissue")
DimPlot(subset(x=Mix_an, subset = Pb_Status == "Infected"), reduction = "umap", label = TRUE, pt.size = 0.1, group.by = "Tissue")
DimPlot(subset(x=Mix_an, subset = Pb_Status == "Infected"), reduction = "umap", pt.size = 0.1, split.by ="Tissue")
DimPlot(subset(x=Mix_an, subset = Pb_Status == "Infected"), reduction = "umap", label = TRUE, pt.size = 0.1, repel = TRUE)
DimPlot(subset(x=Mix_an, subset = Pb_Status == "Infected"), reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()



#Plot cell counts per tissue

counts.bm <- as.data.frame(table(Idents(subset(x = Seurats_Analysis_BM, 
                                               idents = c("Monocytes", "B-Cells","Macrophages", "T-Cells"))), 
                                 subset(x = Seurats_Analysis_BM, 
                                        idents = c("Monocytes", "B-Cells","Macrophages", "T-Cells"))$Pb_Status))
ggplot(counts.bm, aes(x=Var1, y=Freq, fill=Var2))+
  geom_bar(stat="identity", position=position_dodge())

counts.ss <- as.data.frame(table(Idents(subset(x = Seurats_Analysis_SS, 
                                               idents = c("Monocytes", "B-Cells","Macrophages", "T-Cells"))), 
                                 subset(x = Seurats_Analysis_SS, 
                                        idents = c("Monocytes", "B-Cells","Macrophages", "T-Cells"))$Pb_Status))
ggplot(counts.ss, aes(x=Var1, y=Freq, fill=Var2))+
  geom_bar(stat="identity", position=position_dodge())















