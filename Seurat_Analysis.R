---
#title: "Seurats analysis for spleen data"
#note: "To run this analysis for the bone marrow data, import the bone marrow 10X files at code line 15-17
---

rm(list=ls())

# Load neccessary packages

library(dplyr)
library(Seurat)
library(patchwork)



# Loading the datasets

Run1_V7_SS.data <- Read10X(data.dir = "../../Malaria_data/Malaria_Mouse_v7/Run1_V7_SS")
Run2_V7_SS2.data <- Read10X(data.dir = "../../Malaria_data/Malaria_Mouse_v7/Run2_V7_SS2")
Run2_V7_SS3.data <- Read10X(data.dir = "../../Malaria_data/Malaria_Mouse_v7/Run2_V7_SS3")

SS1 <- CreateSeuratObject(counts = Run1_V7_SS.data, project = "Malaria_data", min.cells = 3, min.features = 200)
SS2 <- CreateSeuratObject(counts = Run2_V7_SS2.data, project = "Malaria_data", min.cells = 3, min.features = 200)
SS3 <- CreateSeuratObject(counts = Run2_V7_SS3.data, project = "Malaria_data", min.cells = 3, min.features = 200)

#Calculating the percentage mt genes

SS1[["MousePercent.mt"]] <- PercentageFeatureSet(SS1, pattern = "mt-")
SS1[["PbPercent.mt"]] <- PercentageFeatureSet(SS1, pattern = "*PBANKA-MIT")
SS2[["MousePercent.mt"]] <- PercentageFeatureSet(SS2, pattern = "mt-")
SS2[["PbPercent.mt"]] <- PercentageFeatureSet(SS2, pattern = "*PBANKA-MIT")
SS3[["MousePercent.mt"]] <- PercentageFeatureSet(SS3, pattern = "mt-")
SS3[["PbPercent.mt"]] <- PercentageFeatureSet(SS3, pattern = "*PBANKA-MIT")

# Adding the variable that identifies the condition

SS1[["Condition"]] <- "Mouse-1"
SS2[["Condition"]] <- "Mouse-2"
SS3[["Condition"]] <- "Mouse-3"

SS1[["Pb_Status"]] <- "Infected"
SS2[["Pb_Status"]] <- "Infected"
SS3[["Pb_Status"]] <- "Control"

SS.list <- list(SS1, SS2, SS3)

FeatureScatter(SS.list[[1]], feature1 = "nCount_RNA", feature2 = "MousePercent.mt")
FeatureScatter(SS.list[[1]], feature1 = "nCount_RNA", feature2 = "PbPercent.mt")
FeatureScatter(SS.list[[1]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

FeatureScatter(SS.list[[2]], feature1 = "nCount_RNA", feature2 = "MousePercent.mt")
FeatureScatter(SS.list[[2]], feature1 = "nCount_RNA", feature2 = "PbPercent.mt")
FeatureScatter(SS.list[[2]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

FeatureScatter(SS.list[[3]], feature1 = "nCount_RNA", feature2 = "MousePercent.mt")
FeatureScatter(SS.list[[3]], feature1 = "nCount_RNA", feature2 = "PbPercent.mt")
FeatureScatter(SS.list[[3]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(SS.list[[1]], features = c("nFeature_RNA","nCount_RNA", "PbPercent.mt","MousePercent.mt"), ncol = 4)
VlnPlot(SS.list[[2]], features = c("nFeature_RNA","nCount_RNA", "PbPercent.mt","MousePercent.mt"), ncol = 4)
VlnPlot(SS.list[[3]], features = c("nFeature_RNA","nCount_RNA", "PbPercent.mt","MousePercent.mt"), ncol = 4)


# Filtering features

SS1 <- subset(SS.list[[1]], subset = nFeature_RNA > 150 & nFeature_RNA < 5000 & PbPercent.mt < 5 & MousePercent.mt < 5)
SS2 <- subset(SS.list[[2]], subset = nFeature_RNA > 150 & nFeature_RNA < 5000 & PbPercent.mt < 5 & MousePercent.mt < 5)
SS3 <- subset(SS.list[[3]], subset = nFeature_RNA > 150 & nFeature_RNA < 6000 & PbPercent.mt < 5 & MousePercent.mt < 5)

SS1
SS2
SS3
```

# Normalize and identify variable features for each dataset independently

SS.list <- lapply(X = list(SS1, SS2, SS3), FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})
```

# Identify the 10 most highly variable genes in the datasets


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SS.list[[1]]), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(SS.list[[1]])
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
plot2

### SS2

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SS.list[[2]]), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(SS.list[[2]])
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

### SS3

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SS.list[[3]]), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(SS.list[[3]])
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2



# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = SS.list)


# Perform integration

SS.anchors <- FindIntegrationAnchors(object.list = SS.list, anchor.features = features)
SS_anchors.combined <- IntegrateData(anchorset = SS.anchors) # this command creates an 'integrated' data assay  of all common features


# Perform an integrated analysis

DefaultAssay(SS_anchors.combined) <- "integrated" # specify that we will perform downstream analysis on the corrected data note that the

# Run the standard workflow for visualization and clustering
SS_anchors.combined <- ScaleData(SS_anchors.combined, verbose = FALSE)
SS_anchors.combined <- RunPCA(SS_anchors.combined, npcs = 50 , verbose = FALSE)
SS_anchors.combined <- RunUMAP(SS_anchors.combined, reduction = "pca", dims = 1:32)
SS_anchors.combined <- FindNeighbors(SS_anchors.combined, reduction = "pca", dims = 1:32)
SS_anchors.combined <- FindClusters(SS_anchors.combined, resolution = 0.18)

ElbowPlot(SS_anchors.combined, ndims = 50)


# Visualization
DimPlot(SS_anchors.combined, reduction = "umap", group.by = "Condition")
DimPlot(SS_anchors.combined, reduction = "umap", group.by = "Pb_Status")
DimPlot(SS_anchors.combined, reduction = "umap", label = TRUE, repel = TRUE)

DimPlot(SS_anchors.combined, reduction = "umap", split.by = "Condition")

DimPlot(SS_anchors.combined, reduction = "umap", split.by = "Condition",label = TRUE, repel = TRUE)

DimPlot(SS_anchors.combined, reduction = "umap", split.by = "Pb_Status",label = TRUE, repel = TRUE)



# Multiple cluster exploration

# This loop just runs the FindMarkers function on all of the clusters
lapply(
  levels(SS_anchors.combined[["seurat_clusters"]][[1]]),
  function(x)FindConservedMarkers(SS_anchors.combined,ident.1 = x,min.pct = 0.25, grouping.var = "Pb_Status")
) -> SS_cluster.markers

# This simply adds the cluster number to the results of FindMarkers
sapply(0:(length(SS_cluster.markers)-1),function(x) {
  SS_cluster.markers[[x+1]]$gene <<- rownames(SS_cluster.markers[[x+1]])
  SS_cluster.markers[[x+1]]$cluster <<- x
})

# Finally we collapse the list of hits down to a single table and sort it by FDR to put the most significant ones first
rbindlist(SS_cluster.markers,fill=TRUE) -> test
test %>% arrange(cluster,"Mouse-2_p_val_adj") -> SS_cluster.markers.test


# Labelling UMAP clusters
SS_anchors.combined <- RenameIdents(SS_anchors.combined, 
                                    `0` = "Pb_gene", `1` = "Erynthroblasts",`2` ="Erynthroblasts", `3` = "Monocytes", 
                                    `4` = "B-Cells", `5` = "Macrophages", `6` = "Proerynthroblasts", `7` = "T-Cells",
                                    `6` = "HSPCs", `7` = "DCs", `8` = "Late Pro-Bcells")

### Subseting to ignore clusters 9 and 10. These were small clusters and not easily identifiable
sub_obj <- subset(x = SS_anchors.combined, idents = c("Pb_gene", "Erynthroblasts", "Erynthroblasts", "Monocytes", "B-Cells", "Macrophages", "Proerynthroblasts", "T-Cells", "Late Pro-Bcells"))

DimPlot(sub_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

SS_anchors.combined <- RenameIdents(SS_anchors.combined, 
                                    `0` = "Pb_gene", `1` = "Erynthroblasts",`2` ="Erynthroblasts", `3` = "Monocytes", 
                                    `4` = "B-Cells", `5` = "Macrophages", `6` = "Proerynthroblasts", `7` = "T-Cells",
                                    `6` = "HSPCs", `7` = "DCs", `8` = "Late Pro-Bcells")

### Subseting to ignore clusters 9 and 10
sub_obj <- subset(x = SS_anchors.combined, idents = c("Pb_gene", "Erynthroblasts", "Erynthroblasts", "Monocytes", "B-Cells", "Macrophages", "Proerynthroblasts", "T-Cells", "Late Pro-Bcells"))

DimPlot(sub_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(sub_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


## Integrated data
DimPlot(sub_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Splits
DimPlot(sub_obj, reduction = "umap", split.by = "Condition",label = TRUE, repel = TRUE, pt.size = 0.5)+ NoLegend()

DimPlot(sub_obj, reduction = "umap", split.by = "Pb_Status",label = TRUE, repel = TRUE, pt.size = 0.5)+ NoLegend()


# Saving the final seurats object
saveRDS(SS_anchors.combined, file = "Seurats_Analysis_SS.rds") 
# readRDS("Seurats_Analysis_SS.rds") # to read the object in for later









