---
  title: "Cell-Cell Interactions, Bone marrow"
  # Note: To perform this same analysis for the spleen, import the spleen instead of the bonemarrow
          #seurats object at code line 26.
---

  
rm(list=ls())

# Loading the necessary libraries

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(dplyr)
library(Seurat)

# The working directory
data.dir <- './comparison'
dir.create(data.dir)
setwd(data.dir)


# Converting the Seurat object into a CellChat object

BM_anchors.combined <- readRDS("Seurats_Analysis_BM.rds")

data.input <- GetAssayData(BM_anchors.combined, assay = "RNA", slot = "data") # normalized data matri
labels <- Idents(BM_anchors.combined)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
names(meta)[names(meta) == "group"] <- "labels"
 ## Getting the extra variables for the meta data
extra.meta = BM_anchors.combined@meta.data
meta = dplyr::bind_cols(meta, extra.meta)
 ## Before creating the cell chat object
cell.use1 = rownames(meta)[meta$Pb_Status == "Control"] # extract the cell names from control data  
cell.use2 = rownames(meta)[meta$Pb_Status == "Infected"] # extract the cell names from disease data
 ## Prepare input data for CellChat analysis
data.input1 = data.input[, cell.use1]
meta1 = meta[cell.use1, ]

data.input2 = data.input[, cell.use2]
meta2 = meta[cell.use2, ]

## Creating the cell chat objects

# For the infected mice
cellchat1 <- createCellChat(object = data.input1)
cellchat1 <- addMeta(cellchat1, meta = meta1)
cellchat1 <- setIdent(cellchat1, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat1@idents) # show factor levels of the cell labels
cellchat1@idents <- droplevels(cellchat1@idents, exclude = c("Late Pro-Bcells"))
groupSize <- as.numeric(table(cellchat1@idents)) # number of cells in each cell group

# For the control mice.
cellchat2 <- createCellChat(object = data.input2)
cellchat2 <- addMeta(cellchat2, meta = meta2)
cellchat2 <- setIdent(cellchat2, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat2@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat2@idents)) # number of cells in each cell group


# Processing the cellchat objects for each condition

```{r}
for (i in c(cellchat1,cellchat2)) {
  
  cellchat = i
  
  # 1. Set the ligand-receptor interaction database
  
  CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
  showDatabaseCategory(CellChatDB)
  
  #Deleting the two unfound interactions from the  "subset" function in previous code and then runing the upper chunk again
  
  which(CellChatDB[["interaction"]]$ligand == "H2-BI") # 1887
  CellChatDB[["interaction"]] <- CellChatDB[["interaction"]][-1887,]
  which(CellChatDB[["interaction"]]$ligand == "H2-Ea-ps") #1900
  CellChatDB[["interaction"]] <- CellChatDB[["interaction"]][-1900,]
  
  # use all CellChatDB for cell-cell communication analysis
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  
  # set the used database in the object
  cellchat@DB <- CellChatDB.use
  
  # Editing the gene names to remove the extra string on them
  cellchat@data@Dimnames[[1]] <-gsub("Mouse93------","",as.character(cellchat@data@Dimnames[[1]]))
  
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multiprocess", workers = 4) # do parallel
  
  #Part I: Identify over expressed genes and interactions
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  #Part II: Inference of cell-cell communication network
  cellchat <- computeCommunProb(cellchat, population.size = TRUE)
  #cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10) # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  
  #Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  
  #Calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  
  i = cellchat
  
}


Merging the cellchat objects

cellchat.NL <- cellchat ##### from cellchat1 processing. CMD completed...... waiting for merging
cellchat.LS <- cellchat ##### from cellchat2 processing. CMD completed...... waiting for merging
cellchat.NL <- updateCellChat(cellchat.NL)
cellchat.LS <- updateCellChat(cellchat.LS)

# netAnalysis_computeCentrality

#Compute and visualize the network centrality scores
# Compute the network centrality scores
cellchat.NL <- netAnalysis_computeCentrality(cellchat.NL, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat.NL, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)# Not working!

cellchat.LS <- netAnalysis_computeCentrality(cellchat.LS, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat.LS, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)# Not working!

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ggCtr <- netAnalysis_signalingRole_scatter(cellchat.NL)
ggInf <- netAnalysis_signalingRole_scatter(cellchat.LS)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
ggCtr1 <- netAnalysis_signalingRole_scatter(cellchat.NL, signaling = c("COMPLEMENT", "MHC-I"))
ggInf1 <- netAnalysis_signalingRole_scatter(cellchat.LS, signaling = c("COMPLEMENT", "MHC-I"))

#> Signaling role analysis on the cell-cell communication network from user's input
ggCtr + ggCtr1
ggInf + ggInf1

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Control dataset
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.NL, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.NL, pattern = "incoming")
ht1 
ht2
# Infected dataset
ht3 <- netAnalysis_signalingRole_heatmap(cellchat.LS, pattern = "outgoing")
ht4 <- netAnalysis_signalingRole_heatmap(cellchat.LS, pattern = "incoming")
ht3 
ht4

# Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
library(NMF)
library(ggalluvial)

selectK(cellchat.NL, pattern = "outgoing") # Drops after 3
selectK(cellchat.LS, pattern = "outgoing") # Drops after 5

nPatterns = 3 # for control
cellchat.NL <- identifyCommunicationPatterns(cellchat.NL, pattern = "outgoing", k = nPatterns)
# river plot
netAnalysis_river(cellchat.NL, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function

nPatterns = 5 # for Infected
cellchat.LS <- identifyCommunicationPatterns(cellchat.LS, pattern = "outgoing", k = nPatterns)
# river plot
netAnalysis_river(cellchat.LS, pattern = "outgoing")


# Identify and visualize incoming communication pattern of target cells
selectK(cellchat.NL, pattern = "incoming") # Drops after 3
selectK(cellchat.LS, pattern = "incoming") # Drops after 3

nPatterns = 3 # for control
cellchat.NL <- identifyCommunicationPatterns(cellchat.NL, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat.NL, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function

nPatterns = 3 # for Infected
cellchat.LS <- identifyCommunicationPatterns(cellchat.LS, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat.LS, pattern = "incoming")

#Manifold and classification learning analysis of signaling networks
### Identify signaling groups based on their functional similarity
cellchat.NL <- computeNetSimilarity(cellchat.NL, type = "functional")
cellchat.NL <- netEmbedding(cellchat.NL, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat.NL <- netClustering(cellchat.NL, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat.NL, type = "functional", label.size = 3.5)

cellchat.LS <- computeNetSimilarity(cellchat.LS, type = "functional")
cellchat.LS <- netEmbedding(cellchat.LS, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat.LS <- netClustering(cellchat.LS, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat.LS, type = "functional", label.size = 3.5)

### Identify signaling groups based on structure similarity
cellchat.NL <- computeNetSimilarity(cellchat.NL, type = "structural")
cellchat.NL <- netEmbedding(cellchat.NL, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat.NL <- netClustering(cellchat.NL, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat.NL, type = "structural", label.size = 3.5)

cellchat.LS <- computeNetSimilarity(cellchat.LS, type = "structural")
cellchat.LS <- netEmbedding(cellchat.LS, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat.LS <- netClustering(cellchat.LS, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat.LS, type = "structural", label.size = 3.5)

netVisual_embeddingZoomIn(cellchat.NL, type = "structural", nCol = 2)
netVisual_embeddingZoomIn(cellchat.LS, type = "structural", nCol = 2)

########  Save Cellchat Objects   ####
saveRDS(cellchat.NL, file = "cellchat_bm_control.rds")
saveRDS(cellchat.LS, file = "cellchat_bm_Infected.rds")


    #NEXT STEPS #      

#Lift up CellChat object and merge together
#Since there are additional two populations (i.e., dermal DC and pericytes) specific to E14.5 compared to E13.5, we lift up cellchat.E13 by lifting up the cell groups to the same cell labels as E14.5. liftCellChat will only update the slot related to the cell-cell communication network, including slots object@net, object@netP and object@idents.

# Define the cell labels to lift up
group.new = levels(cellchat.LS@idents)
cellchat.NL <- liftCellChat(cellchat.NL, group.new)

#### Merge the cellchats
object.list <- list(WT= cellchat.NL, MAL = cellchat.LS)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

cellchat
saveRDS(cellchat, file = "cellchat_bm_comparisons.rds") # Saving the final cell chat object


#Visualize the inferred signaling network using the lifted object

# Hierarchy plot
pathways.show <- c("COMPLEMENT") 
pathways.show <- c("MHC-I") 
pathways.show <- c("IL16") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
vertex.receiver = seq(1,10) # Left portion of hierarchy plot the shows signaling to dermal cells and right portion shows signaling to epidermal cells
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

# Chord diagram
pathways.show <- c("COMPLEMENT", "MHC-I","IL16") 
par(mfrow = c(1,2), xpd=TRUE)
while (!is.null(dev.list()))  dev.off()
pdf(file ="ComparisonChordDiagram", width = 20, height =16)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]), font.size = 20)
}

dev.off()

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2



#Differential number of interactions or interaction strength among different cell populations

# Circle Plots (Red==Increases, Blue== decreased. Comaprison is btn second (Infected) as compared to first(Uninfected) dataset)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# Using a heat map(Top== Incoming signal, Right== Outgoing signal)
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 
gg2

mat <- cellchat@net[["WT"]]$weight
#par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = as.numeric(table(cellchat.NL@idents)), weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

mat <- cellchat@net[["MAL"]]$weight
#par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = as.numeric(table(cellchat.LS@idents)), weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

##### Number of  interactions

mat <- cellchat@net[["WT"]]$count
#par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = as.numeric(table(cellchat.NL@idents)), weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

mat <- cellchat@net[["MAL"]]$count
#par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = as.numeric(table(cellchat.LS@idents)), weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Stregth of interactions - ", names(object.list)[i]))
}


#Compare the major sources and targets in 2D space

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)


#we can identify the specific signaling changes of the cells that have changed between WT and MAL. 
  ## Identify signaling changes associated with one cell group
```{r}
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "DCs", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Macrophages", signaling.exclude = "MIF")
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Monocytes", signaling.exclude = c("MIF"))
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "T-Cells", signaling.exclude = c("MIF"))
gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Granulocytes", signaling.exclude = c("MIF"))

gg1
gg2
gg3
gg4
gg5


#Part II: Identify the conserved and context-specific signaling pathways

#Identify signaling networks with larger (or less) difference as well as signaling groups based on their functional/structure similarity

#Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)

#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)


#Compare the overall information flow of each signaling pathway
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(2, 1))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(2, 1))
#rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, color.use = cellchat@idents[["WT"]])
gg1 + gg2


# combining all the identified signaling pathways from different datasets 
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 7)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 7)
draw(ht1 + ht2, ht_gap = unit(1.0, "cm"))

#### Incoming signals
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 7, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 7, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(1.0, "cm"))

##### Overall signalling
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 7, color.heatmap = "OrRd", font.size = 5)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 7, color.heatmap = "OrRd", font.size = 5)
draw(ht1 + ht2, ht_gap = unit(1.0, "cm"))


# Visualize the upgulated and down-regulated signaling ligand-receptor pairs using Chord diagram
   ### Chord diagram

saveRDS(cellchat, file = "cellchat_bm_comparison.rds")

pdf(file  =  "UpregandDownreg_plot.pdf",
    width = 12, # The width of the plot in inches
    height = 12)

netVisual_chord_gene(object.list[[2]], sources.use = 2, targets.use = c(4:11), slot.name = 'net', net = net.up, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 2, targets.use = c(4:11), slot.name = 'net', net = net.down, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.up, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 5, targets.use = c(2,4,6:11), slot.name = 'net', net = net.up, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 5, targets.use = c(2,4,6:11), slot.name = 'net', net = net.down, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 6, targets.use = c(2,4,5,7:11), slot.name = 'net', net = net.up, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 6, targets.use = c(2,4,5,7:11), slot.name = 'net', net = net.down, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 7, targets.use = c(2,4:6,8:11), slot.name = 'net', net = net.up, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 7, targets.use = c(2,4:6,8:11), slot.name = 'net', net = net.down, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 8, targets.use = c(2,4:7,9:11), slot.name = 'net', net = net.up, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 8, targets.use = c(2,4:7,9:11), slot.name = 'net', net = net.down, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 9, targets.use = c(2,4:8,10:11), slot.name = 'net', net = net.up, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#netVisual_chord_gene(object.list[[1]], sources.use = 9, targets.use = c(2,4:8,10:11), slot.name = 'net', net = net.down, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 10, targets.use = c(2,4:9,11), slot.name = 'net', net = net.up, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#netVisual_chord_gene(object.list[[1]], sources.use = 10, targets.use = c(2,4:9,11), slot.name = 'net', net = net.down, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

#netVisual_chord_gene(object.list[[2]], sources.use = 11, targets.use = c(2,4:10), slot.name = 'net', net = net.up, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#netVisual_chord_gene(object.list[[1]], sources.use = 11, targets.use = c(2,4:10), slot.name = 'net', net = net.down, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

dev.off()


# Visualize the upgulated and down-regulated signaling ligand-receptor pairs using Chord diagram
#.....Only selected cells

pdf(file  =  "Selected_cells.pdf",
    width = 12, # The width of the plot in inches
    height = 12)

netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5,6,8,10), slot.name = 'net', net = net.up, lab.cex = 1.3, small.gap = 3.5)
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5,6,8,10), slot.name = 'net', net = net.down, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 5, targets.use = c(4,6,8,10), slot.name = 'net', net = net.up, lab.cex = 1.3, small.gap = 3.5)
netVisual_chord_gene(object.list[[1]], sources.use = 5, targets.use = c(4,6,8,10), slot.name = 'net', net = net.down, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 6, targets.use = c(4,5,8,10), slot.name = 'net', net = net.up, lab.cex = 1.3, small.gap = 3.5)
netVisual_chord_gene(object.list[[1]], sources.use = 6, targets.use = c(4,5,8,10), slot.name = 'net', net = net.down, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 8, targets.use = c(4,5,6,10), slot.name = 'net', net = net.up, lab.cex = 1.3, small.gap = 3.5)
netVisual_chord_gene(object.list[[1]], sources.use = 8, targets.use = c(4,5,6,10), slot.name = 'net', net = net.down, lab.cex = 1.3, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 10, targets.use = c(4,5,6,8), slot.name = 'net', net = net.up, lab.cex = 1.3, small.gap = 3.5)

dev.off()



