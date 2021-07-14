# R script:




# Read in spata object CellChat analysis



# Level 1 Exploration -----------------------------------------------------
library(CellChat)


# load data 
object <- readRDS("~/Desktop/T-Cell Project/Revisions2/Analysis/seurat.obj_MNN_ref.RDS")
samplr.read <- "275_T"
spata.obj <-  readRDS(paste0("~/Desktop/SpatialTranscriptomics/Visium/Visium/All_SPATA_Revisions/", samplr.read, "_SPATA_CNV_Pred.RDS"))


Seurat::DimPlot(object, reduction="umap")



# For Level 1, cell will be screened by CellChat for interaction across clusters
# Remove clusters which not interesting for cell cell interaction
bc.to.remove <- NFCN2::getPoints(object, reduction="umap")
cells.keep <- object@meta.data[!rownames(object@meta.data) %in% bc.to.remove, ] %>% rownames()
object <- subset(object, cells=cells.keep)
Seurat::DimPlot(object, reduction="umap")

# Subset Doublets

obj.list <- NFCN2::findDouplets(object, subset = T)


object <- obj.list[[1]]
Seurat::DimPlot(object, reduction="umap")

# Re-Run basic analysis 
object <- Seurat::RunUMAP(object, reduction="mnn", dims=1:30)
object <- Seurat::FindNeighbors(object, reduction="mnn")
object <- Seurat::FindClusters(object, reduction="mnn",resolution=0.7)
Seurat::DimPlot(object, reduction="umap",label=T)
Seurat::DimPlot(object, reduction="umap", group.by = "predicted.celltype.l2")


#save object
saveRDS(object, "~/Desktop/T-Cell Project/Revisions2/Analysis/seurat.obj_cleaned_inter.RDS")
object <- readRDS("~/Desktop/T-Cell Project/Revisions2/Analysis/seurat.obj_cleaned_inter.RDS")

# CellChat


# Prepare Data ------------------------------------------------------------



library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

object@meta.data$Clusters_Cell_Chat <- paste0("cluster_",object@meta.data$seurat_clusters)

cellchat <- CellChat::createCellChat(object, group.by = "Clusters_Cell_Chat")
cellchat@DB <- CellChat::CellChatDB.human
cellchat <- CellChat::subsetData(cellchat)

future::plan("multiprocess", workers = 8)

cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- CellChat::projectData(cellchat, PPI.human)

cellchat <- CellChat::computeCommunProb(cellchat)
cellchat <- CellChat::filterCommunication(cellchat, min.cells = 1)
df.net <- CellChat::subsetCommunication(cellchat, thresh = 1)


cellchat@DB$interaction %>% filter(ligand=="IL10")
cellchat@LR$LRsig %>% filter(ligand=="IL6")




# Visualisation -----------------------------------------------------------


cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


pathways.show <- c("TGFb") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair

# Hierarchy plot
vertex.receiver = seq(1,4) 
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver, layout = "hierarchy")

cellchat@meta
netVisual_bubble(cellchat, sources.use = c(2:9), targets.use = c(1,11,10), remove.isolate = FALSE)












# NFCN2

# Select pathway and source-Target population -----------------------------

object <- Seurat::FindNeighbors(object, reduction="mnn")
object <- Seurat::FindClusters(object, reduction="mnn",resolution=0.6)
Seurat::DimPlot(object, reduction="umap",label=T)
Seurat::DimPlot(object, reduction="umap", group.by = "predicted.celltype.l2", label=T)

target <- c(2,9,8)
source <- c(0:10)[!0:10 %in% target]


# NFCN Analysis -----------------------------
object <- NFCN2::createSlotNFCN(object,target=target, source=source)
object <- NFCN2::inferConnectedPathways(object, pathway = "IL10")
object <- NFCN2::inferNFCN(object, model.span = 0.3)

NFCN2::plotNFCN(object)
NFCN2::plotNFCNDimRed(object)
NFCN2::plotNFCNDotplot(object)
NFCN2::plotNFCNDimRed3D(object, pt.size = 1.5, l.alpha = 0.5, l.color = "black", l.width = 0.4)

future::plan("multiprocess", workers = 1)
ls.objects <- NFCN2::inferSpatialPosition(object, spata.obj)

# Interactive Surface plot
SPATA2::plotSurfaceInteractive(ls.objects[[2]])

#3Dplot
NFCN2::plotNFCNSpatial3D(object, ls.objects, l.width = 0.5, pt.size = 2,l.alpha = 0.2)




# Get some correlation plots

df.cor <- SPATA2::getFeatureDf(ls.objects[[2]]) %>% NFCN2::getCleaned(., feat="top_inter_connected", q=0.01)

ggplot(df.cor, aes(x=top_inter_connected, y=AC))+geom_point()+theme_classic()+
  geom_smooth(se=F, model="lm")

ggplot(df.cor, aes(x=top_inter_connected, y=Chr17))+geom_point()+theme_classic()+
  geom_smooth(se=F, model="lm")

ggplot(df.cor, aes(x=top_inter_connected, y=MES))+geom_point()+theme_classic()+
  geom_smooth(se=F, model="lm")


# Functional enrichment

functional <- NFCN2::findFunctionalEnrichment(ls.objects$marker$connected)

NFCN2::plotFunctionalDotplot(functional[[1]], top=30)  
NFCN2::plotFunctionalDotplot(functional[[2]], top=30)

clusterProfiler::cnetplot(functional[[1]], showCategory = 10, colorEdge = TRUE)


































