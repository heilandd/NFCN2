# R script:




# Read in spata object CellChat analysis



# Level 1 Exploration -----------------------------------------------------
library(CellChat)
library(tidyverse)

# load data 
object.all <- readRDS("~/Desktop/T-Cell Project/Revisions2/Analysis/seurat.obj_MNN_ref.RDS")

#Spata
samplr.read <- "275_T"
spata.obj <-  readRDS(paste0("~/Desktop/SpatialTranscriptomics/Visium/Visium/All_SPATA_Revisions/", samplr.read, "_SPATA_CNV_Pred.RDS"))



# For Level 1, cell will be screened by CellChat for interaction across clusters
# Remove clusters which not interesting for cell cell interaction
bc.to.remove <- NFCN2::getPoints(object.all, reduction="umap")
length(bc.to.remove)
cells.keep <- object.all@meta.data[!rownames(object.all@meta.data) %in% bc.to.remove, ] %>% rownames()
object.all <- subset(object.all, cells=cells.keep)

# Subset Doublets
obj.list <- NFCN2::findDouplets(object.all, subset = T)



# Save and process non doublets -------------------------------------------



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


#Import CD4/8 cluster annotations
setwd("~/Desktop/T-Cell Project/Revisions2/Analysis")
CD4 <- readRDS("CD4.Tcells.RDS")
CD8 <- readRDS("CD8.Tcells.RDS")
Seurat::DimPlot(CD4)
Seurat::DimPlot(CD8)

inter <- intersect(object@meta.data %>% rownames() , CD4@meta.data %>% rownames())
CD4 <- subset(CD4, cells=inter)

inter <- intersect(object@meta.data %>% rownames() , CD8@meta.data %>% rownames())
CD8 <- subset(CD8, cells=inter)

object@meta.data$seurat_clusters <- as.character(object@meta.data$seurat_clusters)
object@meta.data[CD4@meta.data %>% rownames(), ]$seurat_clusters <- paste0("CD4_", CD4@meta.data$seurat_clusters)
object@meta.data[CD8@meta.data %>% rownames(), ]$seurat_clusters <- paste0("CD8_", CD8@meta.data$seurat_clusters)




# Myeloid cell types ------------------------------------------------------


# Define Myeloid cells


myeloid <- subset(object, cells=rownames(object@meta.data %>% filter(seurat_clusters %in% source)))
DimPlot(myeloid, reduction = "ref.umap")

keep <- myeloid@meta.data[myeloid@meta.data$predicted.celltype.l2 %in% c("CD16 Mono","CD14 Mono",  "cDC2"  ,"cDC1"  ,"pDC"  ), ] %>% rownames()
myeloid <- subset(myeloid, cells=keep)

DimPlot(myeloid, reduction = "ref.umap")

myeloid@meta.data

#NMF Cluster
myeloid <- STutility::RunNMF(myeloid)
myeloid <- FindNeighbors(myeloid, dims = 1:10, reduction = "NMF")
myeloid <- FindClusters(myeloid, resolution = 0.2)
myeloid <- Seurat::RunUMAP(myeloid, reduction="NMF", dims=1:20)
DimPlot(myeloid, reduction = "umap", label = T)

FeaturePlot(object = myeloid, 
            reduction="umap", 
            features = "CD163", 
            pt.size = 1, order=T)+
  scale_color_viridis_c(option="C", limit=c(2,3), oob = scales::squish)

FeaturePlot(object = myeloid, 
            reduction="umap", 
            features = "HMOX1", 
            pt.size = 1, order=T)+
  scale_color_viridis_c(option="C", limit=c(2.5,4), oob = scales::squish)

Compare_Barplot(myeloid, Feat="predicted.celltype.l2", compare_to = "seurat_clusters")
Compare_Barplot(myeloid, Feat="sample", compare_to = "seurat_clusters")



#Marker gene and Heatmap
myeloid.top.nmf <- Seurat::FindAllMarkers(myeloid)




#MNN Cluster

myeloid <- SeuratWrappers::RunFastMNN(Seurat::SplitObject(myeloid, "sample"))
myeloid <- SCTransform(myeloid, method = "glmGamPoi", vars.to.regress = c("percent.mt", "percent.RB"))
myeloid <- RunPCA(myeloid)
myeloid <- RunUMAP(myeloid, reduction = "mnn", dims = 1:30)

myeloid <- FindNeighbors(myeloid, dims = 1:30, reduction = "mnn")
myeloid <- FindClusters(myeloid, resolution = 0.3)
DimPlot(myeloid, reduction = "umap")
DimPlot(myeloid, reduction = "umap", group.by = "sample")

#Marker gene and Heatmap
Diff.myeloid.mnn <- Seurat::FindAllMarkers(myeloid)

#ref.spca Cluster

myeloid <- RunUMAP(myeloid, reduction = "ref.spca", dims = 1:30)
myeloid <- FindNeighbors(myeloid, dims = 1:30, reduction = "ref.spca")
myeloid <- FindClusters(myeloid, resolution = 0.3)
DimPlot(myeloid, reduction = "umap")
DimPlot(myeloid, reduction = "umap", group.by = "sample")


Diff.myeloid.ref <- Seurat::FindAllMarkers(myeloid)

Compare_Barplot(myeloid, Feat="sample", compare_to = "seurat_clusters")
Compare_Barplot(myeloid, Feat="predicted.celltype.l2", compare_to = "seurat_clusters")







#Marker gene and Heatmap
#Diff.myeloid <- Seurat::FindAllMarkers(myeloid)

#Get teh top 10 genes of each cluster
myeloid.top <- Diff.myeloid.ref
myeloid.top <- myeloid.top %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)


mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(100)
DoHeatmap(myeloid, features = myeloid.top$gene) + NoLegend() + scale_fill_gradientn(colours = rev(mapal))


# Estimate the mean gene expression in each cluster

genes.go <- intersect(myeloid.top$gene, myeloid@assays$SCT@scale.data %>% rownames())
df.plot.heatmap <- purrr::map(.x=genes.go, function(gene){
  
  clust.arrange <- unique(myeloid@meta.data$seurat_clusters) %>% sort() %>% as.character()
  
  mean.per.cluster <- purrr::map_df(.x=clust.arrange, function(clusters){
    
    cells <- 
      myeloid@meta.data %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column("barcodes") %>% 
      filter(seurat_clusters=={{clusters}}) %>% 
      pull(barcodes)
    
    out <- data.frame(C1 =myeloid@assays$SCT@scale.data[gene, cells] %>% mean())
    
  }) %>% t() %>% as.data.frame()
  names(mean.per.cluster) <- paste0("C_", clust.arrange)
  rownames(mean.per.cluster) <- gene
  
  return(mean.per.cluster)
}) %>% do.call(rbind, .) %>% as.data.frame()

#Scale and remove outliers
for(i in names(df.plot.heatmap)){
  
  quantil <- c(stats::quantile(df.plot.heatmap[,i], 0.2),
               stats::quantile(df.plot.heatmap[,i], 0.8)
  )
  max <- df.plot.heatmap[df.plot.heatmap[i]>quantil[2] ,i] %>% max()
  min <- df.plot.heatmap[df.plot.heatmap[i]<quantil[1] ,i] %>% min()
  
  v.max <- df.plot.heatmap[,i]>quantil[2]
  v.min <- df.plot.heatmap[,i]<quantil[1]
  
  #oob limits
  
  df.plot.heatmap[v.max, i] <- max
  df.plot.heatmap[v.min, i] <- min
  
  #rescale
  df.plot.heatmap[,i] <- scales::rescale(df.plot.heatmap[,i], c(0,1))
  
}


pheatmap::pheatmap(as.matrix(df.plot.heatmap),
                   scale="none",
                   cluster_rows = F,
                   cluster_cols = F,
                   color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "PRGn")))(100))


cell.types <- c("aM1", "aM2", "microM", "moM", "Mono", "cDC2", "prolif", "cDC2")
cell.types <- myeloid@meta.data %>% mutate(cell.type = case_when(
  seurat_clusters==0~ cell.types[1],
  seurat_clusters==1~ cell.types[2],
  seurat_clusters==2~ cell.types[3],
  seurat_clusters==3~ cell.types[4],
  seurat_clusters==4~ cell.types[5],
  seurat_clusters==5~ cell.types[6],
  seurat_clusters==6~ cell.types[7],
  seurat_clusters==7~ cell.types[8],
)) %>% pull(cell.type)
myeloid@meta.data$seurat_clusters <- cell.types


# Annotate Cell Types ------------------------------------------------------

object@meta.data[myeloid@meta.data %>% rownames(), ]$seurat_clusters <- myeloid@meta.data$seurat_clusters

object@meta.data$seurat_clusters <- as.factor(object@meta.data$seurat_clusters)

object <- Seurat::SetIdent(object, value="seurat_clusters")
DimPlot(object,reduction="ref.umap")

object@reductions$umap <- object@reductions$ref.umap
colnames(object@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")

target <- c(paste0("CD4_", CD4@meta.data$seurat_clusters %>% unique()),
            paste0("CD8_", CD8@meta.data$seurat_clusters%>% unique()))

source <- c("aM1", "aM2", "microM", "moM", "Mono", "cDC2", "prolif", "cDC2")

saveRDS(object, "Seurat_all_processed_28_08_NFCN.RDS")

object <- NFCN2::createSlotNFCN(object,target=target, source=source)
object <- NFCN2::inferConnectedPathways(object, pathway = "IL10")
object <- NFCN2::inferNFCN(object, model.span = 0.3)

NFCN2::plotNFCN(object)
NFCN2::plotNFCNDimRed(object)
NFCN2::plotNFCNDotplot(object)

NFCN2::plotNFCNDimRed3D(object, pt.size = 1.5, l.alpha = 0.5, l.color = "black", l.width = 0.4)

future::plan("multiprocess", workers = 1)
inferSpatialPosition.object <- NFCN2::inferSpatialPosition(object, spata.obj)



# Interactive Surface plot
SPATA2::plotSurfaceInteractive(inferSpatialPosition.object[[2]])

#3Dplot
NFCN2::plotNFCNSpatial3D(object, inferSpatialPosition.object, l.width = 0.5, pt.size = 2,l.alpha = 0.2)

#Correct for patial position
object <- NFCN2::inferSpatialCorrection(object, inferSpatialPosition.object, assay="SCT", estimator=1, nr.genes=20)
NFCN2::plotNFCN(object, correction.distance = 150, pt.size = 2, pt.alpha = 0.9)+theme_classic()
NFCN2::plotNFCNDimRed(object, correction.distance = 150, pt.con=c(1,0.1))+theme_classic()

FeaturePlot(object = object, 
            reduction="umap", 
            features = "HMOX1", 
            pt.size = 1, order=T)+
  scale_color_viridis_c(option="C", limit=c(3,5), oob = scales::squish)

FeaturePlot(object = object, 
            reduction="umap", 
            features = "IL10", 
            pt.size = 0.5, order=T)+
  scale_color_viridis_c(option="C", limit=c(1,2), oob = scales::squish)


inferSpatialPosition.object <- NFCN2::inferSpatialPosition(object, spata.obj, correction.distance = 100)






# Get some correlation plots

df.cor <- SPATA2::getFeatureDf(inferSpatialPosition.object[[2]]) %>% NFCN2::getCleaned(., feat="top_inter_connected", q=0.01)

ggplot(df.cor, aes(x=top_inter_connected, y=AC))+geom_point()+theme_classic()+
  geom_smooth(se=F, model="lm")

ggplot(df.cor, aes(x=top_inter_connected, y=Chr17))+geom_point()+theme_classic()+
  geom_smooth(se=F, model="lm")

ggplot(df.cor, aes(x=top_inter_connected, y=MES))+geom_point()+theme_classic()+
  geom_smooth(se=F, model="lm")


# Functional enrichment

functional <- NFCN2::findFunctionalEnrichment(inferSpatialPosition.object$marker$connected)

NFCN2::plotFunctionalDotplot(functional[[1]], top=30)  
NFCN2::plotFunctionalDotplot(functional[[2]], top=30)

clusterProfiler::cnetplot(functional[[1]], showCategory = 10, colorEdge = TRUE)



# Doublet analysis --------------------------------------------------------

library(Seurat)

doublet.seurat <- obj.list[[2]]

doublet.seurat %>% DimPlot()

object <- readRDS("~/Desktop/T-Cell Project/Revisions2/Analysis/Seurat_all_processed_28_08_NFCN.RDS")



# Isolate CD3 - AIF1 pos doublets

doublet.seurat <- 
  doublet.seurat %>% 
  FindNeighbors(reduction="mnn") %>% 
  FindClusters(resolution=0.4)

doublet.seurat %>% DimPlot(label = F)



doublet.seurat <- doublet.seurat %>% RunUMAP(reduction="mnn",dims=1:10)
doublet.seurat %>% FeaturePlot(feature="AIF1", reduction = "umap")
doublet.seurat %>% FeaturePlot(feature="CD3D", reduction = "umap")
FeaturePlot(object = doublet.seurat, 
            reduction="umap", 
            features = "AIF1", 
            pt.size = 3, order=T)+
  scale_color_viridis_c(option="C", limit=c(1,3), oob = scales::squish)
FeaturePlot(object = doublet.seurat, 
            reduction="umap", 
            features = "CD3D", 
            pt.size = 3, order=T)+
  scale_color_viridis_c(option="C", limit=c(0.5,2), oob = scales::squish)


FeaturePlot(object = doublet.seurat, 
            reduction="umap", 
            features = "HMOX1", 
            pt.size = 3, order=T)+
  scale_color_viridis_c(option="C", limit=c(1,3), oob = scales::squish)

FeaturePlot(object = doublet.seurat, 
            reduction="umap", 
            features = "CD163", 
            pt.size = 3, order=T)+
  scale_color_viridis_c(option="C", limit=c(1,3), oob = scales::squish)

marker.doub <- FindAllMarkers(doublet.seurat)
top.doub <- marker.doub %>% group_by(cluster) %>% top_n(10, wt=avg_log2FC)

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(100)
DoHeatmap(doublet.seurat, features = top.doub$gene) + scale_fill_gradientn(colours = rev(mapal))



# Doublets cluster 1 

keep <- doublet.seurat@meta.data %>% as.data.frame() %>% filter(seurat_clusters %in% c(1)) %>% rownames()
doublet.seurat.select <- subset(doublet.seurat, cells=keep)
doublet.seurat.select@meta.data$type="D"

#The most connected cells
ligand <- subset(object, cells=object@assays$NFCN$model$top.partner %>% filter(distance<150) %>% pull(Ligand))
ligand@meta.data$type="L"

receptor <- subset(object, cells=object@assays$NFCN$model$top.partner %>% filter(distance<150) %>% pull(Receptor))
receptor@meta.data$type="R"

# combibe 

doublet.seurat.comb <-  merge(doublet.seurat.select, ligand)
doublet.seurat.comb <- merge(doublet.seurat.comb, receptor)

# Heatmap
doublet.seurat.comb <- 
doublet.seurat.comb %>% Seurat::SCTransform() %>% RunPCA() %>% RunUMAP(dims=1:10)
DimPlot(doublet.seurat.comb, group.by = "type", reduction = "pca", pt.size = 3)

doublet.seurat.comb <- Seurat::SetIdent(doublet.seurat.comb, value="type")

FeaturePlot(object = doublet.seurat.comb, 
            reduction="pca", 
            features = "HMOX1", 
            pt.size = 3, order=T)+
  scale_color_viridis_c(option="C", limit=c(1,3), oob = scales::squish)

FeaturePlot(object = doublet.seurat.comb, 
            reduction="pca", 
            features = "CD163", 
            pt.size = 3, order=T)+
  scale_color_viridis_c(option="C", limit=c(1,2), oob = scales::squish)



marker.doub <- FindAllMarkers(doublet.seurat.comb)
top.doub <- marker.doub %>% filter(cluster!="D") %>% group_by(cluster) %>% top_n(30, wt=avg_log2FC)

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(100)
DoHeatmap(doublet.seurat.comb, features = top.doub$gene) + scale_fill_gradientn(colours = rev(mapal))

plot.df <- doublet.seurat.comb@meta.data
ggplot(data=plot.df, mapping = aes(x=1:nrow(doublet.seurat.comb@meta.data), y=nCount_RNA, color=type))+theme_classic()+geom_point()




#  Integartion of the spatial data sets -----------------------------------

library(SPATA2)
library(tidyverse)
spata.obj <- SPATA::loadSpataObject("~/Desktop/T-Cell Project/Revisions/RNAVelo/Spatial/spata-obj-S4.RDS")
spata.obj@fdata$segmentation <- ""

create.SPATA2 <- function(spata.obj){
  spata.obj@fdata$segmentation <- ""

  spata.obj <- SPATA2::initiateSpataObject_CountMtr(spata.obj@data@counts, 
                                                    coords_df=spata.obj@coordinates,
                                                    sample_name = spata.obj@samples)
  spata.obj <-  SPATA2::runAutoencoderDenoising(spata.obj, activation="relu", bottleneck=8)

  #Add gene sets
  
  spata.obj@used_genesets <- readRDS("~/Desktop/SpatialTranscriptomics/Visium/Visium/GS_new.RDS")
  
  return(spata.obj)
}

saveRDS(spata.ob.S4, "~/Desktop/T-Cell Project/Revisions/RNAVelo/Spatial/spata-obj-S4.RDS")

spata.ob.S5 <- create.SPATA2(SPATA::loadSpataObject("~/Desktop/T-Cell Project/Revisions/RNAVelo/Spatial/spata-obj-S5.RDS"))
spata.ob.S4 <- create.SPATA2(SPATA::loadSpataObject("~/Desktop/T-Cell Project/Revisions/RNAVelo/Spatial/spata-obj-S4.RDS"))
spata.ob.S3 <- create.SPATA2(SPATA::loadSpataObject("~/Desktop/T-Cell Project/Revisions/RNAVelo/Spatial/spata-obj-S3.RDS"))

spata.ob.S5 <- NFCN2::inferSpatialPosition(object, spata.ob.S5)[[2]]
spata.ob.S4 <- NFCN2::inferSpatialPosition(object, spata.ob.S4)[[2]]
spata.ob.S3 <- NFCN2::inferSpatialPosition(object, spata.ob.S3)[[2]]

SPATA2::plotSurfaceInteractive(spata.ob.S5)


#Function
spatial.plot <- function(object,feature, plot=T){
  
  
  #Run analysis 
  message(paste0(Sys.time(), "  Run MC Simulation for Spatial  Correlations analysis"))
  
  cor.mat <- matrix(NA, length(feature),length(feature));colnames(cor.mat) = rownames(cor.mat) <- feature
  cor.mat <- reshape2::melt(cor.mat)
  
  # fill data 
  cor.mat$value <- pbmcapply::pbmclapply(1:nrow(cor.mat), function(i){
    cor.out <- spatial.mc(P1=SPATA2::getFeatureDf(object) %>% pull(cor.mat[i,1]),
                          P2=SPATA2::getFeatureDf(object) %>% pull(cor.mat[i,2]),
                          n=200)
    return(cor.out$Cor)
    
  }, mc.cores = 8) %>% unlist()
  
  message(paste0(Sys.time(), "  Run Auto Corelation analysis"))
  
  cor.Auto <- spatial.ac(object, feature)
  cor.Auto <- as.data.frame(do.call(rbind,cor.Auto)) %>% 
    tibble::rownames_to_column("Var2") %>% 
    mutate(Var1="Autocorelation",
           value=V1) %>%
    dplyr::select(Var1,Var2,value)
  
  cor.mat <- rbind(cor.mat, cor.Auto)
  
  
  message(paste0(Sys.time(), " Plotting "))
  
  if(plot==T){
    corrplot::corrplot.mixed(t(reshape2::acast(cor.mat, Var1~Var2, value.var="value")))
  }
  
  
  return(cor.mat)
  
}
spatial.mc <- function(P1,P2, n = 599){
  # Monte-Carlo Simulation of spatial correlation
  if(length(P2)!=length(P1)) stop("Unequal Inputs")
  message(paste0(Sys.time(), "Start Model"))
  M <- glm(P1 ~ P2)
  #coef(M)[2]
  message(paste0(Sys.time(), "Start MC Simulation"))
  #MC
  I.r <- purrr::map_dbl(.x=1:n, .f=function(i){
    x <- sample(P1, replace=FALSE)
    y <- P2
    # Compute new set of lagged values
    #x.lag <- lag.listw(lw, x)
    # Compute the regression slope and store its value
    M.r    <- glm(y ~ x)
    I.r <- coef(M.r)[2]
    return(I.r)
  })
  
  #hist(I.r, main=NULL, xlab="Spatial-Cor-MC", las=1, xlim=c(-0.5,0.5))
  #abline(v=coef(M)[2], col="red")
  
  #p-value
  message(paste0(Sys.time(), "P-Value Extraction"))
  N.greater <- sum(coef(M)[2] > I.r)
  p <- min(N.greater + 1, n + 1 - N.greater) / (n + 1)
  p
  out <- list(coef(M)[2], p)
  names(out) <- c("Cor", "p")
  return(out)
}
spatial.ac <- function(object,feature){
  
  message(paste0(Sys.time(), "  Start Autocorrelation"))
  
  tissue.pos <- SPATA2::getCoordsDf(object)
  plot.new()
  message(paste0(Sys.time(), "  Transforme in spatial dataset"))
  segments <- pbmcapply::pbmclapply(1:nrow(tissue.pos), function(xx){
    
    segment_c <-plotrix::draw.circle(x=tissue.pos$x[xx], 
                                     y=tissue.pos$y[xx],
                                     radius=5*2) 
    
    segment <- as.data.frame(do.call(cbind, segment_c))
    
    segment <- sp::Polygon(cbind(segment$x, segment$y))
    segment <- sp::Polygons(list(segment), tissue.pos$barcodes[xx])
    return(segment)
  }, mc.cores = 8)
  
  SpP = sp::SpatialPolygons(segments, 1:length(segments))
  
  attr <- SPATA2::getFeatureDf(object) %>% as.data.frame()
  rownames(attr) <- attr$barcodes
  attr <- attr[,feature]
  
  SrDf = sp::SpatialPolygonsDataFrame(SpP, attr)
  
  message(paste0(Sys.time(), "  Define neighboring Spots"))
  
  nb <- spdep::poly2nb(SrDf, queen=F, snap=25)
  lw <- spdep::nb2listw(nb, style="W", zero.policy=TRUE)
  
  message(paste0(Sys.time(), "  Computing the Moran’s I statistic"))
  
  stat <- purrr::map(.x=feature, .f=function(x){
    print(x)
    model <- spdep::moran.test(SrDf@data[,x],lw, zero.policy=T)
    model$estimate[[1]]
  })
  
  names(stat) <- feature
  
  return(stat)
  
  
}



# Load data and estimate spatial correlation: Modules  ------------------------------

spata.ob <- spata.ob.S3


features <- spata.ob %>% SPATA2::getFeatureDf() %>% dplyr::select(10:11) %>% names()
gs <- c("Neftel_AClike","Neftel_OPClike", "Neftel_Mes_Comb","Neftel_NPC_Comb")

spata.ob <- 
  spata.ob %>% 
  SPATA2::joinWithGeneSets(., gene_sets = gs) %>% 
  dplyr::select(barcodes, {{gs}}) %>% 
  SPATA2::addFeatures(spata.ob, .)

export <- spatial.plot(object=spata.ob, feature = c(features,gs) , plot=F)


export <- export %>% filter(Var1 %in% gs) %>% filter(!Var2 %in% gs) #%>% mutate(value=scales::rescale(value, c(0,1)))

col2 = colorRampPalette(c('#67001F', '#B2182B', '#D6604D', '#F4A582',
                          '#FDDBC7', '#FFFFFF', '#D1E5F0', '#92C5DE',
                          '#4393C3', '#2166AC', '#053061'))
corrplot::corrplot(t(reshape2::acast(export, Var1~Var2, value.var="value")), is.corr=F, col=rev(col2(200)))


#Interactive Surface plot
SPATA2::plotSurfaceInteractive(inferSpatialPosition.object[[2]])




















