#' @title inferConnectedPathways
#' @author Dieter Henrik Heiland
#' @description inferConnectedPathways
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

inferConnectedPathways <- function(object, 
                                   pathway, 
                                   q=0.1, 
                                   group.by="seurat_clusters", 
                                   smooth = T,
                                   smooth_span = 0.02){
  
  # Collect information
  
  select <- CellChat::CellChatDB.human$interaction %>% filter(pathway_name %in% pathway)
  
  ligand <- purrr::map(.x=select$interaction_name_2, .f=function(x){
    ligand <- stringr::str_split(x, pattern=" - ") %>% unlist()
    return(ligand[1])
  } ) %>% unlist() %>% unique()
  
  receptor <- purrr::map(.x=select$interaction_name_2, .f=function(x){
    receptor <- c(stringr::str_split(x, pattern=" - ") %>% unlist())[2]
    receptor <- stringr::str_remove_all(receptor, pattern="[)]") %>% unlist()
    receptor <- stringr::str_remove_all(receptor, pattern="[(]") %>% unlist()
    receptor <- stringr::str_split(receptor, pattern="[+]") %>% unlist()
    
    return(receptor)
  } ) %>% unlist() %>% unique()
  

  
  # Infer pathway information
  
  pw.ligand <- SPATA2::gsdf %>% filter(gene %in% ligand) %>% pull(ont) %>% unique()
  pw.ligand <- SPATA2::gsdf %>% filter(ont %in% pw.ligand) %>% pull(gene) %>% unique()
  
  pw.receptor <- SPATA2::gsdf %>% filter(gene %in% receptor) %>% pull(ont) %>% unique()
  pw.receptor <- SPATA2::gsdf %>% filter(ont %in% pw.receptor) %>% pull(gene) %>% unique()
  
  
  
  #filter genes by target/source expression
  
  genes <- Seurat::GetAssayData(object) %>% rownames()
  
  message("## 1/4 Set up soure up-stream")
  source.mat <- Seurat::GetAssayData(object)[intersect(genes, pw.ligand), NFCN2::getSource(object)]
  #print(dim(source.mat))
  upper <- base::rowMeans(source.mat %>% as.matrix()) %>% stats::quantile(1-q) %>% as.numeric()
  pw.ligand <- base::rowMeans(source.mat %>% as.matrix()) %>% as.data.frame() %>% dplyr::filter(.>upper) %>% rownames()
  pw.ligand <- data.frame(ont="NFCN_Ligand", gene=pw.ligand)
  
  message("## 2/4 Set up soure down-stream")
  target.mat <- Seurat::GetAssayData(object)[intersect(genes, pw.receptor), NFCN2::getTargets(object)]
  #dim(target.mat)
  upper <- rowMeans(target.mat %>% as.matrix()) %>% stats::quantile(1-q) %>% as.numeric()
  pw.receptor <- rowMeans(target.mat %>% as.matrix()) %>% as.data.frame() %>% dplyr::filter(.>upper) %>% rownames()
  pw.receptor <- data.frame(ont="NFCN_Receptor", gene=pw.receptor)
  
  message("## 3/4 Infer scores for up- and downstream activation")
  spata.obj <- SPATA2::transformSeuratToSpata(object, sample_name="NFCN", method='single_cell',assay_name="RNA", coords_from = "umap", gene_set_path=SPATA2::gsdf, verbose=F)
  
  
  #Add new pathways
  spata.obj@used_genesets <- SPATA2::gsdf
  spata.obj@used_genesets <- rbind(spata.obj@used_genesets, rbind(pw.ligand, pw.receptor))
  
  pw <- c("NFCN_Ligand", "NFCN_Receptor")
  
  receptor <- intersect(SPATA2::getGenes(spata.obj), receptor)
  
  message("## 4/4 Prepare data")
  data.nfcn <- 
    SPATA2::joinWithGeneSets(spata.obj, gene_sets = pw, verbose = F, smooth=smooth, smooth_span=smooth_span) %>% 
    dplyr::left_join(., SPATA2::joinWithFeatures(spata.obj, feature=group.by, verbose = F) %>% dplyr::select(-sample, x, y)) %>% 
    dplyr::left_join(., SPATA2::joinWithGenes(spata.obj, genes=ligand, verbose = F,smooth=smooth, smooth_span=smooth_span,average_genes=T) %>% dplyr::select(-sample, x, y)) %>% 
    dplyr::rename("ligand":=mean_genes) %>% 
    dplyr::left_join(., SPATA2::joinWithGenes(spata.obj, genes=receptor, verbose = F,smooth=smooth, smooth_span=smooth_span, average_genes=T) %>% dplyr::select(-sample, x, y)) %>% 
    dplyr::rename("receptor":=mean_genes) %>% 
    as.data.frame()
  
  
  object@assays$NFCN$data <- data.nfcn
  
  return(object)
  
  
  } 
 
#' @title inferNFCN
#' @author Dieter Henrik Heiland
#' @description inferNFCN
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

inferNFCN <- function(object,
                      group.by="seurat_clusters",
                      model.span=0.5,
                      q=0.1,
                      NN=500){
  
  message("Start Fitting the Model by using sc_RNAseq Data")
  #Fit model
  norm_s=function(x){(x-min(x))/(max(x)-min(x))}
  
  data <- object@assays$NFCN$data
  
  
  
  message("## 1/3  Fitting Receptor")
  score.r <- 
    data %>% 
    dplyr::filter(barcodes %in% NFCN2::getTargets(object)) %>% 
    dplyr::filter(NFCN_Receptor<quantile(NFCN_Receptor, 1-q)) %>% 
    dplyr::filter(NFCN_Receptor>quantile(NFCN_Receptor, q)) %>% 
    dplyr::mutate(x=2-scales::rescale(NFCN_Receptor, c(0,1) )) %>% 
    dplyr::mutate(y=2-scales::rescale(receptor, c(0,1) )) %>% 
    dplyr::select(barcodes,x,y, {{group.by}})%>% 
    dplyr::mutate(type="R")

  message("## 2/3  Fitting Ligand ")
  score.l <- 
    data %>% 
    dplyr::filter(barcodes %in% NFCN2::getSource(object)) %>% 
    dplyr::filter(NFCN_Ligand<quantile(NFCN_Ligand, 1-q)) %>% 
    dplyr::filter(NFCN_Ligand>quantile(NFCN_Ligand, q)) %>% 
    dplyr::mutate(x=scales::rescale(NFCN_Ligand, c(0,1) )) %>% 
    dplyr::mutate(y=scales::rescale(ligand, c(0,1) )) %>% 
    dplyr::select(barcodes,x,y, {{group.by}}) %>% 
    dplyr::mutate(type="L")
  
  
  df.plot <- rbind(score.r, score.l)
  
  
  
  message("## 3/3  Fitting Model ")
  model <- loess(y~x, df.plot, span = model.span)
  f <- c(seq(from=7, to=0, length.out = nrow(df.plot)/2),
         seq(to=0, from=7, length.out = nrow(df.plot)/2)) %>% exp()
  
  summary(model)
  

  df.plot <- 
    df.plot %>% 
    dplyr::arrange(desc(x)) %>% 
    dplyr::mutate(x.p=df.plot$x %>% jitter(., factor=f)) %>% 
    dplyr::mutate(y.p=predict(model, df.plot$x) %>% jitter(., factor=f))

  message("## Quantify NFCN's ")
  
  #Post Process and cleaning
  receptor.top <- 
    df.plot %>% 
    dplyr::filter(type=="R") %>% 
    dplyr::mutate(sum.f = x.p+y.p) %>% 
    dplyr::arrange((sum.f)) %>% 
    utils::head(., NN)
  
  ligand.top <- 
    df.plot %>% 
    dplyr::filter(type=="L") %>% 
    dplyr::mutate(sum.f = x.p+y.p) %>% 
    dplyr::arrange(desc(sum.f)) %>% 
    utils::head(., NN)
  
  top.partner <- data.frame(Ligand=ligand.top$barcodes, 
                            Receptor=receptor.top$barcodes)
  
  
  df.plot$top=0;df.plot[df.plot$barcodes %in% c(top.partner$Ligand, top.partner$Receptor),  ]$top=1
  
  
  R2 <- cor(df.plot$x.p, df.plot$y.p)^2
  p.vale <- cor.test(df.plot$x.p, df.plot$y.p)$p.value
  
  
  out.list <- list(model, df.plot, receptor.top, ligand.top, top.partner, R2, p.vale)
  names(out.list) <- c("model", "df.plot", "receptor.top", "ligand.top", "top.partner", "R2", "p.value")
  
  object@assays$NFCN$model <- out.list
  
  return(object)
}

#' @title inferSpatialPosition
#' @author Dieter Henrik Heiland
#' @description inferSpatialPosition
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

inferSpatialPosition <- function(object, 
                                 spata.obj,
                                 group.by="seurat_clusters",
                                 NN=500,
                                 intergrate.group.by=F,
                                 RL_interaction=F,
                                 correction.distance=NULL){
  
  future::plan("multiprocess", workers = 1)
  
  message("## Preprocess data ")
  
  marker <- list()

# R-L interaction ---------------------------------------------------------

  if(RL_interaction==T){
  
  #Create a new variable for connected cells
  top <- object@assays$NFCN$model$top.partner
  
  if(!is.null(correction.distance)){
    
    if(!is.null(object@assays$NFCN$model$top.partner$distance)){
      
      top <- object@assays$NFCN$model$top.partner %>% filter(distance<correction.distance)
      
    }else{ message("Distance is not quantified, please run the inferSpatialCorrection function ")}
  }
  
  
  object@meta.data$top <- "no"
  object@meta.data[rownames(object@meta.data) %in% c(top$Ligand, top$Receptor), ]$top <- "both"
  object@meta.data[rownames(object@meta.data) %in% c(top$Ligand), ]$top <- "Ligand"
  object@meta.data[rownames(object@meta.data) %in% c(top$Receptor), ]$top <- "Receptor"
  
  message("## Screen for marker genes: L-R ")
  Seurat::Idents(object) <-  object@meta.data$top
  marker.df.1 <- Seurat::FindAllMarkers(object, verbose=T)
  marker$R_L <- marker.df.1
  
  message("## Start with L-R interactions")
  
  set.seed(123)
  spotlight_ls <- SPOTlight::spotlight_deconvolution(se_sc = object,
                                                     counts_spatial = spata.obj %>% SPATA2::getCountMatrix(),
                                                     clust_vr = "top",
                                                     cluster_markers = marker.df.1,
                                                     cl_n = 100,
                                                     hvg = 3000,
                                                     ntop = NULL,
                                                     transf = "uv",
                                                     method = "nsNMF",
                                                     min_cont = 0.09)
  
  decon_mtrx <- 
    spotlight_ls[[2]] %>% 
    as.data.frame()
  
  names(decon_mtrx) <- paste0("top_inter_", names(decon_mtrx))
  
  decon_mtrx <- 
    decon_mtrx %>% 
    dplyr::mutate(barcodes=SPATA2::getBarcodes(spata.obj),
                  inter_score=c(top_inter_Ligand+top_inter_Receptor) %>% scales::rescale(c(0,1)))
  
  
  spata.obj <- spata.obj %>% SPATA2::addFeatures(feature_df = decon_mtrx, overwrite = T)
  
  
  
  }
  
  
# Connected - Non-Connected -----------------------------------------------

  
  message("## Screen for marker genes: connected-nonconnected ")
  
  get.upper <- 
    object@assays$NFCN$model$df.plot %>% 
    dplyr::mutate(score=x.p+y.p) %>% 
    dplyr::arrange((score)) %>% 
    utils::head(NN) %>% 
    dplyr::pull(barcodes)
   
  get.bottom <- 
    object@assays$NFCN$model$df.plot %>% 
    dplyr::mutate(score=x.p+y.p) %>% 
    dplyr::arrange((score)) %>% 
    utils::tail(NN) %>% 
    dplyr::pull(barcodes)
  
  object.con <- subset(object, cells=c(get.upper, get.bottom))
  object.con@meta.data$conNN <- "no"
  object.con@meta.data[rownames(object.con@meta.data) %in% get.upper, ]$conNN <- "connected"
  object.con@meta.data[rownames(object.con@meta.data) %in% get.bottom, ]$conNN <- "non_connected"

  
  
  Seurat::Idents(object.con) <-  object.con@meta.data$conNN
  marker.df.con <- Seurat::FindAllMarkers(object.con, verbose=T)
  marker$connected <- marker.df.con
  
  
  message("## Infer connected and non connected positions")
  
  set.seed(123)
  spotlight_ls <- SPOTlight::spotlight_deconvolution(se_sc = object.con,
                                                     counts_spatial = spata.obj %>% SPATA2::getCountMatrix(),
                                                     clust_vr = "conNN",
                                                     cluster_markers = marker.df.con,
                                                     cl_n = 100,
                                                     hvg = 3000,
                                                     ntop = NULL,
                                                     transf = "uv",
                                                     method = "nsNMF",
                                                     min_cont = 0.09)
  
  decon_mtrx <- 
    spotlight_ls[[2]] %>% 
    as.data.frame()
  
  names(decon_mtrx) <- paste0("top_inter_", names(decon_mtrx))
  
  decon_mtrx <- 
    decon_mtrx %>% 
    dplyr::mutate(barcodes=SPATA2::getBarcodes(spata.obj))
  
  
  spata.obj <- spata.obj %>% SPATA2::addFeatures(feature_df = decon_mtrx, overwrite = T)
  
  
  
  
  
  
  
  
  
  

# Group by Parameter ------------------------------------------------------

  
  
  if(intergrate.group.by==T){
  message("## Screen for marker genes: group.by ")
  
  Seurat::Idents(object) <-  object@meta.data[,group.by]
  marker.df.2 <- Seurat::FindAllMarkers(object, verbose=T)
  
  message("## Reference in spatial dataset: Using SPOTlight. Please cite: doi: https://doi.org/10.1093/nar/gkab043 ")
  message("## Start with group.by")
  set.seed(123)
  spotlight_ls <- SPOTlight::spotlight_deconvolution(se_sc = object,
                                                     counts_spatial = spata.obj %>% SPATA2::getCountMatrix(),
                                                     clust_vr = group.by,
                                                     cluster_markers = marker.df.2,
                                                     cl_n = 100,
                                                     hvg = 3000,
                                                     ntop = NULL,
                                                     transf = "uv",
                                                     method = "nsNMF",
                                                     min_cont = 0.09)
  
  decon_mtrx <- 
    spotlight_ls[[2]] %>% 
    as.data.frame()
  
  names(decon_mtrx) <- paste0(group.by,"_", names(decon_mtrx))
  
  decon_mtrx <- 
    decon_mtrx %>% 
    dplyr::mutate(barcodes=SPATA2::getBarcodes(spata.obj))
  
  
  spata.obj <- spata.obj %>% SPATA2::addFeatures(feature_df = decon_mtrx, overwrite = T)
  
  }
  
  
  
  

# Output ------------------------------------------------------------------

  
  
  message("## Prepare output ")
  out <- list(object, spata.obj, marker)
  names(out) <- list("object", "spata.obj", "marker")
  return(out)
  
  message("## Done ")
}


#' @title inferSpatialCorrection
#' @author Dieter Henrik Heiland
#' @description inferSpatialCorrection
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

inferSpatialCorrection<- function(object, inferSpatialPosition.object, assay="SCT", estimator=1, nr.genes=20){
  
  
  # The most estimated cell paires need to be corrected by the likelihood of physical interaction
  
  cell.paires <- object@assays$NFCN$model$top.partner
  
  #Find most matched position
  
  message("## 1/3  Estimation of juxtaposition probability")
  
  p <- dplyr::progress_estimated(nrow(cell.paires))
  
  map.dist <- purrr::map(.x=1:nrow(cell.paires),
                         .f=function(i){
  
  p$tick()$print()
  
  Ligand.sig <- 
    object@assays[[assay]]@scale.data[, cell.paires[i,1]] %>% 
    as.data.frame() %>% 
    dplyr::rename("gene":=.) %>% 
    dplyr::arrange(desc(gene)) %>% 
    head(nr.genes) %>% 
    rownames()
  Receptor.sig <- 
    object@assays[[assay]]@scale.data[, cell.paires[i,2]] %>% 
    as.data.frame() %>% 
    dplyr::rename("gene":=.) %>% 
    dplyr::arrange(desc(gene)) %>% 
    head(nr.genes) %>% 
    rownames()
  
  Ligand.pos <- 
    inferSpatialPosition.object[[2]] %>% 
    SPATA2::joinWithGenes(genes=Ligand.sig, average_genes = T, verbose = F) %>% 
    dplyr::rename("Ligant":=mean_genes) %>% 
    dplyr::arrange(desc(Ligant))
  
  Receptor.pos <- 
    inferSpatialPosition.object[[2]] %>% 
    SPATA2::joinWithGenes(genes=Receptor.sig, average_genes = T, verbose = F)%>% 
    dplyr::rename("Receptor":=mean_genes) %>% 
    dplyr::arrange(desc(Receptor))
  
  relative.dist <- NFCN2::getDistance(Ligand.pos, Receptor.pos)[1:estimator] %>% mean()
  
  return(relative.dist)
                         }) %>% unlist()
  
  cell.paires$distance <- map.dist
  
  
  cell.paires <-
    object@assays$NFCN$model$receptor.top %>% 
    dplyr::select(barcodes, sum.f) %>% 
    dplyr::rename("Receptor":=barcodes) %>% 
    dplyr::left_join(cell.paires, ., by="Receptor") %>% 
    dplyr::rename("Receptor.f":=sum.f)
  
  cell.paires <-
    object@assays$NFCN$model$ligand.top %>% 
    dplyr::select(barcodes, sum.f) %>% 
    dplyr::rename("Ligand":=barcodes) %>% 
    dplyr::left_join(cell.paires, ., by="Ligand") %>% 
    dplyr::rename("Ligand.f":=sum.f)
  
  #ggplot(data=cell.paires, mapping=aes(x=distance, y=Ligand.f+Receptor.f))+geom_point()+theme_classic()+geom_smooth(se=F)
  
  
  object@assays$NFCN$model$top.partner <- cell.paires
  
  return(object)
  
}









