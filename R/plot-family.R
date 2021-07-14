#' @title plotNFCN
#' @author Dieter Henrik Heiland
#' @description plotNFCN
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

plotNFCN <- function(object,
                     group.by="seurat_clusters",
                     pt.size=0.2,
                     pt.alpha=0.5,
                     color.top=F,
                     plot.connection=T){
  
  
  # Get data
  plot.df <- object@assays$NFCN$model$df.plot
  
  
  #Create plot
  
  x.m <- max(plot.df$x.p)
  y.m <- max(plot.df$y.p)
  
  x.i <- min(plot.df$x.p)
  y.i <- min(plot.df$y.p)
  
  
  p <- 
  ggplot2::ggplot()+
    theme_void()+
    geom_point(data=plot.df, mapping=aes(x.p,y.p,color=!!sym(group.by)), size=pt.size+plot.df$top, alpha=pt.alpha)
    #geom_point(data=plot.df, mapping=aes(x.p,y.p,color=y), size=pt.size, alpha=pt.alpha)
  
  p <- p+
    geom_segment(mapping=aes(x = x.i, y = y.i, xend = x.i+0.2, yend = y.i), color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"))+
    geom_segment(mapping=aes(x = x.m, y = y.m, xend = x.m-0.2, yend = y.m), color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"))+
    geom_segment(mapping=aes(x = x.i, y = y.i, xend = x.i, yend = y.i+0.2), color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"))+
    geom_segment(mapping=aes(x = x.m, y = y.m, xend = x.m, yend = y.m-0.2), color="black",arrow = arrow(length = unit(0.2, "cm"), type = "closed"))

  
  if(plot.connection==T){
    
    segment.data <- data.frame(x=object@assays$NFCN$model$receptor.top$x.p,
               y=object@assays$NFCN$model$receptor.top$y.p,
               xend=object@assays$NFCN$model$ligand.top$x.p,
               yend=object@assays$NFCN$model$ligand.top$y.p)
    
    p <- p+
      geom_segment(data=segment.data, mapping=aes(x=xend,y=yend,xend=x,yend=y),size=seq(from=0.1, to=0.001, length.out = nrow(segment.data)), color="black", arrow = arrow(length = unit(0.1, "cm"), type = "closed"))
    
    
  }
  
  #p <- p+
  #  geom_segment(mapping=aes(x = 1, y = 1.1, xend = x.m, yend = 1.1), color="lightgrey", linetype=2, alpha=0.5)+
  #  geom_segment(mapping=aes(x = 1, y = 0.9, xend = x.i, yend = 0.9), color="lightgrey", linetype=2, alpha=0.5)+
  #  geom_segment(mapping=aes(x = 1, y = 1.1, xend = 1, yend = y.m), color="lightgrey", linetype=2, alpha=0.5)+
  #  geom_segment(mapping=aes(x = 1, y = 0.9, xend = 1, yend = y.i), color="lightgrey", linetype=2, alpha=0.5)
  
  
  p <- p+
    geom_text(mapping=aes(x = x.i+0.1, y = y.i-0.1, label="Expression Ligand"), color="black")+
    geom_text(mapping=aes(x = x.m-0.1, y = y.m+0.1, label="Expression Receptor"), color="black")+
    geom_text(mapping=aes(x = x.i-0.1, y = y.i+0.1, label="Score Induction"), color="black",angle=90)+
    geom_text(mapping=aes(x = x.m+0.1, y = y.m-0.1, label="Score Activation"), color="black",angle=90)
  
  p <- p+theme(legend.position="none")
    
  return(p)
    

}




#' @title plotNFCNDimRed
#' @author Dieter Henrik Heiland
#' @description plotNFCNDimRed
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

plotNFCNDimRed <- function(object,
                           reduction="umap",
                           group.by="seurat_clusters",
                           pt.size=0.2,
                           pt.alpha=0.5,
                           plot.connection=T){
  
  
  # Get data
  corrdinates <- 
    object@reductions[[reduction]]@cell.embeddings %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("barcodes") %>% 
    dplyr::left_join(.,object@meta.data %>% 
                as.data.frame() %>% 
                tibble::rownames_to_column("barcodes") %>% 
                dplyr::select(barcodes, {{group.by}}),
              by="barcodes")
  
  p <- 
    ggplot2::ggplot()+
    ggplot2::geom_point(data=corrdinates, mapping=aes(x=UMAP_1, y=UMAP_2, color=!!sym(group.by)), size=pt.size, alpha=pt.alpha)+
    theme_void()+
    theme(legend.position="none")
  
  # Add connections
  if(plot.connection==T){
    
    r.umap <- corrdinates %>% dplyr::filter(barcodes %in% object@assays$NFCN$model$receptor.top$barcodes)
    l.umap <- corrdinates %>% dplyr::filter(barcodes %in% object@assays$NFCN$model$ligand.top$barcodes)
    
    segment.data <- data.frame(x=l.umap$UMAP_1,
                               y=l.umap$UMAP_2,
                               xend=r.umap$UMAP_1,
                               yend=r.umap$UMAP_2)
    
    p <- p+
      geom_segment(data=segment.data, mapping=aes(x=x,y=y,xend=xend,yend=yend, color=l.umap[,group.by] ), alpha=pt.alpha, size=(seq(from=0.4, to=0.000001, length.out = nrow(segment.data))), arrow = arrow(length = unit(0.15, "cm"), type = "closed"))
    
    
  }
  
  
  return(p)
  
  
}


#' @title plotNFCNDimRed
#' @author Dieter Henrik Heiland
#' @description plotNFCNDimRed
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

plotNFCNDimRed3D <- function(object,
                           reduction="umap",
                           group.by="seurat_clusters",
                           pt.size=0.2,
                           pt.alpha=1,
                           l.width=0.1,
                           l.alpha=0.5,
                           l.color="lightgrey"){
  
  
  # Get data
  corrdinates <- 
    object@reductions[[reduction]]@cell.embeddings %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("barcodes") %>% 
    dplyr::left_join(.,object@meta.data %>% 
                       as.data.frame() %>% 
                       tibble::rownames_to_column("barcodes") %>% 
                       dplyr::select(barcodes, {{group.by}}),
                     by="barcodes")
  
  corrdinates3d <- rbind(
    corrdinates %>% dplyr::mutate(z=0),
    corrdinates %>% dplyr::mutate(z=1)
  )
  color.ext <- corrdinates3d[,group.by]
  
  p <- plotly::plot_ly() 
  
  p <- p %>% plotly::add_trace(data=corrdinates3d, 
                                x= ~UMAP_1, 
                                y= ~UMAP_2,
                                z= ~z, 
                                color = color.ext,
                                marker=list(#colorscale = "Plasma",
                                  reversescale =T,
                                  size=pt.size,
                                  alpha=pt.alpha),
                                type = "scatter3d") 
  
  
  
  p <- p %>% plotly::layout(   scene = list(xaxis = list(title = "UMAP 1"), 
                                    yaxis = list(title = "UMAP 2"),
                                    zaxis = list(title = " Ligand-Receptor ", 
                                           range = c(-1,2) 
                                           ) 
                                    )
  )
  
  r.umap <- corrdinates %>% dplyr::filter(barcodes %in% object@assays$NFCN$model$receptor.top$barcodes)
  l.umap <- corrdinates %>% dplyr::filter(barcodes %in% object@assays$NFCN$model$ligand.top$barcodes)
  
  segment.data.all <- data.frame(x=l.umap$UMAP_1,
                             y=l.umap$UMAP_2,
                             xend=r.umap$UMAP_1,
                             yend=r.umap$UMAP_2,
                             z=1,
                             zend=0)

  for(i in 1:nrow(segment.data.all)){
    segment.data <- segment.data.all[i, ]
    p <- p %>% plotly::add_paths(
      x = c(segment.data$x,segment.data$xend), 
      y = c(segment.data$y,segment.data$yend), 
      z = c(segment.data$z, segment.data$zend), 
      line = list(width = l.width, alpha=l.alpha, color = l.color, reverscale = FALSE))
  }
  
  return(p)
  
  
}


#' @title plotNFCNDimRed
#' @author Dieter Henrik Heiland
#' @description plotNFCNDimRed
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

plotNFCNSpatial3D <- function(object,
                              ls.objects,
                              reduction="umap",
                              group.by="seurat_clusters",
                              pt.size=0.2,
                              pt.alpha=1,
                              l.width=0.1,
                              l.alpha=0.5){
  
  
  spata.obj <- ls.objects[[2]]
  # Get data Seurat
  corrdinates.sc <- 
    object@reductions[[reduction]]@cell.embeddings %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("barcodes") %>% 
    dplyr::left_join(.,object@meta.data %>% 
                       as.data.frame() %>% 
                       tibble::rownames_to_column("barcodes") %>% 
                       dplyr::select(barcodes, {{group.by}}),
                     by="barcodes")
  
  corrdinates.spatial <- 
    SPATA2::joinWithFeatures(spata.obj, features = group.by) %>% 
    as.data.frame() %>% 
    dplyr::select(barcodes, x,y,{{group.by}}) %>% 
    dplyr::mutate(x=scales::rescale(x, c(min(corrdinates.sc$UMAP_1), max(corrdinates.sc$UMAP_1))),
                  y=scales::rescale(y, c(min(corrdinates.sc$UMAP_2), max(corrdinates.sc$UMAP_2)))
                  )

  names(corrdinates.spatial) <- names(corrdinates.sc)
  
  corrdinates.sc <- 
  corrdinates.sc %>% dplyr::mutate(z=1)
  
  corrdinates.spatial <- 
  corrdinates.spatial %>% dplyr::mutate(z=0)
  

  
  
  
  
  
  
  p <- plotly::plot_ly() 
  
  p <- p %>% plotly::add_trace(data=corrdinates.sc, 
                               x= ~UMAP_1, 
                               y= ~UMAP_2,
                               z= ~z, 
                               color = corrdinates.sc[,group.by],
                               marker=list(#colorscale = "Plasma",
                                 reversescale =T,
                                 size=pt.size,
                                 alpha=pt.alpha),
                               type = "scatter3d") 
  
  p <- p %>% plotly::add_trace(data=corrdinates.spatial, 
                               x= ~UMAP_1, 
                               y= ~UMAP_2,
                               z= ~z, 
                               color = corrdinates.spatial[,group.by],
                               marker=list(#colorscale = "Plasma",
                                 reversescale =T,
                                 size=pt.size+1,
                                 alpha=pt.alpha),
                               type = "scatter3d") 
  
  
  
  p <- p %>% plotly::layout(   scene = list(xaxis = list(title = "UMAP 1"), 
                                            yaxis = list(title = "UMAP 2"),
                                            zaxis = list(title = " Spatial-SingleCell ", 
                                                         range = c(-1,2) 
                                            ) 
  )
  )
  
  
  #Draw the connection into space
  
  
  
  r.umap <- corrdinates.sc %>% dplyr::filter(barcodes %in% object@assays$NFCN$model$receptor.top$barcodes)
  l.umap <- corrdinates.sc %>% dplyr::filter(barcodes %in% object@assays$NFCN$model$ligand.top$barcodes)
  spatial.pos <- 
    SPATA2::joinWithFeatures(spata.obj, features="top_inter_connected") %>% 
    dplyr::mutate(x=scales::rescale(x, c(min(corrdinates.sc$UMAP_1), max(corrdinates.sc$UMAP_1))),
                  y=scales::rescale(y, c(min(corrdinates.sc$UMAP_2), max(corrdinates.sc$UMAP_2)))
    ) %>% 
    dplyr::arrange(desc(top_inter_connected)) %>% 
    utils::head(nrow(l.umap))
  
  
  #Ligand
  
  segment.data.all <- data.frame(x=l.umap$UMAP_1,
                                 y=l.umap$UMAP_2,
                                 xend=spatial.pos$x,
                                 yend=spatial.pos$y,
                                 z=1,
                                 zend=0)
  
  
  for(i in 1:nrow(segment.data.all)){
    segment.data <- segment.data.all[i, ]
    p <- p %>% plotly::add_paths(
      x = c(segment.data$x,segment.data$xend), 
      y = c(segment.data$y,segment.data$yend), 
      z = c(segment.data$z, segment.data$zend), 
      line = list(width = l.width, alpha=l.alpha, color = "green", reverscale = FALSE))
  }
  
  #Receptor
  
  segment.data.all <- data.frame(x=r.umap$UMAP_1,
                                 y=r.umap$UMAP_2,
                                 xend=spatial.pos$x,
                                 yend=spatial.pos$y,
                                 z=1,
                                 zend=0)
  
  
  for(i in 1:nrow(segment.data.all)){
    segment.data <- segment.data.all[i, ]
    p <- p %>% plotly::add_paths(
      x = c(segment.data$x,segment.data$xend), 
      y = c(segment.data$y,segment.data$yend), 
      z = c(segment.data$z, segment.data$zend), 
      line = list(width = l.width, alpha=l.alpha, color = "red", reverscale = FALSE))
  }
  
  
  return(p)
  
  
}




#' @title plotNFCNDotplot
#' @author Dieter Henrik Heiland
#' @description plotNFCNDotplot
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

plotNFCNDotplot <- function(object,
                           group.by="seurat_clusters"){
  
  r.stat <- object@assays$NFCN$model$receptor.top %>% dplyr::count(!!sym(group.by))
  l.stat <- object@assays$NFCN$model$ligand.top %>% dplyr::count(!!sym(group.by))

  object@assays$NFCN$model$top.partner$rank <- 1:nrow(object@assays$NFCN$model$top.partner)
  partner <- object@assays$NFCN$model$top.partner
  partner <- 
    partner %>% 
    dplyr::mutate(ligand.p= object@meta.data[Ligand,{{group.by}}],
                  receptor.p= object@meta.data[Receptor,{{group.by}}]) %>% 
    dplyr::count(ligand.p,receptor.p) %>% 
    dplyr::mutate(n=scales::rescale(n, c(0,1) )) %>% 
    reshape2::acast(ligand.p~receptor.p)
  partner[is.na(partner)] <- 0
  
  corrplot::corrplot(partner %>% t(), is.corr = F, col=viridis::viridis(50))
  
  
}


#' @title plotFunctionalDotplot
#' @author Dieter Henrik Heiland
#' @description plotFunctionalDotplot
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

plotFunctionalDotplot <- function(functional,top=30){
  
  functional %>% clusterProfiler::dotplot(showCategory=top)+
    ggplot2::theme_classic()+
    ggplot2::theme(legend.position="none")+
    SPATA2::scale_color_add_on(clrsp="SunsetDark")
  
}






