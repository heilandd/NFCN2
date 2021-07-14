#' @title findDouplets
#' @author Dieter Henrik Heiland
#' @description findDouplets
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export


findDouplets <- function(object, subset=F, return.object=T, doublet.rate=0.075){
  
  if(!class(object)=="Seurat") stop("object os not of class 'Seurat")
  
  message("The package 'DoubletFinder' is used please cite: https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30073-0 ")
  nExp_poi <- round(doublet.rate*nrow(object@meta.data))  
  object <- DoubletFinder::doubletFinder_v3(object, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  
  DF <- object@meta.data %>% tibble::rownames_to_column("barcodes") %>%  dplyr::select(barcodes,DF.classifications)
  
  
  
  if(subset==T){
    
  message(" The output will be a list of two Seurat obj, the first contained 'Singlet' and the second 'Doublet' cells ")
   out <- list(
     subset(object, cells=DF %>% dplyr::filter(DF.classifications=="Singlet") %>% pull(barcodes)),
     subset(object, cells=DF %>% dplyr::filter(DF.classifications=="Doublet") %>% pull(barcodes))
     )
    
  }else{
    message(" the classification will be in the variable  'DF.classifications' ")
    out <- object
    
  }
  
  if(return.object==F){
    
    message(" The object will not be exported from the function !!! ")
    out <- DF
    
  }
  
  
  return(out)
  
} 




#'@title findFunctionalEnrichment
#' @author Dieter Henrik Heiland
#' @description findFunctionalEnrichment
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

findFunctionalEnrichment <- function(marker,top=100){
  
  out <- 
  purrr::map(.x=unique(marker$cluster), .f=function(x){
    genes <- 
      marker %>% 
      dplyr::filter(cluster=={{x}}) %>% 
      dplyr::top_n(n = top, wt = avg_log2FC) %>% 
      dplyr::pull(gene)
  
    #fE
    gse <- clusterProfiler::enrichGO(gene=genes,
                                     ont ="ALL", 
                                     keyType = "SYMBOL", 
                                     OrgDb = "org.Hs.eg.db")
    })
  
  return(out)
  
}
