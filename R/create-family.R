#' @title getPoints
#' @author Dieter Henrik Heiland
#' @description getPoints
#' @inherit 
#' @param coordinates A data.frame will rownames as identifier and coordinates names as x and y
#' @return 
#' @examples 
#' 
#' @export

createSlotNFCN <- function(object, source, target, group.by="seurat_clusters"){
  
  object@assays$NFCN <- NULL  
  object@assays$NFCN$target <- NULL 
  object@assays$NFCN$target$cells <- object@meta.data %>% dplyr::filter(!!sym(group.by) %in% target) %>% rownames()
  
  object@assays$NFCN$source <- NULL 
  object@assays$NFCN$source$cells <- object@meta.data %>% dplyr::filter(!!sym(group.by) %in% source) %>% rownames()
  
  return(object)
  
}