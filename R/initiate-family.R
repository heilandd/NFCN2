#' @title initiate NFCN2 
#' @author NFCN2
#' @description initiate NFCN2 
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
#' 


initiate.NFCN2 <- function(object){
  
  input.data <- NULL
  
  if(class(object)=="spata") input.data <- "spata"
  if(class(object)=="seurat") input.data <- "seurat"
  
  
  # check if the NFCN2 slot is free
  
  return(list(input.data <- input.data))

  
}