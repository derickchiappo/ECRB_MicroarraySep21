#' @title emperical_p
#' 
#' @description calculates an emperical p-value
#' 
#' @param cormatrix correlation matrix between gene signature of interest and that of the cell type(s) of interest
#' @param randmatrices 
#' correlation matrix between the randomly generated gene signature of interest and gene signature of the cell type(s) of interest
#' @param celltype names of the cel type(s) of interest

#calculates p-value of correlations by cell-type in the correlation matrix generated 
emperical_p <- function(cormatrix,randmatrices,celltype) { 
  
  comparisons <- table(randmatrices[as.character(celltype),] >= cormatrix[as.character(celltype),] )
  
  p <-  ( comparisons[["TRUE"]] + 1 ) / (ncol(randmatrices) + 1)
  
  return(p)  
  
}


}