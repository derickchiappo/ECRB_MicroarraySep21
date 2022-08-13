#' @title random_distribution
#' 
#' @description Generates n-number of correlations between each cell signature and 
#' a random set of genes the same size as that of the gene query
#' 
#' @param exprsVals gene expression values (expressed as fold changes)
#' @param run_number (number of loops to run to generate gene signature of the same size)
#' as that of interest
#' @param celltypes name(s) of the cell types that the gene signatures correspond to
#' @param cell_enrichemnt gene signature of the cell(s) of interest
#' @param sig genes of interest
#' 
#' 
#' @return returns a correlation matrix
#' 
#' @export
random_distribution <- function(exprsVals,run_number,celltypes,cellenrichment,sig) {
  
  set.seed(120)
  
  matrices_names <- c()
  
  random_SPECs <- c()
  
  for(i in 1:run_number)  { 
    assign(x = paste0("Sym_matrix_",i),
           value = sample(x  = Symbols,size = length(sig)) )
    matrices_names[i] <- paste0("Sym_matrix_",i) }
  
  random_SPECs <- 
    sapply(matrices_names,FUN = function(x){SPEC(exprsVals,queries = get(x),cellenrichment) }) 
  
  rownames(random_SPECs) <- names(celltypes)
  
  return(random_SPECs)
}
  