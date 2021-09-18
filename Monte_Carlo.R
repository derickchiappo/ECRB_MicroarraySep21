#calculates p-value of correlations by cell-type in the correlation matrix generated 
P_value_mc  <- function(cormatrix,randmatrices,celltype) { 
  
  comparisons <- table(randmatrices[as.character(celltype),] >= cormatrix[as.character(celltype),] )
  
  p <-  ( comparisons[["TRUE"]] + 1 ) / (ncol(randmatrices) + 1)
  
  return(p)  
  
}

#Generates n-number of correlations between each cell signature and 
#a random set of genes the same size as that of the gene query
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