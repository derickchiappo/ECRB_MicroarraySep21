####################################################################
## This code implements the SPEC method as described in 
## C. Bolen et al., BMC Bioinformatics 2011
## Author - Christopher Bolen   


#' Signature Enrichment
#' @description find enrichment score of gene set
#' @param exprsVals - a matrix of expression values (genes x patients)
#' @param geneSet - a vector representing genes to calculate the enrichment for 
#'           can either be a list of gene names (matching the rownames of exprsVals)
#'           or a numeric vector representing the locations of the genes in exprsVals
#' @param dontInclude - genes to be removed from exprsVals before calculating enrichment
#' @param weighted - the weighting of the enrichment score calculation
#'            (refer to Subramanian 2005 for in depth discussion of this)
#' @param fast - if set to false, will do an in-depth search of the gene names to find matches
#'        (if a probe has multiple gene names, you can label it using "ABC1 /// ABC2" and 
#'         set fast=F to search for all instances of "ABC1")
#' @return numeric vector; the per-sample enrichment of the genes in the gene set
signatureEnrichment <- function(exprsVals, geneSet, dontInclude=NULL, weighted=1, fast=T){
   #find which of the genes from the gene set is in our array
  genesFromArray = rownames(exprsVals)
  takeOut = (genesFromArray %in% dontInclude)
  genesFromArray=genesFromArray[!takeOut]
  
  #just in case I forget elsewhere
  geneSet = geneSet[!is.na(geneSet)]
  

  if(is.numeric(geneSet)){
    locList = geneSet
  ##check for the " /// " syntax (that was common back on old Affy microarrays)
  }else if(!fast && any(grepl("///", genesFromArray))){
    genesFromArray = strsplit(genesFromArray, " /// ")
    allLocs = rep(1:length(genesFromArray), times=sapply(genesFromArray,length))
    genesFromArray = unlist(genesFromArray)
    
    locList = allLocs[which(genesFromArray %in% geneSet)]
  }else{
    locList = which(genesFromArray %in% geneSet)
  }
  if(length(locList)==0){
    warning("No genes in gene set match rows of expression matrix; returning NA\n")
    return(rep(NA, ncol(exprsVals)))
  }
  
  patients =  1:ncol(exprsVals)

  patientGenes = numeric(length(patients))
  for(p in patients){
    #ordered genes(probes)
  	ord = order(exprsVals[!takeOut,p],decreasing=T)
  	map = match(1:length(ord), ord)
  	sortedExpr = exprsVals[!takeOut,p][ord]

    if(!is.list(locList)){numSet=map[locList]
    }else{
      numSet = c()
      for(i in 1:length(locList))
        numSet[i] = min(map[locList[[i]]])
    }
    #patientGenes[p] = GSEA.EnrichmentScore2(1:length(sortedExpr), numSet, weighted, sortedExpr)[[1]]
    patientGenes[p] = GSEA.EnrichmentScore2(1:length(sortedExpr), numSet, weighted, (length(sortedExpr):1)/length(sortedExpr))[[1]]
  }
  patientGenes
}
