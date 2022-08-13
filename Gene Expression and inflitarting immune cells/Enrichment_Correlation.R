#'@Title: Enrichment_Correlation
#'
#'@description 
#'This function takes gene expression data and two gene signature of a specific size 
#'to calculate the correlation between the gene signature of interest and that of specific cell type. 
#'The significance of the correlation is determined by the calculation of an emperical p-value. 
#'The null condition is simulated by correlating a random selection of genes with the gene signature of each cell type
#'The dependencies include the SPEC function (See http://clip.med.yale.edu/SPEC/)
#'Two custom functions, random_distribution and empirical_p
#'
#'@param exprsset gene_expresion data calculated as fold changes
#'@param enrichment_genes vector of genes in the signature of the cell types(s) of interest
#'@param query vector of gene(s) you want to correlate with the enrichment genes
#'
#' @return a nested list with the:
#' 1) cell_enrichment_scores: gene enrichment score of the genes associated with each cell type
#' 2)gene_signature_scores:gene enrichment score of the gene(s) of interest in each sample
#' 3)correlation_matrix: reports the R-values of correlation
#' 4)sig_cor: list of correlation that reach statistical significance
#'
#'@export
enrichment_correlation <- function(exprsset,enrichment_genes,query) {
  
  cell_enrichment <- spec::cellSubsetEnrich(exprsset,subsets = enrichment_genes)
  
  query_enrichment <- signatureEnrichment(exprsset,geneSet = query)
  
  SPEC_cor_matrix <- SPEC(exprsset, query, cell_enrichment)
  
  #generation of correlation values between enrichment scores and randomly chosen genes,looped 1000 times
  random_SPECs  <- random_distribution(exprsset,1000,enrichment_genes,cell_enrichment,query)
  
  p_values <- 
    sapply(names(enrichment_genes),FUN = function(x) { empirical_p(SPEC_cor_matrix,random_SPECs,celltype = x) })
  
  significant_correlations <- SPEC_cor_matrix %>% cbind(.,data.frame(p_values))
  
  output <- list(cell_enrichemnt_scores = cell_enrichment,
                 gene_signature_scores = query_enrichment,
                 correlation_matrix = SPEC_cor_matrix,
                 sig_cor = significant_correlations)
  
}
