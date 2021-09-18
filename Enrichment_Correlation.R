#Function that calculates gene signature enrichment and determines significant correlation 
#Need to have loaded the P_value_mc and random_distribution functions to work
#P-value_mc and random_distribution function are found in Monte_Carlo.R
enrichment_correlation <- function(exprsset,enrichment_genes,query) {
  
  cell_enrichment <- spec::cellSubsetEnrich(exprsset,subsets = enrichment_genes)
  
  query_enrichment <- signatureEnrichment(exprsset,geneSet = query)
  
  SPEC_cor_matrix <- SPEC(exprsset, query, cell_enrichment)
  
  #generation of correlation values between enrichment scores and randomly chosen genes,looped 1000 times
  random_SPECs  <- random_distribution(exprsset,1000,enrichment_genes,cell_enrichment,query)
  
  p_values <- 
    sapply(names(enrichment_genes),FUN = function(x) { P_value_mc(SPEC_cor_matrix,random_SPECs,celltype = x) })
  
  significant_correlations <- SPEC_cor_matrix %>% cbind(.,data.frame(p_values))
  
  output <- list(cell_enrichemnt_scores = cell_enrichment,
                 gene_signature_scores = query_enrichment,
                 correlation_matrix = SPEC_cor_matrix,
                 sig_cor = significant_correlations)
  
}