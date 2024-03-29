#loading of packages 
install.packages(c("BiocManager","pheatmap","ggpubr","ggrepel","RSQLite","tidyverse"))

library(pheatmap)
library(ggpubr)
library(ggrepel)
library(RSQLite)
library(tidyverse)
library(BiocManager)
BiocManager::install(c("affy","limma","annotate",
                     "AnnotationDbi","GEOquery",
                    "org.Hs.eg.db","hgu133plus2.db"))
library(affy)
library(limma)
library(annotate)
library(AnnotationDbi)
library(GEOquery)
library(org.Hs.eg.db)
library(hgu133plus2.db)

#Graph variables, and aesthetics and functions for ggplot
Graph_Settings <- 
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 14, color = "black",face = "bold"), 
        axis.title = element_text(size = 14, color = "black", face = "bold"))
    
    
#signatureEnrichment and GSEA function of the SPEC package needs to be loaded seperately
install.packages("spec_0.5.0.tar.gz")
library(spec)
source("signatureEnrichment.R")
source("GSEA.R")
source("Monte_Carlo.R")
source("Enrichment_Correlation.R")


#platform custom CDF brain Array
GSE104948 <- getGEO("GSE104948",AnnotGPL = T)[[1]]


#Customization of Expression set
pData(GSE104948)$'diagnosis:ch1' <- 
str_replace(string = pData(GSE104948)$'diagnosis:ch1', pattern = "ANCA Associated Vasculitis", 
            replacement = "AAV")

pData(GSE104948)$'diagnosis:ch1' <- 
  str_replace(string = pData(GSE104948)$'diagnosis:ch1', pattern = "Diabetic Nephropathy", 
              replacement = "DN")

pData(GSE104948)$'diagnosis:ch1' <- 
  str_replace(string = pData(GSE104948)$'diagnosis:ch1', pattern = "Focal Segmental Glomerular Sclerosis", 
              replacement = "FSGS")

pData(GSE104948)$'diagnosis:ch1' <- 
  str_replace(string = pData(GSE104948)$'diagnosis:ch1', pattern = "FSGS/Minimal Change Disease", 
              replacement = "FSGS/MCD")

pData(GSE104948)$'diagnosis:ch1' <- 
  str_replace(string = pData(GSE104948)$'diagnosis:ch1', pattern = "Minimal Change Disease", 
              replacement = "MCD")

pData(GSE104948)$'diagnosis:ch1' <- 
  str_replace_na(string = pData(GSE104948)$'diagnosis:ch1',
                 replacement = "HC")

pData(GSE104948)$'diagnosis:ch1' <- 
  str_replace(string = pData(GSE104948)$'diagnosis:ch1', pattern = "Tumor Nephrectomy", 
              replacement = "TN")

exprs(GSE104948) <- normalizeBetweenArrays(exprs(GSE104948))

GSM_AAV_HC <- 
pData(GSE104948) %>% 
filter(.$'diagnosis:ch1' == "AAV" | .$'diagnosis:ch1' == "HC") %>% 
rownames(.)

GSE104948_AAV <- GSE104948[,GSM_AAV_HC]

info <- pData(GSE104948_AAV) %>% dplyr::select(.,"diagnosis:ch1")

corMatrix <- cor(exprs(GSE104948_AAV),use = "c")

pheatmap(corMatrix,annotation_col = info)

#Selection of genes
Axis_goi<- fData(GSE104948) %>% 
           filter(Symbol %in% c("S100A8","S100A9","AGER","TLR4","S100A12","HMGB1"))

RAGE_goi <- fData(GSE104948) %>% 
            filter(Symbol %in% c("IL10","AGER","DIAPH1","IRF7","ARG1","NOS2","CCL2",
                                 "IL1B","VCAM1","ICAM1","TNF","RELA"))

disease_glom <- pData(GSE104948_AAV)[,"diagnosis:ch1"] %>% as.factor(.)

#limma pipeline
design_glom <- model.matrix(~0+disease_glom,data = GSE104948_AAV)
colnames(design_glom) <- levels(disease_glom)
contrastmatrix_glom <- makeContrasts(AAVvsHC = AAV - HC,
                                levels = design_glom)

fit_glom <- lmFit(GSE104948_AAV,design_glom)
fit2_glom <- contrasts.fit(fit_glom,contrasts = contrastmatrix_glom)
fit2_glom <- eBayes(fit2_glom)
pvaluehistogram_glom <- hist(fit2_glom$p.value)
results_glom <- decideTests(fit2_glom,adjust.method = "BH", p.value = 0.05)
summary(results_glom)

resultstable_AAVvsHC_glom <- topTable(fit2_glom,number = nrow(fit2_glom),coef = "AAVvsHC")

sig_goi_AAVvsHC_glom <- resultstable_AAVvsHC_glom[c(Axis_goi$ID,RAGE_goi$ID),] %>% 
                        filter(logFC >= log2(1.5) | logFC <= log2(0.5), adj.P.Val <= 0.05) %>%
                        mutate(Group = ifelse(.$ID %in% Axis_goi$ID,"Axis","RAGE")) 
                        

#Glomerualr_SPEC_Analysis............................................................................
#Requirements: locally downloaded spec tar file (source: http://clip.med.yale.edu/SPEC.)

#list of genes to for cell enrichment (from Defining cell type specificity at the transcriptional level in human disease. Genome Research. 23:1862-1873. 2013),
#PBMC cell gene signatures from Immune response in silico (IRIS): immune-specific genes identified from a compendium of microarray expression data
#http://nano.princeton.edu/standard/view/HPRD--Neutrophil/

glom_genes <- read.table("HPRD_Glomerulus.txt", sep = "/")

glom_genes <- glom_genes %>% 
              as.character(.)

Neutrophils <- irisSignatures$IrisNeutro[!irisSignatures$IrisNeutro %in% c(sig_goi_AAVvsHC_glom_Axis$Symbol)]

glomerulus <- list(mesangial_cells = c("ITGA8","COL4A1", "SMAD1",  "HMOX1", "CD46", "ACVRL1", "GHR",
                                       "SCUBE2","FCAMR",  "KNG1", "PTGIR","AQP1"), 
                   podocytes = c("PODN","PTGS2", "DBN1",  "LRRC7",  "CD59","WTIP","KNG1", "PTGIR",
                                 "CXCR3",  "CXCR1", "KIRREL3"), 
                   glomerular_endothelial_cells = c("ITGA5", "CLDN5","VWF","PECAM1",  "PDGFB","TEK","THBD","THBD",
                                                    "PCDH12", "ICAM2",  "EHD3", "CD34","ENG","FLT1",  "EDN1","KDR","VWF"),
                   glomerlular_cells = glom_genes,
                   neutrophils = Neutrophils,
                   monocytes = irisSignatures$IrisMono[!irisSignatures$IrisMono %in% c(sig_goi_AAVvsHC_glom_Axis$Symbol,sig_goi_AAVvsHC_glom_RAGE$Symbol)],
                   B_cells = irisSignatures$IrisB ,
                   T_cells = irisSignatures$IrisT,
                   NK_cells = irisSignatures$IrisNK,
                   DC = irisSignatures$IrisDC,
                   palmer_granulocytes = palmerSignatures$PalmerGranulocyte,
                   palmer_Tcells = palmerSignatures$PalmerT,
                   palmer_Bcells = palmerSignatures$PalmerB
)

query_Axis =  c(sig_goi_AAVvsHC_glom_Axis$Symbol)

querySig_RAGE = c(sig_goi_AAVvsHC_glom_RAGE$Symbol)

querySig_Tcells = c(sig_goi_AAVvsHC_glom_Tcells$Symbol)


#AAV_glomerulus_Dataset
GSM_AAV <- pData(GSE104948) %>% 
           filter(.$'diagnosis:ch1' == "AAV" ) %>% 
           rownames(.)

GSE104948_AAV <- GSE104948[,GSM_AAV] %>% 
                 exprs(.) %>% 
                 data.frame(.)

rownames(GSE104948_AAV) <- fData(GSE104948)[,"Symbol"]

#LD_glomerulus_Dataset
GSM_HC <- pData(GSE104948) %>% 
          filter(.$'diagnosis:ch1' == "HC" ) %>% 
          rownames(.)

GSE104948_HC <- GSE104948[,GSM_HC] %>% 
                exprs(.) %>% 
                data.frame(.)

rownames(GSE104948_HC) <- fData(GSE104948)[,"Symbol"]


#Glom_AAV
glom_AAV <- enrichment_correlation(GSE104948_AAV,glomerulus,query_Axis) 

glom_AAV_RAGE <- enrichment_correlation(GSE104948_AAV,glomerulus,querySig_RAGE)


#Glom_HC
glom_HC <- enrichment_correlation(GSE104948_HC,glomerulus,query_Axis) 

glom_HC_RAGE <- enrichment_correlation(GSE104948_HC,glomerulus,querySig_RAGE)  

#Graphing of correlations 
Neutrophil_Axis_graph <- 
ggplot(,aes(x = glom_AAV$cell_enrichemnt_scores["neutrophils",],y = glom_AAV$gene_signature_scores)) + 

#First Dataset: AAV
geom_point(aes(color = "glom_AAV$cell_enrichemnt_scores"))  + 
geom_smooth(method = "lm",aes(color = "glom_AAV$cell_enrichemnt_scores"),se = FALSE) + 

#Second Dataset: HC
geom_point(aes(x = glom_HC$cell_enrichemnt_scores["neutrophils",], y = glom_HC$gene_signature_scores,
            color = "glom_HC$cell_enrichemnt_scores")) +
geom_smooth(method = "lm",aes(x = glom_HC$cell_enrichemnt_scores["neutrophils",], y = glom_HC$gene_signature_scores,
            color = "glom_HC$cell_enrichemnt_scores"), se = FALSE) +

#Graph customization
labs(color = "Legend") + 
scale_color_discrete(labels = c("AAV","HC")) + 
Graph_Settings + 
xlab("Cell Signature Enrichment Score") + 
ylab("S100A12/Calprotectin - RAGE/TLR4 axis enrichment score") +
geom_text(aes(x = 0.35, y = 0.75),label = paste0("HC: r =",round(glom_HC$sig_cor["neutrophils",1],digits = 2),",",
         "p =",round(glom_HC$sig_cor["neutrophils",2],digits = 2))) +
geom_text(aes(x = 0.35, y = 0.65),label = paste0("AAV: r =",round(glom_AAV$sig_cor["neutrophils",1],digits = 2),",",
         "p =",round(glom_AAV$sig_cor["neutrophils",2],digits = 2)))  + 
ggtitle("Neutrophil") + theme(plot.title = element_text(face = "bold", size =14, hjust = 0.5)) 



Monocyte_Axis_graph <-  
ggplot(,aes(x = glom_AAV$cell_enrichemnt_scores["monocytes",],y = glom_AAV$gene_signature_scores)) + 

#First Dataset: AAV
geom_point(aes(color = "AAV Monocytes"))  + 
geom_smooth(method = "lm",aes(color = "AAV Monocytes"),se = FALSE) + 

#Second Dataset: HC
geom_point(aes(x = glom_HC$cell_enrichemnt_scores["monocytes",], y = glom_HC$gene_signature_scores,
            color = "HC Monocytes")) +
geom_smooth(method = "lm",aes(x = glom_HC$cell_enrichemnt_scores["monocytes",], y = glom_HC$gene_signature_scores,
            color = "HC Monocytes"), se = FALSE) + 

#Graph customization
labs(color = "Legend") + 
scale_color_discrete(labels = c("AAV","HC")) + 
Graph_Settings + 
xlab("Cell Signature Enrichment Score") + 
ylab("S100A12/Calprotectin - RAGE/TLR4 axis enrichment score") +
geom_text(aes(x = 0.32, y = 0.3),label = paste0("HC: r =",round(glom_HC$sig_cor["monocytes",1],digits = 2),",",
        "p =",round(glom_HC$sig_cor["monocytes",2],digits = 2))) +
geom_text(aes(x = 0.32, y = 0.35),label = paste0("AAV: r =",round(glom_AAV$sig_cor["monocytes",1],digits = 2),",",
        "p =",round(glom_AAV$sig_cor["monocytes",2],digits = 3))) + 
ggtitle("Monocytes") + theme(plot.title = element_text(face = "bold", size =14,hjust = 0.5))                        
