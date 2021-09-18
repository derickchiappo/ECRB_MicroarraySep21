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
plot.mean <- function(x) { 
  m <- mean(x,na.rm = TRUE) 
  c(y = m, ymin = m, ymax =m) }

plot.median <- function(x) { 
  m <- median(x, na.rm = TRUE) 
  c(y = m, ymin = m, ymax =m) }

Graph_Settings <- 
  theme(panel.background = element_rect(fill = "white"),axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black", face = "bold"))

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
                        
AAVvsHC_volcano_plot <- 
  ggplot(sig_goi_AAVvsHC_glom,aes(x = 2^logFC, y = -log10(adj.P.Val),label = Symbol, color = Group)) + 
  geom_point() + 
  geom_hline(aes(yintercept = log10(0.05)),linetype = "dashed",size = 1) + 
  geom_vline(xintercept = 1.5,color = "black",size = 1,linetype = "dashed") +
  geom_label_repel(nudge_x = 0.1,nudge_y = -0.1,size = 3) + 

#Graph Customization
Graph_Settings + 
  theme(plot.title = element_text(face = "bold", size =14, hjust = 0.5)) + 
  scale_x_continuous(limits = c(-2,4),breaks = seq(-2,4,by = 1)) + 
  xlab("Fold Change")
                        
                        
                        
                        
