# ECRB_MicroarraySep21
Analysis of publicly available microarray data from GEO involving differential expression analysis and correlation.

This piece of work was a part of my Ph.D. where I used publically avaliable microarray data from GEO database (https://www.ncbi.nlm.nih.gov/geo/) to determine if 
my pathway of interest was upregulated in patients with the rare systemic autoimmune disease, ANCA-associateed vascultis, was upregulated or not. sis 

This involved multiple steps including remotely accesing the raw data files from GEO,cleaning and processing the raw data files, and running the differental expression analysis using the limma package (https://www.bioconductor.org/packages/release/bioc/html/limma.html). 


I also performed a partial deconvulotion technique, via Bolen et al.'s SPEC package(http://clip.med.yale.edu/SPEC/), and custom function to relate my upregulated genes of interest to different types of immune cells. The gene signatures of the immune cells were obtained from Abbas et al.(https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0006098#s4). 

I changed the calculation of the p-value in the SPEC package to calculate an empirical p-value. Similarly, I generated a null distribution but unlike the original SPEC package I don't fit a normal distrubition to the generated null distribution. 


