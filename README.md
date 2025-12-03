# Pinfish-omics
This repository contains scripts used to analyzing transcriptomic, metabolomic, and DNA methylation data of the pinfish (_Lagodon rhomboides_).

## Transcriptomics
Tag-seq at UT Austin was performed with extracted RNA from muscle tissue samples. 


## Metabolomics
Untargeted metabolomics by Creative Proteomics was performed on flash frozen muscle tissue samples. returned an excel file for positive and negative ion modes. These files were analyzed separately. These files also already contained annotated metabolites. We then performed statistical analysis to find significantly differentially expression metabolites. 

Metaboanalyst R package was used for analysis. Downloading_Metaboanalyst.R script was used to download this package. Many package dependencies were necessary. Then, we performed statistical analysis on negative ion mode and positive ion mode annotated metabolites using neg_stat_analysis.R and pos_stat_analysis.R scripts. Functional and enrichment analysis was performed manually. 

## DNA methylation
Whole Genome Bisulfite Sequencing (WGBS) was performed on extracted DNA from muscle tissue samples. We followed to CpG_Me pipeline for pre-processing analysis (citation). We then used methylKit package in R for downstream analysis.
