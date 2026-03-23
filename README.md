# Pinfish-omics
This repository contains all of the scripts used to analyze transcriptomic, metabolomic, and DNA methylation data of the pinfish (_Lagodon rhomboides_).

## Transcriptomics
TagSeq at UT Austin was performed with extracted RNA from muscle tissue samples. See README within the Transcriptomics folder for more info.

## Metabolomics
Untargeted metabolomics by Creative Proteomics was performed on flash frozen muscle tissue samples. Creative Proteomics returned one Excel file for positive and negative ion modes, respectively, with annotated metabolites. We then performed statistical analysis to find differentially expressed metabolites. See the README within the Metabolomics folder for more info.

## DNA methylation
Whole Genome Bisulfite Sequencing (WGBS) by Novogene was performed on extracted DNA from flash frozen muscle tissue samples. We followed to CpG_Me pipeline for pre-processing (https://github.com/ben-laufer/CpG_Me), and the methylKit package in R for downstream analysis.

## Buoy data
To determine ecological significance of temperatures, data from the National Data Buoy Center was downloaded and analyzed. A guide was provided on how to download these data (GuidetoNDCBdata.docx) and the script that was used to determine summary statistics is also included (buoy_temp_data_Chp1.R).
