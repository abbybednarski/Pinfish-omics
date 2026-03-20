# Pipeline for analyzing untargeted metabolomics data from Creative Proteomics
1. First, follow the Downloading_Metaboanalyst.R script to be able to download the software
   - This script contains all of the package dependencies that are needed before you can install MetaboAnalyst.
2. Next, follow neg_stat_analysis.R and/or pos_stat_analysis.R for data analysis
    - Untargeted metabolomics utilized UPLC MS/MS which is run in negative ion mode and positive ion mode. These datasets need to be analyzed separately since they have their own internal standards.
    - These scripts are essentially the same thing except the different input files that came from Creative Proteomics.
    - Note it does not matter if you start with positive or negative data analysis - the process is the same regardless of which ion mode the mass spec was run in.
