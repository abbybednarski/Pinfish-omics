# Pipeline for analyzing WGBS data from Novogene
1. First, follow the GuidetoCpG_Me.docx file.
   - This will take you thorugh all of the necessary pre-processing steps (e.g., trimming, aligning, mapping).
   - These were run on the Alabama supercomputer (ASC) since the files were so big!
2. Next, you can run methylKit analysis. This will use the DMR_*.R scripts in this repository
   - DMR_MvS.R, DMRs_CvM.R, DMRs_CvS.R are all running a pipeline for methylKit analysis but for the different pairwise comparisons in my analysis. The MvS script is for Multiple v. Single, the CvM script is for Multiple v. Control, and the CvS script is for Single v. Control.
   - These do not need to be run in a specific order
   - Since these scripts are still using big input files, they were run on the ASC. I created individual scripts to submit each individual pairwise R script to the ASC. The bash script was submitted to ASC to run the R script through ASC. The script combos are as follows:
        - DMRs_CvS.R and methylation_analysis.sh
        - DMRs_CvM.R and methylation_analysis_CvM.sh
        - DMR_MvS.R and MvS_methylation_analysis.sh
   - Once you have completed differential methylation analysis, files can be transferred to your local computer and analyzed further
   - The R scripts also contain code for manhattan plots, violin plots, ECDFs, methylation summary statistics, and more.
