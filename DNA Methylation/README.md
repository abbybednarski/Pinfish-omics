# Script order:
1. bismark_prep.sh
    a. This script is preparing the reference genome for bismark and generating bowtie2 indices. 
2. merge_lanes.sh
    a. This script merges sequences that were run across different lanes. This should generate 2 files per sample (forward and reverse files).
3. fastqc.sh
    a. This script quality checks the files before running any additional analysis.
4. methyl_bias.sh
    a. This script has the rest of the pre-processing steps in blocks. 
    b. submit_job.sh was used to submit one block at a time to the Alabama SuperComputer. 
