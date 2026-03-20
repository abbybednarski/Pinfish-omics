# Pipeline for analyzing Tag-Seq data from UT Austin GSAF 
### heavily based on https://github.com/z0on/tag-based_RNAseq
1. Begin by downloading data from BaseSpace and verify downloaded files with md5sum checks
    - You can do this by running md5sum *.gz and then comparing the numbers of your downloaded files and the files still on BaseSpace 
2. Concatenate sequence files for the same sample run on different lanes
    - You may not need to do this!
    - Tag-Seq generates 100-bp single-end reads (not paired-end!)
    - Each sample should have at least 1 fastq file containing raw read data. In some cases, you will have more than one file because sequencing facilities will split up a single sample, generting multiple files.
    - For example, I had C-T1-03-L_S3_L001_R1_001.fastq.gz and C-T1-03-L_S3_L002_R1_001.fastq.gz that needed to be concatenated to make one single file containing all reads for that sample
    - Start by generating a .txt file with all sample names, one per line (this is the sample-IDs.txt file in this repository)
    - Then run the following loop: <br/>
for sample in `cat sample-IDs.txt` <br/>
    do <br/>
cat ${sample}*.fastq.gz > ${sample}.fq.gz <br/>
    done <br/>
    - Move the concatenated fastq files to a new directory to stay organized!
2. Remove adapters with the run_cutadapt.sh script in this repository
    - First, you will need to create a new conda environment with the following code: <br/>
     module load python/anaconda <br/>
     export CONDA_PKGS_DIRS=~/.conda/pkgs <br/>
     conda create -n env1 <br/>
    - The script will activate this conda environment and run cutadapt on the concatenated files
4. Quality check with FastQC by loading and running FastQC on your trimmed files
    - Download the FastQC report files onto your local computer to inspect the trimmed data
    - You are hoping to see: limited number of sequences flagged, ~100 sequence length, above 30 quality score, no N's, all sequence lengths around 100, no adapters
6. Map to reference genome
    - Now that you have verified that you have good, high-quality trimmed sequences, you can map your reads to the genes in the pinfish transcriptome!
    - Using the bowtie.sh script in this repository, you will build a reference index, and map your reads per sample to the reference
    - Make sure to check the output file for mapping rate at this stage! We want the mapping rate to be no lower than 75%.
8. Generate count data
    - First, download samcount.pl and expression_compiler.pl and transfer them to your working directory.
    - Then run the following code if your transcriptome is in a file called Transcriptome-Sequences-Trimmed-Names.fa and it has been built by Trinity: <br/>
      grep ">" Transcriptome-Sequences-Trimmed-Names.fa > Transcript_IDs.txt <br/>
      cat Transcript_IDs.txt | sed 's/>\(TRINITY_DN[0-9]*_c[0-9]*_g[0-9]*\)\(_i[0-9]*\)/\1\2\t\1/' > transcriptome_seq2iso.tab <br/>
    - Following the samcount.sh script in this repository, generate count data per sample, making sure the *.pl and transcriptome_seq2iso.tab files are in your working directory.
    - You now have count data for each sample!
10. Run differential expression analysis with DESeq2 in R
