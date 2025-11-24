#!/bin/bash
#PBS -N bismark_coverage
#PBS -o bismark_coverage.out
#PBS -e bismark_coverage.err
#PBS -l select=1:ncpus=5:mem=100gb
#PBS -l walltime=21:00:00
#PBS -q large
#PBS -M alb0293@auburn.edu
#PBS -m abe

#module load samtools/1.18

#module load fastq-screen/0.15.3
#module load bowtie2/2.4.2

#module load picard/2.26.2
#module load trimgalore/0.6.7

#load anaconda to be able to use cutadapt - cutadapt is included in anaconda so if anaconda is loaded, cutadapt will be too
#module load anaconda/3-2020.02

module load bismark

#load in R to generate histogram for insert size metrics step
#source /apps/profiles/modules_asax.sh.dyn
#module load R/4.4.0

cd /home/aubaxb001/methylation_data/trimmed

#run big script, one block at a time

#trim
#/home/aubaxb001/methyl_bias.sh trim
#echo "Trimming complete"

#align
#/home/aubaxb001/methyl_bias.sh align
#echo "Aligning complete"

#deduplicate
#/home/aubaxb001/methyl_bias.sh deduplicate
#echo "Deduplicate complete"

#insert
#/home/aubaxb001/methyl_bias.sh insert
#echo "Insert complete"

#coverage
#/home/aubaxb001/methyl_bias.sh coverage
#echo "coverage complete"

#extract
/home/aubaxb001/methyl_bias.sh extract
echo "extraction complete"

#cytosinereport with QC (bismark2report)
#/home/aubaxb001/methyl_bias.sh cytosineReport
#echo "Cytosine report and QC complete"
