#!/bin/bash
#PBS -N fastqc
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=2:mem=16gb
#PBS -q medium
#PBS -o fastqc.out
#PBS -e fastqc.err
#PBS -M alb0293@auburn.edu
#PBS -m abe

#load in fastqc
module load fastqc/0.12.1

#path to fastq files
trimmed_dir="/home/aubaxb001/methylation_data/trimmed"
cd "$trimmed_dir"

#make directory for results
mkdir -p "${trimmed_dir}/fastqc_results"

#since theyre all in the same directory, have to grab prefixes for R1 and R2 files
#loop through R1 and R2 files and run fastqc on both
for fastq1 in *_1.fq.gz; do
	prefix="${fastq1%_1.fq.gz}"
	fastq2="${prefix}_2.fq.gz"

	if [[ ! -f "$fastq2" ]]; then
		echo "Warning: $fastq2 not found, skipping sample $prefix"
		continue
	fi

	echo "Running FastQ Screen on sample ${prefix}"
	echo "Read1: ${fastq1}"
	echo "Read2: ${fastq2}"

#run fastqc
fastqc -o "${trimmed_dir}/fastqc_results" "${fastq1}" "${fastq2}"

done
