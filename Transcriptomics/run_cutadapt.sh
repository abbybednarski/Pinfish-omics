#!/bin/bash
#SBATCH --job-name=cutadapt
#SBATCH --ntasks=1
#SBATCH --partition=mab0205_bg4
#SBATCH --mem=20G
#SBATCH --time=08:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alb0293@auburn.edu

module load python/anaconda

source activate env1

cutadapt -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o C-T1-03-M-trim.fq.gz C-T1-03-M.fastq.gz
cutadapt -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o C-T1-14-M-trim.fq.gz C-T1-14-M.fastq.gz
cutadapt -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o C-T1-16-M-trim.fq.gz C-T1-16-M.fastq.gz
cutadapt -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o C-T5-12-M-trim.fq.gz C-T5-12-M.fastq.gz
cutadapt -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o M-T2-09-M-trim.fq.gz M-T2-09-M.fastq.gz
cutadapt -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o M-T2-16-M-trim.fq.gz M-T2-16-M.fastq.gz
cutadapt -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o M-T6-06-M-trim.fq.gz M-T6-06-M.fastq.gz
cutadapt -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o M-T6-10-M-trim.fq.gz M-T6-10-M.fastq.gz
cutadapt -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o M-T6-12-M-trim.fq.gz M-T6-12-M.fastq.gz
cutadapt -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o S-T1-01-M-trim.fq.gz S-T1-01-M.fastq.gz
cutadapt -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o S-T1-07-M-trim.fq.gz S-T1-07-M.fastq.gz
cutadapt -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o S-T1-12-M-trim.fq.gz S-T1-12-M.fastq.gz
cutadapt -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o S-T3-09-M-trim.fq.gz S-T3-09-M.fastq.gz
cutadapt -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o S-T3-14-M-trim.fq.gz S-T3-14-M.fastq.gz
