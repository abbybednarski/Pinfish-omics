#!/bin/bash
#SBATCH --job-name=bowtie
#SBATCH --ntasks=8
#SBATCH --partition=mab0205_bg4
#SBATCH --mem=20G
#SBATCH --time=08:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alb0293@auburn.edu

module load bowtie2

bowtie2 --local -x Transcriptome-Sequences-Trimmed-Names.fa -U C-T1-03-M-trim.fq.gz -S C-T1-03-M-t.sam --no-hd --no-sq --no-unal -k 5 --threads 8
bowtie2 --local -x Transcriptome-Sequences-Trimmed-Names.fa -U C-T1-14-M-trim.fq.gz -S C-T1-14-M-t.sam --no-hd --no-sq --no-unal -k 5 --threads 8
bowtie2 --local -x Transcriptome-Sequences-Trimmed-Names.fa -U C-T1-16-M-trim.fq.gz -S C-T1-16-M-t.sam --no-hd --no-sq --no-unal -k 5 --threads 8
bowtie2 --local -x Transcriptome-Sequences-Trimmed-Names.fa -U C-T5-12-M-trim.fq.gz -S C-T5-12-M-t.sam --no-hd --no-sq --no-unal -k 5 --threads 8
bowtie2 --local -x Transcriptome-Sequences-Trimmed-Names.fa -U M-T2-09-M-trim.fq.gz -S M-T2-09-M-t.sam --no-hd --no-sq --no-unal -k 5 --threads 8
bowtie2 --local -x Transcriptome-Sequences-Trimmed-Names.fa -U M-T2-16-M-trim.fq.gz -S M-T2-16-M-t.sam --no-hd --no-sq --no-unal -k 5 --threads 8
bowtie2 --local -x Transcriptome-Sequences-Trimmed-Names.fa -U M-T6-06-M-trim.fq.gz -S M-T6-06-M-t.sam --no-hd --no-sq --no-unal -k 5 --threads 8
bowtie2 --local -x Transcriptome-Sequences-Trimmed-Names.fa -U M-T6-10-M-trim.fq.gz -S M-T6-10-M-t.sam --no-hd --no-sq --no-unal -k 5 --threads 8
bowtie2 --local -x Transcriptome-Sequences-Trimmed-Names.fa -U M-T6-12-M-trim.fq.gz -S M-T6-12-M-t.sam --no-hd --no-sq --no-unal -k 5 --threads 8
bowtie2 --local -x Transcriptome-Sequences-Trimmed-Names.fa -U S-T1-01-M-trim.fq.gz -S S-T1-01-M-t.sam --no-hd --no-sq --no-unal -k 5 --threads 8
bowtie2 --local -x Transcriptome-Sequences-Trimmed-Names.fa -U S-T1-07-M-trim.fq.gz -S S-T1-07-M-t.sam --no-hd --no-sq --no-unal -k 5 --threads 8
bowtie2 --local -x Transcriptome-Sequences-Trimmed-Names.fa -U S-T1-12-M-trim.fq.gz -S S-T1-12-M-t.sam --no-hd --no-sq --no-unal -k 5 --threads 8
bowtie2 --local -x Transcriptome-Sequences-Trimmed-Names.fa -U S-T3-09-M-trim.fq.gz -S S-T3-09-M-t.sam --no-hd --no-sq --no-unal -k 5 --threads 8
bowtie2 --local -x Transcriptome-Sequences-Trimmed-Names.fa -U S-T3-14-M-trim.fq.gz -S S-T3-14-M-t.sam --no-hd --no-sq --no-unal -k 5 --threads 8
