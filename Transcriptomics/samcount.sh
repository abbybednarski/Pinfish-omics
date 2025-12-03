#!/bin/bash
#SBATCH --job-name=counts
#SBATCH --ntasks=1
#SBATCH --partition=mab0205_bg4
#SBATCH --mem=20G
#SBATCH --time=08:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alb0293@auburn.edu

./samcount.pl C-T1-03-M-t.sam transcriptome_seq2iso.tab aligner=bowtie2 >C-T1-03-M.counts
./samcount.pl C-T1-14-M-t.sam transcriptome_seq2iso.tab aligner=bowtie2 >C-T1-14-M.counts
./samcount.pl C-T1-16-M-t.sam transcriptome_seq2iso.tab aligner=bowtie2 >C-T1-16-M.counts
./samcount.pl C-T5-12-M-t.sam transcriptome_seq2iso.tab aligner=bowtie2 >C-T5-12-M.counts
./samcount.pl M-T2-09-M-t.sam transcriptome_seq2iso.tab aligner=bowtie2 >M-T2-09-M.counts
./samcount.pl M-T2-16-M-t.sam transcriptome_seq2iso.tab aligner=bowtie2 >M-T2-16-M.counts
./samcount.pl M-T6-06-M-t.sam transcriptome_seq2iso.tab aligner=bowtie2 >M-T6-06-M.counts
./samcount.pl M-T6-10-M-t.sam transcriptome_seq2iso.tab aligner=bowtie2 >M-T6-10-M.counts
./samcount.pl M-T6-12-M-t.sam transcriptome_seq2iso.tab aligner=bowtie2 >M-T6-12-M.counts
./samcount.pl S-T1-01-M-t.sam transcriptome_seq2iso.tab aligner=bowtie2 >S-T1-01-M.counts
./samcount.pl S-T1-07-M-t.sam transcriptome_seq2iso.tab aligner=bowtie2 >S-T1-07-M.counts
./samcount.pl S-T1-12-M-t.sam transcriptome_seq2iso.tab aligner=bowtie2 >S-T1-12-M.counts
./samcount.pl S-T3-09-M-t.sam transcriptome_seq2iso.tab aligner=bowtie2 >S-T3-09-M.counts
./samcount.pl S-T3-14-M-t.sam transcriptome_seq2iso.tab aligner=bowtie2 >S-T3-14-M.counts

./expression_compiler.pl *.counts > allcounts.txt
