#!/bin/bash
#PBS -N bismark_prep
#PBS -l walltime=04:00:00
#PBS -l select=1:mem=50gb
#PBS -q large
#PBS -o bismark_prep.out
#PBS -e bismark_prep.err
#PBS -M alb0293@auburn.edu
#PBS -m abe

module load bismark/0.23.0
module load bowtie2

echo "Starting Bowtie2 index build..."

#make Bowtie2 directory to store the index files
mkdir -p /home/aubaxb001/genomes/Lrhomboides/Bowtie2

#move to bowtie2 directory
cd /home/aubaxb001/genomes/Lrhomboides/Bowtie2

#build bowtie2 indexes from the genome fasta in the parent directory
bowtie2-build --threads 12 --verbose /home/aubaxb001/genomes/Lrhomboides/Lrhomboides_genome.fa Lrhomboides

echo "Starting Bismark genome prep (creating bisulfite indexes)..."

#move back to Lrhomboides directory so it creates Bisulfite_Genome directory in that directory
cd /home/aubaxb001/genomes/Lrhomboides

#run bismark genome prep with indexed files (will output into Bismark_Genome directory that it makes automatically when running bismark_genome_prep)
bismark_genome_preparation --bowtie2 --parallel 6 --verbose .

echo "Bismark genome preparation complete"
