#!/bin/bash
#PBS -N merge_lanes
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=6:mem=300gb
#PBS -q bigmem
#PBS -o merge_lanes.out
#PBS -e merge_lanes.err
#PBS -M alb0293@auburn.edu
#PBS -m abe

#have paired end reads from novogene WGBS
#check for the right number of unique sample IDs for both R1 and R2
#script modified from: https://www.benlaufer.com/CpG_Me/#merging-lanes

#want to make it loop through all of the directories per sample inside methylation_data directory
base_dir="/home/aubaxb001/methylation_data"

#loop through each sample directory
	 # Extract the directory name/sample name
for sample_dir in "$base_dir"/*/; do
	echo "Found directory: $sample_dir"
    sample_name=$(basename "$sample_dir")

    echo "Processing sample: $sample_name"

    # Change into the sample directory
    cd "$sample_dir" || { echo "Cannot enter directory $sample_dir"; continue; }


if ! ls *_1.fq.gz 1> /dev/null 2>&1 || ! ls *_2.fq.gz 1> /dev/null 2>&1; then
    echo "No fq.gz files found in $sample_dir, skipping..."
    cd "$base_dir" || exit
    continue
fi

#extract the sample ID from the file names, sort IDs, remove duplicates, count unique sample IDs
	#want to include all of the sample ID info after the first 4 _ (this will include C_T1_03_FF as the naming conventions for samples
countFASTQ(){
    awk -F '_' '{print $1"_"$2"_"$3"_"$4}' | \
    sort -u | \
    wc -l
}
export -f countFASTQ

#make sure script cd into each sample directory
echo "Current directory: $(pwd)"
ls -l
ls -1 *_1.fq.gz

#make variables for the counts of R1 and R2 - counting unique sample IDs for both separately and stores them into variables
  #list all files that match the patterns, pipe to countFASTQ to now have sample IDs, count of nuimber of unique IDs
R1=`ls -1 *_1.fq.gz | countFASTQ`
R2=`ls -1 *_2.fq.gz | countFASTQ`

#make sure the numbner of unique forward (R1) and reverse (R2) read samples match before merging

if [ ${R1} = ${R2} ]
then

	#count number of "lanes" (sequencing runs) the samples were sequenced across, list all forward read files
                 #print the first parts of the file name
		#count how many times each unique ID appears (how many lanes per sample)
		#extract those counts
		#get unique counts
		#says how many lanes and samples that were found
		#prints error if the counts of forward and reverse dont match

	lanes=`ls -1 *_1.fq.gz | \
        awk -F '_' '{print $1"_"$2"_"$3"_"$4}' | \
        sort | \
        uniq -c | \
        awk -F ' ' '{print $1}' | \
        sort -u`
        echo "${R1} samples sequenced across ${lanes} lanes identified for merging"
else
        echo "ERROR: There are ${R1} R1 files and ${R2} R2 files"
        exit 1
fi

#create file of unique IDs based on _ delimiter and first string
	#list all files ending with fq.gz
	#extract unique sample ID
	#sort and redirect the output into new txt file

ls -1 *fq.gz | \
awk -F '_' '{print $1"_"$2"_"$3"_"$4}' | \
sort -u > \
task_samples.txt

#test merge commands for each read

#prints the commands that would concatenate all R1 fastq files for wach sample into a single file
#read sample IDs from txt file, run mergeLanesTest function on each sample in parallel, 6 jobs at once, output will be the list of commands printed for each sample
mergeLanesTest(){
    i=$1
    echo cat ${i}\_*_1.fq.gz \> ${i}\_1.fq.gz
    echo cat ${i}\_*_2.fq.gz \> ${i}\_2.fq.gz
}
while read sample; do
        mergeLanesTest "$sample"
done < task_samples.txt


#use merge commands for each read by removing echo and actually running
mergeLanes(){
    i=$1
    cat ${i}\_*_1.fq.gz > ${i}\_1.fq.gz
    cat ${i}\_*_2.fq.gz > ${i}\_2.fq.gz
}
while read sample; do
	mergeLanes "$sample"
done < task_samples.txt

cd "$base_dir" || exit
done

echo "All samples processed!"

#if the script runs, the samples have been merged across lanes
#also have a task_samples.txt file for the next steps (ex correcting for methylation bias)
