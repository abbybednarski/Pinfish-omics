#First, find overlapping methylation sites across pairwise comparisons

#read in all csv files and load in all libraries

library(dplyr)
library(purrr)

#input files are DMC files
files <- list.files("/Users/abbybednarski/Desktop/data/", full.names = TRUE)
dfs <- map(files, read.csv)

#name each list element by the file name
names(dfs) <- basename(files)

#standardize and create a kind of key - should be fairly similar format already
dfs <- map(dfs, ~ .x %>%
             rename(chr = chr,       # already called chr
             pos = start) %>%      # start = end = position of the CpG
             mutate(location = paste(chr, pos, sep = ":")))

#find CpGs that are common across all files
common_locations <- reduce(map(dfs, ~ .x$location), #extract location column from each file
  intersect) #find value that exists in every file

#this returned character (empty) meaning there is no overlapping locations
  #ORRRR something went wrong (ex the namings are slightly off)

#run the following code to see if there is overlap between the first 2 files
head(dfs[[1]]$location)
#output some contig names

head(dfs[[2]]$location)
#output some contig names

length(intersect(dfs[[1]]$location, dfs[[2]]$location))
#128

#this shows there IS overlap for at least the first 2 files indicating overlapping sites for pairwise comparisons
    #SOOOOO no single overlapping CpG sites across all 3 treatments

#lets try pairwise overlaps instead now

# find pairwise overlaps
find_pairwise_overlap <- function(df1, df2, name1, name2) { #dataframe of dmc files
  common <- intersect(df1$location, df2$location) #find all locations that are in both files
  
  if(length(common) == 0) return(NULL)  # if no overlaps, stop
  
  map_dfr(list(df1, df2), #loop over the dataframe you made
    ~ filter(.x, location %in% common), #keep only the locations in the common vector you made above
    .id = "source_file" #add a column called source file that tells you where the row came from initially
  ) %>%
    mutate(pair = paste(name1, name2, sep = "_vs_")) #make a new column to show pwc vs pwc
}

# make list of all potential possibilities of pairwise combinations
file_pairs <- combn(names(dfs), 2, simplify = FALSE)

# loop over each combo to find overlaps for each file/comparison
pairwise_overlaps <- map(file_pairs, function(pair) {
  df1 <- dfs[[pair[1]]]
  df2 <- dfs[[pair[2]]]
  
  overlap <- find_pairwise_overlap(df1, df2, pair[1], pair[2]) #call function from earlier
  
  if(!is.null(overlap)) {
    #write to csv
    outname <- paste0("overlap_", pair[1], "_vs_", pair[2], ".csv") #make files for overlapping dmcs one at a time
    write.csv(overlap, outname, row.names = FALSE)
  }
  
  overlap  # return overlap to list
})

#only output files for CvM vs. CvS and CvS vs. MvS
  #so no overlapping DMCs between CvM vs. MvS

#write csvs separately for files of overlapping DMCs

#list 2 element (CvM vs CvS)
write.csv(pairwise_overlaps[[2]], "MHvC_vs_SHvC_overlap_DMCs.csv", row.names = FALSE)

#list 3 element (CvS vs MvS)
write.csv(pairwise_overlaps[[3]], "SHvC_vs_MHvSH_overlap_DMCs.csv", row.names = FALSE)


#look at the overlapping sites

#read in files
#setwd("~/Desktop/Methylation")
CMCS <- read.csv("MHvC_vs_SHvC_overlap_DMCs.csv")
CSMS <- read.csv("SHvC_vs_MHvSH_overlap_DMCs.csv")

############################ CvM vs CvS ##############################
#quick summary stats
CMCS %>%
  group_by(pair) %>%
  summarize(
    n_overlaps = n(),
    mean_methdiff = mean(meth.diff),
    median_methdiff = median(meth.diff))

#mean = 29.9, median = 49.3
  #average is hypermethylated then

#try to match to ref genome to get locations
str(CMCS)

#convert to GRanges
library(GenomicRanges)

dmc_gr <- GRanges(
  seqnames = CMCS$chr,
  ranges   = IRanges(start = CMCS$pos, end = CMCS$end),
  strand   = CMCS$strand)

dmc_gr
#looks good
#256 DMCs

#look through gff annotation file for one specific PGA scaffold you want to see if this is going to be possible
#run in command line: grep PGA_scaffold0__49_contigs__length_39332093 Lrhomboides_complete_MAKER_annotation.gff
  #where the scaffold is just one of the DMCs you found and the file at the end is the annotation file you want to search through
  #when i did this, grep output HUNDREDS of lines - which is a good thing!

#see if all DMCs are found in gff
all(seqlevels(dmc_gr) %in% seqlevels(gff))
#TRUE
  #so yes!!!

#test all DMCs against ref genome
hits <- findOverlaps(dmc_gr, gff)
hits
#R found 1000 hits but i only have 256 DMCs so that means many DMCs overlapped with multiple annotations

#queryHits - index of the DMC
  #example: queryHits = 2 refers to 2nd CpG in DMC list

#subjectHits - index of annotation feature
  #example subjectHits = 8424516 refers to 8424516th row in GFF

#SOOO output is showing DMCs overlapping with multiple annotation features like

#Hits object with 1000 hits and 0 metadata columns:
#queryHits subjectHits
#<integer>   <integer>
#  [1]         1     8424516
#[2]         2     8424516 --> here 
#[3]         2     8430828 --> here 
#[4]         2     8430829 --> here


############# find unique DMCs rather than overlapping ########

#first, read in all DMCs files and overlapping DMC files

#MvS
setwd("~/Desktop/Methylation/***output files with T2***/MvS/3L")
MvS_DMCs <- read.csv("DMRs_MvS.csv")
#CvM
setwd("~/Desktop/Methylation/***output files with T2***/CvM/3L")
CvM_DMCs <- read.csv("DMRs_CvM.csv")
#CvS
setwd("~/Desktop/Methylation/***output files with T2***/CvS/3L")
CvS_DMCs <- read.csv("DMRs_CvS.csv")

#overlapping DMCs
setwd("~/Desktop/Methylation")
overlapping_CMCS <- read.csv("MHvC_vs_SHvC_overlap_DMCs.csv")
overlapping_CSMS <- read.csv("SHvC_vs_MHvSH_overlap_DMCs.csv")

#make sure they all have the same headers/structure

#make sure start column nrow is the same number of obs in envi for each file
length(unique(MvS_DMCs$start))
#260 --> good
length(unique(CvS_DMCs$start))
#2044 --> good
length(unique(CvM_DMCs$start))
#2217 --> good
length(unique(overlapping_CMCS$pos)) #pos instead of start column name
#128 --> half bc overlapping
length(unique(overlapping_CSMS$pos))
#1 --> half bc overlapping

#filter for nonoverlapping DMCs for each file

#CvM non overlapping
CvM_nonoverlap <- CvM_DMCs[
  !CvM_DMCs$start %in% overlapping_CMCS$pos,
]

#MvS
MvS_nonoverlap <- MvS_DMCs[
  !MvS_DMCs$start %in% overlapping_CSMS$pos,
]

#CvS
CvS_nonoverlap <- CvS_DMCs[
  !CvS_DMCs$start %in% c(overlapping_CMCS$pos, overlapping_CSMS$pos),
]

#make sure everything matches up math wise by checking # of obs in environment on the right!

#Ive got
  #CvM 2089 unique DMCs
  #CvS 1915 unique DMCs
  #MvS 259 unique DMCs

#write csvs for nonoverlapping DMCs
write.csv(CvS_nonoverlap, "SHvC-nonoverlapping-DMCs.csv", row.names = FALSE)
write.csv(MvS_nonoverlap, "MHvSH-nonoverlapping-DMCs.csv", row.names = FALSE)
write.csv(CvM_nonoverlap, "MHvC-nonoverlapping-DMCs.csv", row.names = FALSE)



###################### annotate intron, exon, UTRs, genes #####################
library(GenomicRanges)
library(rtracklayer)

#read in nonoverlapping and overlapping DMC files
#MvS
setwd("~/Desktop/Methylation/***output files with T2***/MvS/3L")
MvS_nonoverlap <- read.csv("DMRs_MvS.csv")
#CvM
setwd("~/Desktop/Methylation/***output files with T2***/CvM/3L")
CvM_nonoverlap <- read.csv("DMRs_CvM.csv")
#CvS
setwd("~/Desktop/Methylation/***output files with T2***/CvS/3L")
CvS_nonoverlap <- read.csv("DMRs_CvS.csv")

#overlapping DMCs
setwd("~/Desktop/Methylation")
overlapping_CMCS <- read.csv("MHvC_vs_SHvC_overlap_DMCs.csv")
overlapping_CSMS <- read.csv("SHvC_vs_MHvSH_overlap_DMCs.csv")

#load in gff file
setwd("~/Desktop")
gff <- import("Lrhomboides_complete_MAKER_annotation.gff")

#make function to convert DMCs to GRanges
make_dmc_gr <- function(DMC_file) {
  GRanges(
    seqnames = DMC_file$chr,
    ranges = IRanges(start = DMC_file$pos, end = DMC_file$end),
    strand = DMC_file$strand
  )
}

#apply to all DMC comparisons
CvS_gr <- make_dmc_gr(CvS_nonoverlap)
CvM_gr <- make_dmc_gr(CvM_nonoverlap)
MvS_gr <- make_dmc_gr(MvS_nonoverlap)
overlapping_CMCS_gr <- make_dmc_gr(overlapping_CMCS) #change ranges in make_gr function to pos and end rather than start and end
overlapping_CSMS_gr <- make_dmc_gr(overlapping_CSMS)

#look at the gene features in the annotation file
unique(gff$type)
# [1] contig                   match                    match_part               expressed_sequence_match
#[5] gene                     mRNA                     exon                     CDS                     
#[9] protein_match            five_prime_UTR           three_prime_UTR        

#lets annotate exons, 5' UTR and 3' UTR

#extract gene features from annotation file
exons <- gff[gff$type == "exon"]
five_utr <- gff[gff$type == "five_prime_UTR"]
three_utr <- gff[gff$type == "three_prime_UTR"]
genes <- gff[gff$type == "gene"]

#annotate DMCs with feature types you extracted above
annotate_DMC <- function(dmc_gr, DMC_file) {
  
  DMC_file$feature <- "NA"  #add new column called feature, default NA
  
  overlaps_exon  <- findOverlaps(dmc_gr, exons)
  overlaps_5utr  <- findOverlaps(dmc_gr, five_utr)
  overlaps_3utr  <- findOverlaps(dmc_gr, three_utr)
  overlaps_gene  <- findOverlaps(dmc_gr, genes)
  
  #find gene overlaps and add to table
  overlaps_gene <- findOverlaps(dmc_gr, genes)
  DMC_file$feature[queryHits(overlaps_gene)] <- "gene"
  
  #find introns overlaps
  #this is including only DMCs in genes not already labeled as gene
  exon_hits <- findOverlaps(dmc_gr, exons)
  five_hits <- findOverlaps(dmc_gr, five_utr)
  three_hits <- findOverlaps(dmc_gr, three_utr)
  
  for(i in seq_len(nrow(DMC_file))) {
    feats <- c()
    
    if(i %in% queryHits(overlaps_3utr)) feats <- c(feats, "three_prime_UTR")
    if(i %in% queryHits(overlaps_5utr)) feats <- c(feats, "five_prime_UTR")
    if(i %in% queryHits(overlaps_exon)) feats <- c(feats, "exon")
    
    # intron = inside gene but not exon/UTR
    if(i %in% queryHits(overlaps_gene) & length(feats) == 0) feats <- c(feats, "intron")
    
    # gene = highest priority, always append
    if(i %in% queryHits(overlaps_gene)) feats <- c(feats, "gene")
    
    if(length(feats) > 0) DMC_file$feature[i] <- paste(feats, collapse=";")
  }
  
  return(DMC_file)
}

#run function on all DMC sets
CvS_anno <- annotate_DMC(CvS_gr, CvS_nonoverlap)
CvM_anno <- annotate_DMC(CvM_gr, CvM_nonoverlap)
MvS_anno <- annotate_DMC(MvS_gr, MvS_nonoverlap)
overlapping_CMCS_anno <- annotate_DMC(overlapping_CMCS_gr, overlapping_CMCS)
overlapping_CSMS_anno <- annotate_DMC(overlapping_CSMS_gr, overlapping_CSMS)

library(dplyr)
library(tidyr)

table(CvS_anno$feature)
#exon;gene               intron;gene                        NA 
#59                       325                      1530 
#three_prime_UTR;exon;gene 
#1 
#so going to categorize as 325 + 60 = 385 genes / 1915 DMCs

table(CvM_anno$feature)
#exon;gene               intron;gene                        NA 
#34                       383                      1671 
#three_prime_UTR;exon;gene 
#1 
#so going to categorize as 383 + 34 + 1 = 417 genes / 2089 DMCs

table(MvS_anno$feature)
#exon;gene intron;gene          NA 
#4          48         207 
#so going to categorize as 4 + 48 = 52 genes / 259 DMCs

table(overlapping_CMCS_anno$feature)
#exon;gene intron;gene          NA 
#8          46         202 
#so going to categorize as 8 + 46 = 54 genes / 256 DMCs

table(overlapping_CSMS_anno$feature)
#NA = 2
#so nothing! going to stop doing analysis on this dataset since nothing here**

######################## match the Lrhomboides gene names to descriptive gene names ###########################

#match with annotations.csv
setwd("~/Desktop")
annotation_csv <- read.csv("all_annotations.csv")
colnames(annotation_csv)
#[1] "Lrhomboides.Gene.Name"  "Percent.Identity"       "E.value"               
#[4] "Bit.Score"              "Query.Coverage.Percent" "Match.Ensembl.ID"      
#[7] "Gene.Name"              "Gene.Description"  

#function to add gene name and gene description
annotate_dmcs_full <- function(dmc_table, annotation) {
  #make columns with default NA
  dmc_table$Gene.Name <- NA
  dmc_table$Gene.Description <- NA
  
  #annotate rows where gene ID exists
  has_gene <- !is.na(dmc_table$gene)
  
  #add gene name
  dmc_table$Gene.Name[has_gene] <- annotation$Gene.Name[
    match(dmc_table$gene[has_gene], annotation$Lrhomboides.Gene.Name)
  ]
  
  #add gene description
  dmc_table$Gene.Description[has_gene] <- annotation$Gene.Description[
    match(dmc_table$gene[has_gene], annotation$Lrhomboides.Gene.Name)
  ]
  
  return(dmc_table)
}

#run through each dataframe
CvM_nonoverlap <- annotate_dmcs_full(CvM_anno, annotation_csv)
MvS_nonoverlap <- annotate_dmcs_full(MvS_anno, annotation_csv)
CvS_nonoverlap <- annotate_dmcs_full(CvS_anno, annotation_csv)
overlapping <- annotate_dmcs_full(overlapping_CMCS_anno, annotation_csv)

#write to csv
write.csv(CvS_nonoverlap, "CvS_nonoverlap_annotated.csv", row.names = FALSE)
write.csv(CvM_nonoverlap, "CvM_nonoverlap_annotated.csv", row.names = FALSE)
write.csv(MvS_nonoverlap, "MvS_nonoverlap_annotated.csv", row.names = FALSE)
write.csv(overlapping, "CvS_CvM_overlapping_annotated.csv", row.names = FALSE)

#extract only chr lines that have a gene in gene column rather than NA
head(CvS_nonoverlap)
CvS_nonoverlap_gene <- CvS_nonoverlap[!is.na(CvS_nonoverlap$gene), ]
head(CvS_nonoverlap_gene)
#385 - good!

CvM_nonoverlap_gene <- CvM_nonoverlap[!is.na(CvM_nonoverlap$gene), ]
#418 - good!

MvS_nonoverlap_gene <- MvS_nonoverlap[!is.na(MvS_nonoverlap$gene), ]
#52 - good!

overlapping_gene <- overlapping_CMCS[!is.na(overlapping_CMCS$gene), ]
#54 - good!

#then write csvs for those files - those have all of the DMCs that matched to a gene but might not have a real gene name
write.csv(CvS_nonoverlap_gene, "SHvC-DMGs.csv", row.names = FALSE)
write.csv(CvM_nonoverlap_gene, "MHvC-DMGs.csv", row.names = FALSE)
write.csv(MvS_nonoverlap_gene, "MHvSH-DMGs.csv", row.names = FALSE)
write.csv(overlapping_gene, "overlapping-DMGs.csv", row.names = FALSE)


################ hypothesis driven promoter analysis ###################

library(GenomicRanges)
library(rtracklayer)
library(dplyr)

#read in unannotated DMCs
#MvS
setwd("~/Desktop/Methylation/***output files with T2***/MvS/3L")
MvS_nonoverlap <- read.csv("DMRs_MvS.csv")
#CvM
setwd("~/Desktop/Methylation/***output files with T2***/CvM/3L")
CvM_nonoverlap <- read.csv("DMRs_CvM.csv")
#CvS
setwd("~/Desktop/Methylation/***output files with T2***/CvS/3L")
CvS_nonoverlap <- read.csv("DMRs_CvS.csv")

#overlapping DMCs
setwd("~/Desktop/Methylation")
overlapping <- read.csv("MHvC_vs_SHvC_overlap_DMCs.csv")

#read in annotation file
setwd("~/Downloads")
gff <- import("Lrhomboides_complete_MAKER_annotation.gff")

#extract gene features
genes <- gff[gff$type == "gene"]
exons <- gff[gff$type == "exon"]
five_utr <- gff[gff$type == "five_prime_UTR"]
three_utr <- gff[gff$type == "three_prime_UTR"]

#make GRanges for each comparison

####### CvS #########
CvS_gr <- GRanges(
  seqnames = CvS_nonoverlap$chr,
  ranges = IRanges(start = CvS_nonoverlap$start, end = CvS_nonoverlap$end),
  strand = CvS_nonoverlap$strand
)

########## CvM ##########
CvM_gr <- GRanges(
  seqnames = CvM_nonoverlap$chr,
  ranges = IRanges(start = CvM_nonoverlap$start, end = CvM_nonoverlap$end),
  strand = CvM_nonoverlap$strand
)

########### MvS #########
MvS_gr <- GRanges(
  seqnames = MvS_nonoverlap$chr,
  ranges = IRanges(start = MvS_nonoverlap$start, end = MvS_nonoverlap$end),
  strand = MvS_nonoverlap$strand
)

########### overlapping #########
overlapping_gr <- GRanges(
  seqnames = overlapping$chr,
  ranges = IRanges(start = overlapping$pos, end = overlapping$end),
  strand = overlapping$strand
)

#make function to add hypothesized promoter information to DMC dataframes
    #following: https://www.tandfonline.com/doi/full/10.1080/15592294.2020.1795597#d1e480
        #their promoter region was 500bp upstream so going with that
        #also ran with window of 300bp but didnt change numbers much so leaving at 500bp window
annotate_promoter <- function(dmc_file, dmc_gr, genes, upstream = 500) { #pick the window you want
  
  #actually make promoter window around each DMC
  dmc_window <- GRanges(
    seqnames = seqnames(dmc_gr),
    ranges = IRanges(
      start = pmax(start(dmc_gr) - upstream, 1), #start DMC - upstream value (500 as indicated above) but make sure start is never <1
      end   = start(dmc_gr) -1 #end DMC 1 bp before DMC - not considering anything downstream for promoter
    )
  )
  
  #find genes that overlap each DMC window you just made
  overlaps <- findOverlaps(dmc_window, genes)
  #store DMCs that overlap a gene
  dmc_indices <- queryHits(overlaps)
  #store genes
  gene_indices <- subjectHits(overlaps)
  
  #make a new column in DMC dataframes with default NA
  dmc_file$potential_promoter_gene <- NA
  
  #be sure youre able to add multiple genes if you need to - temporary dataframe where each row is one DMC-gene overlap only if there are multiple overlaps
  if(length(dmc_indices) > 0) {
    
    # Make temporary dataframe of DMC-gene overlaps
    promoter_df <- data.frame(
      DMC_index = dmc_indices,
      GeneID = as.character(genes$ID[gene_indices]),
      stringsAsFactors = FALSE
    )
    
    # Keep only the first gene annotation per DMC
    promoter_df <- promoter_df[!duplicated(promoter_df$DMC_index), ]
    
    # Add gene annotation back to original DMC dataframe
    dmc_file$potential_promoter_gene[promoter_df$DMC_index] <- 
      promoter_df$GeneID
  }
  return(dmc_file)
}

#apply function to each comparison
CvS_nonA_promoters <- annotate_promoter(CvS_nonoverlap, CvS_gr, genes)
CvM_nonA_promoters <- annotate_promoter(CvM_nonoverlap, CvM_gr, genes)
MvS_nonA_promoters <- annotate_promoter(MvS_nonoverlap, MvS_gr, genes)
overlapping_nonA_promoters <- annotate_promoter(overlapping, overlapping_gr, genes)

#extract gene counts

sum(!is.na(CvS_nonA_promoters$potential_promoter_gene))
#782 gene hits - 782 DMCs overlap promoter window defined above
length(unique(CvS_nonA_promoters$potential_promoter_gene))
#706 - map to 706 unique genes

sum(!is.na(CvM_nonA_promoters$potential_promoter_gene))
#842
length(unique(CvM_nonA_promoters$potential_promoter_gene))
#774 - map to 778 unique genes

sum(!is.na(MvS_nonA_promoters$potential_promoter_gene))
#128
length(unique(MvS_nonA_promoters$potential_promoter_gene))
#124 - map to 124 unique genes

sum(!is.na(overlapping_nonA_promoters$potential_promoter_gene))
#114
length(unique(overlapping_nonA_promoters$potential_promoter_gene))
#55 - map to 55 unique genes

sum(!is.na(overlapping_nonA_promoters$potential_promoter_gene))
#114
length(unique(overlapping_nonA_promoters$potential_promoter_gene))
#56

#annotate like you did above

#add gene name and desc
#match with annotations.csv
setwd("~/Desktop")
annotation_csv <- read.csv("all_annotations.csv")
colnames(annotation_csv)
#[1] "Lrhomboides.Gene.Name"  "Percent.Identity"       "E.value"               
#[4] "Bit.Score"              "Query.Coverage.Percent" "Match.Ensembl.ID"      
#[7] "Gene.Name"              "Gene.Description"  

#function to add gene name and gene description
annotate_promoters_full <- function(promoter_table, annotation) {
  #make columns with default NA
  promoter_table$Gene.Name <- NA
  promoter_table$Gene.Description <- NA
  
  #annotate rows where gene ID exists
  has_gene <- !is.na(promoter_table$potential_promoter_gene)
  
  #add gene name
  promoter_table$Gene.Name[has_gene] <- annotation$Gene.Name[
    match(promoter_table$potential_promoter_gene[has_gene], annotation$Lrhomboides.Gene.Name)
  ]
  
  #add gene description
  promoter_table$Gene.Description[has_gene] <- annotation$Gene.Description[
    match(promoter_table$potential_promoter_gene[has_gene], annotation$Lrhomboides.Gene.Name)
  ]
  
  return(promoter_table)
}

#run through each dataframe
CvM_nonoverlap_promoters_anno <- annotate_promoters_full(CvM_promoters, annotation_csv)
MvS_nonoverlap_promoters_anno <- annotate_promoters_full(MvS_promoters, annotation_csv)
CvS_nonoverlap_promoters_anno <- annotate_promoters_full(CvS_promoters, annotation_csv)
overlapping_promoters_anno <- annotate_promoters_full(overlapping_promoters, annotation_csv)

#SO HIGH!

#most DMCs close to annotated genes already so playing role in regulation?
#SOO most DMCs overlapping the upstream 300-500bp promoter region I made
#these have been annotated to show these are genes
#also known (kinda) before since a lot of DMCs annotated to genes - the DMC is close to the gene annotation

#changing window from 500bp to 300bp BARELY changed the numbers showing DMCs are pretty much right at where the gene starts

#the DMC is reinforced to be close to the gene annotation by the fact that the shrinking window didnt change the numbers
#the DMC is right where the gene starts

#filter to only include rows with a potential promoter gene
CvM_promoters_filtered <- CvM_nonoverlap_promoters_anno[!is.na(CvM_nonoverlap_promoters_anno$potential_promoter_gene), ]
CvS_promoters_filtered <- CvS_nonoverlap_promoters_anno[!is.na(CvS_nonoverlap_promoters_anno$potential_promoter_gene), ]
MvS_promoters_filtered <- MvS_nonoverlap_promoters_anno[!is.na(MvS_nonoverlap_promoters_anno$potential_promoter_gene), ]
overlapping_promoters_filtered <- overlapping_promoters_anno[!is.na(overlapping_promoters_anno$potential_promoter_gene), ]

#write to csv
write.csv(CvS_promoters_filtered, "CvS_annotated_promoters.csv", row.names = FALSE)
write.csv(CvM_promoters_filtered, "CvM_annotated_promoters.csv", row.names = FALSE)
write.csv(MvS_promoters_filtered, "MvS_annotated_promoters.csv", row.names = FALSE)
write.csv(overlapping_promoters_filtered, "overlapping_annotated_promoters.csv", row.names = FALSE)


############################## mean and median methylation for pairwise comparisons ####################################
  #note that mean for this is taking hypo and hyper together meaning they kinda cancel each other out and means will be semi close to 0 because of that
  #also why you want to find median
library(rtracklayer)

#read in files
setwd("~/Downloads")
CvM <- read.csv("CvM_nonoverlap_annotated.csv")
CvS <- read.csv("CvS_nonoverlap_annotated.csv")
MvS <- read.csv("MvS_nonoverlap_annotated.csv")

#CvM
CvM %>%
  summarize(
    n_overlaps = n(),
    mean_methdiff = mean(meth.diff),
    median_methdiff = median(meth.diff))
#overlaps = 2089, mean = -3.38, median = -31.58

#CvS
CvS %>%
  summarize(
    n_overlaps = n(),
    mean_methdiff = mean(meth.diff),
    median_methdiff = median(meth.diff))
#overlaps = 1915, mean = -4.28, median = -33.33

#MvS
MvS %>%
  summarize(
    n_overlaps = n(),
    mean_methdiff = mean(meth.diff),
    median_methdiff = median(meth.diff))
#overlaps = 259, mean = -2.26, median = -29.27

#out of unique genes, how many were hypo vs hypermethylated
library(dplyr)

######## find unique gene counts ########
# CvM
CvM_unique_genes <- unique(CvM$Gene.Name)
length(CvM_unique_genes)
#288 - 1 for NA entry so 287 unique genes

# CvS
CvS_unique_genes <- unique(CvS$Gene.Name)
length(CvS_unique_genes)
#243 - 1 so 242 unique genes

# MvS
MvS_unique_genes <- unique(MvS$Gene.Name)
length(MvS_unique_genes)
#42 - 1 so 41 unique genes

# overlapping
overlapping_unique_genes <- unique(overlapping$Gene.Name)
length(overlapping_unique_genes)
#20 - 1 so 19 unique

###### CvM #######
gene_direction <- CvM %>%
  filter(!is.na(Gene.Name)) %>%
  mutate(direction = case_when(
    meth.diff > 24 ~ "hyper",
    meth.diff < -24 ~ "hypo"
  )) %>%
  distinct(Gene.Name, direction)

gene_direction %>%
  count(direction)
#hyper 139
#hypo 157
#caveat: some genes might be hyper AND hypomethylated

######## CvS ########
gene_direction <- CvS %>%
  filter(!is.na(Gene.Name)) %>%
  mutate(direction = case_when(
    meth.diff > 24 ~ "hyper",
    meth.diff < -24 ~ "hypo"
  )) %>%
  distinct(Gene.Name, direction)

gene_direction %>%
  count(direction)
#hyper 127
#hypo 126

####### MvS #########
gene_direction <- MvS %>%
  filter(!is.na(Gene.Name)) %>%
  mutate(direction = case_when(
    meth.diff > 24 ~ "hyper",
    meth.diff < -24 ~ "hypo"
  )) %>%
  distinct(Gene.Name, direction)

gene_direction %>%
  count(direction)
#hyper 19
#hypo 23

###### overlapping #########
gene_direction <- overlapping %>%
  filter(!is.na(Gene.Name) & Gene.Name != "") %>%
  mutate(direction = case_when(
    meth.diff > 24 ~ "hyper",
    meth.diff < -24 ~ "hypo"
  )) %>%
  distinct(Gene.Name, direction)

gene_direction %>%
  count(direction)
#hyper 13
#hypo  5

######### same thing as above but for promoter associated ! #########

setwd("~/Downloads")
CvM <- read.csv("CvM_annotated_promoters.csv")
CvS <- read.csv("CvS_annotated_promoters.csv")
MvS <- read.csv("MvS_annotated_promoters.csv")
overlapping <- read.csv("overlapping_annotated_promoters.csv")

######## find unique gene counts ########
# CvM
CvM_unique_genes <- unique(CvM$Promoter.Gene.Name)
length(CvM_unique_genes)
#565 - 1 for NA entry so 564 unique genes

# CvS
CvS_unique_genes <- unique(CvS$Promoter.Gene.Name)
length(CvS_unique_genes)
#511 - 1 so 510 unique genes

# MvS
MvS_unique_genes <- unique(MvS$Promoter.Gene.Name)
length(MvS_unique_genes)
#89 - 1 so 88 unique genes

# overlapping
overlapping_unique_genes <- unique(overlapping$Promoter.Gene.Name)
length(overlapping_unique_genes)
#44 - 1 so 43 unique

###### CvM #######
gene_direction <- CvM %>%
  filter(!is.na(Promoter.Gene.Name)) %>%
  mutate(direction = case_when(
    meth.diff > 24 ~ "hyper",
    meth.diff < -24 ~ "hypo"
  )) %>%
  distinct(Promoter.Gene.Name, direction)

gene_direction %>%
  count(direction)
#hyper 278
#hypo 306
#caveat: some genes might be hyper AND hypomethylated

######## CvS ########
gene_direction <- CvS %>%
  filter(!is.na(Promoter.Gene.Name)) %>%
  mutate(direction = case_when(
    meth.diff > 24 ~ "hyper",
    meth.diff < -24 ~ "hypo"
  )) %>%
  distinct(Promoter.Gene.Name, direction)

gene_direction %>%
  count(direction)
#hyper 249
#hypo 287

####### MvS #########
gene_direction <- MvS %>%
  filter(!is.na(Promoter.Gene.Name)) %>%
  mutate(direction = case_when(
    meth.diff > 24 ~ "hyper",
    meth.diff < -24 ~ "hypo"
  )) %>%
  distinct(Promoter.Gene.Name, direction)

gene_direction %>%
  count(direction)
#hyper 43
#hypo 47

###### overlapping #########
gene_direction <- overlapping %>%
  filter(!is.na(Promoter.Gene.Name) & Promoter.Gene.Name != "") %>%
  mutate(direction = case_when(
    meth.diff > 24 ~ "hyper",
    meth.diff < -24 ~ "hypo"
  )) %>%
  distinct(Promoter.Gene.Name, direction)

gene_direction %>%
  count(direction)
#hyper 32
#hypo  10



#---------------------------------------------------- DEG MATCHING AND GO ANALYSIS -----------------------------------------------------

#now try to match DMGs with DEGs

#read in DMG and DEG files
setwd("~/Downloads")
CvS_nonoverlap <- read.csv("SHvC-DMGs.csv")
CvM_nonoverlap <- read.csv("MHvC-DMGs.csv")
MvS_nonoverlap <- read.csv("MHvSH-DMGs.csv")
overlapping_CMCS <- read.csv("overlapping-DMGs.csv")

#DEGs
setwd("~/Desktop/Tagseq")
library(readxl)
DEGs_CvM <- read_excel("ALL-sig-annotated-DEGs.xlsx", sheet = 1)
DEGs_CvS <- read_excel("ALL-sig-annotated-DEGs.xlsx", sheet = 2)
DEGs_MvS <- read_excel("ALL-sig-annotated-DEGs.xlsx", sheet = 3)

#see if column names match up 
colnames(DEGs_CvM)
colnames(DEGs_CvS)
colnames(DEGs_MvS)
colnames(CvS_nonoverlap)
colnames(CvM_nonoverlap)
colnames(MvS_nonoverlap)
colnames(overlapping_CMCS)

#rework column names and remove empty columns specifically for CvS sheet in DEGs file
colnames(DEGs_CvS) <- gsub("\\.\\.\\..*", "", colnames(DEGs_CvS))
DEGs_CvS <- DEGs_CvS[, colnames(DEGs_CvS) != ""]

#make sure it worked
colnames(DEGs_CvS)

#want to search through nonoverlapping files for DEGs

#match by gene name OR gene description because could potentially be differing gene names

#make a smaller DEG table without trin ID, GO term..., 
#trim white space R read in when read in excel file and white space when reading in gene names and descriptions in cells
prepare_DEGs <- function(DEGs) {
  DEGs_text <- DEGs[, c("GeneName", "GeneDesc", "log2FoldChange", "padj")]
  DEGs_text$GeneName <- trimws(DEGs_text$GeneName)
  DEGs_text$GeneDesc <- trimws(DEGs_text$GeneDesc)
  return(DEGs_text)
}

#run function on each sheet
DEGs_CvM_text <- prepare_DEGs(DEGs_CvM)
DEGs_CvS_text <- prepare_DEGs(DEGs_CvS)
DEGs_MvS_text <- prepare_DEGs(DEGs_MvS)

#find matching rows in nonoverlap DMC file - match gene name
match_DMC_DEGs_by_gene <- function(DMC_file, DEGs_text) {
  
  #match DMC gene names to DEGs
  matches <- DMC_file[
    !is.na(DMC_file$Gene.Name) & 
      DMC_file$Gene.Name %in% DEGs_text$GeneName,
  ]
  
  #merge tables
  DMC_DEG_matches <- merge(
    matches,
    DEGs_text,
    by.x = "Gene.Name",
    by.y = "GeneName",
    all = FALSE
  )
  
  return(DMC_DEG_matches)
}

#run function for each comparison

#CvS DMCs against all DEG comparisons
CvS_DMC_CvS_DEG_matches <- match_DMC_DEGs_by_gene(CvS_nonoverlap, DEGs_CvS_text)
CvS_DMC_CvM_DEG_matches <- match_DMC_DEGs_by_gene(CvS_nonoverlap, DEGs_CvM_text)
CvS_DMC_MvS_DEG_matches <- match_DMC_DEGs_by_gene(CvS_nonoverlap, DEGs_MvS_text)

#CvM
CvM_DMC_CvM_DEG_matches <- match_DMC_DEGs_by_gene(CvM_nonoverlap, DEGs_CvM_text)
CvM_DMC_CvS_DEG_matches <- match_DMC_DEGs_by_gene(CvM_nonoverlap, DEGs_CvS_text)
CvM_DMC_MvS_DEG_matches <- match_DMC_DEGs_by_gene(CvM_nonoverlap, DEGs_MvS_text)

#MvS
MvS_DMC_MvS_DEG_matches <- match_DMC_DEGs_by_gene(MvS_nonoverlap, DEGs_MvS_text)
MvS_DMC_CvS_DEG_matches <- match_DMC_DEGs_by_gene(MvS_nonoverlap, DEGs_CvS_text)
MvS_DMC_CvM_DEG_matches <- match_DMC_DEGs_by_gene(MvS_nonoverlap, DEGs_CvM_text)

#overlapping
CMCS_DMC_CvM_DEG_matches <- match_DMC_DEGs_by_gene(overlapping_CMCS, DEGs_CvM_text)
CMCS_DMC_CvS_DEG_matches <- match_DMC_DEGs_by_gene(overlapping_CMCS, DEGs_CvS_text)
CMCS_DMC_MvS_DEG_matches <- match_DMC_DEGs_by_gene(overlapping_CMCS, DEGs_MvS_text)

#look at matches
#CvS DMCs against all DEG comparisons
CvS_DMC_CvS_DEG_matches
#tram1
CvS_DMC_CvM_DEG_matches
#fbxo9
CvS_DMC_MvS_DEG_matches
#none

#CvM
CvM_DMC_CvM_DEG_matches
#none
CvM_DMC_CvS_DEG_matches
#grpel1
CvM_DMC_MvS_DEG_matches
#none

#MvS
MvS_DMC_MvS_DEG_matches
#none
MvS_DMC_CvS_DEG_matches
#none
MvS_DMC_CvM_DEG_matches
#none

#overlapping
CMCS_DMC_CvM_DEG_matches
CMCS_DMC_CvS_DEG_matches
CMCS_DMC_MvS_DEG_matches
#none for all


#now try to match overlapping DMCs with DEGs

#DEGs
setwd("~/Desktop/Tagseq")
library(readxl)
DEGs_CvM <- read_excel("ALL-sig-annotated-DEGs.xlsx", sheet = 1)
DEGs_CvS <- read_excel("ALL-sig-annotated-DEGs.xlsx", sheet = 2)
DEGs_MvS <- read_excel("ALL-sig-annotated-DEGs.xlsx", sheet = 3)

#see if column names match up 
colnames(DEGs_CvM)
colnames(DEGs_CvS)
colnames(DEGs_MvS)
colnames(CMCS)

#rework column names and remove empty columns specifically for CvS sheet in DEGs file
colnames(DEGs_CvS) <- gsub("\\.\\.\\..*", "", colnames(DEGs_CvS))
DEGs_CvS <- DEGs_CvS[, colnames(DEGs_CvS) != ""]

#make sure it worked
colnames(DEGs_CvS)

#want to search through overlapping files for DEGs

#match by gene name

#make a smaller DEG table without trin ID, GO term..., 
#trim white space R read in when read in excel file and white space when reading in gene names and descriptions in cells
prepare_DEGs <- function(DEGs) {
  DEGs_text <- DEGs[, c("GeneName", "GeneDesc", "log2FoldChange", "padj")]
  DEGs_text$GeneName <- trimws(DEGs_text$GeneName)
  DEGs_text$GeneDesc <- trimws(DEGs_text$GeneDesc)
  return(DEGs_text)
}

#run function on each sheet
DEGs_CvM_text <- prepare_DEGs(DEGs_CvM)
DEGs_CvS_text <- prepare_DEGs(DEGs_CvS)
DEGs_MvS_text <- prepare_DEGs(DEGs_MvS)

#make list to loop through easier
DEG_sets <- list(
  CvS = DEGs_CvS_text,
  CvM = DEGs_CvM_text,
  MvS = DEGs_MvS_text)

#find matching rows in overlap DMC file - match gene name
match_DMC_DEGs_by_gene <- function(DMC_file, DEGs_text) {
  
  #match DMC gene names to DEGs
  matches <- DMC_file[
    !is.na(DMC_file$gene) & 
      DMC_file$gene %in% DEGs_text$GeneName,
  ]
  
  #merge tables
  DMC_DEG_matches <- merge(
    matches,
    DEGs_text,
    by.x = "gene",
    by.y = "GeneName",
    all = FALSE)
  
  return(DMC_DEG_matches)
}

#run function for each DEG set
CMCS_DEG_matches <- lapply(DEG_sets, function(deg_df) {
  match_DMC_DEGs_by_gene(CMCS, deg_df)
})

#check for matches
sapply(CMCS_DEG_matches, nrow)
#none :(

#now try to match promoter potential genes to DEGs following the same workflow above

#read in DEG files
setwd("~/Desktop/Tagseq")
library(readxl)
DEGs_CvM <- read_excel("ALL-sig-annotated-DEGs.xlsx", sheet = 1)
DEGs_CvS <- read_excel("ALL-sig-annotated-DEGs.xlsx", sheet = 2)
DEGs_MvS <- read_excel("ALL-sig-annotated-DEGs.xlsx", sheet = 3)


#see if column names match up 
colnames(DEGs_CvM)
colnames(DEGs_CvS)
colnames(DEGs_MvS)
colnames(CvS_FULLY_annotated)
colnames(CvM_FULLY_annotated)
colnames(MvS_FULLY_annotated)
colnames(CvS_CvM_overlapping)


#rework column names and remove empty columns specifically for CvS sheet in DEGs file
colnames(DEGs_CvS) <- gsub("\\.\\.\\..*", "", colnames(DEGs_CvS))
DEGs_CvS <- DEGs_CvS[, colnames(DEGs_CvS) != ""]

#make sure it worked
colnames(DEGs_CvS)

#want to search through overlapping files for DEGs

#match by gene name

#make a smaller DEG table without trin ID, GO term..., 
#trim white space R read in when read in excel file and white space when reading in gene names and descriptions in cells
prepare_DEGs <- function(DEGs) {
  DEGs_text <- DEGs[, c("GeneName", "GeneDesc", "log2FoldChange", "padj")]
  DEGs_text$GeneName <- trimws(DEGs_text$GeneName)
  DEGs_text$GeneDesc <- trimws(DEGs_text$GeneDesc)
  return(DEGs_text)
}

#run function on each sheet
DEGs_CvM_text <- prepare_DEGs(DEGs_CvM)
DEGs_CvS_text <- prepare_DEGs(DEGs_CvS)
DEGs_MvS_text <- prepare_DEGs(DEGs_MvS)

#make list to loop through easier
DEG_sets <- list(
  CvS = DEGs_CvS_text,
  CvM = DEGs_CvM_text,
  MvS = DEGs_MvS_text
)

#find matching rows in overlap DMC file - match gene name
match_DMC_DEGs_by_promoter <- function(DMC_file, DEGs_text) {
  
  # keep only rows that have promoter gene names
  has_promoter <- !is.na(DMC_file$Promoter.Gene.Name)
  DMC_sub <- DMC_file[has_promoter, ]
  
  # identify rows where ANY promoter gene is a DEG
  keep <- sapply(DMC_sub$Promoter.Gene.Name, function(x) {
    genes <- unlist(strsplit(x, ";"))
    genes <- trimws(genes)
    any(genes %in% DEGs_text$GeneName)
  })
  
  DMC_sub <- DMC_sub[keep, ]
  
  return(DMC_sub)
}

#run function for each DEG set

#### CvS DMCs against CvS then CvM then MvS DEGs
CvS_nonoverlap_CvS_DEG_promoter <- match_DMC_DEGs_by_promoter(CvS_FULLY_annotated, DEGs_CvS_text)
CvS_nonoverlap_CvM_DEG_promoter <- match_DMC_DEGs_by_promoter(CvS_FULLY_annotated, DEGs_CvM_text)
CvS_nonoverlap_MvS_DEG_promoter <- match_DMC_DEGs_by_promoter(CvS_FULLY_annotated, DEGs_MvS_text)

#look at results
CvS_nonoverlap_CvS_DEG_promoter
#3 DEGs
CvS_nonoverlap_CvM_DEG_promoter
#2 DEGs
CvS_nonoverlap_MvS_DEG_promoter
#none

###### CvM DMCs against CvM then CvS then MvS DEGs
CvM_nonoverlap_CvM_DEG_promoter <- match_DMC_DEGs_by_promoter(CvM_FULLY_annotated, DEGs_CvM_text)
CvM_nonoverlap_CvS_DEG_promoter <- match_DMC_DEGs_by_promoter(CvM_FULLY_annotated, DEGs_CvS_text)
CvM_nonoverlap_MvS_DEG_promoter <- match_DMC_DEGs_by_promoter(CvM_FULLY_annotated, DEGs_MvS_text)

#look at results
CvM_nonoverlap_CvM_DEG_promoter
#none
CvM_nonoverlap_CvS_DEG_promoter
#1 DEG
CvM_nonoverlap_MvS_DEG_promoter
#none

####### MvS DMCs against MvS then CvS then CvM DEGs
MvS_nonoverlap_MvS_DEG_promoter <- match_DMC_DEGs_by_promoter(MvS_FULLY_annotated, DEGs_MvS_text)
MvS_nonoverlap_CvS_DEG_promoter <- match_DMC_DEGs_by_promoter(MvS_FULLY_annotated, DEGs_CvS_text)
MvS_nonoverlap_CvM_DEG_promoter <- match_DMC_DEGs_by_promoter(MvS_FULLY_annotated, DEGs_CvM_text)

#look at results
MvS_nonoverlap_MvS_DEG_promoter
#none
MvS_nonoverlap_CvS_DEG_promoter
#1 DEG
MvS_nonoverlap_CvM_DEG_promoter
#1 DEG
#same one


####### CvS vs CvM overlapping DMCs against CvS DEGs then CvM DEGs
CvS_CvM_overlapping_DEG_promoter_CvS <- match_DMC_DEGs_by_promoter(CvS_CvM_overlapping, DEGs_CvS_text)
CvS_CvM_overlapping_DEG_promoter_CvM <- match_DMC_DEGs_by_promoter(CvS_CvM_overlapping, DEGs_CvM_text)
CvS_CvM_overlapping_DEG_promoter_CvM <- match_DMC_DEGs_by_promoter(CvS_CvM_overlapping, DEGs_MvS_text)

#look at results
CvS_CvM_overlapping_DEG_promoter_CvS
CvS_CvM_overlapping_DEG_promoter_CvM
#none for either


######## running GO analysis on DMGs ##########
  #https://royalsocietypublishing.org/rspb/article/289/1974/20220670/79449/Thermal-regime-during-parental-sexual-maturation
  #they followed this pipeline: https://github.com/enormandeau/go_enrichment and it seems just like normal GO analysis so just going to try it like how we usually do it in lab

#going to copy and paste here to keep everything in one file
#need to do some "preprocessing" of DMC table to make this run properly
#need to develop DMG to GO match table 

#make file to match gene name to GO term with FINAL_ANNOTATION_CORRECT_by_gene.xlsx
setwd("~/Desktop")
annotations <- read.csv("FINAL_ANNOTATION_CORRECT_by_gene.csv")

#first keep only GeneName and GOTermList columns
annotations <- annotations[, c("GeneName", "GOTermList")]

#now remove an NAs
annotations <- annotations[!is.na(annotations$GOTermList), ]
sum(is.na(annotations$GOTermList))
annotations <- annotations[!is.na(annotations$GeneName), ]
sum(is.na(annotations$GeneName))

library(dplyr)

#read in DMGs
setwd("~/Desktop/Chp1/Methylation/Working files/")
CvS_nonoverlap <- read.csv("SHvC-DMGs.csv")
CvM_nonoverlap <- read.csv("MHvC-DMGs.csv")
MvS_nonoverlap <- read.csv("MHvSH-DMGs.csv")
overlapping_CMCS <- read.csv("overlapping-DMGs.csv")

#make sure they all have Gene.Name as column name for gene
colnames(CvS_nonoverlap)
colnames(CvM_nonoverlap)
colnames(MvS_nonoverlap)
colnames(overlapping_CMCS)

# keep only DMGs that have names in Gene.Name column
keep_genes <- function(df) {
  df %>% filter(!is.na(Gene.Name), Gene.Name != "")
}

CvS <- keep_genes(CvS_nonoverlap)
#266
CvM <- keep_genes(CvM_nonoverlap)
#309
MvS <- keep_genes(MvS_nonoverlap)
#40
overlapping <- keep_genes(overlapping_CMCS)
#40

#clean gene names the same way you cleaned the names in the file above
clean_gene_names <- function(df) {
  df %>%
    filter(!is.na(Gene.Name), Gene.Name != "") %>%
    mutate(Gene.Name = sub(" .*", "", Gene.Name))
}

CvS <- clean_gene_names(CvS)
CvM <- clean_gene_names(CvM)
MvS <- clean_gene_names(MvS)
overlapping <- clean_gene_names(overlapping)

#see how many DMG genes match GO table
sum(CvS$Gene.Name %in% annotations$GeneName)
#104

#find missing genes
missing_genes <- CvS$Gene.Name[!CvS$Gene.Name %in% annotations$GeneName]
head(missing_genes, 30)
length(missing_genes)
#162

#ignore capitalization
sum(toupper(CvS$Gene.Name) %in% toupper(annotations$GeneName))
#124 so kinda helped

CvS$Gene.Name <- toupper(CvS$Gene.Name)
CvM$Gene.Name <- toupper(CvM$Gene.Name)
MvS$Gene.Name <- toupper(MvS$Gene.Name)
overlapping$Gene.Name <- toupper(overlapping$Gene.Name)

#Clean gene names and try again
annotations$GeneName <- sub(" .*", "", annotations$GeneName)

#Merge duplicate genes created by cleaning
annotations <- annotations %>%
  group_by(GeneName) %>%
  summarise(GOTermList = paste(unique(unlist(strsplit(GOTermList, ";"))), collapse = ";")) %>%
  ungroup()

#Check it worked
head(annotations$GeneName, 20)
table(table(annotations$GeneName))

#export file
write.table(annotations, "GO_annotations.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#nothing to be able to match why there are missing genes so running on genes that have GO terms

#write csv for GO input
write.csv(CvS, "CvS_DMG_GO.csv", row.names = FALSE)
write.csv(CvM, "CvM_DMG_GO.csv", row.names = FALSE)
write.csv(MvS, "MvS_DMG_GO.csv", row.names = FALSE)
write.csv(overlapping, "FINAL_overlapping_all_DMC_GO.csv", row.names = FALSE)



###### run GO #####

#set to where all your files are (the input files below)
setwd("~/Desktop/Chp1/Methylation/New-GO")

##### had to go into all DMC GO input files and change Gene.Name column name to GeneName to get it to work ######
# Edit these to match your data file names: 
input="CvS_DMG_GO.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="GO_annotations.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="MF" # either MF, or BP, or CC
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25) # threshold for merging similar (gene-sharing) terms. See README for details.

#nothing significant for CvS
#nothing significant for CvM
#not enough genes for MvS or overlapping


##### now GO for potential promoter genes ####

##### Go into those files below before reading them in and change Promoter.Gene.Name to just GeneName!!! #########
#tried to do it in R but it didnt work so just doing it in excel

#read in DMC files that have just gene annotations (with promoter annotations)
setwd("~/Desktop/Chp1/Methylation/Working files/")
CvS_nonoverlap <- read.csv("CvS_annotated_promoters.csv")
CvM_nonoverlap <- read.csv("CvM_annotated_promoters.csv")
MvS_nonoverlap <- read.csv("MvS_annotated_promoters.csv")
overlapping_CvS_CvM <- read.csv("overlapping_annotated_promoters.csv")

#make sure they all have Gene.Name as column name for gene
colnames(CvS_nonoverlap)
colnames(CvM_nonoverlap)
colnames(MvS_nonoverlap)
colnames(overlapping_CvS_CvM)

#pick only the first gene in semicolon separated list - Promoter.gene.name is a list of gene names
first_gene <- function(x) {
  sapply(strsplit(x, ";"), `[`, 1)
}

# keep only DMCs that mapped to genes - want Gene.Name column
#also choose the first gene in the list of promoter potential genes
keep_genes <- function(df) {
  df %>%
    filter(!is.na(Gene.Name), Gene.Name != "") %>%  # remove NA or empty
    mutate(Gene.Name = first_gene(Gene.Name)) %>%  # take first gene
    filter(Gene.Name != "")  # remove any leftover empty strings
}

#apply function
CvS <- keep_genes(CvS_nonoverlap)
CvM <- keep_genes(CvM_nonoverlap)
MvS <- keep_genes(MvS_nonoverlap)
overlapping <- keep_genes(overlapping_CvS_CvM)

#clean gene names
clean_gene_names <- function(df) {
  df %>%
    filter(!is.na(Gene.Name), Gene.Name != "") %>%
    mutate(Gene.Name = sub(" .*", "", Gene.Name))
}

#apply function
CvS <- clean_gene_names(CvS)
CvM <- clean_gene_names(CvM)
MvS <- clean_gene_names(MvS)
overlapping <- clean_gene_names(overlapping)

#see how many DMG genes match GO table
sum(CvS$Gene.Name %in% annotations$GeneName)
#240

#find missing genes
missing_genes <- CvS$Gene.Name[!CvS$Gene.Name %in% annotations$GeneName]
head(missing_genes, 30)
length(missing_genes)
#331

#ignore capitalization
CvS$Gene.Name <- toupper(CvS$Gene.Name)
CvM$Gene.Name <- toupper(CvM$Gene.Name)
MvS$Gene.Name <- toupper(MvS$Gene.Name)
overlapping$Gene.Name <- toupper(overlapping$Gene.Name)

#write csv for GO input
write.csv(CvS, "CvS_promoter_DMG_GO.csv", row.names = FALSE)
write.csv(CvM, "CvM_promoter_DMG_GO.csv", row.names = FALSE)
write.csv(MvS, "MvS_promoter_DMG_GO.csv", row.names = FALSE)
write.csv(overlapping, "overlapping_promoter_DMG_GO.csv", row.names = FALSE)


###### run GO #####

#set to where all your files are (the input files below)
setwd("~/Desktop/Chp1/Methylation/New-GO")

input="CvM_promoter_DMG_GO.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="GO_annotations.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="MF" # either MF, or BP, or CC
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25) # threshold for merging similar (gene-sharing) terms. See README for details.

#did not run for any comparison

