#based off: https://www.sciencedirect.com/science/article/pii/S0044848625007240#bb0010
  #Unveiling the epigenetic signatures: DNA methylation patterns in the social hierarchies of male giant freshwater prawns (Macrobrachium rosenbergii)
    #Liping Li, Jiongying Yu, Zhenglong Xia, Quanxin Gao, Qiongying Tang, Shaokui Yi 

#example code here: https://code.google.com/archive/p/methylkit/
#here too: https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationSeq/Seq_Tutorial.html

#Bioconduction methylKit: https://www.bioconductor.org/packages/devel/bioc/vignettes/methylKit/inst/doc/methylKit.html

setwd("/home/aubaxb001/methylation_data/trimmed/")

#BiocManager::install("methylKit")

#methylkit specifically designed for RBBS but can use it on WGBS too
library(methylKit)

#read in output files from methyl_input.sh script with the headers already added

# List of methylKit-formatted files
files_CvS <- list(
   ("/home/aubaxb001/methylation_data/trimmed/C_T1_03_FF_1_val_1_bismark_bt2_pe.deduplicated_methylKit.txt"),
   ("/home/aubaxb001/methylation_data/trimmed/C_T1_14_FF_1_val_1_bismark_bt2_pe.deduplicated_methylKit.txt"),
   ("/home/aubaxb001/methylation_data/trimmed/C_T1_16_FF_1_val_1_bismark_bt2_pe.deduplicated_methylKit.txt"),
   ("/home/aubaxb001/methylation_data/trimmed/C_T5_10_FF_1_val_1_bismark_bt2_pe.deduplicated_methylKit.txt"),
   ("/home/aubaxb001/methylation_data/trimmed/C_T5_12_FF_1_val_1_bismark_bt2_pe.deduplicated_methylKit.txt"),
   ("/home/aubaxb001/methylation_data/trimmed/S_T1_01_FF_1_val_1_bismark_bt2_pe.deduplicated_methylKit.txt"),
   ("/home/aubaxb001/methylation_data/trimmed/S_T1_07_FF_1_val_1_bismark_bt2_pe.deduplicated_methylKit.txt"),
   ("/home/aubaxb001/methylation_data/trimmed/S_T1_12_FF_1_val_1_bismark_bt2_pe.deduplicated_methylKit.txt"),
   ("/home/aubaxb001/methylation_data/trimmed/S_T3_07_FF_1_val_1_bismark_bt2_pe.deduplicated_methylKit.txt"),
   ("/home/aubaxb001/methylation_data/trimmed/S_T3_11_FF_1_val_1_bismark_bt2_pe.deduplicated_methylKit.txt")
)

#pipeline generation
  #read help page for methRead and pipeline thing for more info
  #this is just saying what information is in what column - ex the coverage information is in column 5
#fraction is from 0-100 in input files so FALSE
#chr.col = 2, start.col = 3, end.col = 3, coverage.col = 5, freqC.col = 6, 
  #start and end at the same spot because per base methylation not regional
#strand column = 4 but just put + as placeholder since dont have strand information
pipeline <- list(fraction = FALSE, chr.col = 2, start.col = 3, end.col = 3, coverage.col = 5, freqC.col = 6, strand.col = 4)
  
#read in files into methylRawList object making sure the other parameters are filled in correctly
#can only do pairwise comparisons / cant look at more than 2 groups at once
#using files you manually added headers to so have to manually enter the "pipeline" to follow
#this is just saying what information is in what column - ex the coverage information is in column 5
#name assembly whatever you want doesnt matter
#0 for control, 1 for treatment
#default says there is a header in the input files, default tab delimited format
myobjCvS <- methRead(files_CvS,
                  sample.id=list("C_T1_03","C_T1_14","C_T1_16","C_T5_10", "C_T5_12", "S_T1_01", "S_T1_07", "S_T1_12", "S_T3_07", "S_T3_11"),
                  assembly="Lrhomboides",
                  treatment=c(0,0,0,0,0,1,1,1,1,1),
                  pipeline = pipeline,
                  mincov = 10,
                  dbtype = "tabix"
                  )

#ckeck number of samples
myobjCvS
#should sat 10 methylRaw objects - one entry per sample


##### get histogram of methylation % per sample   ######
  #going to have to loop through all samples since getMethylationstats works per sample


#set folder to save pdfs and stats
pdf_dir <- "/home/aubaxb001/methylation_data/trimmed/methylation_pdfs"
stats_dir <- "/home/aubaxb001/methylation_data/trimmed/methylation_stats"

#create those directories for output pdfs and output stats to organize yourself
dir.create(pdf_dir)
dir.create(stats_dir)

#loop through samples
#both.strands = true IF you want stats for both strands separately - not the case for us since start and end are the same since per base methylation not regional
for (i in seq_along(myobjCvS)) {
  
  sample_name <- myobjCvS[[i]]@sample.id #get sample name for output file names
  
  #PDF
  pdf_file <- file.path(pdf_dir, paste0("methylation_", sample_name, ".pdf")) #output file name will be methylation_samplename.pdf
  pdf(pdf_file, width = 8, height = 6) #open new pdf device
  getMethylationStats(myobjCvS[[i]], plot = TRUE, both.strands = FALSE) #histogram of methyl % for each sample
  dev.off() #close plotting device and saves file to pdf_dir
  
  #numeric stats
  stats_file <- file.path(stats_dir, paste0("methylation_stats_", sample_name, ".txt")) #put stats in file called methylation_stats_samplename.txt
  
  #capture output stats into file per sample
  stats <- capture.output(
    getMethylationStats(myobjCvS[[i]], plot = FALSE, both.strands = FALSE) #plot=FALSE to just get the stats
  )
 writeLines(stats, stats_file)
  
  cat("Finished processing sample:", sample_name, "\n") #check progress
}


#coverage stats and histograms per sample

#make directories for this like before to stay organized
coverage_pdfs <- "/home/aubaxb001/methylation_data/trimmed/coverage_pdfs"
coverage_stats <- "/home/aubaxb001/methylation_data/trimmed/coverage_stats"

#make the directories
dir.create(coverage_pdfs, showWarnings = FALSE)
dir.create(coverage_stats, showWarnings = FALSE)

#loop through samples
for (i in seq_along(myobjCvS)) {
  sample_name <- myobjCvS[[i]]@sample.id
  
  #PDF
  pdf_file <- file.path(coverage_pdfs, paste0("coverage_", sample_name, ".pdf"))
  pdf(pdf_file, width = 8, height = 6)
  getCoverageStats(myobjCvS[[i]], plot = TRUE, both.strands = FALSE)
  dev.off()
  
  #stats
  stats_file <- file.path(coverage_stats, paste0("coverage_stats_", sample_name, ".txt"))
  stats_out <- capture.output(getCoverageStats(myobjCvS[[i]], plot = FALSE, both.strands = FALSE))
  writeLines(stats_out, stats_file)
  
  cat("Finished coverage stats for sample:", sample_name, "\n")
}


#### FILTERING #####

#did lo.perc = 5 instead of lo.count = 10 to try to retain more sites and decreased hi.perc to 99.5 instead of 99.9
myobjCvS.filt <- filterByCoverage(myobjCvS, lo.perc = 5, hi.perc = 99.5)

#######         Normalizing     ########
#normalize coverage values bw samples
  #uses a scaling factor derived from differences bw the median of the coverage distributions
myobjCvS.filt.norm <- normalizeCoverage(myobjCvS.filt, method = "median")

########        Make single table for more analysis     ########
#merge all samples with base pair locations that are covered in all samples
#destrand = FALSE is default
#destrand = TRUE - use when just want to focus on CpG methylation
#min.per.group - default says that only CpGs covered in ALL samples will be output, want to relax this to have some sites still be considered
#save.db = TRUE to save files while uniting to working directory instead of in R - will take up too much storage in R if you dont do this and will fail
  #also needs ~18gb of data from your laptop so make sure theres enough open storage
methCvS <- unite(myobjCvS.filt.norm, min.per.group = 5L, destrand = FALSE, save.db = TRUE)

methCvS

#get percent methylation data
meth_perc <- percMethylation(methCvM)

#save as csv file
write.csv(meth_perc, "CvS_meth_perc_all_samples.csv", row.names = FALSE)


####### PLOTTING ###########

###### check correlation bw samples #########

pdf("/home/aubaxb001/methylation_data/trimmed/methylKit_db/CvS_correlation.pdf", width = 8, height = 6)
getCorrelation(methCvS, plot = TRUE, method = "pearson")
dev.off()

#get correlation matrix to make pretty in R later if we want
correlation_matrix <-getCorrelation(methCvS, plot = FALSE)
write.csv(correlation_matrix, file = "/home/aubaxb001/methylation_data/trimmed/methylKit_db/CvS_correlation_matrix.csv", row.names = TRUE)


####### visualize data in a dendrogram using hierarchical clustering of distance measures from each samples percenage methylation #######
#clustering used to group dta points by similarity

pdf("/home/aubaxb001/methylation_data/trimmed/methylKit_db/CvS_clustered.pdf", width = 8, height = 6)
clusterSamples(methCvS, dist = "correlation", method = "ward", plot = TRUE)
dev.off()


####### PCA screeplot #######

pdf("/home/aubaxb001/methylation_data/trimmed/methylKit_db/CvS_PCA.pdf", width = 8, height = 6)
PCASamples(methCvS, screeplot = TRUE)
dev.off()

###### normal PCA #######

pdf("/home/aubaxb001/methylation_data/trimmed/methylKit_db/CvS_normalPCA.pdf", width = 8, height = 6)
PCASamples(methCvS, screeplot = FALSE, adj.lim = c(0.5,0.1), obj.return = TRUE)
dev.off()

#save PC scores to make better PCA
pca_obj <- PCASamples(methCvS, screeplot = FALSE, adj.lim = c(0.5,0.1), obj.return = TRUE)

#extract sample scores
pca_scores <- as.data.frame(pca_obj$x)

#print scores to copy and paste into text file
pca_scores

#levels it outputs are $x $sdev $rotation $center $scale

#extract standard deviation score and square it
eigenvalues = variance of each PC
eigs <- pca_obj$sdev^2

#proportion of variance explained by each PC
var_explained <- eigs / sum(eigs)

#cummulative variance
cum_var_explained <- cumsum(var_explained)

#export to txt file to make better PCA
variance_df <- data.frame(
  PC = paste0("PC", seq_along(eigs)),
  Variance_Explained = var_explained,
  Cumulative_Variance = cum_var_explained)

#save to csv
write.csv(variance_df, "PCA_variance_summary_CvS.csv", row.names = FALSE)

####### calculate differential methylation stats bw 2 groups ########
  #this is comparing methylation levels bw treatment and control

#overdispersion option will make sure no false positives - basic overdispersion correction
#benjamini hochberg correction for p-values adjustments
#want to use logistic regression (not fishers)
#McCullagh and Nelder (MN) overdispersion correction applied
#Benjamini hochberg correction to find adjusted p values
#Chisq = default but want to use this (not fishers)

myDiffCvS <- calculateDiffMeth(methCvS, overdispersion = "MN", adjust = "BH", test = "Chisq")
head(myDiffCvS)

#save output to csv file
myDiffCvS_df <- getData(myDiffCvS)
write.csv(myDiffCvS_df, "DiffMeth_CvS.csv", row.names = FALSE)

#get differnetially methylated regions with 25% difference and q-value < 0.05
  #lots of other literature uses 25% and 0.01 but want to stay consistent with the same significance threshold as other datasets in chapter
#type = all indicates I want ALL differentially methylated bases - not just hypo or hyper methylated regions
myDiff25q1 <- getMethylDiff(myDiffCvS, difference = 25, qvalue = 0.05, type = "all")


#save output
myDiff25q1_df <- getData(myDiff25q1)
write.csv(myDiff25q1_df, "DMRs_CvS.csv", row.names = FALSE)

#go by hypo
myDiff25q1_hypo <- getMethylDiff(myDiffCvS, difference = 25, qvalue = 0.05, type = "hypo")

#save output
myDiff25q1_df <- getData(myDiff25q1_hypo)
write.csv(myDiff25q1_df, "DMRs_CvS_hypo.csv", row.names = FALSE)

#go by hyper
myDiff25q1_hyper <- getMethylDiff(myDiffCvS, difference = 25, qvalue = 0.05, type = "hyper")

#save output
myDiff25q1_df <- getData(myDiff25q1_hyper)
write.csv(myDiff25q1_df, "DMRs_CvS_hyper.csv", row.names = FALSE)

##### tiling window analysis just in case #######
  #groups methylation into tiling windows - contiguous regions of the genome of a fixed size
  #within each window, the function sums up the counts of methyl and unmethyl reds to get summary of methyl profile per window instead of per base
  #window size = size of tiling window
  #step size = if the same as window size, start of each window = previous window end, no overlapping windows
  #large windows (like 1-2000 give lower res but smoother data, detect broad changes)
  #small windows (like 500-1000bp common for genome wide tiling in verts)
  #SMALL windows (like 100-500bp) captude fine scale methylation changes - work best for CpG dense regions like promoters or CpG islands

#The tilling function adds up C and T counts from each covered cytosine and returns a total C and T count for each tile.
#might want to set the initial per base coverage threshold to a lower value and then filter based on the number of bases (cytosines) per region.

#try different window sizes to see differences - range of all
  #github vignette for methylKit uses 1000 but lots of lit using ~200
  #should output info per window:
    #chromosome: which chromosome the window is on
    #start and end position of the window
    #total methyl counts in that window
    #total unmethyl counts in that window
    #coverage info (how many bases in the window had how many reads)
  #overall should have each row correspond to a tiling window not a single base
#tiles1000 = tileMethylCounts(methCvS, win.size = 1000, step.size = 1000, cov.bases = 3)

#tiles500 = tileMethylCounts(methCvS, win.size = 500, step.size = 500, cov.bases = 3)

#tiles200 = tileMethylCounts(methCvS, win.size = 200, step.size = 200, cov.bases = 3)

#tiles25 = tileMethylCounts(methCvS, win.size = 25, step.size = 25, cov.bases = 3)


#run differential methylation analysis on tiled data
#myDiffTiles1000 <- calculateDiffMeth(tiles1000, overdispersion = "MN", adjust = "BH", test = "Chisq")

#myDiffTiles500 <- calculateDiffMeth(tiles500, overdispersion = "MN", adjust = "BH", test = "Chisq")

#myDiffTiles200 <- calculateDiffMeth(tiles200, overdispersion = "MN", adjust = "BH", test = "Chisq")

#myDiffTiles25 <- calculateDiffMeth(tiles25, overdispersion = "MN", adjust = "BH", test = "Chisq")

#save output to csv file
#myDiffTiles1000_df <- getData(myDiffTiles1000)
#write.csv(myDiffTiles1000_df, "DiffMeth_Tiled1000_CvS.csv", row.names = FALSE)

#myDiffTiles500_df <- getData(myDiffTiles500)
#write.csv(myDiffTiles500_df, "DiffMeth_Tiled500_CvS.csv", row.names = FALSE)

#myDiffTiles200_df <- getData(myDiffTiles200)
#write.csv(myDiffTiles200_df, "DiffMeth_Tiled200_CvS.csv", row.names = FALSE)

#myDiffTiles25_df <- getData(myDiffTiles25)
#write.csv(myDiffTiles25_df, "DiffMeth_Tiled25_CvS.csv", row.names = FALSE)


#get differnetially methylated regions with 25% difference and q-value < 0.05
#lots of other literature uses 25% and 0.01 but want to stay consistent with the same significance threshold>
#type = all indicates I want ALL differentially methylated bases - not just hypo or hyper methylated regions
############Tiled1000############
#myDiff25q5_Tiles1000 <- getMethylDiff(myDiffTiles1000, difference = 25, qvalue = 0.05, type = "all")

#less strict thresholds, qvalue = 0.01, diff = 10
#myDiff10q1_Tiles1000 <- getMethylDiff(myDiffTiles1000, difference = 10, qvalue = 0.01, type = "all")

#save output
#myDiff25q5_Tiles1000_df <- getData(myDiff25q5_Tiles1000)
#write.csv(myDiff25q5_Tiles1000_df, "DMRs_Tiled1000_q5d25_CvS.csv", row.names = FALSE)

#save output
#myDiff10q1_Tiles1000_df <- getData(myDiff10q1_Tiles1000)
#write.csv(myDiff10q1_Tiles1000_df, "DMRs_Tiled1000_q1d10_CvS.csv", row.names = FALSE)


###########Tiles500###############
#myDiff25q5_Tiles500 <- getMethylDiff(myDiffTiles500, difference = 25, qvalue = 0.05, type = "all")

#less strict thresholds, qvalue = 0.01, diff = 10
#myDiff10q1_Tiles500 <- getMethylDiff(myDiffTiles500, difference = 10, qvalue = 0.01, type = "all")

#save output
#myDiff25q5_Tiles500_df <- getData(myDiff25q5_Tiles500)
#write.csv(myDiff25q5_Tiles500_df, "DMRs_Tiled500_q5d25_CvS.csv", row.names = FALSE)

#save output
#myDiff10q1_Tiles500_df <- getData(myDiff10q1_Tiles500)
#write.csv(myDiff10q1_Tiles500_df, "DMRs_Tiled500_q1d10_CvS.csv", row.names = FALSE)

#############Tiles200###############
#myDiff25q5_Tiles200 <- getMethylDiff(myDiffTiles200, difference = 25, qvalue = 0.05, type = "all")

#less strict thresholds, qvalue = 0.01, diff = 10
#myDiff10q1_Tiles200 <- getMethylDiff(myDiffTiles200, difference = 10, qvalue = 0.01, type = "all")

#save output
#myDiff25q5_Tiles200_df <- getData(myDiff25q5_Tiles200)
#write.csv(myDiff25q5_Tiles200_df, "DMRs_Tiled200_q5d25_CvS.csv", row.names = FALSE)

#save output
#myDiff10q1_Tiles200_df <- getData(myDiff10q1_Tiles200)
#write.csv(myDiff10q1_Tiles200_df, "DMRs_Tiled200_q1d10_CvS.csv", row.names = FALSE)


############Tiles25################
#myDiff25q5_Tiles25 <- getMethylDiff(myDiffTiles25, difference = 25, qvalue = 0.05, type = "all")

#less strict thresholds, qvalue = 0.01, diff = 10
#myDiff10q1_Tiles25 <- getMethylDiff(myDiffTiles25, difference = 10, qvalue = 0.01, type = "all")

#save output
#myDiff25q5_Tiles25_df <- getData(myDiff25q5_Tiles25)
#write.csv(myDiff25q5_Tiles25_df, "DMRs_Tiled25_q5d25_CvS.csv", row.names = FALSE)

#save output
#myDiff10q1_Tiles25_df <- getData(myDiff10q1_Tiles25)
#write.csv(myDiff10q1_Tiles25_df, "DMRs_Tiled25_q1d10_CvS.csv", row.names = FALSE)


##### CpG annotation #######
#this is using the genomation package: https://bioconductor.org/packages/release/bioc/html/genomation.html

#want to do this with the tiling window you made, not per base methylation (so not with the myDiff variables above)

#load genomation
#library(genomation)

#read in bed file
#absolute path to bed file
#set remove.unsual to false so it doesnt get rid of "weirdly named" contigs
#have to pick how upstream/downstream from TSS you want to look
  #https://www.mdpi.com/1422-0067/26/7/3224 used as references (chose 2500bp)
#unique.prom = TRUE to just get one promoter per gene

#SOOOO
  #read bed file, add 2500bp up- and downstream of TSS as part of promoter region, do not overlap promoters
  #should output gene.obj = list of genomic features split into: promoters, exon,s introns, gene bodies
  #this can then be used to map methylation data onto gene features
#gene.obj = readTranscriptFeatures("/home/aubaxb001/Lrhomboides_genes_BED12.bed", remove.unusual = FALSE, up.flank = 2500, down.flank = 2500, unique.prom - TRUE)

#saveRDS(gene.obj, file = "gene-obj.rds")

#this annotates methylation regions with respect to gene parts
  #for each methylation region, the function will tell you which part of a gene it overlaps with (ex promoter, exon, intron...)
#should output dataframe or GRanges object where each row corresponds to a methylation region
  #columns = chr, overlapping gene name or ID, which gene part it overlaps (exon, intron...)
  #so pretty much saying this methylated region falls in the promoter of GeneX or this region is in an intron of GeneY
#annotateWithGeneParts(as(myDiff25q1_Tiles1000, "GRanges"), gene.obj)

#this annotates methylation regions with respect to gene parts
#for each methylation region, the function will tell you which part of a gene it overlaps with (ex promoter, exon, intron...)
#should output dataframe or GRanges object where each row corresponds to a methylation region
#columns = chr, overlapping gene name or ID, which gene part it overlaps (exon, intron...)
#so pretty much saying this methylated region falls in the promoter of GeneX or this region is in an intron of GeneY

##############Tiles1000################
#annotate tiles based on 25, 0.05
#annotate_Tiles_q5d25 <- annotateWithGeneParts(as(myDiff25q5_Tiles1000, "GRanges"), gene.obj)
#save as rds to be able to put into R studio for further analysis
#saveRDS(annotate_Tiles_q5d25, file = "annotate_Tiles1000_q5d25_CvS.rds")

#annotate tiles based on 10, 0.01 - less stringent
#annotate_Tiles_q1d10 <- annotateWithGeneParts(as(myDiff10q1_Tiles1000, "GRanges"), gene.obj)
#save as rds to be able to put into R studio for further analysis
#saveRDS(annotate_Tiles_q1d10, file = "annotate_Tiles1000_q1d10_CvS.rds")

#############Tiles500###################
#annotate tiles based on 25, 0.05
#annotate_Tiles500_q5d25 <- annotateWithGeneParts(as(myDiff25q5_Tiles500, "GRanges"), gene.obj)
#save as rds to be able to put into R studio for further analysis
#saveRDS(annotate_Tiles500_q5d25, file = "annotate_Tiles500_q5d25_CvS.rds")

#annotate tiles based on 10, 0.01 - less stringent
#annotate_Tiles500_q1d10 <- annotateWithGeneParts(as(myDiff10q1_Tiles500, "GRanges"), gene.obj)
#save as rds to be able to put into R studio for further analysis
#saveRDS(annotate_Tiles500_q1d10, file = "annotate_Tiles500_q1d10_CvS.rds")

############Tiles200###################
#annotate_Tiles200_q5d25 <- annotateWithGeneParts(as(myDiff25q5_Tiles200, "GRanges"), gene.obj)
#save as rds to be able to put into R studio for further analysis
#saveRDS(annotate_Tiles200_q5d25, file = "annotate_Tiles200_q5d25_CvS.rds")

#annotate tiles based on 10, 0.01 - less stringent
#annotate_Tiles200_q1d10 <- annotateWithGeneParts(as(myDiff10q1_Tiles200, "GRanges"), gene.obj)
#save as rds to be able to put into R studio for further analysis
#saveRDS(annotate_Tiles200_q1d10, file = "annotate_Tiles200_q1d10_CvS.rds")

############Tiles25###################
#annotate_Tiles25_q5d25 <- annotateWithGeneParts(as(myDiff25q5_Tiles25, "GRanges"), gene.obj)
#save as rds to be able to put into R studio for further analysis
#saveRDS(annotate_Tiles25_q5d25, file = "annotate_Tiles25_q5d25_CvS.rds")

#annotate tiles based on 10, 0.01 - less stringent
#annotate_Tiles25_q1d10 <- annotateWithGeneParts(as(myDiff10q1_Tiles25, "GRanges"), gene.obj)
#save as rds to be able to put into R studio for further analysis
#saveRDS(annotate_Tiles25_q1d10, file = "annotate_Tiles25_q1d10_CvS.rds")

#read in annotation files
#library(qs)
#setwd("~/Downloads")
#annotations <- qread("all_df_annot_Lrhomboides.qs")
#colnames(annotations)

#setwd("~/Desktop/Methylation/***output files with T2***/CvS/3L")

#tiles1000 with cov bases 3
#annotate_Tiles1000_q1d10 <- readRDS("annotate_Tiles1000_q1d10_CvS.rds")
#annotate_Tiles1000_q5d25 <- readRDS("annotate_Tiles1000_q5d25_CvS.rds")

#tiles 500 with cov bases 3
#annotate_Tiles500_q1d10 <- readRDS("annotate_Tiles500_q1d10_CvS.rds")
#annotate_Tiles500_q5d25 <- readRDS("annotate_Tiles500_q5d25_CvS.rds")

#tiles 200 with cov bases 3
#annotate_Tiles200_q1d10 <- readRDS("annotate_Tiles200_q1d10_CvS.rds")
#annotate_Tiles200_q5d25 <- readRDS("annotate_Tiles200_q5d25_CvS.rds")

#tiles 25 with cov bases 3
#annotate_Tiles25_q1d10 <- readRDS("annotate_Tiles25_q1d10_CvS.rds")
#annotate_Tiles25_q5d25 <- readRDS("annotate_Tiles25_q5d25_CvS.rds")

### annotate tiles ###
#library(dplyr)

############## annotate tiles 1000, q1, d10 ##############
#assoc_Tiles1000_q1d10 <- getAssociationWithTSS(annotate_Tiles1000_q1d10)
#assoc_Tiles1000_q1d10

#read feature name as character not number
#assoc_Tiles1000_q1d10$feature.name <- as.character(assoc_Tiles1000_q1d10$feature.name)

#read rowname in annotations file as character not number
#annotations$rowname <- as.character(annotations$rowname)

#run annotation
#assoc_annot_Tiles1000_q1d10 <- assoc_Tiles1000_q1d10 %>%
#  left_join(annotations, by = c("feature.name" = "rowname"))

#look at it in environment not console or terminal because lots of columns

#save as new file to just read in annotations
#saveRDS(assoc_annot_Tiles1000_q1d10,
#        "CvS_assoc_annot_Tiles1000_q1d10.rds")

############# annotate tiles 1000, q5, d25 ###########
#assoc_Tiles1000_q5d25 <- getAssociationWithTSS(annotate_Tiles1000_q5d25)
#assoc_Tiles1000_q5d25

#read feature name as character not number
#assoc_Tiles1000_q5d25$feature.name <- as.character(assoc_Tiles1000_q5d25$feature.name)

#read rowname in annotations file as character not number
#annotations$rowname <- as.character(annotations$rowname)

#run annotation
#assoc_annot_Tiles1000_q5d25 <- assoc_Tiles1000_q5d25 %>%
#  left_join(annotations, by = c("feature.name" = "rowname"))

#look at it in environment not console or terminal because lots of columns

#save as new file to just read in annotations
#saveRDS(assoc_annot_Tiles1000_q5d25,
 #       "CvS_assoc_annot_Tiles1000_q5d25.rds")

############# annotate tiles 500, q1, d10 ############
#assoc_Tiles500_q1d10 <- getAssociationWithTSS(annotate_Tiles500_q1d10)
#assoc_Tiles500_q1d10

#read feature name as character not number
#assoc_Tiles500_q1d10$feature.name <- as.character(assoc_Tiles500_q1d10$feature.name)

#read rowname in annotations file as character not number
#annotations$rowname <- as.character(annotations$rowname)

#run annotation
#assoc_annot_Tiles500_q1d10 <- assoc_Tiles500_q1d10 %>%
#  left_join(annotations, by = c("feature.name" = "rowname"))

#look at it in environment not console or terminal because lots of columns

#save as new file to just read in annotations
#saveRDS(assoc_annot_Tiles500_q1d10,
 #       "CvS_assoc_annot_Tiles500_q1d10.rds")

########### annotate tiles 500, q5, d25 ##############
#assoc_Tiles500_q5d25 <- getAssociationWithTSS(annotate_Tiles500_q5d25)
#assoc_Tiles500_q5d25

#read feature name as character not number
#assoc_Tiles500_q5d25$feature.name <- as.character(assoc_Tiles500_q5d25$feature.name)

#read rowname in annotations file as character not number
#annotations$rowname <- as.character(annotations$rowname)

#run annotation
#assoc_annot_Tiles500_q5d25 <- assoc_Tiles500_q5d25 %>%
#  left_join(annotations, by = c("feature.name" = "rowname"))

#look at it in environment not console or terminal because lots of columns

#save as new file to just read in annotations
#saveRDS(assoc_annot_Tiles500_q5d25,
 #       "CvS_assoc_annot_Tiles500_q5d25.rds")

########### annotate tiles 200, q5, d25 ##############
#assoc_Tiles200_q5d25 <- getAssociationWithTSS(annotate_Tiles200_q5d25)
#assoc_Tiles200_q5d25

#read feature name as character not number
#assoc_Tiles200_q5d25$feature.name <- as.character(assoc_Tiles200_q5d25$feature.name)

#read rowname in annotations file as character not number
#annotations$rowname <- as.character(annotations$rowname)

#run annotation
#assoc_annot_Tiles200_q5d25 <- assoc_Tiles200_q5d25 %>%
#  left_join(annotations, by = c("feature.name" = "rowname"))

#look at it in environment not console or terminal because lots of columns

#save as new file to just read in annotations
#saveRDS(assoc_annot_Tiles200_q5d25,
 #       "CvS_assoc_annot_Tiles200_q5d25.rds")

############# annotate tiles 200, q1, d10 ############
#assoc_Tiles200_q1d10 <- getAssociationWithTSS(annotate_Tiles200_q1d10)
#assoc_Tiles200_q1d10

#read feature name as character not number
#assoc_Tiles200_q1d10$feature.name <- as.character(assoc_Tiles200_q1d10$feature.name)

#read rowname in annotations file as character not number
#annotations$rowname <- as.character(annotations$rowname)

#run annotation
#assoc_annot_Tiles200_q1d10 <- assoc_Tiles200_q1d10 %>%
#  left_join(annotations, by = c("feature.name" = "rowname"))

#look at it in environment not console or terminal because lots of columns

#save as new file to just read in annotations
#saveRDS(assoc_annot_Tiles200_q1d10,
 #       "CvS_assoc_annot_Tiles200_q1d10.rds")

############# annotate tiles 25, q1, d10 ############
#assoc_Tiles25_q1d10 <- getAssociationWithTSS(annotate_Tiles25_q1d10)
#assoc_Tiles25_q1d10

#read feature name as character not number
#assoc_Tiles25_q1d10$feature.name <- as.character(assoc_Tiles25_q1d10$feature.name)

#read rowname in annotations file as character not number
#annotations$rowname <- as.character(annotations$rowname)

#run annotation
#assoc_annot_Tiles25_q1d10 <- assoc_Tiles25_q1d10 %>%
 # left_join(annotations, by = c("feature.name" = "rowname"))

#look at it in environment not console or terminal because lots of columns

#save as new file to just read in annotations
#saveRDS(assoc_annot_Tiles25_q1d10,
 #       "CvS_assoc_annot_Tiles25_q1d10.rds")

########### annotate tiles 25, q5, d25 ##############
#assoc_Tiles25_q5d25 <- getAssociationWithTSS(annotate_Tiles25_q5d25)
#assoc_Tiles25_q5d25

#read feature name as character not number
#assoc_Tiles25_q5d25$feature.name <- as.character(assoc_Tiles25_q5d25$feature.name)

#read rowname in annotations file as character not number
#annotations$rowname <- as.character(annotations$rowname)

#run annotation
#assoc_annot_Tiles25_q5d25 <- assoc_Tiles25_q5d25 %>%
#  left_join(annotations, by = c("feature.name" = "rowname"))

#look at it in environment not console or terminal because lots of columns

#save as new file to just read in annotations
#saveRDS(assoc_annot_Tiles25_q5d25,
#        "CvS_assoc_annot_Tiles25_q5d25.rds")




#make manhattan plot of methylated CpG sites that are sig differentially methylated
#based off https://royalsocietypublishing.org/rspb/article/284/1864/20171667/78684/Persistent-and-plastic-effects-of-temperature-on

#load libraries
library(dplyr)
library(ggplot2)

#read in DMRs
setwd("~/Desktop/Methylation/***output files with T2***/CvS/3L")
dmr <- read.csv("DMRs_CvS.csv")
head(dmr)
colnames(dmr)
str(dmr)

#since these are methylation sites not regions dont have to find midpoint to make one genomic coordinate for the manhattan plot
#can just use the start coordinate as the coordinate
#SOOOO call them DMC (differentially methylated cytosines)

#Order the contigs
dmc <- dmr
dmc$chr <- as.factor(dmc$chr)
dmc <- dmc[order(dmc$chr, dmc$start), ]

# Step 2: Calculate max position for each contig
chr_lengths <- dmc %>%
  group_by(chr) %>%
  summarise(max_pos = max(start)) %>%
  arrange(chr)

# Step 3: Calculate cumulative offsets
chr_lengths <- chr_lengths %>%
  mutate(cumstart = lag(cumsum(max_pos), default = 0))

# Step 4: Add cumulative position to DMC table
dmc <- dmc %>%
  left_join(chr_lengths[, c("chr", "cumstart")], by = "chr") %>%
  mutate(pos = start + cumstart)

# Step 5: Classify direction
dmc$direction <- ifelse(dmc$meth.diff > 0, "Hyper", "Hypo")

# Now `pos` can go on x-axis and `meth.diff` on y-axis for the Manhattan plot
head(dmc)


#plot
ggplot(dmc, aes(x = pos, y = meth.diff, color = direction)) +
  geom_point(alpha = 0.6, size = 1) +
  theme_bw() +
  labs(
    x = "Genomic position (cumulative)",
    y = "Methylation difference (%)",
    color = "Direction"
  ) +
  scale_color_manual(values = c("Hyper" = "red", "Hypo" = "blue"))

#add vertical lines on plot
# Make a line for each chromosome to separate visually - line will be at beginning of each
#24 scaffolds in ref genome so 24 vertical lines
chr_starts <- dmc %>%
  group_by(chr) %>%
  summarise(start_pos = min(pos)) %>%  # start of each chromosome
  arrange(start_pos) %>%
  mutate(label = paste0("chr", row_number()))  # chr1, chr2, ...

# Plot
ggplot(dmc, aes(x = pos, y = meth.diff, color = direction)) + #genomic position on x axis, meth perc on y axis
  geom_point(alpha = 0.6, size = 1) +
  geom_vline(data = chr_starts, aes(xintercept = start_pos), linetype = "dashed", color = "grey70") + #add vertical lines
  theme_classic() +
  labs(
    x = "Genomic position (chromosome)",
    y = "Methylation difference (%)",
    color = "Direction"
  ) +
  scale_color_manual(values = c("Hyper" = "red", "Hypo" = "blue")) +
  scale_x_continuous(breaks = NULL)



######### bar plot to show DMC% per chromosome ###########

#### decided not to do this since methylation % are so low looks quite wonky ######

#library(dplyr)
#library(ggplot2)
#library(tidyr)

#read in all CpGs file
#setwd("~/Desktop")
#all_CpGs <- read.csv("Control_all_CpGs.csv")
#head(all_CpGs)

#prepare DMCs - make sure there is hyper adn hypo information and genomic location
#dmc <- dmc %>%
  # make sure you have a direction column
 # mutate(direction = ifelse(meth.diff > 0, "Hyper", "Hypo"))

# Count DMCs per chromosome and direction
#chrom_freq <- dmc %>%
 # group_by(chr, direction) %>%
 # summarise(n_DMCs = n(), .groups = "drop") %>%
  # make sure every chr has both Hyper and Hypo rows
 # complete(chr, direction = c("Hyper", "Hypo"), fill = list(n_DMCs = 0))

# Count total CpGs per chromosome
#total_CpGs <- all_CpGs %>%
#  group_by(chr) %>%
 # summarise(n_CpGs = n(), .groups = "drop")

# Join DMCs with total CpGs and calculate percent
#chrom_freq <- chrom_freq %>%
#  left_join(total_CpGs, by = "chr") %>%
#  mutate(percent = n_DMCs / n_CpGs * 100)

# Keep only the chromosomes with non-zero percent DMCs
#chrom_freq_filtered <- chrom_freq %>%
#  filter(percent > 0)

# Ensure Hyper is bottom, Hypo is top
#chrom_freq_filtered$direction <- factor(chrom_freq_filtered$direction, levels = c("Hyper", "Hypo"))

# Plot
#ggplot(chrom_freq_filtered, aes(x = chr, y = n_DMCs, fill = direction)) +
#  geom_bar(stat = "identity") +
 # scale_fill_manual(values = c("Hyper" = "darkgrey", "Hypo" = "lightgrey")) +
 # theme_classic() +
 # labs(x = "Chromosome", y = "Number of DMCs", fill = "Direction") +
 # theme(axis.text.x = element_text(angle = 45, hjust = 1))



#########  violin plot  ##########

#need percent methylation file

setwd("~/Desktop/Methylation/***output files with T2***/CvS/3L")
perc_methyl <- read.csv("meth_perc_CvS_with_coords.csv")
head(perc_methyl)

#change formatting of data table
treatment_info <- data.frame(
  sample_id = c("C_T1_03","C_T1_14","C_T1_16","C_T5_10","C_T5_12",
                "S_T1_01","S_T1_07","S_T1_12","S_T3_07","S_T3_11"),
  treatment = c("Control","Control","Control","Control","Control",
                "Single","Single","Single","Single","Single"))

# Convert data table from wide → long
perc_methyl_long <- perc_methyl %>%
  pivot_longer(cols = starts_with("C_") | starts_with("S_"),
               names_to = "sample_id",
               values_to = "perc_meth") %>%
  left_join(treatment_info, by = "sample_id")

#make sure it worked
head(perc_methyl_long)

# Filter out any missing values
perc_methyl_long_clean <- perc_methyl_long %>%
  filter(!is.na(perc_meth))

#make sure it worked
head(perc_methyl_long_clean)

# Violin plot
ggplot(perc_methyl_long_clean, aes(x = treatment, y = perc_meth, fill = treatment)) +
  geom_violin(trim = FALSE, scale = "width") +       # violin showing density
  geom_boxplot(width = 0.1, outlier.size = 0.5) +   # embedded boxplot
  scale_fill_manual(values = c(
    "Control" = "grey",
    "Single"= "red"
  )) +
  theme_classic() +
  labs(
    x = "Treatment",
    y = "% methylation",
    fill = "Treatment"
  )



######### run stats #######

#want to get summary of methylation per individual (so mean, median, and ranges)

#then want to see if any treatments are different on average

#then want to look at distribution
#distributional test (K-S)
#are methylation distributions different even if means arent
#code based on: https://onlinelibrary.wiley.com/doi/full/10.1111/mec.15764

#first think of what data I have situated
#I have methylation summaries per CpG per sample
#I want mean and median % methylation per sample

#load libraries
library(dplyr)
library(ggplot2)
library(tidyr)

#read in perc methylation data file
setwd("~/Desktop/Methylation/***output files with T2***/CvS/3L")
perc_methyl <- read.csv("meth_perc_CvS_with_coords.csv")
head(perc_methyl)

#do the same thing to the table you did before when making the violin plot

#change formatting of data table
treatment_info <- data.frame(
  sample_id = c("C_T1_03","C_T1_14","C_T1_16","C_T5_10","C_T5_12",
                "S_T1_01","S_T1_07","S_T1_12","S_T3_07","S_T3_11"),
  treatment = c("Control","Control","Control","Control","Control",
                "Single","Single","Single","Single","Single"))

# Convert data table from wide → long
perc_methyl_long <- perc_methyl %>%
  pivot_longer(cols = starts_with("C_") | starts_with("S_"),
               names_to = "sample_id",
               values_to = "perc_meth") %>%
  left_join(treatment_info, by = "sample_id")

#make sure it worked
head(perc_methyl_long)
#this table is saying: at this genomic position, this individual has X% methylation

#make new dataframe to calc average % methylation across all CpGs in each individual (mean_meth) AND midpoint of methyl dist for that individual (median_meth)
sample_summary <- perc_methyl_long %>%
  filter(!is.na(perc_meth)) %>%       #filter out NAs
  group_by(sample_id, treatment) %>%      #group by sample ID and treatment so one row per sample
  summarise(
    mean_meth = mean(perc_meth),
    median_meth = median(perc_meth),
    .groups = "drop"
  )

#so this will allow us to determine what fraction or % of CpGs are methylated on average
#making genome wide methylation summaries per individual

sample_summary

#looks good!

#sample_id treatment mean_meth median_meth
# C_T1_03   Control        72.6        83.3
# C_T1_14   Control        74.0        84.6
# C_T1_16   Control        72.5        83.3
# C_T5_10   Control        73.0        83.3
# C_T5_12   Control        72.6        83.3
# S_T1_01   Single         72.8        83.3
# S_T1_07   Single         71.5        81.8
# S_T1_12   Single         72.5        83.3
# S_T3_07   Single         71.7        83.3
# S_T3_11   Single         72.3        83.3

#now convert treatment to a factor
sample_summary$treatment <- factor(
  sample_summary$treatment,
  levels = c("Control", "Single"))


#### run stats now that you have the mean, median ###

wilcox.test(mean_meth ~ treatment, data = sample_summary)
#W = 22, p-value = 0.05556
#n = 5 per group
#no sig difference in average genome-wide methylation per individual SO global methylation doesnt shift much with treatment
  #but close!!

wilcox.test(
  perc_meth ~ treatment,
  data = perc_methyl_long,
  exact = FALSE
)
#n = millions of CpGs
#W = 6.3411e+13, p-value < 2.2e-16
#stat sig difference in overall CpG distribution between treatments
#caveat = sample size is huge and effect size is small so might not be biologically significant


#### cumulative distribution plot of per-CpG methylation ####
#Empirical Cumulative Distribution Function (ECDF) for each treatment group
#plotting the fraction of CpGs with methylation that is less that or equ2al to that value
#overall, showing how CpG methylation is distributed across the genome
#potentially able to see shifts in distrbution

#remove NAs
perc_methyl_long_clean <- perc_methyl_long %>%
  filter(!is.na(perc_meth))

# Make treatment a factor
perc_methyl_long_clean$treatment <- factor(
  perc_methyl_long_clean$treatment, 
  levels = c("Control", "Single")
)

# Cumulative distribution plot
ggplot(perc_methyl_long_clean, aes(x = perc_meth, color = treatment)) +
  stat_ecdf(linewidth = 1) +                                # empirical cumulative distribution function
  scale_color_manual(values = c("Control" = "grey", "Single" = "red")) +
  theme_classic() +
  labs(
    x = "CpG methylation (%)",
    y = "Cumulative fraction of CpGs",
    title = "Cumulative distribution of per-CpG methylation"
  )

#plot again but with vertical line at the methylation median

#plot
ggplot(perc_methyl_long_clean, aes(x = perc_meth, color = treatment)) +
  stat_ecdf(linewidth = 1) +
  geom_vline(data = sample_summary, aes(xintercept = median_meth, color = treatment),
             linetype = "dashed") +
  scale_color_manual(values = c("Control" = "grey", "Single" = "red")) +
  theme_classic() +
  labs(
    x = "CpG methylation (%)",
    y = "Cumulative fraction of CpGs",
    title = "Cumulative distribution of per-CpG methylation with median"
  )

#median is the same for both so cant see control because single line over it





