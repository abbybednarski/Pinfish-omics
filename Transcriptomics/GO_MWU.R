
# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) to identify GO categories that are significantly enriched with either up- or down-regulated genes. The advantage - no need to impose arbitrary significance cutoff.
# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based Fisher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure. 
# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.
# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.
# The fraction next to GO category name indicates the fracton of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValue cutoff (option in gomwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.
# Stretch the plot manually to match tree to text
# Based from Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu

################################################################

#since MWU, looking at log2fold changes of ALL genes found
  #this is independent of the measure of significance used for DEGs
setwd("~/Desktop/Tagseq")

# Edit these to match your data file names: 
input="Multiple-v-Control-ALL-GO.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="GO_List_for_MWU_really_correct.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")

#change input file depending on pairwise comparison
#switch up BP, MF, CC for all of the above too (so BP, MF, and CC for each pwc so 9 output files)

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
gomwuStats(input, goDatabase, goAnnotations, goDivision,
	perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
	largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
	smallest=5,   # a GO category should contain at least this many genes to be considered
	clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
#	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
#	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
#	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.

#27 GO terms for MF for Single v Control
#23 GO term for CC for single v control
#43 GO terms for BP for single v control

#6 GO terms for CC for multiple v control
#22 GO terms for BP for multiple v control
#7 GO terms for MF for multiple v control

#17 GO terms for MF for single v multiple
#0 GO terms for BP for single v multiple
#15 GO terms for CC for single v multiple
