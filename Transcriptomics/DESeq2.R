#set to wherever coldata.txt (sample info) and allcounts.txt (gene count data) are to run properly
setwd("~/Desktop/Tagseq")

#load in packages
#install.packages("tidyverse")
library(tidyverse)
library(magrittr)

#install.packages("WGCNA")
#library(WGCNA)
#library(BiocManager)
#BiocManager::install('GO.db')
library(GO.db)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

# Read in gene expression count data
countdata <- as.matrix(read.table("allcounts.txt", sep = "\t", row.names = 1, header = T))
#remove the X.1 column that has no information in it, may not need this line for other datasets
countdata <- countdata[, colnames(countdata) != "X.1"]

# Remove samples taken at incorrect temp (from tower 2)
countdata <- countdata[, colnames(countdata) != "M.T2.09.M.counts"]
countdata <- countdata[, colnames(countdata) != "M.T2.16.M.counts"]

# Read in sample metadata
coldata <- read.table("coldata.txt")
#make sure sampling temp and tower are read in as a factor
coldata$SamplingTemp <- as.factor(coldata$SamplingTemp)
coldata$Tower <- as.factor(coldata$Tower)

# Remove samples from tower 2
coldata <- coldata[rownames(coldata) != "M.T2.09.M.counts",]
coldata <- coldata[rownames(coldata) != "M.T2.16.M.counts",]

#need for DESeq, should say true, seeing if the column names match countdata and colnames
all(rownames(coldata) == colnames(countdata))

#make big matrix to be able to check for outliers, design=~ Treatment to compare TREATMENT rather than anything else, name the titles of parts of matrix, factor in weight, front/back of tower
#had to center weight to make mean smaller, just did in excel in coldata.txt
dds_noT2 <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design =~ FrontBack + CenterWeight + Treatment)

#make sure the reference is control and not another treatment, will go alphabetically so will most likely be control anyways but run just in case
dds_noT2$Treatment <- relevel(dds_noT2$Treatment, ref = "Control")

#run deseq on matrix to actually have a matrix
dds_noT2 <- DESeq(dds_noT2)

#can only look at pairwise comparison so comparing only 2 different treatments to each other, not 3
#for the following code, this structure is followed:
  #contrast = c("Treatment", "InterestGroup", "ReferenceGroup")
  #this means control is the ref for SvC and MvC but multiple is the ref for SvM *******

#compare multiple to control, alpha should match p value
res_control_multiple <- results(dds_noT2, alpha = 0.05, contrast = c("Treatment", "Multiple", "Control"))
#subset to only show significantly changed genes, will be written to csv to be annotated
res_control_multiple_sigs <- subset(res_control_multiple, padj < 0.05)
#show amount of sig changed genes
summary(res_control_multiple_sigs)
#87 up
#137 down
#total: 224

#same thing as above but now comparing single vs control
res_control_single <- results(dds_noT2, alpha = 0.05, contrast = c("Treatment", "Single", "Control"))
res_control_single_sigs <- subset(res_control_single, padj < 0.05)
summary(res_control_single_sigs)
#55 up
#112 down
#total: 167

#same as above but now comparing single vs multiple
  #multiple is reference
res_single_multiple <- results(dds_noT2, alpha = 0.05, contrast = c("Treatment", "Single", "Multiple"))
res_single_multiple_sigs <- subset(res_single_multiple, padj < 0.05)
summary(res_single_multiple_sigs)
#least difference here
#10 up
#11 down
#total: 21

#compare front/back of towers to see if effect
res_front_back <- results(dds_noT2, alpha = 0.05, contrast = c("FrontBack", "Front", "Back"))
res_front_back_sigs <- subset(res_front_back, padj < 0.05)
summary(res_front_back_sigs)
#4 down, nothing else

# Generate sample-by-sample distance matrix to check for outliers, will be really light if outlier
vsd <- varianceStabilizingTransformation(dds_noT2)
sampleDists <- dist(t(assay(vsd)))
#make as a matrix
sampleDistMatrix <- as.matrix(sampleDists)
#name rows to TREATMENT (what youre looking at)
rownames(sampleDistMatrix) <- paste(vsd$Treatment)
#doesnt matter but run
colnames(sampleDistMatrix) <- NULL
#blue distribution of color, white will show sample is far away from others
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#make PCA to look at outliers another way, change TREATMENT AND TOWER to see reason for weird single samples
DESeq2::plotPCA(vsd, intgroup = "Tower", ntop = 100000)
DESeq2::plotPCA(vsd, intgroup = "Treatment", ntop = 100000)

#export sig genes to csv files
write.csv(res_control_multiple_sigs, "Updated-Multiple-vs-Control.csv")
write.csv(res_control_single_sigs, "Updated-Single-vs-Control.csv")
write.csv(res_single_multiple_sigs, "Updated-Single-vs-Multiple.csv")

#export all genes to csv files
write.csv(res_control_multiple, "Multiple-v-Control-ALL.csv")
write.csv(res_control_single, "Single-v-Control-ALL.csv")
write.csv(res_single_multiple, "Single-v-Multiple-ALL.csv")
write.csv(res_front_back, "Front-v-Back-ALL.csv")
#these used eventually for GO analysis with some changes
#make separate file that just has Trinity ID and log2foldchange, use this for GO analysis

# Write significant output files containing only significant genes
write.csv(res_control_multiple_sigs, "Multiple-v-Control-sigs-full.csv")
write.csv(res_control_single_sigs, "Single-v-Control-sigs-full.csv")
write.csv(res_single_multiple_sigs, "Single-v-Multiple-sigs-full.csv")
write.csv(res_front_back_sigs, "Front-v-Back-sigs-full.csv")

# Read in significant gene files
MvC_sig <- read.csv("Multiple-v-Control-sigs-full.csv")
SvC_sig <- read.csv("Single-v-Control-sigs-full.csv")
SvM_sig <- read.csv("Single-v-Multiple-sigs-full.csv")
FvB_sig <- read.csv("Front-v-Back-sigs-full.csv")

# Read in annotation file
#install.packages("openxlsx")
library(openxlsx)
txm_anots <- read.xlsx("FINAL_ANNOTATION_CORRECT_by_gene.xlsx") # This should be the pinfish transcriptome annotation file I emailed you
# Create new column in annotation file that removes the _i isoform designation on the Trinity ID
txm_anots$TrinIDGene <- gsub("_i\\d+", "", txm_anots$TrinityID)

# Merge transcriptome annotations with significant gene tables
MvC_sig_anot <- merge.data.frame(MvC_sig, txm_anots, by.x = "X", by.y = "TrinIDGene", all.x = T)
SvC_sig_anot <- merge.data.frame(SvC_sig, txm_anots, by.x = "X", by.y = "TrinIDGene", all.x = T)
SvM_sig_anot <- merge.data.frame(SvM_sig, txm_anots, by.x = "X", by.y = "TrinIDGene", all.x = T)
FvB_sig_anot <- merge.data.frame(FvB_sig, txm_anots, by.x = "X", by.y = "TrinIDGene", all.x = T)

# Output significant gene tables with annotations
write.csv(MvC_sig_anot, "Multiple-v-Control-significant-annotated.csv")
write.csv(SvC_sig_anot, "Single-v-Control-significant-annotated.csv")
write.csv(SvM_sig_anot, "Single-v-Multiple-significant-annotated.csv")
write.csv(FvB_sig_anot, "Front-v-Back-significant-annotated.csv")


### generate pretty PCA ###
# Extract the transformed data from vsd
pca_data <- prcomp(t(assay(vsd)))

# Create a data frame of PCA results
pca_df <- as.data.frame(pca_data$x)
pca_df$Treatment <- vsd$Treatment  # Add treatment to data frame to be able to color and shape by treatment

# Plot the PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment, shape = Treatment)) +
  geom_point(size = 5) +
  #stat_ellipse(linewidth = 0.5, level = 0.95) + #add ellipses
  labs(title = NULL,
    x = paste0("PC1 (", round(summary(pca_data)$importance[2, 1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca_data)$importance[2, 2] * 100, 1), "% variance)")) +
  theme_classic() +
  scale_color_manual(values = c("grey", "blue", "red")) + 
  theme(legend.background = element_rect(fill = "white", color = "black"))


#PCA with ellipses
library(ggforce)

x_buffer <- (max(pca_df$PC1) - min(pca_df$PC1)) * 0.5
y_buffer <- (max(pca_df$PC2) - min(pca_df$PC2)) * 0.5

#plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment, shape = Treatment)) +
  geom_point(size = 5) +
  ggforce::geom_mark_ellipse(aes(color = Treatment)) +
  labs(title = NULL,
       x = paste0("PC1 (", round(summary(pca_data)$importance[2, 1] * 100, 1), "% variance)"),
       y = paste0("PC2 (", round(summary(pca_data)$importance[2, 2] * 100, 1), "% variance)")) +
  theme_classic() +
  scale_color_manual(values = c("grey", "blue", "red")) + 
  theme(legend.background = element_rect(fill = "white", color = "black")) +
  coord_cartesian(
    xlim = c(min(pca_df$PC1) - x_buffer, max(pca_df$PC1) + x_buffer),
    ylim = c(min(pca_df$PC2) - y_buffer, max(pca_df$PC2) + y_buffer),
    expand = TRUE,
    clip = "off")
