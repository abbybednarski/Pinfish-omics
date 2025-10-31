#set to wherever coldata and allcounts are to run properly
setwd("~/Desktop/Tagseq")

#load in packages, maybe have to install with BiocManager in the future
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

#check for outliers on PCA, single HW tower 3

# Remove outlier samples (detected on heatmap of sample-to-sample distance matrices), specific for Katie
#countdata <- countdata[, colnames(countdata) != "C.T5.10.counts"]
# Remove samples taken at incorrect temp (from tower 2)
countdata <- countdata[, colnames(countdata) != "M.T2.09.M.counts"]
countdata <- countdata[, colnames(countdata) != "M.T2.16.M.counts"]

# Read in sample metadata
coldata <- read.table("coldata.txt")
#make sure sampling temp is a factor, not going to be changed bw treatments/over the dataset
coldata$SamplingTemp <- as.factor(coldata$SamplingTemp)
#make sure the towers are read as factors instead of numerical data
coldata$Tower <- as.factor(coldata$Tower)

# Remove outlier samples (if any) and samples from tower 2
coldata <- coldata[rownames(coldata) != "M.T2.09.M.counts",]
coldata <- coldata[rownames(coldata) != "M.T2.16.M.counts",]

#need for DESeq, should say true, seeing if the column names match countdata and colnames, check link for more explanation
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
  #this means control is the ref for SvC and MvC but multiple is the ref for SvM

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
  scale_color_manual(values = c("red", "blue", "green")) + 
  theme(legend.background = element_rect(fill = "white", color = "black"))







#stopped hereeeeee - nothing else past this but keep in case
#after this, have a lot to change

# Run WGCNA on your dataset
coldata$fishID <- rownames(coldata)
normalized_counts <- assay(vsd) %>%
  t()

sft <- pickSoftThreshold(normalized_counts, dataIsExpr = T, networkType = "signed")
sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) + 
  geom_point() + geom_text(nudge_y = 0.1) + 
  geom_hline(yintercept = 0.80, col = "red") + 
  ylim(c(min(sft_df$model_fit), 1.05)) + 
  xlab("Soft Threshold (power)") + 
  ylab("Scale Free Topology Model Fit, signed R^2") + 
  ggtitle("Scale Independence") + theme_classic()

bwnet <- blockwiseModules(normalized_counts, maxBlockSize = 5000, 
                          TOMType = "signed", power = 12, numericLabels = T, 
                          randomSeed = 1234)

readr::write_rds(bwnet, file = "WGCNA_results.RDS")

module_eigengenes <- bwnet$MEs
head(module_eigengenes)

all.equal(coldata$fishID, rownames(module_eigengenes))

des_mat <- model.matrix(~ coldata$Treatment)
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)
fit <- limma::eBayes(fit)
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")
head(stats_df)
# Three modules significantly different across groups
# ME167, ME140, ME285

module_df <- module_eigengenes %>%
  tibble::rownames_to_column("FishID") %>%
  dplyr::inner_join(coldata %>%
                      dplyr::select(fishID, Treatment), by = c("FishID" = "fishID"))

ggplot(module_df, aes(x = Treatment, y = ME285, color = Treatment)) + 
  geom_boxplot(width = 0.2, outlier.shape = NA) + ggforce::geom_sina(maxwidth = 0.3) + 
  theme_classic()

gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
  dplyr::mutate(module = paste0("ME", module))

gene_module_key_ME167 <- gene_module_key %>%
  dplyr::filter(module == "ME167")

make_module_heatmap <- function(module_name, expression_mat = normalized_counts, # These next few lines specify the function name (make_module_heatmap) and any arguments the function takes
                                metadata_df = coldata, gene_module_key_df = gene_module_key, # as well as any default values of those arguments. For example, it requires a module_name, as well as an expression matrix
                                module_eigengenes_df = module_eigengenes) { # There's no default module name, but the default expression matrix name is normalized_counts, etc.
  
  module_eigengene <- module_eigengenes_df %>% # This block of code creates a little table with the fish ID and the module eigengene value associated with each fish
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("fishID")
  
  col_annot_df <- metadata_df %>% # This block of code adds the treatment information (i.e., Dec or Feb) to that table
    dplyr::select(Treatment, fishID) %>%
    dplyr::inner_join(module_eigengene, by = "fishID") %>%
    dplyr::arrange(Treatment, fishID) %>%
    tibble::column_to_rownames("fishID")
  
  col_annot <- ComplexHeatmap::HeatmapAnnotation( # This block of code creates a small barplot below the actual heatmap, which will contain the eigengene (i.e., summary) expression values for each sample. The samples will be blocked and colored by treatment.
    treatment = col_annot_df$Treatment,
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    col = list(Treatment = c("Control" = "blue2", "Single" = "coral2", "Multiple" = "green"))
  )
  
  module_genes <- gene_module_key_df %>% # This creates an object called module_genes, which selects all of the gene names in our module of interest
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)
  
  mod_mat <- expression_mat %>% # This grabs the expression values for each gene in our module
    t() %>%
    as.data.frame() %>%
    dplyr::filter(rownames(.) %in% module_genes) %>%
    dplyr::select(rownames(col_annot_df)) %>%
    as.matrix()
  
  mod_mat <- mod_mat %>% # This normalizes the expression values
    t() %>%
    scale() %>%
    t()
  
  color_func <- circlize::colorRamp2( # And this creates our colorization scale
    c(-2, 0, 2),
    c("blue2", "white", "coral2")
  )
  
  heatmap <- ComplexHeatmap::Heatmap(mod_mat, # Finally, this plots the expression values of each gene in each sample on a heatmap
                                     name = module_name,
                                     col = color_func,
                                     bottom_annotation = col_annot,
                                     cluster_columns = F,
                                     show_row_names = F,
                                     show_column_names = F)
  
  return(heatmap) # And the function will output the heatmap
  
}

module_167_heatmap <- make_module_heatmap(module_name="ME167")
