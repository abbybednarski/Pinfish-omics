#2.3 Statistical analysis (one-factor)
#code from https://www.metaboanalyst.ca/resources/vignettes/Statistical_Analysis_Module.html

setwd("~/Desktop/Metabolomics/neg_stats_files")

# Clean global environment if needed
rm(list = ls())

# Load MetaboAnalystR
library(MetaboAnalystR)


####### Univariate methods ########

##Initialize the data, conc = input data represents concentration values, stat = intended analysis is statistical analysis
mSet<-InitDataObjects("conc", "stat", FALSE);

#Load in dataset, rowu = rows represent metabolites, disc = discontinuous data for grouping
mSet <- Read.TextData(mSet, "~/Desktop/Metabolomics/neg_stats_files/neg_all_samples.csv", "rowu", "disc");

#check for missing values
sum(is.na(mSet$dataSet$norm))
#no missing values

#check for missing values, invalid or formatting issues
mSet<-SanityCheckData(mSet);
#should say 3 groups here since big dataset includes control, single, multiple

#replace minimal values (missing or extremely low values)
mSet<-ReplaceMin(mSet);

#prepare data for normalization, overwriting old prenorm object / making a new one
mSet<-PreparePrenormData(mSet);

#normalize data, NULL = no scaling, LogNorm = log transformation, meancenter = shifts data so each feature has mean of 0,
#ref = NULL because dont want to focus on one particular sample over another, could input sample information if normalizing to a specific reference
#ratioNum=20 = use top 20 most stable features for normalization
mSet<-Normalization(mSet, "NULL", "LogNorm", "MeanCenter", ref=NULL, ratio=FALSE, ratioNum=20);

#check dimensions of normalized data
print(dim(mSet$dataSet$norm))
#15 456

#plot normalized metabolites
mSet<-PlotNormSummary(mSet, "norm_0_", format = "png", dpi=100, width=NA);
#plot normalized samples overall instead of by metabolites
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "png", dpi=100, width=NA);

############### Fold change analysis ###############

# Perform fold-change analysis - doing pairwise comparisons eventually
# tried to subset the data in R MULTIPLEEEE different ways and nothing worked correctly 
# So I went into the csv file, removed samples I didnt want (ex doing Control vs Multiple so removed single samples), initialized a new mSet, then loaded in the edited csv file to run each pairwise comparison
# because you are uploading a new dataset for each comparison, make sure to normalize, check for missing info, etc. before running FC analysis

# Perform fold-change analysis on uploaded data for the WHOLE dataset
#indicated fc threshold to be 2.0, cmp.type = 0 indicates group 1 - group 2 (not ideal since 3 groups here but not using fc anal for any analysis on the big dataset, eg will use for control vs. multiple where this makes sense)
mSet<-FC.Anal(mSet, 2.0, 0, FALSE)

# Plot fold-change analysis
mSet<-PlotFC(mSet, "Total_FC_", "png", 100, width=NA)

# To view fold-change 
head(mSet$analSet$fc$fc.log)

######### ANOVA for big mSet (which contains control, single, multiple) ########

#no using non-parametric test so F because not assuming any distribution (parametric used when doing ANOVA so use the option nonpar=F), p-value threshold = 0.05, using Fisher's LSD test for post-hoc analysis
#options are fisher or tukey but those are for NOT using a non-parametric test (nonpar = F)
mSet <- ANOVA.Anal(mSet, nonpar=F, thresh=0.05, "fisher")
#18 significant features

### make list of metabolites and output as csv file ###
# want to include fc, pvalue, padj, metabolite name

#extract metabolites
anova_df <- mSet$analSet$aov$sig.mat
anova_df$Metabolite <- rownames(anova_df)
rownames(anova_df) <- NULL

#extract fold change and p-value data
fc_data <- mSet$analSet$fc$fc.log
p_values <- mSet$analSet$aov$p.value

#correction test here to add padj to output file, benjamini hochberg correction
#tried bonferroni correction but nothing was significant due to the high number of metabolites in the dataset so using BH method
#install.packages("stats")
library(stats)

padj_values <- p.adjust(p_values, method = "BH")

#filter metabolites
all_metabolites <- data.frame(
  Metabolite = names(fc_data),
  FoldChange = fc_data,
  PValue = p_values[names(fc_data)],
  AdjustedPValue = padj_values
)

print(all_metabolites)

#save to output file
write.csv(all_metabolites, "all_metabolites_TOTAL.csv", row.names = FALSE)


################ CvS ###############
##Initialize the data, conc = input data represents concentration values, stat = intended analysis is statistical analysis
mSet1<-InitDataObjects("conc", "stat", FALSE);

#Load in dataset, rowu = rows represent metabolites, disc = discontinuous data for grouping
mSet1 <- Read.TextData(mSet1, "~/Desktop/Metabolomics/neg_stats_files/neg_CvS.csv", "rowu", "disc");

#check for missing values
sum(is.na(mSet1$dataSet$norm))
#no missing values

#check for missing values, invalid or formatting issues
mSet1<-SanityCheckData(mSet1);
#should output that it found 2 groups, NOT 3!

#replace minimal values (missing or extremely low values)
mSet1<-ReplaceMin(mSet1);

#prepare data for normalization, making new prenorm object for this specific pairwise comparison
mSet1<-PreparePrenormData(mSet1);

#normalize data, NULL = no scaling, LogNorm = log transformation, meancenter = shifts data so each feature has mean of 0,
#ref = NULL because dont want to focus on one particular sample over another, could input sample information if normalizing to a specific reference
#ratioNum=20 = use top 20 most stable features for normalization
mSet1<-Normalization(mSet1, "NULL", "LogNorm", "MeanCenter", ref=NULL, ratio=FALSE, ratioNum=20);

#check dimensions of normalized data
print(dim(mSet1$dataSet$norm))
#output 10 456

#plot normalized metabolites
mSet1<-PlotNormSummary(mSet1, "CvS_norm_0_", format = "png", dpi=100, width=NA);
#plot normalized samples overall instead of by metabolites
mSet1<-PlotSampleNormSummary(mSet1, "CvS_snorm_0_", format = "png", dpi=100, width=NA);

# Perform fold-change analysis on uploaded data, unpaired, cmp.type = 0 meaning group 1 - group 2 
# group 1 will be the first treatment in the csv file which was the control group
mSet1<-FC.Anal(mSet1, 2.0, 0, FALSE)

# Plot fold-change analysis
mSet1<-PlotFC(mSet1, "CvS_FC_", "png", 100, width=NA)

# To view fold-change 
head(mSet1$analSet$fc$fc.log)

####### Perform ANOVA for mSet1 (CvS) #######
mSet1 <- ANOVA.Anal(mSet1, nonpar=F, thresh = 0.05, "fisher")
#18 significant features

#find a list of all metabolites, output to file with log2foldchange, pvalue, padj
#extract fold change and p-value data
fc_data_CvS <- mSet1$analSet$fc$fc.log
p_values_CvS <- mSet1$analSet$aov$p.value

#benjamini hochberg correction to find adjusted p-values
padj_values_CvS <- p.adjust(p_values_CvS, method = "BH")

#filter metabolites
all_metabolites_CvS <- data.frame(
  Metabolite = names(fc_data_CvS),
  FoldChange = fc_data_CvS,
  PValue = p_values_CvS,
  AdjustPValue = padj_values_CvS
)

print(all_metabolites_CvS)

#save to output file
write.csv(all_metabolites_CvS, "all_metabolites_CvS.csv", row.names = FALSE)

############### CvM ###############
##Initialize the data, conc = input data represents concentration values, stat = intended analysis is statistical analysis
mSet2<-InitDataObjects("conc", "stat", FALSE);

#Load in dataset, rowu = rows represent metabolites, disc = discontinuous data for grouping
mSet2 <- Read.TextData(mSet2, "~/Desktop/Metabolomics/neg_stats_files/neg_CvM.csv", "rowu", "disc");

#check for missing values
sum(is.na(mSet2$dataSet$norm))
#no missing values

#check for missing values, invalid or formatting issues
mSet2<-SanityCheckData(mSet2);
#should output 2 groups here

#replace minimal values (missing or extremely low values)
mSet2<-ReplaceMin(mSet2);

#prepare data for normalization, making new prenorm object for mSet2
mSet2<-PreparePrenormData(mSet2);

#normalize data, NULL = no scaling, LogNorm = log transformation, meancenter = shifts data so each feature has mean of 0,
#ref = NULL because dont want to focus on one particular sample over another, could input sample information if normalizing to a specific reference
#ratioNum=20 = use top 20 most stable features for normalization
mSet2<-Normalization(mSet2, "NULL", "LogNorm", "MeanCenter", ref=NULL, ratio=FALSE, ratioNum=20);

#check dimensions of normalized data
print(dim(mSet2$dataSet$norm))
#10 456

#plot normalized metabolites
mSet2 <-PlotNormSummary(mSet2, "CvM_norm_0_", format = "png", dpi=100, width=NA);
#plot normalized samples overall instead of by metabolites
mSet2 <-PlotSampleNormSummary(mSet2, "CvM_snorm_0_", format = "png", dpi=100, width=NA);

# Perform fold-change analysis on uploaded data, unpaired, control = group 1 so basis for fc analysis
mSet2<-FC.Anal(mSet2, 2.0, 0, FALSE)

# Plot fold-change analysis
mSet2 <-PlotFC(mSet2, "CvM_FC_", "png", 100, width=NA)

# To view fold-change 
head(mSet2$analSet$fc$fc.log)

# Perform ANOVA for mSet2 (CvM)
mSet2 <- ANOVA.Anal(mSet2, nonpar=F, thresh = 0.05, "fisher")
#output 0 significant features

#find all metabolites and output to file, include pvalue, log2foldchange, padj

#extract fold change and p-value data
fc_data_CvM <- mSet2$analSet$fc$fc.log
p_values_CvM <- mSet2$analSet$aov$p.value

#add benjamini hochberg correction test to be able to include padj values
padj_values_CvM <- p.adjust(p_values_CvM, method = "BH")

#filter metabolites
all_metabolites_CvM <- data.frame(
  Metabolite = names(fc_data_CvM),
  FoldChange = fc_data_CvM,
  PValue = p_values_CvM,
  AdjustedPValue = padj_values_CvM
)

#save to output file
write.csv(all_metabolites_CvM, "all_metabolites_CvM.csv", row.names = FALSE)
#no significant metabolites (all above 0.05 p-value)


############## SvM ##############
##Initialize the data, conc = input data represents concentration values, stat = intended analysis is statistical analysis
mSet3 <-InitDataObjects("conc", "stat", FALSE);

#Load in dataset, rowu = rows represent metabolites, disc = discontinuous data for grouping
mSet3 <- Read.TextData(mSet3, "~/Desktop/Metabolomics/neg_stats_files/neg_SvM.csv", "rowu", "disc");

#check for missing values
sum(is.na(mSet3$dataSet$norm))
#no missing values

#check for missing values, invalid or formatting issues
mSet3<-SanityCheckData(mSet3);
#should output 2 groups here

#replace minimal values (missing or extremely low values)
mSet3<-ReplaceMin(mSet3);

#prepare data for normalization, making new prenorm object for mSet3
mSet3<-PreparePrenormData(mSet3);

#normalize data, NULL = no scaling, LogNorm = log transformation, meancenter = shifts data so each feature has mean of 0,
#ref = NULL because dont want to focus on one particular sample over another, could input sample information if normalizing to a specific reference
#ratioNum=20 = use top 20 most stable features for normalization
mSet3<-Normalization(mSet3, "NULL", "LogNorm", "MeanCenter", ref=NULL, ratio=FALSE, ratioNum=20);

#check dimensions of normalized data
print(dim(mSet3$dataSet$norm))
#10 456 - good

#plot normalized metabolites
mSet3 <-PlotNormSummary(mSet3, "SvM_norm_0_", format = "png", dpi=100, width=NA);
#plot normalized samples overall instead of by metabolites
mSet3 <-PlotSampleNormSummary(mSet3, "SvM_snorm_0_", format = "png", dpi=100, width=NA);


# Perform fold-change analysis on uploaded data, unpaired, using multiple as the baseline for analysis since CvM showed no significant metabolites and therefore are more similar than CvS metabolites
mSet3 <-FC.Anal(mSet3, 2.0, 0, FALSE)

# Plot fold-change analysis
mSet3 <-PlotFC(mSet3, "SvM_FC_", "png", 100, width=NA)

# To view fold-change 
head(mSet3$analSet$fc$fc.log)

# Perform ANOVA for mSet3 (SvM)
mSet3 <- ANOVA.Anal(mSet3, nonpar=F, thresh = 0.05, "fisher")
#output 13 significant features

#create list of metabolites to make a file of all metabolites with fc, p-value, padj, metabolite name

#extract fold change and p-value data
fc_data_SvM <- mSet3$analSet$fc$fc.log
p_values_SvM <- mSet3$analSet$aov$p.value

#benjamini hochberg correction with stats package to include padj values in output file
padj_values_SvM <- p.adjust(p_values_SvM, method = "BH")

#filter metabolites
all_metabolites_SvM <- data.frame(
  Metabolite = names(fc_data_SvM),
  FoldChange = fc_data_SvM,
  PValue = p_values_SvM,
  AdjustedPValue = padj_values_SvM
)

print(all_metabolites_SvM)

#save to output file
write.csv(all_metabolites_SvM, "all_metabolites_SvM.csv", row.names = FALSE)

###### make individual volcano plots for pairwise comparisons and then combine them #####

library(ggplot2)

# Volcano plot for Control vs Single #

#filter all_metabolites vector for metabolites with padj < 0.05
significant_metabolites_CvS <- all_metabolites_CvS[
  all_metabolites_CvS$AdjustPValue < 0.05, 
]

# Volcano plot for Control vs Single
volcano_data_CvS <- significant_metabolites_CvS

ggplot(volcano_data_CvS, aes(x = FoldChange, y = -log10(AdjustPValue))) +
  geom_point(alpha = 0.8, size = 2, color = "blue") +
  theme_minimal() +
  labs(
    title = "Control vs Single",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme(legend.title = element_blank())

# Volcano plot for Single vs Multiple #

#filter all_metabolites vector for metabolites with padj < 0.05
significant_metabolites_SvM <- all_metabolites_SvM[
  all_metabolites_SvM$AdjustedPValue < 0.05, 
]

# Volcano plot for Single vs Multiple
volcano_data_SvM <- significant_metabolites_SvM

ggplot(volcano_data_SvM, aes(x = FoldChange, y = -log10(AdjustedPValue))) +
  geom_point(alpha = 0.8, size = 2, color = "red") +
  theme_minimal() +
  labs(
    title = "Single vs Multiple",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme(legend.title = element_blank())

#combine the volcano plots
#install.packages("patchwork")
library(patchwork)

volcano_CvS_plot <- ggplot(volcano_data_CvS, aes(x = FoldChange, y = -log10(AdjustPValue))) +
  geom_point(alpha = 0.8, size = 2, color = "blue") +
  theme_minimal() +
  labs(
    title = "Control vs Single",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  ylim(0, 8) + #make y-axis limit 8
  xlim(-2, 6) + #set x-axis limits
  theme(legend.title = element_blank())

volcano_SvM_plot <- ggplot(volcano_data_SvM, aes(x = FoldChange, y = -log10(AdjustedPValue))) +
  geom_point(alpha = 0.8, size = 2, color = "red") +
  theme_minimal() +
  labs(
    title = "Single vs Multiple",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  ylim(0,8) + #make y-axis limit 8
  xlim(-2, 6) + #set x-axis limits
  theme(legend.title = element_blank())

# Combine plots
volcano_CvS_plot + volcano_SvM_plot

############# 2.7 PCA ###########

# Perform PCA analysis
mSet<-PCA.Anal(mSet)

PCA.Anal(mSet)
#use above to be able to find where the PC scores are located - output csv file too (metaboanalyst just does this by itself)
#then look for $analSet$pca$variance and $analSet$pca$cum.var to find coordinates, copy and paste into BBedit or Xcode to know the PC scores

library(ggplot2)

#extract PCA scores, not variance yet
pca_scores <- mSet$analSet$pca$x

#extract class info (treatment info)
group_labels <- as.vector(mSet$dataSet$meta.info$Class)

# Create a dataframe for ggplot2
pca_df <- data.frame(PC1 = pca_scores[,1],
                     PC2 = pca_scores[,2],
                     Group = group_labels)

# ggplot2 to make the PCA plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point(size = 3) +
  labs(title = NULL, #no broad title
       x = paste0("PC1 (", round(mSet$analSet$pca$variance[1] * 100, 1), "% variance)"), #extracting variance here for plot instead of x values
       y = paste0("PC2 (", round(mSet$analSet$pca$variance[2] * 100, 1), "% variance)")) +
  scale_color_manual(values = c("red", "blue", "yellow")) +
  theme_classic() #make no background




# the following code is the same as above however rerunnning analysis to exclude T2 samples - changing output files names #






####### Univariate methods ########

##Initialize the data, conc = input data represents concentration values, stat = intended analysis is statistical analysis
mSet<-InitDataObjects("conc", "stat", FALSE);

#Load in dataset, rowu = rows represent metabolites, disc = discontinuous data for grouping
mSet <- Read.TextData(mSet, "~/Desktop/Metabolomics/neg_stats_files/neg_all_but_no_T2.csv", "rowu", "disc");

#check for missing values
sum(is.na(mSet$dataSet$norm))
#no missing values

#check for missing values, invalid or formatting issues
mSet<-SanityCheckData(mSet);
#should say 3 groups here since big dataset includes control, single, multiple

#replace minimal values (missing or extremely low values)
mSet<-ReplaceMin(mSet);

#prepare data for normalization, overwriting old prenorm object / making a new one
mSet<-PreparePrenormData(mSet);

#normalize data, NULL = no scaling, LogNorm = log transformation, meancenter = shifts data so each feature has mean of 0,
#ref = NULL because dont want to focus on one particular sample over another, could input sample information if normalizing to a specific reference
#ratioNum=20 = use top 20 most stable features for normalization
mSet<-Normalization(mSet, "NULL", "LogNorm", "MeanCenter", ref=NULL, ratio=FALSE, ratioNum=20);

#check dimensions of normalized data
print(dim(mSet$dataSet$norm))
#13 456 - removed T2 samples

#plot normalized metabolites
mSet<-PlotNormSummary(mSet, "no_T2_norm_0_", format = "png", dpi=100, width=NA);
#plot normalized samples overall instead of by metabolites
mSet<-PlotSampleNormSummary(mSet, "no_T2_snorm_0_", format = "png", dpi=100, width=NA);

############### Fold change analysis ###############

# Perform fold-change analysis
# doing pairwise comparisons eventually
# tried to subset the data in R MULTIPLEEEE different ways and nothing worked correctly 
# So I went into the csv file, removed samples I didnt want (ex doing Control vs Multiple so removed single samples), initialized a new mSet, then loaded in the edited csv file to run each pairwise comparison
# because you are uploading a new dataset for each comparison, make sure to normalize, check for missing info, etc. before running FC analysis

# Perform fold-change analysis on uploaded data for the WHOLE dataset
#indicated fc threshold to be 2.0, cmp.type = 0 indicates group 1 - group 2 (not ideal since 3 groups here but not using fc anal for any analysis on the big dataset, eg will use for control vs. multiple where this makes sense)
mSet<-FC.Anal(mSet, 2.0, 0, FALSE)

# Plot fold-change analysis
mSet<-PlotFC(mSet, "no_T2_Total_FC_", "png", 100, width=NA)

# To view fold-change 
head(mSet$analSet$fc$fc.log)

######### ANOVA for big mSet (which contains control, single, multiple) ########

#no using non-parametric test so F because not assuming any distribution (parametric used when doing ANOVA so use the option nonpar=F), p-value threshold = 0.05, using Fisher's LSD test for post-hoc analysis
#options are fisher or tukey but those are for NOT using a non-parametric test (nonpar = F)
mSet <- ANOVA.Anal(mSet, nonpar=F, thresh=0.05, "fisher")
#18 significant features

### make list of metabolites and output as csv file ###
# want to include fc, pvalue, padj, metabolite name

#extract metabolites
anova_df <- mSet$analSet$aov$sig.mat
anova_df$Metabolite <- rownames(anova_df)
rownames(anova_df) <- NULL

#extract fold change and p-value data
fc_data <- mSet$analSet$fc$fc.log
p_values <- mSet$analSet$aov$p.value

#correction test here to add padj to output file, benjamini hochberg correction
#tried bonferroni correction but nothing was significant due to the high number of metabolites in the dataset so using BH method
#install.packages("stats")
library(stats)

padj_values <- p.adjust(p_values, method = "BH")

#filter metabolites
all_metabolites <- data.frame(
  Metabolite = names(fc_data),
  FoldChange = fc_data,
  PValue = p_values[names(fc_data)],
  AdjustedPValue = padj_values
)

print(all_metabolites)

#save to output file
write.csv(all_metabolites, "no_T2_all_metabolites_TOTAL.csv", row.names = FALSE)


################ CvS ###############
##Initialize the data, conc = input data represents concentration values, stat = intended analysis is statistical analysis
mSet1<-InitDataObjects("conc", "stat", FALSE);

#Load in dataset, rowu = rows represent metabolites, disc = discontinuous data for grouping
mSet1 <- Read.TextData(mSet1, "~/Desktop/Metabolomics/neg_stats_files/neg_CvS.csv", "rowu", "disc");
# the following is going to be the same as with all samples since T2 only included multiple treatment samples

#check for missing values
sum(is.na(mSet1$dataSet$norm))
#no missing values

#check for missing values, invalid or formatting issues
mSet1<-SanityCheckData(mSet1);
#should output that it found 2 groups, NOT 3!

#replace minimal values (missing or extremely low values)
mSet1<-ReplaceMin(mSet1);

#prepare data for normalization, making new prenorm object for this specific pairwise comparison
mSet1<-PreparePrenormData(mSet1);

#normalize data, NULL = no scaling, LogNorm = log transformation, meancenter = shifts data so each feature has mean of 0,
#ref = NULL because dont want to focus on one particular sample over another, could input sample information if normalizing to a specific reference
#ratioNum=20 = use top 20 most stable features for normalization
mSet1<-Normalization(mSet1, "NULL", "LogNorm", "MeanCenter", ref=NULL, ratio=FALSE, ratioNum=20);

#check dimensions of normalized data
print(dim(mSet1$dataSet$norm))
#output 10 456 

#plot normalized metabolites
mSet1<-PlotNormSummary(mSet1, "CvS_norm_0_", format = "png", dpi=100, width=NA);
#plot normalized samples overall instead of by metabolites
mSet1<-PlotSampleNormSummary(mSet1, "CvS_snorm_0_", format = "png", dpi=100, width=NA);

# Perform fold-change analysis on uploaded data, unpaired, cmp.type = 0 meaning group 1 - group 2 
# group 1 will be the first treatment in the csv file which was the control group
mSet1<-FC.Anal(mSet1, 2.0, 0, FALSE)

# Plot fold-change analysis
mSet1<-PlotFC(mSet1, "CvS_FC_", "png", 100, width=NA)

# To view fold-change 
head(mSet1$analSet$fc$fc.log)

####### Perform ANOVA for mSet1 (CvS) #######
mSet1 <- ANOVA.Anal(mSet1, nonpar=F, thresh = 0.05, "fisher")
#18 significant features

#find a list of all metabolites, output to file with log2foldchange, pvalue, padj
#extract fold change and p-value data
fc_data_CvS <- mSet1$analSet$fc$fc.log
p_values_CvS <- mSet1$analSet$aov$p.value

#benjamini hochberg correction to find adjusted p-values
padj_values_CvS <- p.adjust(p_values_CvS, method = "BH")

#filter metabolites
all_metabolites_CvS <- data.frame(
  Metabolite = names(fc_data_CvS),
  FoldChange = fc_data_CvS,
  PValue = p_values_CvS,
  AdjustPValue = padj_values_CvS
)

print(all_metabolites_CvS)

#save to output file
write.csv(all_metabolites_CvS, "all_metabolites_CvS.csv", row.names = FALSE)


############### CvM ###############
##Initialize the data, conc = input data represents concentration values, stat = intended analysis is statistical analysis
mSet2<-InitDataObjects("conc", "stat", FALSE);

#Load in dataset, rowu = rows represent metabolites, disc = discontinuous data for grouping
mSet2 <- Read.TextData(mSet2, "~/Desktop/Metabolomics/neg_stats_files/neg_CvM_no_T2.csv", "rowu", "disc");

#check for missing values
sum(is.na(mSet2$dataSet$norm))
#no missing values

#check for missing values, invalid or formatting issues
mSet2<-SanityCheckData(mSet2);
#should output 2 groups here

#replace minimal values (missing or extremely low values)
mSet2<-ReplaceMin(mSet2);

#prepare data for normalization, making new prenorm object for mSet2
mSet2<-PreparePrenormData(mSet2);

#normalize data, NULL = no scaling, LogNorm = log transformation, meancenter = shifts data so each feature has mean of 0,
#ref = NULL because dont want to focus on one particular sample over another, could input sample information if normalizing to a specific reference
#ratioNum=20 = use top 20 most stable features for normalization
mSet2<-Normalization(mSet2, "NULL", "LogNorm", "MeanCenter", ref=NULL, ratio=FALSE, ratioNum=20);

#check dimensions of normalized data
print(dim(mSet2$dataSet$norm))
#8 456

#plot normalized metabolites
mSet2 <-PlotNormSummary(mSet2, "no_T2_CvM_norm_0_", format = "png", dpi=100, width=NA);
#plot normalized samples overall instead of by metabolites
mSet2 <-PlotSampleNormSummary(mSet2, "no_T2_CvM_snorm_0_", format = "png", dpi=100, width=NA);

# Perform fold-change analysis on uploaded data, unpaired, control = group 1 so basis for fc analysis
mSet2<-FC.Anal(mSet2, 2.0, 0, FALSE)

# Plot fold-change analysis
mSet2 <-PlotFC(mSet2, "no_T2_CvM_FC_", "png", 100, width=NA)

# To view fold-change 
head(mSet2$analSet$fc$fc.log)

# Perform ANOVA for mSet2 (CvM)
mSet2 <- ANOVA.Anal(mSet2, nonpar=F, thresh = 0.05, "fisher")
#0 significant features

#find all metabolites and output to file, include pvalue, log2foldchange, padj

#extract fold change and p-value data
fc_data_CvM <- mSet2$analSet$fc$fc.log
p_values_CvM <- mSet2$analSet$aov$p.value

#add benjamini hochberg correction test to be able to include padj values
padj_values_CvM <- p.adjust(p_values_CvM, method = "BH")

#filter metabolites
all_metabolites_CvM <- data.frame(
  Metabolite = names(fc_data_CvM),
  FoldChange = fc_data_CvM,
  PValue = p_values_CvM,
  AdjustedPValue = padj_values_CvM
)

#save to output file
write.csv(all_metabolites_CvM, "no_T2_all_metabolites_CvM.csv", row.names = FALSE)
#no significant metabolites (all above 0.05 p-value)


############## SvM ##############
##Initialize the data, conc = input data represents concentration values, stat = intended analysis is statistical analysis
mSet3 <-InitDataObjects("conc", "stat", FALSE);

#Load in dataset, rowu = rows represent metabolites, disc = discontinuous data for grouping
mSet3 <- Read.TextData(mSet3, "~/Desktop/Metabolomics/neg_stats_files/neg_SvM_no_T2.csv", "rowu", "disc");

#check for missing values
sum(is.na(mSet3$dataSet$norm))
#no missing values

#check for missing values, invalid or formatting issues
mSet3<-SanityCheckData(mSet3);
#should output 2 groups here

#replace minimal values (missing or extremely low values)
mSet3<-ReplaceMin(mSet3);

#prepare data for normalization, making new prenorm object for mSet3
mSet3<-PreparePrenormData(mSet3);

#normalize data, NULL = no scaling, LogNorm = log transformation, meancenter = shifts data so each feature has mean of 0,
#ref = NULL because dont want to focus on one particular sample over another, could input sample information if normalizing to a specific reference
#ratioNum=20 = use top 20 most stable features for normalization
mSet3<-Normalization(mSet3, "NULL", "LogNorm", "MeanCenter", ref=NULL, ratio=FALSE, ratioNum=20);

#check dimensions of normalized data
print(dim(mSet3$dataSet$norm))
#8 456 - good

#plot normalized metabolites
mSet3 <-PlotNormSummary(mSet3, "no_T2_SvM_norm_0_", format = "png", dpi=100, width=NA);
#plot normalized samples overall instead of by metabolites
mSet3 <-PlotSampleNormSummary(mSet3, "no_T2_SvM_snorm_0_", format = "png", dpi=100, width=NA);


# Perform fold-change analysis on uploaded data, unpaired, using multiple as the baseline for analysis since CvM showed no significant metabolites and therefore are more similar than CvS metabolites
mSet3 <-FC.Anal(mSet3, 2.0, 0, FALSE)

# Plot fold-change analysis
mSet3 <-PlotFC(mSet3, "no_T2_SvM_FC_", "png", 100, width=NA)

# To view fold-change 
head(mSet3$analSet$fc$fc.log)

# Perform ANOVA for mSet3 (SvM)
mSet3 <- ANOVA.Anal(mSet3, nonpar=F, thresh = 0.05, "fisher")
#output 9 significant features

#create list of metabolites to make a file of all metabolites with fc, p-value, padj, metabolite name

#extract fold change and p-value data
fc_data_SvM <- mSet3$analSet$fc$fc.log
p_values_SvM <- mSet3$analSet$aov$p.value

#benjamini hochberg correction with stats package to include padj values in output file
padj_values_SvM <- p.adjust(p_values_SvM, method = "BH")

#filter metabolites
all_metabolites_SvM <- data.frame(
  Metabolite = names(fc_data_SvM),
  FoldChange = fc_data_SvM,
  PValue = p_values_SvM,
  AdjustedPValue = padj_values_SvM
)

print(all_metabolites_SvM)

#save to output file
write.csv(all_metabolites_SvM, "no_T2_all_metabolites_SvM.csv", row.names = FALSE)

###### make individual volcano plots for pairwise comparisons and then combine them #####

library(ggplot2)

# Volcano plot for Control vs Single #

#filter all_metabolites vector for metabolites with padj < 0.05
significant_metabolites_CvS <- all_metabolites_CvS[
  all_metabolites_CvS$AdjustPValue < 0.05, 
]

# Volcano plot for Control vs Single
volcano_data_CvS <- significant_metabolites_CvS

ggplot(volcano_data_CvS, aes(x = FoldChange, y = -log10(AdjustPValue))) +
  geom_point(alpha = 0.8, size = 2, color = "blue") +
  theme_minimal() +
  labs(
    title = "Control vs Single",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme(legend.title = element_blank())

# Volcano plot for Single vs Multiple #

#filter all_metabolites vector for metabolites with padj < 0.05
significant_metabolites_SvM <- all_metabolites_SvM[
  all_metabolites_SvM$AdjustedPValue < 0.05, 
]

# Volcano plot for Single vs Multiple
volcano_data_SvM <- significant_metabolites_SvM

ggplot(volcano_data_SvM, aes(x = FoldChange, y = -log10(AdjustedPValue))) +
  geom_point(alpha = 0.8, size = 2, color = "red") +
  theme_minimal() +
  labs(
    title = "Single vs Multiple",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme(legend.title = element_blank())

#combine the volcano plots
#install.packages("patchwork")
library(patchwork)

volcano_CvS_plot <- ggplot(volcano_data_CvS, aes(x = FoldChange, y = -log10(AdjustPValue))) +
  geom_point(alpha = 0.8, size = 2, color = "blue") +
  theme_minimal() +
  labs(
    title = "Control vs Single",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  ylim(0, 8) + #make y-axis limit 8
  xlim(-2, 6) + #set x-axis limits
  theme(legend.title = element_blank())

volcano_SvM_plot <- ggplot(volcano_data_SvM, aes(x = FoldChange, y = -log10(AdjustedPValue))) +
  geom_point(alpha = 0.8, size = 2, color = "red") +
  theme_minimal() +
  labs(
    title = "Single vs Multiple",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  ylim(0,8) + #make y-axis limit 8
  xlim(-2, 6) + #set x-axis limits
  theme(legend.title = element_blank())

# Combine plots
volcano_CvS_plot + volcano_SvM_plot

############# 2.7 PCA ###########

# Perform PCA analysis
mSet<-PCA.Anal(mSet)

PCA.Anal(mSet)
#use above to be able to find where the PC scores are located - output csv file too (metaboanalyst just does this by itself)
#then look for $analSet$pca$variance and $analSet$pca$cum.var to find coordinates, copy and paste into BBedit or Xcode to know the PC scores

library(ggplot2)

#extract PCA scores, not variance yet
pca_scores <- mSet$analSet$pca$x

#extract class info (treatment info)
group_labels <- as.vector(mSet$dataSet$meta.info$Class)

# Create a dataframe for ggplot2
pca_df <- data.frame(PC1 = pca_scores[,1],
                     PC2 = pca_scores[,2],
                     Group = group_labels)

# ggplot2 to make the PCA plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point(size = 5) +
  labs(title = NULL, #no broad title
       x = paste0("PC1 (", round(mSet$analSet$pca$variance[1] * 100, 1), "% variance)"), #extracting variance here for plot instead of x values
       y = paste0("PC2 (", round(mSet$analSet$pca$variance[2] * 100, 1), "% variance)"),
       color = "Treatment",
       shape = "Treatment"
       ) +
  scale_color_manual(values = c("red", "blue", "green")) +
  scale_y_continuous(limits = c(-8, 8)) + #set y-axis limits
  theme_classic() + #make no background
  theme(legend.background = element_rect(fill = "white", color = "black"))





