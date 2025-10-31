#pipeline coming from: https://github.com/xia-lab/MetaboAnalystR/blob/master/README.md
#mandatory tools: xcode and GNU Fortran compiler on laptop
#Xcode downloaded from Apple App Store
#use the following code in command line to install GNU fortran:
#step 1: install homebrew
#/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
#step 2: verify homebrew installation
#brew --version
#step 3: install gfortran using homebrew
#brew install gcc
#step 4: verify installation
#gfortran --version

#download all package dependencies below!!! big code didnt work for the most part so had to install separately

#example code for using BiocManager if needed, many of them needed this
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("devtools")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("RBGL", force = TRUE)
#library(RBGL)

#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("BiocParallel", force = TRUE)
#library(BiocParallel)

#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("edgeR", force = TRUE)
#library(edgeR)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("fgsea", force = TRUE)
#library(fgsea)

#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("impute", force = TRUE)
#library(impute)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("pcaMethods", force = TRUE)
#library(pcaMethods)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("MSnbase", force = TRUE)
#library(MSnbase)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("siggenes", force = TRUE)
#library(siggenes)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("qvalue", force = TRUE)
#library(qvalue)


#SSPA was so annoying, download tar from https://www.bioconductor.org/packages//2.10/bioc/html/SSPA.html
#extract tarball on laptop, go into Namespace and remove the lines of code with the @<- (remove this whole line: importFrom(methods, "@<-", new)
#then on command line, cd into directory of unpacked tarball, run code 'R CMD build .'
#now can install SSPA with the modified tarball using example code below
#install.packages("/Users/abbybednarski/Downloads/SSPA/SSPA_1.12.0.tar.gz", repos = NULL, type = "source", dependencies = TRUE)
#library(SSPA)

#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("edgeR")

#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("impute")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("sva")
#library(sva)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("randomForest")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("lars")
#library(lars)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("metap")
#library(metap)

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("GlobalAncova")
#library(GlobalAncova)

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("KEGGgraph")
#library(KEGGgraph)

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("globaltest")
#library(globaltest)

#library(Rgraphviz)

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("genefilter")
#library(genefilter)

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("reshape")
#library(reshape)

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("spls")
#library(spls)

#install using following code: followed chatgpt a lot because have to specifically change paths for where gfortran is on YOUR laptop (sorry that is not helpful
#need specific library links for where gfortran is on your labtop so had to use comman line to find the paths, make R use those specific paths, then download
#Sys.setenv(LDFLAGS="-L/usr/local/gfortran/lib -L/usr/local/Cellar/gcc/14.2.0_1/lib/gcc/14")
#Sys.setenv(PKG_LIBS="-lgfortran -lquadmath")
#pkgbuild::check_build_tools(debug = TRUE)
#remotes::install_github("xia-lab/MetaboAnalystR")
#tutorials required vignettes to be built so should add "build_vignettes = TRUE" to the above command (i didnt know that before i downloaded whoops)

library(devtools)
library(RBGL)
library(fgsea)
library(pcaMethods)
library(MSnbase)
library(siggenes)
library(Biobase)
library(edgeR)
library(BiocParallel)
library(impute)
library(sva)
library(KEGGgraph)
library(genefilter)
library(Rgraphviz)
library(pls)
library(ellipse)
library(scatterplot3d)
library(randomForest)
library(som)
library(RJSONIO)
library(ROCR)
library(fitdistrplus)
library(lars)
library(reshape)
library(spls)
library(metap)
library(entropy)
library(rsm)
library(globaltest)
library(GlobalAncova)
library(qvalue)

library(MetaboAnalystR)

