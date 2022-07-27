# CoNet

#==================================== Install packages ======================================================
#install_github("hallucigenia-sparsa/seqtime")
library(seqtime)
#install_github("ramellose/CoNetinR")
library(CoNetinR)  

library(gdata)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DirichletMultinomial")
#install_github("hallucigenia-sparsa/seqgroup")
library(seqgroup)

#install.packages("remotes")
#remotes::install_github("Tong-Chen/YSX")
library(dplyr)
library(remotes)
library(reshape2)
library(psych)
library(ggplot2)
library (otuSummary)
#remotes::install_github("taowenmicro/ggClusterNet")
library(ggClusterNet)

#==================================== CoNet =================================================================
setwd("~/PostDoc/project1_2021/00-FINAL_DEC102021/Final_Tables/PopGenandPathogens/AllPathogens/")

# Import data file
data <- read.csv("conet_data.csv", header=TRUE, row.names=1) # 13042 rows because 162 pariwise samples used; sample by variable

# Push Column 1 into index
row.names(data) <- data$X  # this was just a column with numbers 1 to 13041
data$X <- NULL
data_matrix <- as.matrix(t(data)) # convert to variable by sample

#Calculate the CoNet scores
coNetScores = getNetwork(data_matrix,
                         method="spearman","bray","pearson","kld",
                         bh=T,
                         iters=100,
                         report.full=T)
coNetScores = coNetScores$scores

#Calculate the CoNet p values
coNetPValues=getNetwork(data_matrix,
                        method="spearman","bray","pearson","kld",
                        pval.cor=T,
                        bh=T,
                        iters=100,
                        report.full=T)
#warnings()
coNetPValues=coNetPValues$pvalues
coNetPValuesdf = as.data.frame(coNetPValues)

#write out the p values in excel
write.csv(coNetPValuesdf,"coNet_pvalues3_samplebyvar.csv")

#Convert the CoNet scores into a matrix, will output into excel
N=4
adjmatrix=matrix(nrow=4,ncol=4)
adjmatrix[lower.tri(adjmatrix)] = coNetScores
adjmatrix = t(adjmatrix)
adjmatrix[lower.tri(adjmatrix)] = coNetScores
for (i in 1:N){
  for (j in 1:N){
    if (is.na(adjmatrix[i,j])){
      adjmatrix[i,j] = 0
    }
    else if (coNetPValues[i,j] > 1){
      adjmatrix[i,j] = 0
    }
  }
}

adjmatrixdf = as.data.frame(adjmatrix)
#outputs the CoNet correlation values
write.csv(adjmatrixdf,"coNet_finalized_4.csv")
