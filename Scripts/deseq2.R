## Negative Binomial. You flip a coin repeatedly and count the # of times the coin lands on heads.
## You keep flipping until it landed 5 times on head. This is a negative binomial experiment because:
## 1. There are repeated trials (e.g. biological replicates/samples)
## 2. Each trial has two possible outcomes (e.g. present/absent, upregulated/downregulated)
## 3. Probability of each outcome is equal
## 4. Each trial is independent (samples are taken independent of each other)
## 5. Experiment continues until fixed number of successes occurred (i.e. 5 heads in our example)
## NB is good for this dataset 

########## Install DEseq2 for R 4.1+ ################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

#----------------------------------------------------------------------------------------#
#-------------------------- Microbiomes by Variables of Interest ------------------------#
#----------------------------------------------------------------------------------------#
library(DESeq2)
library(dplyr)

# Load data and variables
rm(list=ls())
setwd("~/PostDoc/project1_2021/00-FINAL_DEC102021/Final_Tables/")
domain <- "all" # Domain of interest
FoG <- "family" # Family or Genus?
data <- read.csv(paste0(domain,"_",FoG,"_abun_0.1_Dec242021_nozero.csv"), header=TRUE)

# Get rid of sites that have zero individuals
data <- data[apply(data[,c(-1:-9)],1,function(x) !all(x)),]
# Get rid of families that sum up to 0
data <- data %>% select(all_of(names(.)[c(1:9)]),
                        where(~ is.numeric(.) && sum(., na.rm=TRUE) > 0))

#### ------------------- Development Percent ------------------- ####

column_descriptor <- data$binned_devperc # Variable of interest
string <- "Development Percent" # Variable of interest title string
data$GridID <- paste(data$GridID, data$binned_devperc, sep="_") # update Sites names with bin number
row.names(data) <- data$GridID  # Push Sites into index column
data <- data[order(data$binned_devperc),] # Sort dataframe by the binned column

# Abundance dataframe
data_sp <- data[,c(-1:-9)]
data_sp.counts <- as.data.frame(t(data_sp))

# Metadata and binned column renaming
data_env <- data[,c(1:9)]
# data_env$binned[data_env$binned_devperc == 1] <- "Very_Low"
# data_env$binned[data_env$binned_devperc == 2] <- "Low"
# data_env$binned[data_env$binned_devperc == 3] <- "Moderate"
# data_env$binned[data_env$binned_devperc == 4] <- "High"
# data_env$binned[data_env$binned_devperc == 5] <- "Very_High"
data_env.coldata <- data.frame(rows=colnames(data_sp.counts), condition=as.factor(data_env$binned_devperc))
data_env.coldata$rows <- as.character(data_env.coldata$rows)
data_env.coldata$condition <- as.factor(data_env.coldata$condition)

#### ------------------- Annual Precipitation ------------------- ####

column_descriptor <- data$binned_precip # Variable of interest
string <- "Annual Precipitation" # Variable of interest title string
data$GridID <- paste(data$GridID, data$binned_precip, sep="_") # update Sites names with bin number
row.names(data) <- data$GridID  # Push Sites into index column
data <- data[order(data$binned_precip),] # Sort dataframe by the binned column

# Abundance dataframe
data_sp <- data[,c(-1:-9)]
data_sp.counts <- as.data.frame(t(data_sp))

# Metadata and binned column renaming
data_env <- data[,c(1:9)]
# data_env$binned[data_env$binned_precip == 1] <- "Very_Low"
# data_env$binned[data_env$binned_precip == 2] <- "Low"
# data_env$binned[data_env$binned_precip == 3] <- "Moderate"
# data_env$binned[data_env$binned_precip == 4] <- "High"
# data_env$binned[data_env$binned_precip == 5] <- "Very_High"
data_env.coldata <- data.frame(rows=colnames(data_sp.counts), condition=as.factor(data_env$binned_precip))
data_env.coldata$rows <- as.character(data_env.coldata$rows)
data_env.coldata$condition <- as.factor(data_env.coldata$condition)

#### ------------------- Annual Temperature ------------------- ####

column_descriptor <- data$binned_temp # Variable of interest
string <- "Annual Temperature" # Variable of interest title string
data$GridID <- paste(data$GridID, data$binned_temp, sep="_") # update Sites names with bin number
row.names(data) <- data$GridID  # Push Sites into index column
data <- data[order(data$binned_temp),] # Sort dataframe by the binned column

# Abundance dataframe
data_sp <- data[,c(-1:-9)]
data_sp.counts <- as.data.frame(t(data_sp))

# Metadata and binned column renaming
data_env <- data[,c(1:9)]
# data_env$binned[data_env$binned_temp == 1] <- "Very_Low"
# data_env$binned[data_env$binned_temp == 2] <- "Low"
# data_env$binned[data_env$binned_temp == 3] <- "Moderate"
# data_env$binned[data_env$binned_temp == 4] <- "High"
# data_env$binned[data_env$binned_temp == 5] <- "Very_High"
data_env.coldata <- data.frame(rows=colnames(data_sp.counts), condition=as.factor(data_env$binned_temp))
data_env.coldata$rows <- as.character(data_env.coldata$rows)
data_env.coldata$condition <- as.factor(data_env.coldata$condition)


#### ------------------- Agricultural Percent ------------------- ####
column_descriptor <- data$binned_agriperc # Variable of interest
string <- "Agricultural Percent" # Variable of interest title string
data$GridID <- paste(data$GridID, data$binned_agriperc, sep="_") # update Sites names with bin number
row.names(data) <- data$GridID  # Push Sites into index column
data <- data[order(data$binned_agriperc),] # Sort dataframe by the binned column

# Abundance dataframe
data_sp <- data[,c(-1:-9)]
data_sp.counts <- as.data.frame(t(data_sp))

# Metadata and binned column renaming
data_env <- data[,c(1:9)]
# data_env$binned[data_env$binned_temp == 1] <- "Very_Low"
# data_env$binned[data_env$binned_temp == 2] <- "Low"
# data_env$binned[data_env$binned_temp == 3] <- "Moderate"
# data_env$binned[data_env$binned_temp == 4] <- "High"
# data_env$binned[data_env$binned_temp == 5] <- "Very_High"
data_env.coldata <- data.frame(rows=colnames(data_sp.counts), condition=as.factor(data_env$binned_agriperc))
data_env.coldata$rows <- as.character(data_env.coldata$rows)
data_env.coldata$condition <- as.factor(data_env.coldata$condition)


################################### Convert into DEseqDataSet #########################################
dds <- DESeqDataSetFromMatrix(countData = data_sp.counts, 
                              colData = data_env.coldata, 
                              design = ~ condition)
#dds
#as.data.frame(colData(dds)) # just check how the data looks

################################## Running the pipeline ###############################################

dds <- DESeq(dds) # will fail if there isn't at least 1 family with values > 0 in all sites

# library(pheatmap)
# ntd <- normTransform(dds)
# select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing = TRUE)
# df <- as.data.frame(colData(dds)[,c("condition","rows")])
# pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_col=FALSE,
#          annotation_col=df)


colData(dds)
res <- results(dds)
res

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
dds

normalized_counts <- counts(dds, normalized=TRUE)

variable <- gsub(" ","_",tolower(string))
write.csv(normalized_counts, 
          file=paste0("DEseq2/",domain,"_",FoG,"_",variable,"_normalized.csv"),row.names=TRUE)

################################## Comparing the bins ###############################################

# ---- Bin 1 vs Bin 5 ---- Very Low vs. Very High -------#
b1b5 <- as.data.frame(results(dds, contrast=c("condition",1,5)))
write.csv(b1b5, 
          file=paste0("DEseq2/Genus/",domain,"_",FoG,"_",variable,"_b1b5.csv"),row.names=TRUE)

# ---- Bin 1 vs Bin 4 ---- Very Low vs. High -------#
b1b4 <- as.data.frame(results(dds, contrast=c("condition",1,4)))
write.csv(b1b4, 
          file=paste0("DEseq2/Genus/",domain,"_",FoG,"_",variable,"_b1b4.csv"),row.names=TRUE)

# ---- Bin 1 vs Bin 3 ---- Very Low vs. Moderate -------#
b1b3 <- as.data.frame(results(dds, contrast=c("condition",1,3)))
write.csv(b1b3, 
          file=paste0("DEseq2/Genus/",domain,"_",FoG,"_",variable,"_b1b3.csv"),row.names=TRUE)

# ---- Bin 1 vs Bin 2 ---- Very Low vs. Low -------#
b1b2 <- as.data.frame(results(dds, contrast=c("condition",1,2)))
write.csv(b1b2, 
          file=paste0("DEseq2/Genus/",domain,"_",FoG,"_",variable,"_b1b2.csv"),row.names=TRUE)

# ---- Bin 2 vs Bin 5 ---- Low vs.Very High -------#
b2b5 <- as.data.frame(results(dds, contrast=c("condition",2,5)))
write.csv(b2b5, 
          file=paste0("DEseq2/Genus/",domain,"_",FoG,"_",variable,"_b2b5.csv"),row.names=TRUE)

# ---- Bin 2 vs Bin 4 ---- Low vs.High -------#
b2b4 <- as.data.frame(results(dds, contrast=c("condition",2,4)))
write.csv(b2b4, 
          file=paste0("DEseq2/Genus/",domain,"_",FoG,"_",variable,"_b2b4.csv"),row.names=TRUE)

# ---- Bin 2 vs Bin 3 ---- Low vs. Moderate -------#
b2b3 <- as.data.frame(results(dds, contrast=c("condition",2,3)))
write.csv(b2b3, 
          file=paste0("DEseq2/Genus/",domain,"_",FoG,"_",variable,"_b2b3.csv"),row.names=TRUE)

# ---- Bin 3 vs Bin 5 ---- Moderate vs. Very High -------#
b3b5 <- as.data.frame(results(dds, contrast=c("condition",3,5)))
write.csv(b3b5, 
          file=paste0("DEseq2/Genus/",domain,"_",FoG,"_",variable,"_b3b5.csv"),row.names=TRUE)

# ---- Bin 3 vs Bin 4 ---- Moderate vs.High -------#
b3b4 <- as.data.frame(results(dds, contrast=c("condition",3,4)))
write.csv(b3b4, 
          file=paste0("DEseq2/Genus/",domain,"_",FoG,"_",variable,"_b3b4.csv"),row.names=TRUE)

# ---- Bin 4 vs Bin 5 ---- High vs. Very High -------#
b4b5 <- as.data.frame(results(dds, contrast=c("condition",4,5)))
write.csv(b4b5, 
          file=paste0("DEseq2/Genus/",domain,"_",FoG,"_",variable,"_b4b5.csv"),row.names=TRUE)
