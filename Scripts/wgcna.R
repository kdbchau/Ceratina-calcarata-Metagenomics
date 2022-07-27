#===============================================================================
# Load Libraries and settings
#===============================================================================
library(WGCNA)
library(phyloseq)
library(ggplot2)
library(igraph)
library(vegan)
library(xlsx)
library(Hmisc)
library(caret)
options(stringsAsFactors = FALSE) # setting string not as factor
enableWGCNAThreads() # enable multithreading

# Set working directory
setwd("~/PostDoc/project1_2021/00-FINAL_DEC102021/Final_Tables/WGCNA_OneHotCoded/")
#===============================================================================
# Set some variable names here
#===============================================================================
rm(list=ls())
domain <- "all"
FoG <- "family" # "family" or "genus"
string <- "Annual Precipitation"
#===============================================================================
#   Load Data
#===============================================================================
# Load data and create bins
variable1 <- gsub(" ","_",tolower(string))
datExpr0 <- read.csv(paste0(domain,"_",FoG,"_",variable1,"_normalized_transpose.csv"), header=TRUE)

row.names(datExpr0) <- datExpr0$GridID # Push Sites into index
datExpr0$GridID <- NULL

datExpr <- datExpr0[,c(-1:-8)] # species as "genes" normalized counts
datTraits <- datExpr0[,c(1:8)] # the bin variable

# Isolate just variable of interest
if (variable1 == "annual_temperature") {
  datTraits <- datTraits[5] # temp
} else if (variable1 == "annual_precipitation") {
  datTraits <- datTraits[6] # precip
} else if (variable1 == "development_percent") {
  datTraits <- datTraits[7] # devperc
} else {
  datTraits <- datTraits[8] # agriperc
}

# Rename the binned column name into just "binned"
colnames(datTraits) <- "binned"

#===============================================================================
#   Check Data for Missing Values
#===============================================================================

gsg <- goodSamplesGenes(datExpr0)
gsg$allOK # If True, all "genes" passed the cuts. Otherwise, remove offending genes and samples.

# Run this if gsg$allOK is False
# if (!gsg$allOK)
# {
#   # Optionally, print the gene and sample names that were removed:
#   if (sum(!gsg$goodGenes)>0)
#     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
#   if (sum(!gsg$goodSamples)>0)
#     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
#   # Remove the offending genes and samples from the data:
#   datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
# }

#===============================================================================
#   Check for outliers
#===============================================================================

# Check that the "expression" table (datExpr) and "trait" table (datTraits) match up
table(rownames(datExpr)==rownames(datTraits)) # should say TRUE and # of samples

A <- adjacency(t(datExpr), type="distance") # Sample network based on Euclidean distances
k <- as.numeric(apply(A,2,sum))-1 # Define measure of connectivity
Z.k <- scale(k) # Standardize connectivity (make it between 0-1)
threshZ.k <- -5 # Define outliers if the site's Z.k is below a threshold
outlierColor <- ifelse(Z.k<threshZ.k, "red", "black") # Color vector indicates outliers (red)
sampleTree <- hclust(as.dist(1-A), method="average")

# Convert traits to a color rep
traitColors <- data.frame(numbers2colors(datTraits, signed = FALSE))
dimnames(traitColors)[[2]] <- names(datTraits)
datColors <- data.frame(outlierC=outlierColor, traitColors)

# Plot the dendogram
plotDendroAndColors(sampleTree, groupLabels = names(datColors),colors = datColors, 
                    main=paste(capitalize(domain),":",string,"binning"),cex.dendroLabels = 0.6)

# Remove Outliers
remove.samples <- Z.k<threshZ.k|is.na(Z.k)
datExpr <- datExpr[!remove.samples,]
datTraits <- as.data.frame(datTraits[!remove.samples,])
rownames(datTraits) <- rownames(datExpr)
colnames(datTraits) <- "binned" # Rename it back to "binned"
table(rownames(datTraits)==rownames(datExpr)) # should say TRUE
# Recompute the sample network among remaining samples
A <- adjacency(t(datExpr), type="distance")

# ReCalc cluster tree
sampleTree <- hclust(as.dist(1-A), method="average")
# Convert traits to a color rep
traitColors <- data.frame(numbers2colors(datTraits, signed = FALSE))
dimnames(traitColors)[[2]] <- names(datTraits)
datColors <- data.frame(traitColors)
plotDendroAndColors(sampleTree, groupLabels = names(datColors),colors = datColors,
                    main=paste(capitalize(domain),":",string,"binning"),cex.dendroLabels = 0.6)

collectGarbage()


#===============================================================================
#  One-Hot Code
#===============================================================================

# Rename bins into factors  ##
datTraits$binned[datTraits$binned == 1] <- "Very_Low"
datTraits$binned[datTraits$binned == 2] <- "Low"
datTraits$binned[datTraits$binned == 3] <- "Moderate"
datTraits$binned[datTraits$binned == 4] <- "High"
datTraits$binned[datTraits$binned == 5] <- "Very_High"

# One-hot code
datTraits_dmy <- dummyVars(" ~ .", data = datTraits)
trsf <- data.frame(predict(datTraits_dmy, newdata=datTraits))
trsf # the new datTraits
colnames(trsf)[length(trsf)] <- "Very_Low"
colnames(trsf)[length(trsf)-1] <- "Very_High"
colnames(trsf)[length(trsf)-2] <- "Moderate"
colnames(trsf)[length(trsf)-3] <- "Low"
colnames(trsf)[length(trsf)-4] <- "High"

datTraits <- trsf

#===============================================================================
#  Soft Threshold
#===============================================================================
powers <- c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.60,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#===============================================================================
#  Co expression Similarity, Adjacency and TOM
#===============================================================================
best_sft <- 6 
adj <- adjacency(datExpr, power=best_sft)

# Minimize noise and spurious associations, transform adj into Toplogical Overlap Matrix (TOM), get dissimilarity
TOM <- TOMsimilarity(adj)
dissTOM <- 1-TOM

# Now cluster using TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main="Gene clustering on TOM-based dissimilarity",
     labels=FALSE, hang=0.04)

# We like large modules so we set it to a high value like 30
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

# Now plot the module assignment under the gene dendogram
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#===============================================================================
#  Merge modules with similar expression profiles
#===============================================================================
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# To see the merging effects, plot the gene dendogram again with original and merged
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs

#===============================================================================
#  Relating Modules to External Traits (Heat Map)
#===============================================================================
# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


## Below is based on my dataset because I have binned variables so it is just
## repeated for each variable bin for me to get top modules
#===============================================================================
# VERY LOW
#===============================================================================
#===============================================================================
# VERY LOW: Gene Relationship to Trait and Important Modules (GS and MM)
#===============================================================================
# Define variable Very Low containing the Very Low column of datTrait
varofin = as.data.frame(datTraits$Very_Low)
names(varofin) = "Very_Low"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM.", modNames, sep="")
names(MMPvalue) = paste("p.MM.", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr, varofin, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(varofin), sep="")
names(GSPvalue) = paste("p.GS.", names(varofin), sep="")

#===============================================================================
# VERY LOW: Intramodular Analysis - Identify Genes with High GS and MM
#===============================================================================
modNames
sizeGrWindow(7,7)
par(mar=c(4,4,4,4))
par(mfrow = c(4,5))
for (module in modNames) {
  column <- match(module, modNames)
  moduleGenes <- moduleColors==module
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("MM in", capitalize(module), "Module"),
                     ylab = paste("GS for",string),
                     main = paste("Module Membership vs. Gene Significance\n", capitalize(domain)),
                     cex.main = 0.9, cex.lab = 1, cex.axis = 1, col="black")
}

#names(datExpr)
#names(datExpr)[moduleColors=="brown"]

annot <- read.csv(paste0("Taxonomy_Matrix/",domain,"_",FoG,"_taxmat_66sites.csv"))
probes <- names(datExpr)
probes2annot <- match(probes, annot$Domain)
sum(is.na(probes2annot)) # should say number of families/genera
# Need to add annot$Genus when working with genus level
geneInfo0 <- data.frame(annot$Domain, annot$Kingdom, annot$Phylum, annot$Class, annot$Order, 
                        annot$Family, annot$Genus,
                        moduleColor = moduleColors, geneTraitSignificance, GSPvalue)

# Order modules by their significance for varofin
modOrder = order(-abs(cor(MEs, varofin, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order families in the geneINfo variable first by module colour then by geneTraitSignificance
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Very_Low))
geneInfo <- geneInfo0[geneOrder,]
geneInfo <- tibble::rownames_to_column(geneInfo, "FAMILY2")
geneInfo$FAMILY2 <- NULL
names(geneInfo)[1] <- "Domain" # only the "all" dataframe has this category
# Change to 2:6 and remove "Genus" when doing family-level dataframe
#Family
#names(geneInfo)[2:6] <- c("Kingdom","Phylum","Class","Order","Family")
# Genus
names(geneInfo)[2:7] <- c("Kingdom","Phylum","Class","Order","Family","Genus") # rename these

# Write the data out
var1 <- variable1 # variable name and if we used continuous or binned version
var2 <- paste0("_modsize",as.character(minModuleSize)) # what was the module size used in net?
var3 <- paste0("_p",as.character(best_sft)) # what power did we use?
var4 <- paste0("_",sub("_","",names(varofin)))
write.xlsx(geneInfo, file=paste0(domain,"_",FoG,"_normalized_",var1,var2,var3,var4,".xlsx"), 
           row.names = FALSE)


#===============================================================================
# LOW
#===============================================================================
#===============================================================================
# LOW: Gene Relationship to Trait and Important Modules (GS and MM)
#===============================================================================
# Define variable Very Low containing the Very Low column of datTrait
varofin = as.data.frame(datTraits$Low)
names(varofin) = "Low"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM.", modNames, sep="")
names(MMPvalue) = paste("p.MM.", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr, varofin, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(varofin), sep="")
names(GSPvalue) = paste("p.GS.", names(varofin), sep="")

#===============================================================================
# LOW: Intramodular Analysis - Identify Genes with High GS and MM
#===============================================================================
modNames
sizeGrWindow(7,7)
par(mar=c(4,4,4,4))
par(mfrow = c(4,5))
for (module in modNames) {
  column <- match(module, modNames)
  moduleGenes <- moduleColors==module
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("MM in", capitalize(module), "Module"),
                     ylab = paste("GS for",string),
                     main = paste("Module Membership vs. Gene Significance\n", capitalize(domain)),
                     cex.main = 0.9, cex.lab = 1, cex.axis = 1, col="black")
}

#names(datExpr)
#names(datExpr)[moduleColors=="brown"]

annot <- read.csv(paste0("Taxonomy_Matrix/",domain,"_",FoG,"_taxmat_66sites.csv"))
probes <- names(datExpr)
probes2annot <- match(probes, annot$Domain)
sum(is.na(probes2annot)) # should say number of families/genera
# Need to add annot$Genus when working with genus level
geneInfo0 <- data.frame(annot$Domain, annot$Kingdom, annot$Phylum, annot$Class, annot$Order, 
                        annot$Family, annot$Genus,
                        moduleColor = moduleColors, geneTraitSignificance, GSPvalue)

# Order modules by their significance for varofin
modOrder = order(-abs(cor(MEs, varofin, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
# Order families in the geneINfo variable first by module colour then by geneTraitSignificance
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Low))
geneInfo <- geneInfo0[geneOrder,]
geneInfo <- tibble::rownames_to_column(geneInfo, "FAMILY2")
geneInfo$FAMILY2 <- NULL
names(geneInfo)[1] <- "Domain" # only the "all" dataframe has this category
#Family
#names(geneInfo)[2:6] <- c("Kingdom","Phylum","Class","Order","Family")
# Genus
names(geneInfo)[2:7] <- c("Kingdom","Phylum","Class","Order","Family","Genus") # rename these

# Write the data out
var1 <- variable1 # variable name and if we used continuous or binned version
var2 <- paste0("_modsize",as.character(minModuleSize)) # what was the module size used in net?
var3 <- paste0("_p",as.character(best_sft)) # what power did we use?
var4 <- paste0("_",sub("_","",names(varofin)))
write.xlsx(geneInfo, file=paste0(domain,"_",FoG,"_normalized_",var1,var2,var3,var4,".xlsx"), 
           row.names = FALSE)

#===============================================================================
# MODERATE
#===============================================================================
#===============================================================================
# MODERATE: Gene Relationship to Trait and Important Modules (GS and MM)
#===============================================================================
# Define variable Very Low containing the Very Low column of datTrait
varofin = as.data.frame(datTraits$Moderate)
names(varofin) = "Moderate"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM.", modNames, sep="")
names(MMPvalue) = paste("p.MM.", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr, varofin, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(varofin), sep="")
names(GSPvalue) = paste("p.GS.", names(varofin), sep="")

#===============================================================================
# MODERATE: Intramodular Analysis - Identify Genes with High GS and MM
#===============================================================================
modNames
sizeGrWindow(7,7)
par(mar=c(4,4,4,4))
par(mfrow = c(4,5))
for (module in modNames) {
  column <- match(module, modNames)
  moduleGenes <- moduleColors==module
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("MM in", capitalize(module), "Module"),
                     ylab = paste("GS for",string),
                     main = paste("Module Membership vs. Gene Significance\n", capitalize(domain)),
                     cex.main = 0.9, cex.lab = 1, cex.axis = 1, col="black")
}

#names(datExpr)
#names(datExpr)[moduleColors=="brown"]

annot <- read.csv(paste0("Taxonomy_Matrix/",domain,"_",FoG,"_taxmat_66sites.csv"))
probes <- names(datExpr)
probes2annot <- match(probes, annot$Domain)
sum(is.na(probes2annot)) # should say number of families/genera
# Need to add annot$Genus when working with genus level
geneInfo0 <- data.frame(annot$Domain, annot$Kingdom, annot$Phylum, annot$Class, annot$Order, 
                        annot$Family, annot$Genus,
                        moduleColor = moduleColors, geneTraitSignificance, GSPvalue)

# Order modules by their significance for varofin
modOrder = order(-abs(cor(MEs, varofin, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
# Order families in the geneINfo variable first by module colour then by geneTraitSignificance
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Moderate))
geneInfo <- geneInfo0[geneOrder,]
geneInfo <- tibble::rownames_to_column(geneInfo, "FAMILY2")
geneInfo$FAMILY2 <- NULL
names(geneInfo)[1] <- "Domain" # only the "all" dataframe has this category
#Family
#names(geneInfo)[2:6] <- c("Kingdom","Phylum","Class","Order","Family")
# Genus
names(geneInfo)[2:7] <- c("Kingdom","Phylum","Class","Order","Family","Genus") # rename these

# Write the data out
var1 <- variable1 # variable name and if we used continuous or binned version
var2 <- paste0("_modsize",as.character(minModuleSize)) # what was the module size used in net?
var3 <- paste0("_p",as.character(best_sft)) # what power did we use?
var4 <- paste0("_",sub("_","",names(varofin)))
write.xlsx(geneInfo, file=paste0(domain,"_",FoG,"_normalized_",var1,var2,var3,var4,".xlsx"), 
           row.names = FALSE)

#===============================================================================
# HIGH
#===============================================================================
#===============================================================================
# HIGH: Gene Relationship to Trait and Important Modules (GS and MM)
#===============================================================================
# Define variable Very Low containing the Very Low column of datTrait
varofin = as.data.frame(datTraits$High)
names(varofin) = "High"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM.", modNames, sep="")
names(MMPvalue) = paste("p.MM.", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr, varofin, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(varofin), sep="")
names(GSPvalue) = paste("p.GS.", names(varofin), sep="")

#===============================================================================
# HIGH: GIntramodular Analysis - Identify Genes with High GS and MM
#===============================================================================
modNames
sizeGrWindow(7,7)
par(mar=c(4,4,4,4))
par(mfrow = c(4,5))
for (module in modNames) {
  column <- match(module, modNames)
  moduleGenes <- moduleColors==module
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("MM in", capitalize(module), "Module"),
                     ylab = paste("GS for",string),
                     main = paste("Module Membership vs. Gene Significance\n", capitalize(domain)),
                     cex.main = 0.9, cex.lab = 1, cex.axis = 1, col="black")
}

#names(datExpr)
#names(datExpr)[moduleColors=="brown"]

annot <- read.csv(paste0("Taxonomy_Matrix/",domain,"_",FoG,"_taxmat_66sites.csv"))
probes <- names(datExpr)
probes2annot <- match(probes, annot$Domain)
sum(is.na(probes2annot)) # should say number of families/genera
# Need to add annot$Genus when working with genus level
geneInfo0 <- data.frame(annot$Domain, annot$Kingdom, annot$Phylum, annot$Class, annot$Order, 
                        annot$Family, annot$Genus,
                        moduleColor = moduleColors, geneTraitSignificance, GSPvalue)

# Order modules by their significance for varofin
modOrder = order(-abs(cor(MEs, varofin, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
# Order families in the geneINfo variable first by module colour then by geneTraitSignificance
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.High))
geneInfo <- geneInfo0[geneOrder,]
geneInfo <- tibble::rownames_to_column(geneInfo, "FAMILY2")
geneInfo$FAMILY2 <- NULL
names(geneInfo)[1] <- "Domain" # only the "all" dataframe has this category
#Family
#names(geneInfo)[2:6] <- c("Kingdom","Phylum","Class","Order","Family")
# Genus
names(geneInfo)[2:7] <- c("Kingdom","Phylum","Class","Order","Family","Genus") # rename these

# Write the data out
var1 <- variable1 # variable name and if we used continuous or binned version
var2 <- paste0("_modsize",as.character(minModuleSize)) # what was the module size used in net?
var3 <- paste0("_p",as.character(best_sft)) # what power did we use?
var4 <- paste0("_",sub("_","",names(varofin)))
write.xlsx(geneInfo, file=paste0(domain,"_",FoG,"_normalized_",var1,var2,var3,var4,".xlsx"), 
           row.names = FALSE)

#===============================================================================
# VERY HIGH
#===============================================================================
#===============================================================================
# VERY HIGH: Gene Relationship to Trait and Important Modules (GS and MM)
#===============================================================================
# Define variable Very Low containing the Very Low column of datTrait
varofin = as.data.frame(datTraits$Very_High)
names(varofin) = "Very_High"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM.", modNames, sep="")
names(MMPvalue) = paste("p.MM.", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr, varofin, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(varofin), sep="")
names(GSPvalue) = paste("p.GS.", names(varofin), sep="")

#===============================================================================
# VERY HIGH: Intramodular Analysis - Identify Genes with High GS and MM
#===============================================================================
modNames
sizeGrWindow(7,7)
par(mar=c(4,4,4,4))
par(mfrow = c(4,5))
for (module in modNames) {
  column <- match(module, modNames)
  moduleGenes <- moduleColors==module
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("MM in", capitalize(module), "Module"),
                     ylab = paste("GS for",string),
                     main = paste("Module Membership vs. Gene Significance\n", capitalize(domain)),
                     cex.main = 0.9, cex.lab = 1, cex.axis = 1, col="black")
}

#names(datExpr)
#names(datExpr)[moduleColors=="brown"]

annot <- read.csv(paste0("Taxonomy_Matrix/",domain,"_",FoG,"_taxmat_66sites.csv"))
probes <- names(datExpr)
probes2annot <- match(probes, annot$Domain)
sum(is.na(probes2annot)) # should say number of families/genera
# Need to add annot$Genus when working with genus level
geneInfo0 <- data.frame(annot$Domain, annot$Kingdom, annot$Phylum, annot$Class, annot$Order, 
                        annot$Family, annot$Genus,
                        moduleColor = moduleColors, geneTraitSignificance, GSPvalue)

# Order modules by their significance for varofin
modOrder = order(-abs(cor(MEs, varofin, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
# Order families in the geneINfo variable first by module colour then by geneTraitSignificance
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Very_High))
geneInfo <- geneInfo0[geneOrder,]
geneInfo <- tibble::rownames_to_column(geneInfo, "FAMILY2")
geneInfo$FAMILY2 <- NULL
names(geneInfo)[1] <- "Domain" # only the "all" dataframe has this category
#Family
#names(geneInfo)[2:6] <- c("Kingdom","Phylum","Class","Order","Family")
# Genus
names(geneInfo)[2:7] <- c("Kingdom","Phylum","Class","Order","Family","Genus") # rename these

# Write the data out
var1 <- variable1 # variable name and if we used continuous or binned version
var2 <- paste0("_modsize",as.character(minModuleSize)) # what was the module size used in net?
var3 <- paste0("_p",as.character(best_sft)) # what power did we use?
var4 <- paste0("_",sub("_","",names(varofin)))
write.xlsx(geneInfo, file=paste0(domain,"_",FoG,"_normalized_",var1,var2,var3,var4,".xlsx"), 
           row.names = FALSE)
