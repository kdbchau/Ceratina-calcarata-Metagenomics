# Diversity States
# This is for supplementary table S11 - Diversity using Bray-Curtis

library(vegan)
library(ggplot2)
setwd("~/PostDoc/project1_2021/00-FINAL_DEC102021/Final_Tables/")

#====================================================================================================
#   UPLOAD DATA - MAKE DATAFRAMES
#====================================================================================================
rm(list=ls())
FoG = "genus"
string = "Annual Precipitation"  # Variable of interest
variable1 <- gsub(" ","_",tolower(string))

data <- read.csv(paste0("all_",FoG,"_abun_0.1_Dec242021_nozero.csv"), header=TRUE)

# Convert GridID names to include variable of interest bin value
if (variable1 == "annual_temperature") {
  varofin <- data$binned_temp
} else if (variable1 == "annual_precipitation") {
  varofin <- data$binned_precip 
} else if (variable1 == "developmental_percent") {
  varofin <- data$binned_devperc
} else {
  varofin <- data$binned_agriperc
}
data$GridID <- paste(data$GridID, varofin, sep="_") # Rename GridID with binned variable
row.names(data) <- data$GridID # Push Sites into index column
data$GridID <- NULL # Drop the sites column

# Make dataframes
data_sp <- data[,c(-1:-8)] # Dataframe of just species
data_env <- data[,c(1:8)]  # Dataframe of just environmental variables

# Rename the binned column of interest to just bin for easier handling
# Isolate just variable of interest
if (variable1 == "annual_temperature") {
  names(data_env)[names(data_env) == "binned_temp"] <- "binned"
} else if (variable1 == "annual_precipitation") {
  names(data_env)[names(data_env) == "binned_precip"] <- "binned" 
} else if (variable1 == "developmental_percent") {
  names(data_env)[names(data_env) == "binned_devperc"] <- "binned"
} else {
  names(data_env)[names(data_env) == "binned_agriperc"] <- "binned"
}

# Check that the two dataframes are in the same order, should say TRUE
all.equal(rownames(data_sp), rownames(data_env))

#====================================================================================================
#   BRAY CURTIS NMDS
#====================================================================================================
# Break into groups for individual domain analysis; don't run following 3 lines when using "all"
data_sp <- data_sp[,grepl("F_", names(data_sp))]
# Remove rows (i.e. sites) with 0 species
data_sp <- data_sp[apply(data_sp,1,function(x) !all(x==0)),]

# Keep rows in data_env that match to the new data_sp
data_env <- subset(data_env, rownames(data_env) %in% rownames(data_sp))
all.equal(rownames(data_sp), rownames(data_env)) # check again that the two dataframes match in row names

# Convert species abundance dataframe into % of total for each row
cerdata.spp.rel <- decostand(data_sp, method="total")

# Bray Curtis method
cerdata.spp_distmat <- vegdist(cerdata.spp.rel, method="bray")
data_env$binned <- as.factor(data_env$binned)
treat <- c(data_env$binned)

#====================================================================================================
#   PERMANOVA: DISPERSION FOR BRAY-CURTIS
#====================================================================================================
set.seed(425)

# Test if groups are sig diff from each other
test <- adonis2(cerdata.spp_distmat ~ binned, 
               data=data_env,
               method="bray", 
               by="terms",
               permutations=999,
               strata = NULL,
               contr.unordered = "contr.sum",
               contr.ordered = "contr.poly",
               parallel=getOption("mc.cores"))


test
# the larger the F stat, the more sig the results. R2 explains data.
# Assume variation among groups is equal for this test to be valid, test for this below:

# The bray-curtis dist matrix is "cerdata.spp_dist"
# The binned variable groups is in "treat"

# Run betadisper
betadisp <- betadisper(cerdata.spp_distmat, treat)
betadisp
boxplot(betadisp)
plot(betadisp)

# Check the validity of the permanova test
anova_result <- anova(betadisp) # If sig, permanova violated
anova_result
# permutest(betadisp) # this should be similar to anova_result

# If above tests are valid and significant, run Tukey
test.HSD <- TukeyHSD(betadisp)
test.HSD
plot(test.HSD)
