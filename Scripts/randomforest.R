# Random Forest
library(randomForest)
library(caTools)
library(caret)
setwd("~/PostDoc/project1_2021/00-FINAL_DEC102021/Final_Tables")
rm(list=ls())


# Following StatQuest protocol
#====================================================================================================
#   UPLOAD DATA - MAKE DATAFRAMES
#====================================================================================================
FoG = "genus"
string = "Agriculture Percent"  # Variable of interest
variable1 <- gsub(" ","_",tolower(string))

data <- read.csv(paste0("all_",FoG,"_abun_0.1_Dec242021_nozero.csv"), header=TRUE)

# Isolate just variable of interest
if (variable1 == "annual_temperature") {
  names(data)[names(data) == "binned_temp"] <- "binned"
} else if (variable1 == "annual_precipitation") {
  names(data)[names(data) == "binned_precip"] <- "binned" 
} else if (variable1 == "developmental_percent") {
  names(data)[names(data) == "binned_devperc"] <- "binned"
} else {
  names(data)[names(data) == "binned_agriperc"] <- "binned"
}

# Convert GridID names to include variable of interest bin value
data$GridID <- paste(data$GridID, data$binned, sep="_") # Rename GridID with binned variable
row.names(data) <- data$GridID # Push Sites into index column
data$GridID <- NULL # Drop the sites column

# Make dataframes
data_sp <- data[,c(-1:-8)] # Dataframe of just species
data_env <- data[,c(1:8)]  # Dataframe of just environmental variables

#====================================================================================================
#   RANDOM FOREST (FIRST CHECKING DISPERSION OF OUR DATA)
#====================================================================================================
# # Rename the binned column variables
data_env$binned[data_env$binned == "1"] <- "Very_Low"
data_env$binned[data_env$binned == "2"] <- "Low"
data_env$binned[data_env$binned == "3"] <- "Moderate"
data_env$binned[data_env$binned == "4"] <- "High"
data_env$binned[data_env$binned == "5"] <- "Very_High"

#Regression (do not use the binned columns)
data_sp$binned <- data_env$Agriculture_Percent  # for regression, R2 negative = bad model
set.seed(425)
model <- randomForest(binned ~ ., data=data_sp, proximity=TRUE, ntree=1000)
model # if var explained is negative, it is a bad fit

#Classification (use the binned columns)
# Add the binned column to data_sp dataframe
data_sp$binned <- as.factor(data_env$binned)
set.seed(425)
model <- randomForest(binned ~ ., data=data_sp, proximity=TRUE, ntree=1000)
model # if OOB estimate of error rate is above 50%, it is a bad fit
