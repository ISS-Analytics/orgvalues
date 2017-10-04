## Install Packages
# First, install `devtools` once so that you can install SEMinR using it
install.packages("devtools")

# Now `devtools` can install SEMinR from its Github repo
library(devtools)
devtools::install_github("sem-in-r/seminr")

# load seminr library
library(seminr)

## Load the data
orgvalues <- read.csv(file = "BD PLS AGILITY_formatted.csv")

## Replace NA with missing value (Mean Value)

for(i in 1:ncol(orgvalues)){
  orgvalues[is.na(orgvalues[,i]), i] <- mean(orgvalues[,i], na.rm = TRUE)
}

### First estimate the factor scores for the first order constructs

#Create the Measurement Model

orgvalues_mm <- constructs(
  composite("Clan",      c("P3_1","P4_1","P5_1","P6_1","P7_1","P8_1"), weights = "B"),
  composite("Adhocracy", c("P3_2","P4_2","P5_2","P6_2","P7_2","P8_2"), weights = "B"),
  composite("Market",    c("P3_3","P4_3","P5_3","P6_3","P7_3","P8_3"), weights = "B"),
  composite("Hierarchy", c("P3_4","P4_4","P5_4","P6_4","P7_4","P8_4"), weights = "B"),
  composite("AG_Customers", c("P25_4","P25_5","P25_6"), weights = "A"),
  composite("AG_Operational", c("P25_1","P25_2","P25_3"), weights = "A"),
  composite("AG_Partnering", c("P25_7","P25_8","P25_9","P25_10","P25_11"), weights = "A"),
  composite("Age", single_item("Edad"), weights = "A"),
  composite("Size", single_item("Empleados"),weights = "A")
)

#Create the Structural Model

orgvalues_sm <- relationships(
  paths(from = "Clan",      to = c("AG_Customers", "AG_Operational", "AG_Partnering")),
  paths(from = "Adhocracy", to = c("AG_Customers", "AG_Operational", "AG_Partnering")),
  paths(from = "Market",    to = c("AG_Customers", "AG_Operational", "AG_Partnering")),
  paths(from = "Hierarchy", to = c("AG_Customers", "AG_Operational", "AG_Partnering")),
  paths(from = "Size",      to = c("AG_Customers", "AG_Operational", "AG_Partnering")),
  paths(from = "Age",       to = c("AG_Customers", "AG_Operational", "AG_Partnering"))
)

# Estimate the PLS model with Seminr

orgvalues_pls <- estimate_pls(data = orgvalues,
                         measurement_model = orgvalues_mm,
                         structural_model = orgvalues_sm)

# Display R2 and paths of estimated first order model
print_paths(orgvalues_pls)

## store factor scores for First Order Constructs

AG_Operational_1order <- orgvalues_pls$fscores[,"AG_Operational"]
AG_Customers_1order <- orgvalues_pls$fscores[,"AG_Customers"]
AG_Partnering_1order <- orgvalues_pls$fscores[,"AG_Partnering"]

## Add new construct scores to orgvalues dataset

orgvalues2 <- cbind(orgvalues, AG_Operational_1order, AG_Customers_1order, AG_Partnering_1order)

## Write CSV with new construct scores as items
## (already added - don't add multiple columns with same header)
#write.csv(orgvalues2,file = "BD PLS AGILITY_updated_with_constructs.csv")

### Create Second Order Model with First Order constructs as items

## Create the measurement and structural model
#Create the Matrix of the Structural Model

orgvalues_mm2 <- constructs(
  composite("Clan",      c("P3_1","P4_1","P5_1","P6_1","P7_1","P8_1"), weights = "B"),
  composite("Adhocracy", c("P3_2","P4_2","P5_2","P6_2","P7_2","P8_2"), weights = "B"),
  composite("Market",    c("P3_3","P4_3","P5_3","P6_3","P7_3","P8_3"), weights = "B"),
  composite("Hierarchy", c("P3_4","P4_4","P5_4","P6_4","P7_4","P8_4"), weights = "B"),
  composite("Size", single_item("Empleados"),weights = "A"),
  composite("Age", single_item("Edad"), weights = "A"),
  composite("OA", c("AG_Customers_1order","AG_Operational_1order","AG_Partnering_1order"), weights = "A")
)

orgvalues_sm2 <- relationships(
  paths(from = "Clan",      to = c("OA")),
  paths(from = "Adhocracy", to = c("OA")),
  paths(from = "Market",    to = c("OA")),
  paths(from = "Hierarchy", to = c("OA")),
  paths(from = "Size",      to = c("OA")),
  paths(from = "Age",       to = c("OA"))
)

orgvalues2_pls <- estimate_pls(data = orgvalues2,
                         measurement_model = orgvalues_mm2,
                         structural_model = orgvalues_sm2)

print_paths(orgvalues2_pls)

### Perform 10-fold validation of the model
#Load our Algorithm
source("./lib/PLSpredict.R")
source("./lib/validatePredict.R")

#Load Data
set.seed(222)

# Example of a simple partition and prediction on the data:

trainData <- orgvalues2[1:120, ]
testData <- orgvalues2[121:172, ]

orgvalues_predict <- PLSpredict(trainData = trainData,
                                testData = testData,
                                smMatrix = orgvalues_sm2,
                                mmMatrix = orgvalues_mm2)

# Example of running validate predict with 10 fold, 1 rep
orgvalues_metrics <- validatePredict(orgvalues2, orgvalues_sm2, orgvalues_mm2, kfold = 10, reps = 1)

orgvalues_metrics$PLSRMSE
orgvalues_metrics$PLSMAPE
orgvalues_metrics$PLSMAD
orgvalues_metrics$PLSR2

# Example of running validate predict with 10 fold, 20 rep
orgvalues_metrics <- validatePredict(orgvalues2, orgvalues_sm2, orgvalues_mm2, kfold = 10, reps = 20)

orgvalues_metrics$PLSRMSE
orgvalues_metrics$PLSMAPE
orgvalues_metrics$PLSMAD


# Code to create the results for CONSTRUCT LEVEL RMSE 

#Create 10 equally sized folds
folds <- cut(seq(1,nrow(orgvalues2)),breaks=10,labels=FALSE)

# Predicted values
PLS_predicted_outsample_OA <- matrix(,nrow = 172,ncol = 1,dimnames = list(1:172,"OA"))
# In-sample predictions
PLS_predicted_insample_OA <- matrix(,nrow = 172,ncol = 10,dimnames = list(1:172,1:10))

for(x in 1:10) {
  testIndexes <- which(folds==x,arr.ind=TRUE)
  trainIndexes <- which(folds!=x,arr.ind=TRUE)
  testingData <- orgvalues2[testIndexes, ]
  trainingData <- orgvalues2[-testIndexes, ]
  
  #PLS prediction on testset model
  testHolder <- PLSpredict(trainingData, testingData ,orgvalues_sm2, orgvalues_mm2, maxIt = 300, stopCriterion = 7)
  PLS_predicted_outsample_OA[testIndexes,1] <-  testHolder$compositeScores[,7]
  #PLS prediction on trainset model
  trainHolder <- PLSpredict(trainingData, trainingData ,orgvalues_sm2, orgvalues_mm2, maxIt = 300, stopCriterion = 7)
  PLS_predicted_insample_OA[trainIndexes,x] <- trainHolder$compositeScores[,7]
  PLS_predicted_insample_OA[testIndexes,x] <- 0
  ifelse(x == 10, results <- cbind(PLS_predicted_outsample_OA,rowSums(PLS_predicted_insample_OA)/9),0)
         
}

# Add actuals* values
results <- cbind(results,orgvalues2_pls$fscores[,7])

# Name results columns
colnames(results) <- c("Out-of-sample predictions","In-sample predictions","Actuals")

# Write out CSV of results
write.csv(results, file = "results.csv")
OOSRMSE <- sqrt(mean((orgvalues2_pls$fscores[,7] - results[,1])^2))
ISRMSE <- sqrt(mean((orgvalues2_pls$fscores[,7] - results[,2])^2))

